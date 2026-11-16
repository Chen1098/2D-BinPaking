//! Strip packing driver: constructive start, exploration and compression phases.

use crate::layout::*;
use crate::problem::Problem;
use crate::rng::Rng;
use crate::search::{refine, search_placement, Cand, RotCtx, SearchCfg};
use crate::separator::{Limit, Separator};
use std::time::{Duration, Instant};

#[derive(Clone, Debug)]
pub struct Params {
    pub seed: u64,
    pub workers: usize,
    pub time_limit: f64,
    /// deterministic budget in sample evaluations (0 = use the time limit)
    pub max_evals: u64,
    pub explore_frac: f64,
    pub search: SearchCfg,
    pub cmpr_search: SearchCfg,
    pub gls: (f64, f64, f64),
    /// (c0, c1, c2) of the pair loss, see `LossTable`
    pub loss_c: (f64, f64, f64),
    pub loss_mode: u8,
    pub shrink: f64,
    pub ex_iters: usize,
    pub ex_strikes: usize,
    pub pool_std: f64,
    pub large_pct: f64,
    pub cm_shrink: (f64, f64),
    pub cm_iters: usize,
    pub cm_strikes: usize,
    /// slide items into exact contact (left, then down) after every successful separation
    pub compact: bool,
    /// push-relaxation sweeps tried after every shrink before the stochastic search (0 = off)
    pub relax: usize,
    /// push-relaxation sweeps tried on nearly feasible states during separation (0 = off)
    pub finisher: usize,
    /// depth limit (relative to the largest item dimension) for the finisher
    pub finish_pd: f64,
    pub verbose: u8,
}

impl Default for Params {
    fn default() -> Self {
        let search = SearchCfg { n_container: 50, n_focus: 25, n_descents: 3, cd_pre: (0.25, 0.02), cd_fin: (0.01, 0.001) };
        Params {
            seed: 1,
            workers: 3,
            time_limit: 60.0,
            max_evals: 0,
            explore_frac: 0.8,
            search,
            cmpr_search: search,
            gls: (2.0, 1.2, 0.95),
            loss_c: (0.04, 0.5, 2.0),
            loss_mode: 2,
            shrink: 0.001,
            ex_iters: 200,
            ex_strikes: 3,
            pool_std: 0.25,
            large_pct: 0.75,
            cm_shrink: (0.0005, 0.00001),
            cm_iters: 100,
            cm_strikes: 5,
            compact: true,
            relax: 60,
            finisher: 0,
            finish_pd: 0.03,
            verbose: 1,
        }
    }
}

pub struct Outcome {
    pub lay: Layout,
    pub len: f64,
    pub density: f64,
    pub evals: u64,
    pub moves: u64,
    pub iters: u64,
    pub log: Vec<(f64, f64)>,
}

pub fn density(prob: &Problem, len: f64) -> f64 {
    prob.total_area / (prob.height * len)
}

/// Lower bound on the strip length: area bound and the widest unavoidable item.
pub fn len_lower_bound(prob: &Problem) -> f64 {
    let mut lb = prob.total_area / prob.height;
    for t in &prob.types {
        let w = t.shapes.iter().map(|&s| prob.shapes[s].w).fold(f64::INFINITY, f64::min);
        lb = lb.max(w);
    }
    lb
}

/// Bottom-left-fill style constructive heuristic (largest items first).
pub fn construct(prob: &Problem, lt: &LossTable, rng: &mut Rng) -> Layout {
    let n = prob.items.len();
    let mut order: Vec<usize> = (0..n).collect();
    let key = |i: usize| {
        let s = &prob.shapes[prob.types[prob.items[i]].shapes[0]];
        s.hull_area * s.diam
    };
    order.sort_by(|&a, &b| key(b).partial_cmp(&key(a)).unwrap().then(a.cmp(&b)));
    let cfg = SearchCfg { n_container: 1000, n_focus: 0, n_descents: 3, cd_pre: (0.25, 0.02), cd_fin: (0.01, 0.001) };
    let l0 = (prob.total_area / prob.height * 1.25).max(len_lower_bound(prob));
    let mut lay = Layout {
        x: vec![0.0; n],
        y: vec![0.0; n],
        os: prob.items.iter().map(|&t| prob.types[t].shapes[0] as u32).collect(),
        bins: vec![Bin { x0: 0.0, len: l0 }],
    };
    let mut placed: Vec<usize> = Vec::with_capacity(n);
    let mut evals = 0u64;
    for &i in &order {
        let ty = &prob.types[prob.items[i]];
        loop {
            let c = {
                let mut rcs: Vec<RotCtx> = ty
                    .shapes
                    .iter()
                    .map(|&os| RotCtx::new(prob, lt, &lay, os, placed.iter().map(|&j| (j, 1.0)), true))
                    .collect();
                search_placement(prob, &lay, &mut rcs, None, &cfg, rng, &mut evals, &mut [0.0; 8])
            };
            match c {
                Some(c) if c.e.is_finite() => {
                    lay.x[i] = c.x;
                    lay.y[i] = c.y;
                    lay.os[i] = ty.shapes[c.k] as u32;
                    placed.push(i);
                    break;
                }
                _ => lay.bins[0].len *= 1.2,
            }
        }
    }
    lay.bins[0].len = lay.used_len(prob, 0);
    lay
}

/// Exact gravity compaction of a feasible layout: every item (left to right) slides as far
/// left, then down, as it can without colliding, by coordinate descent on its position with
/// the exact collision test.  Returns the resulting used length of bin `b`.
pub fn compact(prob: &Problem, lt: &LossTable, lay: &mut Layout, rng: &mut Rng, sweeps: usize, evals: &mut u64) -> f64 {
    let n = lay.n();
    let mut order: Vec<usize> = (0..n).collect();
    for _ in 0..sweeps {
        order.sort_by(|&a, &b| lay.x[a].partial_cmp(&lay.x[b]).unwrap().then(a.cmp(&b)));
        let mut moved = 0.0;
        for &i in &order {
            let os = lay.os[i] as usize;
            let c = {
                let mut rcs = vec![RotCtx::new(prob, lt, lay, os, (0..n).filter(|&j| j != i).map(|j| (j, 1.0)), true)];
                let md = rcs[0].w.min(rcs[0].h);
                let bin = lay.bin_of(lay.x[i]);
                let e0 = rcs[0].eval(lay.x[i], lay.y[i], f64::INFINITY);
                if !e0.is_finite() {
                    continue;
                }
                refine(&mut rcs, lay, prob, Cand { k: 0, bin, x: lay.x[i], y: lay.y[i], e: e0 }, 0.02 * md, 1e-7 * md, rng, evals)
            };
            moved += (lay.x[i] - c.x).abs() + (lay.y[i] - c.y).abs();
            lay.x[i] = c.x;
            lay.y[i] = c.y;
        }
        if moved < 1e-7 * prob.max_dim {
            break;
        }
    }
    lay.used_len(prob, 0)
}

/// Jam-squeeze a feasible layout: shrink the strip in ever smaller steps, each time letting
/// the push relaxation absorb the overlaps, until nothing moves any more.
/// Leaves the separator at the tightest feasible layout found and returns its length.
pub fn squeeze(sep: &mut Separator, start: &Layout, sweeps: usize, lb: f64) -> (Layout, f64) {
    let prob = sep.prob;
    let mut best = start.clone();
    let mut best_len = best.used_len(prob, 0);
    best.bins[0].len = best_len;
    let mut step = 0.002;
    while step > 2e-6 {
        let next = (best_len * (1.0 - step)).max(lb);
        if next >= best_len {
            break;
        }
        sep.restore(&best);
        // push from the right wall only: nothing is shifted, items are just clamped
        sep.set_bin_len(0, next, f64::INFINITY);
        if sep.relax(sweeps, f64::INFINITY) {
            best = sep.lay.clone();
            best_len = best.used_len(prob, 0);
            best.bins[0].len = best_len;
        } else {
            step *= 0.5;
        }
    }
    sep.restore(&best);
    (best, best_len)
}

/// Swap two large items to escape a local optimum (sparrow's disruption).
fn disrupt(sep: &mut Separator, large_pct: f64) {
    let prob = sep.prob;
    let n = sep.lay.n();
    if n < 2 {
        return;
    }
    let hull = |i: usize| prob.shapes[sep.lay.os[i] as usize].hull_area;
    let total: f64 = (0..n).map(hull).sum();
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&a, &b| hull(b).partial_cmp(&hull(a)).unwrap());
    let mut cum = 0.0;
    let mut cutoff = 0.0;
    for &i in &idx {
        cum += hull(i);
        if cum > total * large_pct {
            cutoff = hull(i);
            break;
        }
    }
    let large: Vec<usize> = (0..n).filter(|&i| hull(i) >= cutoff).collect();
    if large.is_empty() {
        return;
    }
    let a = large[sep.rng.below(large.len())];
    let sa = &prob.shapes[sep.lay.os[a] as usize];
    let diff: Vec<usize> = large
        .iter()
        .cloned()
        .filter(|&j| {
            let sj = &prob.shapes[sep.lay.os[j] as usize];
            j != a && (sj.area - sa.area).abs() > 0.01 * sa.area && (sj.diam - sa.diam).abs() > 0.01 * sa.diam
        })
        .collect();
    let b = if !diff.is_empty() {
        diff[sep.rng.below(diff.len())]
    } else {
        let mut j = sep.rng.below(n - 1);
        if j >= a {
            j += 1;
        }
        j
    };
    let (osa, osb) = (sep.lay.os[a] as usize, sep.lay.os[b] as usize);
    let (sa, sb) = (&prob.shapes[osa], &prob.shapes[osb]);
    // exchange the centroid positions, keep the orientations
    let ca = (sep.lay.x[a] + sa.centroid.x, sep.lay.y[a] + sa.centroid.y);
    let cb = (sep.lay.x[b] + sb.centroid.x, sep.lay.y[b] + sb.centroid.y);
    sep.move_item_to(a, osa, cb.0 - sa.centroid.x, cb.1 - sa.centroid.y);
    sep.move_item_to(b, osb, ca.0 - sb.centroid.x, ca.1 - sb.centroid.y);
}

pub fn solve_strip(prob: &Problem, pr: &Params) -> Outcome {
    let t0 = Instant::now();
    let mut rng = Rng::new(pr.seed);
    let lt0 = LossTable::new(prob, pr.loss_c, 0);
    let init = construct(prob, &lt0, &mut rng);
    let mut log: Vec<(f64, f64)> = Vec::new();
    let lb = len_lower_bound(prob);
    let mut best = init.clone();
    let mut best_len = init.bins[0].len;
    log.push((t0.elapsed().as_secs_f64(), best_len));
    if pr.verbose > 0 {
        eprintln!("[init ] {:7.2}s len {:.4} density {:.3}%", t0.elapsed().as_secs_f64(), best_len, 100.0 * density(prob, best_len));
    }
    let mut sep = Separator::new(prob, init, pr.workers, pr.loss_c, pr.loss_mode, pr.gls, rng.fork());
    sep.finisher = pr.finisher;
    sep.finish_pd = pr.finish_pd * prob.max_dim;
    let far = t0 + Duration::from_secs(86_400 * 365);
    let (lim_explore, lim_end) = if pr.max_evals > 0 {
        (Limit { t: far, evals: (pr.max_evals as f64 * pr.explore_frac) as u64 }, Limit { t: far, evals: pr.max_evals })
    } else {
        (
            Limit { t: t0 + Duration::from_secs_f64(pr.time_limit * pr.explore_frac), evals: u64::MAX },
            Limit { t: t0 + Duration::from_secs_f64(pr.time_limit), evals: u64::MAX },
        )
    };

    // ---------------- exploration ----------------
    let mut cur_len = best_len;
    let mut pool: Vec<(Layout, f64)> = Vec::new();
    let mut cmp_evals = 0u64;
    // first shrink
    let mut next = (cur_len * (1.0 - pr.shrink)).max(lb);
    if next < cur_len {
        sep.set_bin_len(0, next, cur_len * 0.5);
        cur_len = next;
    }
    while !lim_explore.reached(sep.evals()) {
        let it0 = sep.iters;
        let ts = Instant::now();
        let res = sep.separate(&pr.search, pr.ex_iters, pr.ex_strikes, lim_explore);
        if pr.verbose > 2 {
            // diagnostics: remaining colliding pairs and their depths
            let n = res.lay.n();
            let mut npairs = 0;
            let mut maxpd: f64 = 0.0;
            let mut sumpd = 0.0;
            for i in 0..n {
                for j in i + 1..n {
                    let (si, sj) = (res.lay.os[i] as usize, res.lay.os[j] as usize);
                    let (f, sg) = prob.pair(sj, si);
                    let pd = f.pd(sg * (res.lay.x[i] - res.lay.x[j]), sg * (res.lay.y[i] - res.lay.y[j]));
                    if pd > 0.0 {
                        npairs += 1;
                        maxpd = maxpd.max(pd);
                        sumpd += pd;
                    }
                }
            }
            eprintln!(
                "[sep  ] {:7.2}s len {:.4} iters {:5} {:6.2}s loss {:.4e} pairs {} maxpd {:.4e} sumpd {:.4e}",
                t0.elapsed().as_secs_f64(), cur_len, sep.iters - it0, ts.elapsed().as_secs_f64(), res.loss, npairs, maxpd / prob.max_dim, sumpd / prob.max_dim
            );
        }
        if res.loss == 0.0 {
            let mut res = res;
            if pr.compact {
                let l2 = compact(prob, &sep.lt, &mut res.lay, &mut sep.rng, 2, &mut cmp_evals);
                if l2 < cur_len {
                    cur_len = l2;
                    res.lay.bins[0].len = l2;
                }
                sep.restore(&res.lay);
            }
            if pr.relax > 0 {
                let (l, len) = squeeze(&mut sep, &res.lay, pr.relax, lb);
                if len < cur_len {
                    cur_len = len;
                    res.lay = l;
                }
            }
            if cur_len < best_len {
                best_len = cur_len;
                best = res.lay.clone();
                log.push((t0.elapsed().as_secs_f64(), best_len));
                if pr.verbose > 1 {
                    eprintln!("[expl ] {:7.2}s len {:.4} density {:.3}%", t0.elapsed().as_secs_f64(), best_len, 100.0 * density(prob, best_len));
                }
            }
            next = (cur_len * (1.0 - pr.shrink)).max(lb);
            if next >= cur_len {
                break;
            }
            sep.set_bin_len(0, next, cur_len * 0.5);
            cur_len = next;
            pool.clear();
        } else {
            let pos = pool.iter().position(|(_, l)| res.loss < *l).unwrap_or(pool.len());
            pool.insert(pos, (res.lay, res.loss));
            let s = (sep.rng.normal() * pr.pool_std).abs().min(0.999);
            let k = (s * pool.len() as f64) as usize;
            let chosen = pool[k].0.clone();
            sep.restore(&chosen);
            disrupt(&mut sep, pr.large_pct);
        }
    }
    if pr.verbose > 0 {
        eprintln!(
            "[expl ] {:7.2}s len {:.4} density {:.3}%  ({:.2}M evals/s, {:.1}k moves/s)",
            t0.elapsed().as_secs_f64(),
            best_len,
            100.0 * density(prob, best_len),
            sep.evals() as f64 / t0.elapsed().as_secs_f64() / 1e6,
            sep.moves() as f64 / t0.elapsed().as_secs_f64() / 1e3
        );
    }

    if pr.verbose > 1 {
        eprintln!("[fin  ] finisher tried {} succeeded {}", sep.n_finish_try, sep.n_finish_ok);
        let p = sep.prof();
        let tot: f64 = p.iter().sum();
        eprintln!(
            "[prof ] worker-seconds {:.1}: ctx {:.1}% focus {:.1}% container {:.1}% pre-refine {:.1}% final-refine {:.1}% update {:.1}% sync {:.1}%",
            tot, 100.0 * p[0] / tot, 100.0 * p[1] / tot, 100.0 * p[2] / tot, 100.0 * p[3] / tot, 100.0 * p[4] / tot, 100.0 * p[5] / tot, 100.0 * p[6] / tot
        );
    }

    // ---------------- compression ----------------
    let tc0 = Instant::now();
    let ev0 = sep.evals();
    loop {
        if lim_end.reached(sep.evals()) {
            break;
        }
        let ratio = if pr.max_evals > 0 {
            (sep.evals() - ev0) as f64 / (lim_end.evals.saturating_sub(ev0)).max(1) as f64
        } else {
            (Instant::now() - tc0).as_secs_f64() / (lim_end.t - tc0).as_secs_f64().max(1e-9)
        };
        let step = pr.cm_shrink.0 + ratio * (pr.cm_shrink.1 - pr.cm_shrink.0);
        let next = best_len * (1.0 - step);
        if next < lb || next >= best_len {
            break;
        }
        sep.restore(&best);
        let split = sep.rng.f64() * best_len;
        sep.set_bin_len(0, next, split);
        let res = sep.separate(&pr.cmpr_search, pr.cm_iters, pr.cm_strikes, lim_end);
        if res.loss == 0.0 {
            let mut res = res;
            best_len = next;
            if pr.compact {
                let l2 = compact(prob, &sep.lt, &mut res.lay, &mut sep.rng, 2, &mut cmp_evals);
                if l2 < best_len {
                    best_len = l2;
                    res.lay.bins[0].len = l2;
                }
            }
            if pr.relax > 0 {
                let (l, len) = squeeze(&mut sep, &res.lay, pr.relax, lb);
                if len < best_len {
                    best_len = len;
                    res.lay = l;
                }
            }
            best = res.lay;
            log.push((t0.elapsed().as_secs_f64(), best_len));
        }
    }
    if pr.compact {
        let l2 = compact(prob, &sep.lt, &mut best, &mut sep.rng, 20, &mut cmp_evals);
        if l2 < best_len {
            best.bins[0].len = l2;
        }
    }
    // tighten the strip to the items
    best.bins[0].len = best.used_len(prob, 0);
    let best_len = best.bins[0].len;
    log.push((t0.elapsed().as_secs_f64(), best_len));
    if pr.verbose > 0 {
        eprintln!("[final] {:7.2}s len {:.4} density {:.3}%", t0.elapsed().as_secs_f64(), best_len, 100.0 * density(prob, best_len));
    }
    Outcome { density: density(prob, best_len), len: best_len, lay: best, evals: sep.evals(), moves: sep.moves(), iters: sep.iters, log }
}
