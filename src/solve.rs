//! Solver driver: constructive start, exploration (shrink / separate / squeeze) and
//! compression, for strip packing and for bin packing.

use crate::layout::*;
use crate::problem::{Mode, Problem};
use crate::rng::Rng;
use crate::search::{refine, search_placement, snap, Cand, RotCtx, SearchCfg};
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
    /// (c0, c1, c2) of the depth based pair loss, see `LossTable`
    pub loss_c: (f64, f64, f64),
    /// 0: depth loss, 1: pole proxy evaluated directly, 2: pole proxy field (default)
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
    /// sweeps of the push relaxation used to jam-squeeze feasible layouts (0 = off)
    pub relax: usize,
    /// push-relaxation sweeps tried on nearly feasible states during separation (0 = off)
    pub finisher: usize,
    /// depth limit (relative to the largest item dimension) for the finisher
    pub finish_pd: f64,
    /// when two large items are swapped, items covered by them move along
    pub move_covered: bool,
    /// restarts of the strip exploration: a trajectory whose eliminations failed this many
    /// times in a row is abandoned and a new one starts from a fresh construction, as long
    /// as the remaining budget allows it to get as far as the best one did (0 = never)
    pub restart_fails: usize,
    pub verbose: u8,
}

impl Default for Params {
    fn default() -> Self {
        let search = SearchCfg { n_container: 50, n_focus: 25, n_descents: 3, cd_pre: (0.25, 0.02), cd_fin: (0.01, 0.001), snap: true };
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
            finisher: 30,
            finish_pd: 0.01,
            move_covered: true,
            restart_fails: 0,
            verbose: 1,
        }
    }
}

pub struct Outcome {
    pub lay: Layout,
    pub n_bins: usize,
    /// strip length, or the used length of the last (open) bin
    pub len: f64,
    /// item area / container area actually used
    pub density: f64,
    pub evals: u64,
    pub moves: u64,
    pub iters: u64,
    /// (seconds, length of the open bin) of every improvement
    pub log: Vec<(f64, f64)>,
}

/// Item area divided by the used container area (full bins plus the used part of the last).
pub fn density(prob: &Problem, lay: &Layout) -> f64 {
    let k = lay.bins.len();
    let last = lay.bins[k - 1].len;
    let used = match prob.mode {
        Mode::Strip => prob.height * last,
        Mode::Bin => prob.height * (prob.bin_width * (k - 1) as f64 + last),
    };
    prob.total_area / used
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

fn bin_pitch(prob: &Problem) -> f64 {
    prob.bin_width + 2.5 * prob.max_dim
}

/// Bottom-left-fill style constructive heuristic (largest items first).  In bin mode the
/// bins are filled first-fit and a new bin is opened when an item fits nowhere.
pub fn construct(prob: &Problem, lt: &LossTable, rng: &mut Rng, jitter: f64) -> Layout {
    let n = prob.items.len();
    let mut order: Vec<usize> = (0..n).collect();
    // largest first; `jitter` (the sigma of a log-normal factor) perturbs the order
    let keys: Vec<f64> = (0..n)
        .map(|i| {
            let s = &prob.shapes[prob.types[prob.items[i]].shapes[0]];
            s.hull_area * s.diam * if jitter > 0.0 { (jitter * rng.normal()).exp() } else { 1.0 }
        })
        .collect();
    order.sort_by(|&a, &b| keys[b].partial_cmp(&keys[a]).unwrap().then(a.cmp(&b)));
    let cfg = SearchCfg { n_container: 1000, n_focus: 0, n_descents: 3, cd_pre: (0.25, 0.02), cd_fin: (0.01, 0.0001), snap: true };
    let strip = prob.mode == Mode::Strip;
    let l0 = if strip { (prob.total_area / prob.height * 1.25).max(len_lower_bound(prob)) } else { prob.bin_width };
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
        let mut b = 0;
        loop {
            let c = {
                let mut rcs: Vec<RotCtx> = ty
                    .shapes
                    .iter()
                    .map(|&os| RotCtx::new(prob, lt, &lay, os, placed.iter().map(|&j| (j, 1.0)), true))
                    .collect();
                search_placement(prob, &lay, &mut rcs, None, &cfg, Some(b), rng, &mut evals, &mut [0.0; 8])
            };
            match c {
                Some(c) if c.e.is_finite() => {
                    lay.x[i] = c.x;
                    lay.y[i] = c.y;
                    lay.os[i] = ty.shapes[c.k] as u32;
                    placed.push(i);
                    break;
                }
                _ if strip => lay.bins[0].len *= 1.2,
                _ => {
                    b += 1;
                    if b == lay.bins.len() {
                        lay.bins.push(Bin { x0: b as f64 * bin_pitch(prob), len: prob.bin_width });
                    }
                }
            }
        }
    }
    if strip {
        lay.bins[0].len = lay.used_len(prob, 0);
    }
    lay
}

/// Exact gravity compaction of a feasible layout: every item (left to right) slides as far
/// left, then down, as it can without colliding, by coordinate descent on its position with
/// the exact collision test.
pub fn compact(prob: &Problem, lt: &LossTable, lay: &mut Layout, rng: &mut Rng, sweeps: usize, evals: &mut u64) {
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
                let c = refine(&mut rcs, lay, prob, Cand { k: 0, bin, x: lay.x[i], y: lay.y[i], e: e0 }, 0.02 * md, 1e-6 * md, rng, evals);
                snap(&mut rcs, lay, prob, c, evals)
            };
            moved += (lay.x[i] - c.x).abs() + (lay.y[i] - c.y).abs();
            lay.x[i] = c.x;
            lay.y[i] = c.y;
        }
        if moved < 1e-7 * prob.max_dim {
            break;
        }
    }
}

/// Jam-squeeze a feasible layout: shrink bin `b` in ever smaller steps, each time letting the
/// push relaxation absorb the overlaps, until nothing moves any more.
/// Leaves the separator at the tightest feasible layout found and returns it with its length.
pub fn squeeze(sep: &mut Separator, start: &Layout, b: usize, sweeps: usize, lb: f64, min_step: f64) -> (Layout, f64) {
    let prob = sep.prob;
    let mut best = start.clone();
    let mut best_len = best.used_len(prob, b);
    best.bins[b].len = best_len;
    let mut step = 0.002;
    while step > min_step {
        let next = (best_len * (1.0 - step)).max(lb);
        if next >= best_len {
            break;
        }
        sep.restore(&best);
        // push from the right wall only: nothing is shifted, items are just clamped
        sep.set_bin_len(b, next, f64::INFINITY);
        let n_in = |l: &Layout| (0..l.n()).filter(|&i| l.bin_of(l.x[i]) == b).count();
        let before = n_in(&best);
        // (an item evicted from the shortened bin lands on top of others: never "relaxed")
        let ok = sep.relax(sweeps, f64::INFINITY) && n_in(&sep.lay) == before;
        let new_len = if ok { sep.lay.used_len(prob, b) } else { f64::INFINITY };
        if ok && new_len < best_len * (1.0 - 1e-12) {
            best = sep.lay.clone();
            best_len = new_len;
            best.bins[b].len = best_len;
        } else {
            step *= 0.5;
        }
    }
    sep.restore(&best);
    (best, best_len)
}

/// Swap two large items to escape a local optimum (sparrow's disruption).
fn disrupt(sep: &mut Separator, large_pct: f64, move_covered: bool) {
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
    let old = [(sep.lay.x[a], sep.lay.y[a]), (sep.lay.x[b], sep.lay.y[b])];
    sep.move_item_to(a, osa, cb.0 - sa.centroid.x, cb.1 - sa.centroid.y);
    sep.move_item_to(b, osb, ca.0 - sb.centroid.x, ca.1 - sb.centroid.y);
    if !move_covered {
        return;
    }
    // Items that are now practically covered by a swapped item (their centroid lies inside
    // it) move along into the space that item left behind, keeping their relative position.
    for (k, &s) in [a, b].iter().enumerate() {
        let os = sep.lay.os[s] as usize;
        let (nx, ny) = (sep.lay.x[s], sep.lay.y[s]);
        let (dx, dy) = (old[k].0 - nx, old[k].1 - ny);
        let poly = &prob.shapes[os].poly;
        let covered: Vec<usize> = (0..n)
            .filter(|&c| c != a && c != b)
            .filter(|&c| {
                let sc = &prob.shapes[sep.lay.os[c] as usize];
                let q = crate::geom::p(sep.lay.x[c] + sc.centroid.x - nx, sep.lay.y[c] + sc.centroid.y - ny);
                q.x > 0.0 && q.x < prob.shapes[os].w && q.y > 0.0 && q.y < prob.shapes[os].h && crate::geom::point_in_poly(q, poly)
            })
            .collect();
        for c in covered {
            let osc = sep.lay.os[c] as usize;
            sep.move_item_to(c, osc, sep.lay.x[c] + dx, sep.lay.y[c] + dy);
        }
    }
}

/// Keep track of rejected local optima and restart from a disrupted one of the better ones.
fn restart_from_pool(sep: &mut Separator, pool: &mut Vec<(Layout, f64)>, lay: Layout, loss: f64, pr: &Params) {
    let pos = pool.iter().position(|(_, l)| loss < *l).unwrap_or(pool.len());
    pool.insert(pos, (lay, loss));
    if pool.len() > 64 {
        pool.truncate(64);
    }
    let s = (sep.rng.normal() * pr.pool_std).abs().min(0.999);
    let k = (s * pool.len() as f64) as usize;
    let chosen = pool[k].0.clone();
    sep.restore(&chosen);
    disrupt(sep, pr.large_pct, pr.move_covered);
}

struct Ctx<'a> {
    pr: &'a Params,
    t0: Instant,
    log: Vec<(f64, f64)>,
    cmp_evals: u64,
}

impl<'a> Ctx<'a> {
    fn secs(&self) -> f64 {
        self.t0.elapsed().as_secs_f64()
    }
}

/// Tighten a freshly found feasible layout (gravity compaction + jam squeeze of bin `b`).
/// Returns the layout and the length of bin `b`; the separator is left at that layout.
fn tighten(sep: &mut Separator, cx: &mut Ctx, mut lay: Layout, b: usize, lb: f64, min_step: f64, sweeps_mul: usize) -> (Layout, f64) {
    let prob = sep.prob;
    let mut len = lay.bins[b].len;
    if cx.pr.compact {
        compact(prob, &sep.lt, &mut lay, &mut sep.rng, 2, &mut cx.cmp_evals);
        let l2 = lay.used_len(prob, b);
        if l2 < len {
            len = l2;
            lay.bins[b].len = l2;
        }
    }
    if cx.pr.relax > 0 {
        let (l, l2) = squeeze(sep, &lay, b, cx.pr.relax * sweeps_mul, lb, min_step);
        if l2 < len {
            len = l2;
            lay = l;
        }
    }
    sep.restore(&lay);
    (lay, len)
}

/// Exploration phase on bin `b`: repeatedly shrink the bin, eliminate the collisions, tighten.
/// `best` must be feasible; it is replaced whenever a shorter feasible layout is found.
/// Gives up when `max_fails` eliminations in a row have failed (0 = never); returns whether
/// it did, and the evaluation count at its last improvement.
fn explore(sep: &mut Separator, cx: &mut Ctx, b: usize, best: &mut Layout, lim: Limit, lb: f64, max_fails: usize) -> (bool, u64) {
    let prob = sep.prob;
    let pr = cx.pr;
    let mut best_len = best.bins[b].len;
    let mut cur_len = best_len;
    let mut pool: Vec<(Layout, f64)> = Vec::new();
    let mut fails = 0usize;
    let mut ev_last = sep.evals();
    sep.restore(best);
    let mut next = (cur_len * (1.0 - pr.shrink)).max(lb);
    if next >= cur_len {
        return (false, ev_last);
    }
    sep.set_bin_len(b, next, cur_len * 0.5);
    cur_len = next;
    while !lim.reached(sep.evals()) {
        let it0 = sep.iters;
        let res = sep.separate(&pr.search, pr.ex_iters, pr.ex_strikes, lim);
        if pr.verbose > 2 {
            eprintln!(
                "[sep  ] {:7.2}s bin {} len {:.4} iters {:5} loss {:.4e}",
                cx.secs(),
                b,
                cur_len,
                sep.iters - it0,
                res.loss
            );
        }
        if res.loss == 0.0 {
            let (lay, len) = tighten(sep, cx, res.lay, b, lb, 2e-5, 1);
            cur_len = cur_len.min(len);
            if cur_len < best_len {
                best_len = cur_len;
                *best = lay;
                best.bins[b].len = best_len;
                ev_last = sep.evals();
                cx.log.push((cx.secs(), best_len));
                if pr.verbose > 1 {
                    eprintln!("[expl ] {:7.2}s len {:.4} density {:.3}%", cx.secs(), best_len, 100.0 * density(prob, best));
                }
            }
            fails = 0;
            next = (cur_len * (1.0 - pr.shrink)).max(lb);
            if next >= cur_len {
                break;
            }
            sep.set_bin_len(b, next, cur_len * 0.5);
            cur_len = next;
            pool.clear();
        } else {
            fails += 1;
            if max_fails > 0 && fails >= max_fails {
                return (true, ev_last);
            }
            restart_from_pool(sep, &mut pool, res.lay, res.loss, pr);
        }
    }
    (false, ev_last)
}

/// Compression phase on bin `b`: ever smaller shrink steps, always restarting from the best
/// feasible layout.
fn compress(sep: &mut Separator, cx: &mut Ctx, b: usize, best: &mut Layout, lim: Limit, lb: f64) {
    let pr = cx.pr;
    let mut best_len = best.bins[b].len;
    let tc0 = Instant::now();
    let ev0 = sep.evals();
    loop {
        if lim.reached(sep.evals()) {
            break;
        }
        let ratio = if pr.max_evals > 0 {
            (sep.evals() - ev0) as f64 / (lim.evals.saturating_sub(ev0)).max(1) as f64
        } else {
            (Instant::now() - tc0).as_secs_f64() / (lim.t - tc0).as_secs_f64().max(1e-9)
        };
        let step = pr.cm_shrink.0 + ratio * (pr.cm_shrink.1 - pr.cm_shrink.0);
        let next = best_len * (1.0 - step);
        if next < lb || next >= best_len {
            break;
        }
        sep.restore(best);
        let split = sep.rng.f64() * best_len;
        sep.set_bin_len(b, next, split);
        let res = sep.separate(&pr.cmpr_search, pr.cm_iters, pr.cm_strikes, lim);
        if res.loss == 0.0 {
            let (lay, len) = tighten(sep, cx, res.lay, b, lb, 2e-5, 1);
            best_len = next.min(len);
            *best = lay;
            best.bins[b].len = best_len;
            cx.log.push((cx.secs(), best_len));
        }
    }
}

fn limits(pr: &Params, t0: Instant) -> (Limit, Limit) {
    let far = t0 + Duration::from_secs(86_400 * 365);
    if pr.max_evals > 0 {
        (Limit { t: far, evals: (pr.max_evals as f64 * pr.explore_frac) as u64 }, Limit { t: far, evals: pr.max_evals })
    } else {
        (
            Limit { t: t0 + Duration::from_secs_f64(pr.time_limit * pr.explore_frac), evals: u64::MAX },
            Limit { t: t0 + Duration::from_secs_f64(pr.time_limit), evals: u64::MAX },
        )
    }
}

fn new_separator<'a>(prob: &'a Problem, pr: &Params, lay: Layout, rng: &mut Rng) -> Separator<'a> {
    let mut sep = Separator::new(prob, lay, pr.workers, pr.loss_c, pr.loss_mode, pr.gls, rng.fork());
    sep.finisher = pr.finisher;
    sep.finish_pd = pr.finish_pd * prob.max_dim;
    sep
}

fn report_phase(sep: &Separator, cx: &Ctx, prob: &Problem, best: &Layout, tag: &str) {
    if cx.pr.verbose > 0 {
        let k = best.bins.len();
        eprintln!(
            "[{}] {:7.2}s bins {} len {:.4} density {:.3}%  ({:.2}M evals/s, {:.1}k moves/s)",
            tag,
            cx.secs(),
            k,
            best.bins[k - 1].len,
            100.0 * density(prob, best),
            sep.evals() as f64 / cx.secs().max(1e-9) / 1e6,
            sep.moves() as f64 / cx.secs().max(1e-9) / 1e3
        );
    }
    if cx.pr.verbose > 1 {
        let p = sep.prof();
        let tot: f64 = p.iter().sum::<f64>().max(1e-12);
        eprintln!(
            "[prof ] wall: moves {:.2}s bookkeeping {:.2}s relaxation {:.2}s | workers {:.1}s: ctx {:.0}% focus {:.0}% container {:.0}% refine {:.0}%+{:.0}% update {:.0}% sync {:.0}% | finisher {}/{}",
            sep.wall[0],
            sep.wall[1],
            sep.wall[2],
            tot,
            100.0 * p[0] / tot,
            100.0 * p[1] / tot,
            100.0 * p[2] / tot,
            100.0 * p[3] / tot,
            100.0 * p[4] / tot,
            100.0 * p[5] / tot,
            100.0 * p[6] / tot,
            sep.n_finish_ok,
            sep.n_finish_try
        );
    }
}

fn finish(sep: &mut Separator, cx: &mut Ctx, b: usize, best: Layout, lb: f64) -> Outcome {
    let prob = sep.prob;
    let (mut best, _) = tighten(sep, cx, best, b, lb, 1e-7, 4);
    if cx.pr.compact {
        compact(prob, &sep.lt, &mut best, &mut sep.rng, 20, &mut cx.cmp_evals);
    }
    let len = best.used_len(prob, b);
    best.bins[b].len = len;
    cx.log.push((cx.secs(), len));
    if cx.pr.verbose > 0 {
        eprintln!("[final] {:7.2}s bins {} len {:.4} density {:.3}%", cx.secs(), best.bins.len(), len, 100.0 * density(prob, &best));
    }
    Outcome {
        density: density(prob, &best),
        n_bins: best.bins.len(),
        len,
        lay: best,
        evals: sep.evals(),
        moves: sep.moves(),
        iters: sep.iters,
        log: std::mem::take(&mut cx.log),
    }
}

fn solve_strip(prob: &Problem, pr: &Params) -> Outcome {
    let t0 = Instant::now();
    let mut rng = Rng::new(pr.seed);
    let lt0 = LossTable::new(prob, pr.loss_c, 0);
    let mut best = construct(prob, &lt0, &mut rng, 0.0);
    let mut cx = Ctx { pr, t0, log: vec![(t0.elapsed().as_secs_f64(), best.bins[0].len)], cmp_evals: 0 };
    let lb = len_lower_bound(prob);
    if pr.verbose > 0 {
        eprintln!("[init ] {:7.2}s len {:.4} density {:.3}%", cx.secs(), best.bins[0].len, 100.0 * density(prob, &best));
    }
    let mut sep = new_separator(prob, pr, best.clone(), &mut rng);
    let (lim_explore, lim_end) = limits(pr, t0);
    // Restarts.  The descent to a jammed layout takes seconds, the search around it can stay
    // stuck for the rest of the run, and how good a layout a trajectory gets stuck at varies
    // a lot from one trajectory to the next.  So a trajectory that keeps failing is abandoned
    // for a new one from a fresh (perturbed) construction - unless too little budget is left
    // for a new one to get as far as the best one did; then the best one gets the rest.
    let mut n_traj = 0;
    let mut best_age = 0u64;
    loop {
        let ev0 = sep.evals();
        let mut cur = if n_traj == 0 { best.clone() } else { construct(prob, &lt0, &mut rng, 0.2) };
        let (stalled, ev_impr) = explore(&mut sep, &mut cx, 0, &mut cur, lim_explore, lb, pr.restart_fails);
        if pr.verbose > 1 && pr.restart_fails > 0 {
            eprintln!("[traj ] {:7.2}s trajectory {} ended at {:.4} (best before: {:.4})", cx.secs(), n_traj + 1, cur.bins[0].len, best.bins[0].len);
        }
        if n_traj == 0 || cur.bins[0].len < best.bins[0].len {
            best = cur;
            best_age = ev_impr - ev0;
        }
        n_traj += 1;
        if !stalled || lim_explore.reached(sep.evals()) {
            break;
        }
        let left = if pr.max_evals > 0 {
            lim_explore.evals.saturating_sub(sep.evals()) as f64
        } else {
            lim_explore.t.saturating_duration_since(Instant::now()).as_secs_f64() * sep.evals() as f64 / cx.secs().max(1e-9)
        };
        if left < best_age as f64 {
            explore(&mut sep, &mut cx, 0, &mut best, lim_explore, lb, 0);
            break;
        }
    }
    // the history keeps the improvements of the overall best only
    let mut shortest = f64::INFINITY;
    cx.log.retain(|&(_, l)| {
        let keep = l < shortest;
        if keep {
            shortest = l;
        }
        keep
    });
    if pr.verbose > 0 && n_traj > 1 {
        eprintln!("[traj ] {} trajectories", n_traj);
    }
    report_phase(&sep, &cx, prob, &best, "expl ");
    compress(&mut sep, &mut cx, 0, &mut best, lim_end, lb);
    finish(&mut sep, &mut cx, 0, best, lb)
}

/// Item area per bin.
fn bin_fill(prob: &Problem, lay: &Layout) -> Vec<f64> {
    let mut f = vec![0.0; lay.bins.len()];
    for i in 0..lay.n() {
        f[lay.bin_of(lay.x[i])] += prob.shapes[lay.os[i] as usize].area;
    }
    f
}

/// Move the items of bin `from` to position `to` in the bin order (and vice versa).
fn swap_bins(lay: &mut Layout, a: usize, b: usize) {
    if a == b {
        return;
    }
    let bins: Vec<usize> = (0..lay.n()).map(|i| lay.bin_of(lay.x[i])).collect();
    let (xa, xb) = (lay.bins[a].x0, lay.bins[b].x0);
    for i in 0..lay.n() {
        if bins[i] == a {
            lay.x[i] += xb - xa;
        } else if bins[i] == b {
            lay.x[i] += xa - xb;
        }
    }
    let (la, lb) = (lay.bins[a].len, lay.bins[b].len);
    lay.bins[a].len = lb;
    lay.bins[b].len = la;
}

/// Try to pack everything into one bin less: the emptiest bin is dissolved, its items are
/// dropped into the other bins and the collisions are eliminated.
fn try_eliminate(sep: &mut Separator, cx: &mut Ctx, best: &Layout, lim: Limit) -> Option<Layout> {
    let prob = sep.prob;
    let pr = cx.pr;
    let k = best.bins.len();
    if k <= 1 {
        return None;
    }
    let mut fill = bin_fill(prob, best);
    let victim = (0..k).min_by(|&a, &b| fill[a].partial_cmp(&fill[b]).unwrap()).unwrap();
    let mut lay = best.clone();
    swap_bins(&mut lay, victim, k - 1);
    fill.swap(victim, k - 1);
    let cap = prob.bin_width * prob.height;
    let bins: Vec<usize> = (0..lay.n()).map(|i| lay.bin_of(lay.x[i])).collect();
    let mut movers: Vec<usize> = (0..lay.n()).filter(|&i| bins[i] == k - 1).collect();
    movers.sort_by(|&a, &b| {
        let (sa, sb) = (&prob.shapes[lay.os[a] as usize], &prob.shapes[lay.os[b] as usize]);
        sb.area.partial_cmp(&sa.area).unwrap()
    });
    lay.bins.pop();
    for &i in &movers {
        let s = &prob.shapes[lay.os[i] as usize];
        // the bin with the most free area
        let t = (0..k - 1).max_by(|&a, &b| (cap - fill[a]).partial_cmp(&(cap - fill[b])).unwrap()).unwrap();
        fill[t] += s.area;
        let bin = lay.bins[t];
        lay.x[i] = sep.rng.range(bin.x0, bin.x0 + bin.len - s.w);
        lay.y[i] = sep.rng.range(0.0, prob.height - s.h);
    }
    sep.restore(&lay);
    let mut pool: Vec<(Layout, f64)> = Vec::new();
    loop {
        let res = sep.separate(&pr.search, pr.ex_iters, pr.ex_strikes, lim);
        if res.loss == 0.0 {
            return Some(res.lay);
        }
        if lim.reached(sep.evals()) {
            break;
        }
        restart_from_pool(sep, &mut pool, res.lay, res.loss, pr);
    }
    sep.restore(best);
    None
}

fn solve_bins(prob: &Problem, pr: &Params) -> Outcome {
    let t0 = Instant::now();
    let mut rng = Rng::new(pr.seed);
    let lt0 = LossTable::new(prob, pr.loss_c, 0);
    let mut best = construct(prob, &lt0, &mut rng, 0.0);
    let cap = prob.bin_width * prob.height;
    let k_lb = ((prob.total_area / cap - 1e-9).ceil() as usize).max(1);
    let mut cx = Ctx { pr, t0, log: Vec::new(), cmp_evals: 0 };
    if pr.verbose > 0 {
        eprintln!("[init ] {:7.2}s bins {} (area lower bound {})", cx.secs(), best.bins.len(), k_lb);
    }
    let mut sep = new_separator(prob, pr, best.clone(), &mut rng);
    let (lim_explore, lim_end) = limits(pr, t0);

    // ---- phase A: use fewer bins
    while best.bins.len() > k_lb && !lim_explore.reached(sep.evals()) {
        // every attempt may use up to 40% of what is left of the exploration budget
        let lim_a = if pr.max_evals > 0 {
            Limit { t: lim_explore.t, evals: sep.evals() + (lim_explore.evals.saturating_sub(sep.evals()) as f64 * 0.4) as u64 }
        } else {
            let now = Instant::now();
            Limit { t: now + (lim_explore.t.saturating_duration_since(now)).mul_f64(0.4), evals: u64::MAX }
        };
        match try_eliminate(&mut sep, &mut cx, &best, lim_a) {
            Some(lay) => {
                best = lay;
                if pr.verbose > 0 {
                    eprintln!("[bins ] {:7.2}s packed into {} bins", cx.secs(), best.bins.len());
                }
            }
            None => break,
        }
    }

    // ---- phase B: shorten the emptiest ("open") bin; drop it if it runs empty
    let mut b;
    loop {
        let k = best.bins.len();
        let fill = bin_fill(prob, &best);
        let open = (0..k).min_by(|&a, &c| fill[a].partial_cmp(&fill[c]).unwrap()).unwrap();
        swap_bins(&mut best, open, k - 1);
        b = k - 1;
        let lb = 1e-6 * prob.bin_width;
        let (lay, len) = tighten(&mut sep, &mut cx, best.clone(), b, lb, 2e-5, 1);
        best = lay;
        best.bins[b].len = len.min(prob.bin_width);
        cx.log.push((cx.secs(), best.bins[b].len));
        explore(&mut sep, &mut cx, b, &mut best, lim_explore, lb, 0);
        let empty = (0..best.n()).all(|i| best.bin_of(best.x[i]) != b);
        if empty && k > 1 {
            best.bins.pop();
            if pr.verbose > 0 {
                eprintln!("[bins ] {:7.2}s open bin emptied, {} bins left", cx.secs(), best.bins.len());
            }
            if lim_explore.reached(sep.evals()) {
                b = best.bins.len() - 1;
                break;
            }
            continue;
        }
        break;
    }
    report_phase(&sep, &cx, prob, &best, "expl ");
    let lb = 1e-6 * prob.bin_width;
    if best.bins[b].len > lb {
        compress(&mut sep, &mut cx, b, &mut best, lim_end, lb);
    }
    finish(&mut sep, &mut cx, b, best, lb)
}

pub fn solve(prob: &Problem, pr: &Params) -> Outcome {
    match prob.mode {
        Mode::Strip => solve_strip(prob, pr),
        Mode::Bin => solve_bins(prob, pr),
    }
}
