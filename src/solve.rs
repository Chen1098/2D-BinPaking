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
    /// adaptive shrink step of the exploration: after this many consecutive failed
    /// eliminations at a target the step is halved (down to `shrink_min`) and a fresh attempt
    /// is made from the best layout, cut at a random position; every success doubles the
    /// step again (0 = fixed step)
    pub adapt_fails: usize,
    pub shrink_min: f64,
    /// cut the layout at a random position (instead of the middle) when the bin is shrunk
    pub split_random: bool,
    /// probability that a failed elimination is followed by a fresh attempt from the best
    /// layout with a new random cut (instead of a disrupted layout from the pool)
    pub recut: f64,
    /// futility cut-off of the exploration: once `cut_min` eliminations have failed at a
    /// target, a further attempt is given up after `cut_iter` iterations unless its loss is
    /// by then below the `cut_q` quantile of the losses those failures ended with (0 = off)
    /// tighten (compact + squeeze) every feasible layout found by the exploration; with
    /// `false` the exploration works on the layouts as the eliminations leave them and only
    /// its final result is tightened
    pub tighten_explore: bool,
    pub cut_min: usize,
    pub cut_iter: usize,
    pub cut_q: f64,
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
    /// number of random swaps tried per disruption (the best fitting one is applied)
    pub swap_cands: usize,
    /// restart rule of the strip exploration: a trajectory whose eliminations failed this
    /// many times in a row is parked and a new one is started from a fresh construction, as
    /// long as the remaining budget allows it to get as far (0 = never restart)
    pub restart_fails: usize,
    /// basin hopping: after this many consecutive failed eliminations at a target the search
    /// hops to another local optimum (0 = off, i.e. keep trying at the same target)
    pub hop_fails: usize,
    /// a hop loosens the bin by this fraction of its length before it perturbs the layout
    pub hop_loosen: f64,
    /// a hop is accepted as the new starting point when the local optimum it leads to is at
    /// most this fraction longer than the best layout
    pub hop_tau: f64,
    /// number of swaps per hop
    pub hop_swaps: usize,
    /// how a hop loosens the layout: 0 = a gap at a random position, 1 = uniform dilation
    pub hop_mode: u8,
    /// tighten the layouts of the descents that follow a hop (false: only their end points)
    pub hop_tight: bool,
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
            adapt_fails: 0,
            shrink_min: 0.00002,
            split_random: false,
            recut: 0.0,
            tighten_explore: true,
            cut_min: 0,
            cut_iter: 50,
            cut_q: 0.5,
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
            swap_cands: 1,
            restart_fails: 0,
            hop_fails: 0,
            hop_loosen: 0.03,
            hop_tau: 0.01,
            hop_swaps: 1,
            hop_mode: 0,
            hop_tight: true,
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
    // largest first; `jitter` perturbs the order (used to diversify restarts)
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

/// The "large" items: those that make up the top `large_pct` of the total convex hull area.
fn large_items(sep: &Separator, large_pct: f64) -> Vec<usize> {
    let prob = sep.prob;
    let n = sep.lay.n();
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
    (0..n).filter(|&i| hull(i) >= cutoff).collect()
}

/// Two large items of different shape (any second item if there is no such pair).
fn pick_pair(sep: &mut Separator, large: &[usize]) -> (usize, usize) {
    let prob = sep.prob;
    let n = sep.lay.n();
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
    (a, b)
}

/// Exchange the positions of items `a` and `b` (orientations are kept).  Items that end up
/// practically covered by a swapped item (their centroid lies inside it) move along into the
/// space that item left behind, keeping their relative position.
fn swap_items(sep: &mut Separator, a: usize, b: usize, move_covered: bool) {
    let prob = sep.prob;
    let n = sep.lay.n();
    let (osa, osb) = (sep.lay.os[a] as usize, sep.lay.os[b] as usize);
    let (sa, sb) = (&prob.shapes[osa], &prob.shapes[osb]);
    let ca = (sep.lay.x[a] + sa.centroid.x, sep.lay.y[a] + sa.centroid.y);
    let cb = (sep.lay.x[b] + sb.centroid.x, sep.lay.y[b] + sb.centroid.y);
    let old = [(sep.lay.x[a], sep.lay.y[a]), (sep.lay.x[b], sep.lay.y[b])];
    sep.move_item_to(a, osa, cb.0 - sa.centroid.x, cb.1 - sa.centroid.y);
    sep.move_item_to(b, osb, ca.0 - sb.centroid.x, ca.1 - sb.centroid.y);
    if !move_covered {
        return;
    }
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

/// Swap two large items to escape a local optimum (sparrow's disruption).  With
/// `n_cands > 1` several random swaps are tried and the one that fits best right away (the
/// smallest total collision loss) is kept.
fn disrupt(sep: &mut Separator, large_pct: f64, move_covered: bool, n_cands: usize) {
    if sep.lay.n() < 2 {
        return;
    }
    let large = large_items(sep, large_pct);
    if large.is_empty() {
        return;
    }
    if n_cands <= 1 {
        let (a, b) = pick_pair(sep, &large);
        swap_items(sep, a, b, move_covered);
        return;
    }
    let lay0 = sep.lay.clone();
    let mat0 = sep.tr.loss.clone();
    let mut best = (f64::INFINITY, 0, 0);
    for _ in 0..n_cands {
        let (a, b) = pick_pair(sep, &large);
        swap_items(sep, a, b, move_covered);
        let l = sep.tr.total_loss();
        if l < best.0 {
            best = (l, a, b);
        }
        sep.lay.clone_from(&lay0);
        sep.tr.loss.copy_from_slice(&mat0);
    }
    swap_items(sep, best.1, best.2, move_covered);
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
    disrupt(sep, pr.large_pct, pr.move_covered, pr.swap_cands);
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
/// Stops early when the eliminations failed `max_fails` times in a row (0 = never) and
/// returns whether it did so, together with the evaluation count at the last improvement.
///
/// With `Params::hop_fails > 0` the phase is a basin-hopping search over jammed layouts: when
/// the eliminations at a target keep failing, the layout the descent stands on is loosened,
/// perturbed by a swap and tightened again, which leads to a different local optimum; the
/// search moves there if it is not much longer than the best layout (threshold acceptance).
fn explore(sep: &mut Separator, cx: &mut Ctx, b: usize, best: &mut Layout, lim: Limit, lb: f64, max_fails: usize) -> (bool, u64) {
    let prob = sep.prob;
    let pr = cx.pr;
    let mut best_len = best.bins[b].len;
    // length of the bin the separator currently works on
    let mut cur_len = best_len;
    let mut pool: Vec<(Layout, f64)> = Vec::new();
    // the (sorted) losses that the failed eliminations at the current target ended with
    let mut finals: Vec<f64> = Vec::new();
    // current shrink step
    let mut r = pr.shrink;
    // number of consecutive failures at the current target / since the last improvement that
    // was worth a full shrink step (`ref_len` is the best length at that time)
    let mut fails = 0usize;
    let mut idle_fails = 0usize;
    let mut ref_len = best_len;
    // basin hopping: the accepted local optimum hops start from, and the feasible layout of
    // the descent in progress
    let mut home = best.clone();
    let mut home_len = best_len;
    let mut desc: Option<(Layout, f64)> = None;
    let max_len = if prob.mode == Mode::Bin { prob.bin_width } else { f64::INFINITY };
    let mut n_hops = (0usize, 0usize);
    // whether the descent in progress tightens the layouts it finds
    let mut tight = pr.tighten_explore;
    sep.restore(best);
    let mut next = (cur_len * (1.0 - pr.shrink)).max(lb);
    let mut ev_last = sep.evals();
    if next >= cur_len {
        return (false, ev_last);
    }
    let split = if pr.split_random { sep.rng.f64() * cur_len } else { cur_len * 0.5 };
    sep.set_bin_len(b, next, split);
    cur_len = next;
    while !lim.reached(sep.evals()) {
        let it0 = sep.iters;
        let futile = if pr.cut_min > 0 && finals.len() >= pr.cut_min {
            Some((pr.cut_iter, finals[((pr.cut_q * finals.len() as f64) as usize).min(finals.len() - 1)]))
        } else {
            None
        };
        let n_cut0 = sep.n_cut;
        let res = sep.separate(&pr.search, pr.ex_iters, pr.ex_strikes, lim, futile);
        if pr.verbose > 2 {
            eprintln!(
                "[sep  ] {:7.2}s bin {} len {:.4} iters {:5} loss {:.4e} curve {:.3e} {:.3e} {:.3e} {:.3e} {:.3e} {:.3e} {:.3e}",
                cx.secs(),
                b,
                cur_len,
                sep.iters - it0,
                res.loss,
                sep.curve[0],
                sep.curve[1],
                sep.curve[2],
                sep.curve[3],
                sep.curve[4],
                sep.curve[5],
                sep.curve[6]
            );
        }
        if res.loss == 0.0 {
            let (mut lay, len) = if tight { tighten(sep, cx, res.lay, b, lb, 2e-5, 1) } else { (res.lay, cur_len) };
            let len = cur_len.min(len);
            lay.bins[b].len = len;
            if len < best_len {
                best_len = len;
                best.clone_from(&lay);
                if best_len <= ref_len * (1.0 - pr.shrink) {
                    ref_len = best_len;
                    ev_last = sep.evals();
                    idle_fails = 0;
                }
                cx.log.push((cx.secs(), best_len));
                if pr.verbose > 1 {
                    eprintln!("[expl ] {:7.2}s len {:.4} density {:.3}%", cx.secs(), best_len, 100.0 * density(prob, best));
                }
            }
            fails = 0;
            r = (r * 2.0).min(pr.shrink);
            next = (len * (1.0 - r)).max(lb);
            if next >= len {
                return (false, ev_last);
            }
            let split = if pr.split_random { sep.rng.f64() * len } else { len * 0.5 };
            sep.set_bin_len(b, next, split);
            cur_len = next;
            desc = Some((lay, len));
            pool.clear();
            finals.clear();
        } else {
            if sep.n_cut == n_cut0 {
                let pos = finals.partition_point(|&l| l < res.loss);
                finals.insert(pos, res.loss);
            }
            fails += 1;
            idle_fails += 1;
            if max_fails > 0 && idle_fails >= max_fails {
                return (true, ev_last);
            }
            if pr.hop_fails > 0 && fails >= pr.hop_fails {
                // The descent has arrived at a local optimum.  Move there if it is good enough...
                if let Some((mut lay, mut len)) = desc.take() {
                    if !tight {
                        // the descent did not tighten its layouts: do it now for its end point
                        let (l2, len2) = tighten(sep, cx, lay, b, lb, 2e-5, 1);
                        lay = l2;
                        len = len.min(len2);
                        lay.bins[b].len = len;
                        if len < best_len {
                            best_len = len;
                            best.clone_from(&lay);
                            cx.log.push((cx.secs(), best_len));
                            if pr.verbose > 1 {
                                eprintln!("[expl ] {:7.2}s len {:.4} density {:.3}%", cx.secs(), best_len, 100.0 * density(prob, best));
                            }
                        }
                    }
                    n_hops.0 += 1;
                    if len <= best_len * (1.0 + pr.hop_tau) {
                        n_hops.1 += 1;
                        home = lay;
                        home_len = len;
                    }
                }
                // ...and hop: loosen the layout, swap large items, descend again.
                sep.restore(&home);
                let loose = (home_len * (1.0 + pr.hop_loosen)).min(max_len);
                if pr.hop_mode == 1 {
                    sep.dilate_bin(b, loose / home_len);
                } else {
                    let split = sep.rng.f64() * home_len;
                    sep.set_bin_len(b, loose, split);
                }
                for _ in 0..pr.hop_swaps {
                    disrupt(sep, pr.large_pct, pr.move_covered, pr.swap_cands);
                }
                tight = pr.hop_tight;
                cur_len = loose;
                fails = 0;
                pool.clear();
            finals.clear();
            } else if pr.adapt_fails > 0 && fails >= pr.adapt_fails && r * 0.5 >= pr.shrink_min {
                // The target is out of reach of the arrangements at hand: ask for less.  A
                // fresh attempt from the best layout, cut at a random position, with half
                // the step.  Whatever it finds is tightened, which usually gains more than
                // the step that was asked for.
                r *= 0.5;
                fails = 0;
                sep.restore(best);
                next = (best_len * (1.0 - r)).max(lb);
                if next >= best_len {
                    return (false, ev_last);
                }
                let split = sep.rng.f64() * best_len;
                sep.set_bin_len(b, next, split);
                cur_len = next;
                pool.clear();
            finals.clear();
            } else if pr.recut > 0.0 && pr.hop_fails == 0 && sep.rng.f64() < pr.recut {
                // a fresh attempt at the same target from the best layout, cut elsewhere
                sep.restore(best);
                let split = sep.rng.f64() * best_len;
                sep.set_bin_len(b, cur_len, split);
            } else {
                restart_from_pool(sep, &mut pool, res.lay, res.loss, pr);
            }
        }
    }
    if !tight {
        // the end point of an untightened descent that was cut short by the limit
        if let Some((lay, len)) = desc.take() {
            let (mut l2, len2) = tighten(sep, cx, lay, b, lb, 2e-5, 1);
            let len = len.min(len2);
            if len < best_len {
                l2.bins[b].len = len;
                *best = l2;
                cx.log.push((cx.secs(), len));
            }
        }
    }
    if pr.verbose > 1 && pr.cut_min > 0 {
        eprintln!("[cut  ] {} eliminations given up early", sep.n_cut);
    }
    if pr.verbose > 1 && pr.hop_fails > 0 {
        eprintln!("[hops ] {} local optima visited, {} accepted", n_hops.0, n_hops.1);
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
        let res = sep.separate(&pr.cmpr_search, pr.cm_iters, pr.cm_strikes, lim, None);
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
    // stuck for a long time, and how good a jammed layout a trajectory ends up at varies a lot
    // from one trajectory to the next.  So a trajectory that stopped making progress is
    // parked and a new one is started from a fresh (perturbed) construction - unless there
    // is not enough budget left for a new one to get as far as the best one did, in which
    // case the best one gets the rest.
    let mut n_traj = 0;
    let mut best_age = 0u64;
    loop {
        let ev0 = sep.evals();
        let mut cur = if n_traj == 0 { best.clone() } else { construct(prob, &lt0, &mut rng, 0.2) };
        let (stalled, ev_impr) = explore(&mut sep, &mut cx, 0, &mut cur, lim_explore, lb, pr.restart_fails);
        if pr.verbose > 1 && pr.restart_fails > 0 {
            eprintln!(
                "[traj ] {:7.2}s trajectory {} ended at {:.4} (best before {:.4})",
                cx.secs(),
                n_traj + 1,
                cur.bins[0].len,
                best.bins[0].len
            );
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
    // keep only the improvements of the overall best in the history
    let mut m = f64::INFINITY;
    cx.log.retain(|&(_, l)| {
        let keep = l < m;
        if keep {
            m = l;
        }
        keep
    });
    if pr.verbose > 0 && n_traj > 1 {
        eprintln!("[traj ] {} trajectories", n_traj);
    }
    if !pr.tighten_explore {
        let len0 = best.bins[0].len;
        let (lay, len) = tighten(&mut sep, &mut cx, best, 0, lb, 2e-5, 1);
        best = lay;
        best.bins[0].len = len.min(len0);
        cx.log.push((cx.secs(), best.bins[0].len));
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
        let res = sep.separate(&pr.search, pr.ex_iters, pr.ex_strikes, lim, None);
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
