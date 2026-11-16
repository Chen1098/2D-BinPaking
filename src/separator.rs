//! Collision elimination by guided local search ("separation").
//!
//! The overall scheme (move every colliding item to its best sampled position, raise the
//! weights of pairs that keep colliding, several workers race from the same state and the
//! best one wins) follows Gardeyn, Vanden Berghe & Wauters, "An open-source heuristic to
//! reboot 2D nesting research" (EJOR 2026); the evaluation underneath is the exact NFP field.

use crate::layout::*;
use crate::problem::Problem;
use crate::rng::Rng;
use crate::search::{search_placement, RotCtx, SearchCfg};
use crate::pool::SpinPool;
use std::time::Instant;

pub struct Worker {
    pub lay: Layout,
    pub loss: Vec<f64>,
    pub rng: Rng,
    pub evals: u64,
    pub moves: u64,
    /// profile: ctx build, focus, container, pre-refine, final refine, loss update, sync
    pub prof: [f64; 8],
    order: Vec<usize>,
}

impl Worker {
    fn load(&mut self, lay: &Layout, loss: &[f64]) {
        let t = Instant::now();
        self.lay.clone_from(lay);
        self.loss.clear();
        self.loss.extend_from_slice(loss);
        self.prof[6] += t.elapsed().as_secs_f64();
    }

    /// Give every colliding item (in random order) the chance to move to a better placement.
    fn move_items(&mut self, prob: &Problem, lt: &LossTable, wt: &[f64], cfg: &SearchCfg) {
        let n = self.lay.n();
        self.order.clear();
        for i in 0..n {
            if item_loss(&self.loss, n, i) > 0.0 {
                self.order.push(i);
            }
        }
        self.rng.shuffle(&mut self.order);
        for idx in 0..self.order.len() {
            let i = self.order[idx];
            if item_loss(&self.loss, n, i) > 0.0 {
                self.move_item(prob, lt, wt, cfg, i);
            }
        }
    }

    fn move_item(&mut self, prob: &Problem, lt: &LossTable, wt: &[f64], cfg: &SearchCfg, i: usize) {
        let n = self.lay.n();
        let ty = &prob.types[prob.items[i]];
        let cur_os = self.lay.os[i] as usize;
        let k0 = ty.shapes.iter().position(|&s| s == cur_os).unwrap_or(0);
        let t0 = Instant::now();
        let c = {
            let lay = &self.lay;
            let mut rcs: Vec<RotCtx> = ty
                .shapes
                .iter()
                .map(|&os| RotCtx::new(prob, lt, lay, os, (0..n).filter(|&j| j != i).map(|j| (j, wt[i * n + j])), false))
                .collect();
            self.prof[0] += t0.elapsed().as_secs_f64();
            search_placement(prob, lay, &mut rcs, Some((k0, lay.x[i], lay.y[i])), cfg, None, &mut self.rng, &mut self.evals, &mut self.prof)
        };
        let t1 = Instant::now();
        if let Some(c) = c {
            self.lay.x[i] = c.x;
            self.lay.y[i] = c.y;
            self.lay.os[i] = ty.shapes[c.k] as u32;
            recompute_item(&mut self.loss, prob, lt, &self.lay, i);
            self.moves += 1;
        }
        self.prof[5] += t1.elapsed().as_secs_f64();
    }
}

pub struct Separator<'a> {
    pub prob: &'a Problem,
    pub lt: LossTable,
    pub lay: Layout,
    pub tr: Tracker,
    pub workers: Vec<Worker>,
    pub rng: Rng,
    pool: SpinPool,
    pub gls: (f64, f64, f64),
    pub iters: u64,
    /// push-relaxation sweeps tried on promising (nearly feasible) states; 0 = off
    pub finisher: usize,
    /// only states whose deepest overlap is below this depth are worth finishing
    pub finish_pd: f64,
    /// work done outside sample evaluations (relaxation queries), in query units
    pub aux: u64,
    pub n_finish_try: u64,
    pub n_finish_ok: u64,
    /// wall-clock profile: parallel move phase, serial bookkeeping, relaxation
    pub wall: [f64; 4],
}

pub struct SepResult {
    pub lay: Layout,
    pub loss: f64,
}

/// Stopping criterion: wall clock and/or a deterministic evaluation budget.
#[derive(Clone, Copy, Debug)]
pub struct Limit {
    pub t: Instant,
    pub evals: u64,
}

impl Limit {
    #[inline]
    pub fn reached(&self, evals: u64) -> bool {
        evals >= self.evals || Instant::now() >= self.t
    }
}

impl<'a> Separator<'a> {
    pub fn new(prob: &'a Problem, lay: Layout, n_workers: usize, loss_c: (f64, f64, f64), loss_mode: u8, gls: (f64, f64, f64), mut rng: Rng) -> Self {
        let lt = LossTable::new(prob, loss_c, loss_mode);
        let tr = Tracker::new(prob, &lt, &lay);
        let workers = (0..n_workers.max(1))
            .map(|_| Worker { lay: lay.clone(), loss: tr.loss.clone(), rng: rng.fork(), evals: 0, moves: 0, prof: [0.0; 8], order: Vec::new() })
            .collect();
        let pool = SpinPool::new(n_workers.max(1));
        Separator { prob, lt, lay, tr, workers, rng, pool, gls, iters: 0, finisher: 0, finish_pd: 0.0, aux: 0, n_finish_try: 0, n_finish_ok: 0, wall: [0.0; 4] }
    }

    /// Work counter used for budgets: sample evaluations plus relaxation work
    /// (five field queries of the relaxation count as one sample evaluation).
    pub fn evals(&self) -> u64 {
        self.workers.iter().map(|w| w.evals).sum::<u64>() + self.aux / 5
    }
    pub fn moves(&self) -> u64 {
        self.workers.iter().map(|w| w.moves).sum()
    }
    pub fn prof(&self) -> [f64; 8] {
        let mut p = [0.0; 8];
        for w in &self.workers {
            for k in 0..8 {
                p[k] += w.prof[k];
            }
        }
        p
    }

    /// Replace the layout; losses are recomputed and the weights reset.
    pub fn restore(&mut self, lay: &Layout) {
        self.lay.clone_from(lay);
        self.tr.recompute_all(self.prob, &self.lt, &self.lay);
        self.tr.reset_weights();
    }

    pub fn move_item_to(&mut self, i: usize, os: usize, x: f64, y: f64) {
        let mut b = self.lay.bin_of(x);
        if self.prob.shapes[os].w > self.lay.bins[b].len {
            // does not fit into the target bin (a shortened open bin): stay in the old one
            b = self.lay.bin_of(self.lay.x[i]);
        }
        let (x, y) = self.lay.clamp_to_bin(self.prob, os, b, x, y);
        self.lay.os[i] = os as u32;
        self.lay.x[i] = x;
        self.lay.y[i] = y;
        recompute_item(&mut self.tr.loss, self.prob, &self.lt, &self.lay, i);
    }

    /// Change the length of bin `b`; items whose centroid lies right of `split` (relative to
    /// the bin) are shifted along, everything is clamped back into the bin.
    pub fn set_bin_len(&mut self, b: usize, new_len: f64, split: f64) {
        let old = self.lay.bins[b].len;
        let delta = new_len - old;
        let x0 = self.lay.bins[b].x0;
        self.lay.bins[b].len = new_len;
        for i in 0..self.lay.n() {
            if self.lay.bin_of(self.lay.x[i]) != b {
                continue;
            }
            let s = &self.prob.shapes[self.lay.os[i] as usize];
            if self.lay.x[i] + s.centroid.x - x0 > split {
                self.lay.x[i] += delta;
            }
            let (x, y) = self.lay.clamp_to_bin(self.prob, self.lay.os[i] as usize, b, self.lay.x[i], self.lay.y[i]);
            self.lay.x[i] = x;
            self.lay.y[i] = y;
            if s.w > new_len {
                self.evict(i, b);
            }
        }
        self.tr.recompute_all(self.prob, &self.lt, &self.lay);
        self.tr.reset_weights();
    }

    /// Item `i` no longer fits into bin `b` in its current orientation: turn it into an
    /// orientation that does, or drop it at a random position of another bin.
    fn evict(&mut self, i: usize, b: usize) {
        let prob = self.prob;
        let ty = &prob.types[prob.items[i]];
        let len = self.lay.bins[b].len;
        if let Some(&os) = ty.shapes.iter().find(|&&os| prob.shapes[os].w <= len) {
            let (x, y) = self.lay.clamp_to_bin(prob, os, b, self.lay.x[i], self.lay.y[i]);
            self.lay.os[i] = os as u32;
            self.lay.x[i] = x;
            self.lay.y[i] = y;
            return;
        }
        let nb = self.lay.bins.len();
        let start = self.rng.below(nb);
        for d in 0..nb {
            let k = (start + d) % nb;
            if k == b {
                continue;
            }
            if let Some(&os) = ty.shapes.iter().find(|&&os| prob.shapes[os].w <= self.lay.bins[k].len) {
                let s = &prob.shapes[os];
                let bin = self.lay.bins[k];
                self.lay.os[i] = os as u32;
                self.lay.x[i] = self.rng.range(bin.x0, bin.x0 + bin.len - s.w);
                self.lay.y[i] = self.rng.range(0.0, prob.height - s.h);
                return;
            }
        }
    }

    /// Push relaxation ("squeeze"): every colliding pair is pushed apart along its exact
    /// penetration direction (half the depth each), Gauss-Seidel style, with positions
    /// clamped to the container.  Pressure propagates through contact chains, so small
    /// overlaps are absorbed by slack anywhere in the layout without any sampling.
    /// Returns true when the layout became collision free.
    pub fn relax(&mut self, sweeps: usize, max_pd_start: f64) -> bool {
        let t_relax = Instant::now();
        let r = self.relax_inner(sweeps, max_pd_start);
        self.wall[2] += t_relax.elapsed().as_secs_f64();
        r
    }

    fn relax_inner(&mut self, sweeps: usize, max_pd_start: f64) -> bool {
        let prob = self.prob;
        let n = self.lay.n();
        let tol = prob.tol;
        let mut order: Vec<usize> = (0..n).collect();
        let mut ok = false;
        let mut hist = [f64::INFINITY; 6];
        for sweep in 0..sweeps {
            let mut sum_pd = 0.0;
            // alternate the sweep direction so pushes travel both ways
            order.sort_by(|&a, &b| self.lay.x[a].partial_cmp(&self.lay.x[b]).unwrap().then(a.cmp(&b)));
            if sweep % 2 == 0 {
                order.reverse();
            }
            let mut any = false;
            for a in 0..n {
                let i = order[a];
                for &j in order.iter().skip(a + 1) {
                    let (si, sj) = (self.lay.os[i] as usize, self.lay.os[j] as usize);
                    let (shi, shj) = (&prob.shapes[si], &prob.shapes[sj]);
                    let (xi, yi, xj, yj) = (self.lay.x[i], self.lay.y[i], self.lay.x[j], self.lay.y[j]);
                    if xj >= xi + shi.w || xj + shj.w <= xi || yj >= yi + shi.h || yj + shj.h <= yi {
                        continue;
                    }
                    let (f, sg) = prob.pair(sj, si);
                    let (pd, ux, uy) = f.pd_dir(sg * (xi - xj), sg * (yi - yj));
                    self.aux += 1;
                    if pd <= 0.0 {
                        continue;
                    }
                    if sweep == 0 && pd > max_pd_start {
                        // too deep to be absorbed by pushing; not worth trying
                        self.tr.recompute_all(prob, &self.lt, &self.lay);
                        return false;
                    }
                    any = true;
                    sum_pd += pd;
                    let m = 0.5 * (pd + 2.0 * tol);
                    let bi = self.lay.bin_of(xi);
                    let bj = self.lay.bin_of(xj);
                    let (nxi, nyi) = self.lay.clamp_to_bin(prob, si, bi, xi + sg * ux * m, yi + sg * uy * m);
                    let (nxj, nyj) = self.lay.clamp_to_bin(prob, sj, bj, xj - sg * ux * m, yj - sg * uy * m);
                    self.lay.x[i] = nxi;
                    self.lay.y[i] = nyi;
                    self.lay.x[j] = nxj;
                    self.lay.y[j] = nyj;
                }
            }
            self.aux += (n * n / 64) as u64;
            if !any {
                ok = true;
                break;
            }
            // jammed: the total overlap is no longer shrinking geometrically
            if sum_pd > 0.5 * hist[sweep % 6] {
                break;
            }
            hist[sweep % 6] = sum_pd;
        }
        self.tr.recompute_all(prob, &self.lt, &self.lay);
        ok && self.tr.total_loss() == 0.0
    }

    fn move_items_multi(&mut self, cfg: &SearchCfg) {
        let (prob, lt) = (self.prob, &self.lt);
        let (lay, loss, wt) = (&self.lay, &self.tr.loss, &self.tr.wt);
        struct WPtr(*mut Worker);
        unsafe impl Sync for WPtr {}
        impl WPtr {
            fn at(&self, k: usize) -> *mut Worker {
                unsafe { self.0.add(k) }
            }
        }
        let wp = WPtr(self.workers.as_mut_ptr());
        let nw = self.workers.len();
        self.pool.run(|k| {
            if k < nw {
                // SAFETY: every participant index is handed to exactly one thread
                let w = unsafe { &mut *wp.at(k) };
                w.load(lay, loss);
                w.move_items(prob, lt, wt, cfg);
            }
        });
        let n = self.tr.n;
        let mut bi = 0;
        let mut bl = f64::INFINITY;
        for (k, w) in self.workers.iter().enumerate() {
            let l = total_wloss(&w.loss, wt, n);
            if l < bl {
                bl = l;
                bi = k;
            }
        }
        let bw = &self.workers[bi];
        self.lay.clone_from(&bw.lay);
        self.tr.loss.copy_from_slice(&bw.loss);
    }

    /// Try to remove all collisions from the current layout.  Stops when feasible, after
    /// `k_max` consecutive attempts of `n_max` non-improving iterations, or at the deadline.
    /// The separator is left at (and returns) the least infeasible layout found.
    pub fn separate(&mut self, cfg: &SearchCfg, n_max: usize, k_max: usize, limit: Limit) -> SepResult {
        let mut min_loss = self.tr.total_loss();
        let mut min_lay = self.lay.clone();
        let mut min_mat = self.tr.loss.clone();
        let mut strikes = 0;
        'outer: while strikes < k_max && min_loss > 0.0 {
            let mut no_imp = 0;
            let init_loss = self.tr.total_loss();
            while no_imp < n_max {
                let tw = Instant::now();
                self.move_items_multi(cfg);
                self.wall[0] += tw.elapsed().as_secs_f64();
                let tw = Instant::now();
                self.iters += 1;
                let loss = self.tr.total_loss();
                if loss == 0.0 {
                    min_loss = 0.0;
                    min_lay.clone_from(&self.lay);
                    min_mat.copy_from_slice(&self.tr.loss);
                    break 'outer;
                } else if loss < min_loss {
                    if loss < min_loss * 0.98 {
                        no_imp = 0;
                    }
                    min_loss = loss;
                    min_lay.clone_from(&self.lay);
                    min_mat.copy_from_slice(&self.tr.loss);
                    if self.finisher > 0 {
                        // a new best state: if only shallow overlaps remain, try to push them
                        // into whatever slack is left
                        self.n_finish_try += 1;
                        if self.relax(self.finisher, self.finish_pd) {
                            self.n_finish_ok += 1;
                            min_loss = 0.0;
                            min_lay.clone_from(&self.lay);
                            min_mat.copy_from_slice(&self.tr.loss);
                            break 'outer;
                        }
                        self.lay.clone_from(&min_lay);
                        self.tr.loss.copy_from_slice(&min_mat);
                    }
                } else {
                    no_imp += 1;
                }
                self.tr.update_weights(self.gls.0, self.gls.1, self.gls.2);
                self.wall[1] += tw.elapsed().as_secs_f64();
                if limit.reached(self.evals()) {
                    break 'outer;
                }
            }
            if init_loss * 0.98 <= min_loss {
                strikes += 1;
            } else {
                strikes = 0;
            }
            // roll back to the best layout of this attempt, keeping the weights
            self.lay.clone_from(&min_lay);
            self.tr.loss.copy_from_slice(&min_mat);
        }
        self.lay.clone_from(&min_lay);
        self.tr.loss.copy_from_slice(&min_mat);
        SepResult { lay: min_lay, loss: min_loss }
    }
}
