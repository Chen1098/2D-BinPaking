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
use rayon::prelude::*;
use std::time::Instant;

pub struct Worker {
    pub lay: Layout,
    pub loss: Vec<f64>,
    pub rng: Rng,
    pub evals: u64,
    pub moves: u64,
    order: Vec<usize>,
}

impl Worker {
    fn load(&mut self, lay: &Layout, loss: &[f64]) {
        self.lay.clone_from(lay);
        self.loss.clear();
        self.loss.extend_from_slice(loss);
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
        let c = {
            let lay = &self.lay;
            let mut rcs: Vec<RotCtx> = ty
                .shapes
                .iter()
                .map(|&os| RotCtx::new(prob, lt, lay, os, (0..n).filter(|&j| j != i).map(|j| (j, wt[i * n + j])), false))
                .collect();
            search_placement(prob, lay, &mut rcs, Some((k0, lay.x[i], lay.y[i])), cfg, &mut self.rng, &mut self.evals)
        };
        if let Some(c) = c {
            self.lay.x[i] = c.x;
            self.lay.y[i] = c.y;
            self.lay.os[i] = ty.shapes[c.k] as u32;
            recompute_item(&mut self.loss, prob, lt, &self.lay, i);
            self.moves += 1;
        }
    }
}

pub struct Separator<'a> {
    pub prob: &'a Problem,
    pub lt: LossTable,
    pub lay: Layout,
    pub tr: Tracker,
    pub workers: Vec<Worker>,
    pub rng: Rng,
    pool: Option<rayon::ThreadPool>,
    pub gls: (f64, f64, f64),
    pub iters: u64,
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
    pub fn new(prob: &'a Problem, lay: Layout, n_workers: usize, loss_c: (f64, f64, f64), gls: (f64, f64, f64), mut rng: Rng) -> Self {
        let lt = LossTable::new(prob, loss_c);
        let tr = Tracker::new(prob, &lt, &lay);
        let workers = (0..n_workers.max(1))
            .map(|_| Worker { lay: lay.clone(), loss: tr.loss.clone(), rng: rng.fork(), evals: 0, moves: 0, order: Vec::new() })
            .collect();
        let pool = if n_workers > 1 { rayon::ThreadPoolBuilder::new().num_threads(n_workers).build().ok() } else { None };
        Separator { prob, lt, lay, tr, workers, rng, pool, gls, iters: 0 }
    }

    pub fn evals(&self) -> u64 {
        self.workers.iter().map(|w| w.evals).sum()
    }
    pub fn moves(&self) -> u64 {
        self.workers.iter().map(|w| w.moves).sum()
    }

    /// Replace the layout; losses are recomputed and the weights reset.
    pub fn restore(&mut self, lay: &Layout) {
        self.lay.clone_from(lay);
        self.tr.recompute_all(self.prob, &self.lt, &self.lay);
        self.tr.reset_weights();
    }

    pub fn move_item_to(&mut self, i: usize, os: usize, x: f64, y: f64) {
        let b = self.lay.bin_of(x);
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
        }
        self.tr.recompute_all(self.prob, &self.lt, &self.lay);
        self.tr.reset_weights();
    }

    fn move_items_multi(&mut self, cfg: &SearchCfg) {
        let (prob, lt) = (self.prob, &self.lt);
        let (lay, loss, wt) = (&self.lay, &self.tr.loss, &self.tr.wt);
        let run = |w: &mut Worker| {
            w.load(lay, loss);
            w.move_items(prob, lt, wt, cfg);
        };
        match &self.pool {
            Some(pool) if self.workers.len() > 1 => {
                let ws = &mut self.workers;
                pool.install(|| ws.par_iter_mut().for_each(run));
            }
            _ => self.workers.iter_mut().for_each(run),
        }
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
                self.move_items_multi(cfg);
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
                } else {
                    no_imp += 1;
                }
                self.tr.update_weights(self.gls.0, self.gls.1, self.gls.2);
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
