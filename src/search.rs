//! Searching a new placement for one item: sampling + coordinate descent on the weighted
//! collision loss, evaluated with the exact NFP fields.

use crate::layout::{Layout, LossTable};
use crate::nfp::{Field, Grid, Kind};
use crate::problem::Problem;
use crate::rng::Rng;

#[derive(Clone, Copy, Debug)]
pub struct SearchCfg {
    pub n_container: usize,
    pub n_focus: usize,
    pub n_descents: usize,
    /// (initial, final) step of the first refinement, relative to the item's smaller side
    pub cd_pre: (f64, f64),
    /// (initial, final) step of the final refinement
    pub cd_fin: (f64, f64),
    /// polish the result into exact contact with its neighbours
    pub snap: bool,
}

#[derive(Clone, Copy)]
struct Nb<'a> {
    /// grid coordinates of a position: fx = ax * x + bx, fy = ay * y + by
    ax: f64,
    bx: f64,
    ay: f64,
    by: f64,
    kind: *const u8,
    nx: usize,
    ny: usize,
    grid: Option<&'a Grid>,
    f: &'a Field,
    sg: f64,
    px: f64,
    py: f64,
    w: f64,
    k2: f64,
    k1: f64,
    k0: f64,
    poles: &'a [[f64; 3]],
    pe: f64,
    pl: f64,
}

unsafe impl<'a> Send for Nb<'a> {}
unsafe impl<'a> Sync for Nb<'a> {}

/// Everything needed to evaluate placements of one item in one orientation.
pub struct RotCtx<'a> {
    pub os: usize,
    pub w: f64,
    pub h: f64,
    lox: Vec<f64>,
    hix: Vec<f64>,
    loy: Vec<f64>,
    hiy: Vec<f64>,
    nb: Vec<Nb<'a>>,
    local: Vec<u32>,
    lcx: f64,
    lcy: f64,
    lr: f64,
    /// constructive mode: any collision is infeasible, clear placements are scored by position
    lbf: bool,
    cxo: f64,
    cyo: f64,
    mode: u8,
    poles: &'a [[f64; 3]],
}

impl<'a> RotCtx<'a> {
    /// `others`: (item index, weight) of every item that can collide with the moving one.
    pub fn new(prob: &'a Problem, lt: &LossTable, lay: &Layout, os: usize, others: impl Iterator<Item = (usize, f64)>, lbf: bool) -> Self {
        let s = &prob.shapes[os];
        let cap = lay.n();
        let mut c = RotCtx {
            os,
            w: s.w,
            h: s.h,
            lox: Vec::with_capacity(cap),
            hix: Vec::with_capacity(cap),
            loy: Vec::with_capacity(cap),
            hiy: Vec::with_capacity(cap),
            nb: Vec::with_capacity(cap),
            local: Vec::new(),
            lcx: f64::NAN,
            lcy: f64::NAN,
            lr: 0.0,
            lbf,
            cxo: s.centroid.x,
            cyo: s.centroid.y,
            mode: lt.mode,
            poles: &s.poles,
        };
        for (j, w) in others {
            let sj = lay.os[j] as usize;
            let (f, sg) = prob.pair(sj, os);
            let (px, py) = (lay.x[j], lay.y[j]);
            // positions p of the moving item with sg * (p - pj) inside the NFP bounding box
            let (lx, hx, ly, hy) = if sg > 0.0 {
                (px + f.bb[0], px + f.bb[2], py + f.bb[1], py + f.bb[3])
            } else {
                (px - f.bb[2], px - f.bb[0], py - f.bb[3], py - f.bb[1])
            };
            let k = lt.get(os, sj);
            c.lox.push(lx);
            c.hix.push(hx);
            c.loy.push(ly);
            c.hiy.push(hy);
            let (grid, kind, nx, ny, ax, bx, ay, by) = match &f.kind {
                Kind::Grid(g) => (
                    Some(&**g),
                    g.kind.as_ptr(),
                    g.nx,
                    g.ny,
                    sg * g.inv_h,
                    (-sg * px - g.ox) * g.inv_h,
                    sg * g.inv_h,
                    (-sg * py - g.oy) * g.inv_h,
                ),
                _ => (None, std::ptr::null(), 1, 1, 0.0, 0.0, 0.0, 0.0),
            };
            let e = lt.get_pp(os, sj);
            c.nb.push(Nb { ax, bx, ay, by, kind, nx, ny, grid, f, sg, px, py, w, k2: k[0], k1: k[1], k0: k[2], poles: &prob.shapes[sj].poles, pe: e[0], pl: e[1] });
        }
        c
    }

    #[inline(always)]
    fn contrib(&self, j: usize, x: f64, y: f64) -> f64 {
        unsafe {
            if x > *self.lox.get_unchecked(j) && x < *self.hix.get_unchecked(j) && y > *self.loy.get_unchecked(j) && y < *self.hiy.get_unchecked(j) {
                let nb = self.nb.get_unchecked(j);
                let pd = if let Some(g) = nb.grid {
                    let fx = nb.ax * x + nb.bx;
                    let fy = nb.ay * y + nb.by;
                    let ix = (fx as usize).min(nb.nx - 1);
                    let iy = (fy as usize).min(nb.ny - 1);
                    let c = iy * nb.nx + ix;
                    let k = *nb.kind.add(c);
                    if k == 0 {
                        return 0.0;
                    }
                    if self.mode == 2 {
                        let a = g.proxy_cell(c, k, ix, iy, fx, fy, nb.sg * (x - nb.px), nb.sg * (y - nb.py), nb.f.eps);
                        if a < 0.0 {
                            return 0.0;
                        }
                        return nb.w * (a + nb.pe * nb.pe).sqrt() * nb.pl;
                    }
                    g.pd_cell(c, k, ix, iy, fx, fy, nb.sg * (x - nb.px), nb.sg * (y - nb.py), nb.f.eps)
                } else {
                    nb.f.pd(nb.sg * (x - nb.px), nb.sg * (y - nb.py))
                };
                if pd > 0.0 {
                    if self.mode == 1 {
                        let a = crate::geom::pole_proxy(self.poles, x, y, nb.poles, nb.px, nb.py, nb.pe);
                        return nb.w * (a + nb.pe * nb.pe).sqrt() * nb.pl;
                    }
                    return nb.w * ((nb.k2 * pd + nb.k1) * pd + nb.k0).sqrt();
                }
            }
        }
        0.0
    }

    /// Weighted loss of placing the item at (x, y); returns +inf as soon as `bound` is reached.
    #[inline]
    pub fn eval(&self, x: f64, y: f64, bound: f64) -> f64 {
        if self.lbf {
            for j in 0..self.nb.len() {
                if self.contrib(j, x, y) > 0.0 {
                    return f64::INFINITY;
                }
            }
            return 10.0 * (2.0 * x + self.cxo) + (2.0 * y + self.cyo);
        }
        let mut acc = 0.0;
        for j in 0..self.nb.len() {
            let c = self.contrib(j, x, y);
            if c > 0.0 {
                acc += c;
                if acc >= bound {
                    return f64::INFINITY;
                }
            }
        }
        acc
    }

    pub fn build_local(&mut self, x: f64, y: f64, r: f64) {
        self.local.clear();
        for j in 0..self.nb.len() {
            if x > self.lox[j] - r && x < self.hix[j] + r && y > self.loy[j] - r && y < self.hiy[j] + r {
                self.local.push(j as u32);
            }
        }
        self.lcx = x;
        self.lcy = y;
        self.lr = r;
    }

    /// Same as `eval` but only scans the neighbours near the last local-list centre.
    #[inline]
    pub fn eval_local(&mut self, x: f64, y: f64, bound: f64, radius: f64) -> f64 {
        if !((x - self.lcx).abs() <= self.lr && (y - self.lcy).abs() <= self.lr) {
            self.build_local(x, y, radius);
        }
        if self.lbf {
            for &j in &self.local {
                if self.contrib(j as usize, x, y) > 0.0 {
                    return f64::INFINITY;
                }
            }
            return 10.0 * (2.0 * x + self.cxo) + (2.0 * y + self.cyo);
        }
        let mut acc = 0.0;
        for &j in &self.local {
            let c = self.contrib(j as usize, x, y);
            if c > 0.0 {
                acc += c;
                if acc >= bound {
                    return f64::INFINITY;
                }
            }
        }
        acc
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Cand {
    /// index into the rotation contexts
    pub k: usize,
    pub bin: usize,
    pub x: f64,
    pub y: f64,
    pub e: f64,
}

struct BestSamples {
    cap: usize,
    v: Vec<Cand>,
}

impl BestSamples {
    #[inline]
    fn bound(&self) -> f64 {
        if self.v.len() == self.cap {
            self.v[self.cap - 1].e
        } else {
            f64::INFINITY
        }
    }

    fn report(&mut self, s: Cand, thr: f64) -> bool {
        if !(s.e < self.bound()) {
            return false;
        }
        let sim = |a: &Cand| a.k == s.k && (a.x - s.x).abs() < thr && (a.y - s.y).abs() < thr;
        if self.v.iter().any(sim) {
            if self.v.iter().filter(|a| sim(a)).all(|a| s.e < a.e) {
                self.v.retain(|a| !sim(a));
            } else {
                return false;
            }
        } else if self.v.len() == self.cap {
            self.v.pop();
        }
        let pos = self.v.iter().position(|a| s.e < a.e).unwrap_or(self.v.len());
        self.v.insert(pos, s);
        true
    }
}

/// Adaptive coordinate descent (step grows on success, halves on failure, random axis among
/// horizontal / vertical / both diagonals), after Loshchilov et al. and sparrow.
pub fn refine(rcs: &mut [RotCtx], lay: &Layout, prob: &Problem, start: Cand, step0: f64, step_min: f64, rng: &mut Rng, evals: &mut u64) -> Cand {
    let mut cur = start;
    let rc = &mut rcs[cur.k];
    let (mut sx, mut sy) = (step0, step0);
    let mut axis = rng.below(4);
    let radius = 3.0 * step0;
    let mut guard = 0;
    while (sx >= step_min || sy >= step_min) && guard < 2000 {
        guard += 1;
        let (c1, c2) = match axis {
            0 => ((cur.x + sx, cur.y), (cur.x - sx, cur.y)),
            1 => ((cur.x, cur.y + sy), (cur.x, cur.y - sy)),
            2 => ((cur.x + sx, cur.y + sy), (cur.x - sx, cur.y - sy)),
            _ => ((cur.x - sx, cur.y + sy), (cur.x + sx, cur.y - sy)),
        };
        let c1 = lay.clamp_to_bin(prob, rc.os, cur.bin, c1.0, c1.1);
        let c2 = lay.clamp_to_bin(prob, rc.os, cur.bin, c2.0, c2.1);
        let e1 = rc.eval_local(c1.0, c1.1, cur.e, radius);
        let e2 = rc.eval_local(c2.0, c2.1, cur.e.min(e1), radius);
        *evals += 2;
        let (bc, be) = if e2 < e1 { (c2, e2) } else { (c1, e1) };
        let better = be < cur.e;
        if be <= cur.e {
            cur.x = bc.0;
            cur.y = bc.1;
            cur.e = be;
        }
        let m: f64 = if better { 1.1 } else { 0.5 };
        match axis {
            0 => sx *= m,
            1 => sy *= m,
            _ => {
                sx *= m.sqrt();
                sy *= m.sqrt();
            }
        }
        if !better {
            axis = rng.below(4);
        }
    }
    cur
}

/// Exact-contact polish.  Sampling and coordinate descent only get *close* to a contact;
/// here the placement is projected onto the exact NFP boundary of nearby neighbours
/// (alternating between them for corner positions), which is what makes exact fits reachable.
/// In constructive / compaction mode a snap is kept when it lowers the position score, in
/// separation mode when it lowers the weighted collision loss.
pub fn snap(rcs: &mut [RotCtx], lay: &Layout, prob: &Problem, start: Cand, evals: &mut u64) -> Cand {
    let mut cur = start;
    let rc = &mut rcs[cur.k];
    let md = rc.w.min(rc.h);
    let reach = if rc.lbf { 0.01 * md } else { f64::INFINITY };
    let radius = 0.05 * md;
    let mut cands: Vec<(f64, f64)> = Vec::new();
    for _round in 0..4 {
        if !((cur.x - rc.lcx).abs() <= rc.lr && (cur.y - rc.lcy).abs() <= rc.lr) {
            rc.build_local(cur.x, cur.y, radius);
        }
        cands.clear();
        for &j in &rc.local {
            let nb = &rc.nb[j as usize];
            let (dx, dy) = (nb.sg * (cur.x - nb.px), nb.sg * (cur.y - nb.py));
            if !rc.lbf && nb.f.pd(dx, dy) <= 0.0 {
                continue;
            }
            if let Some((qx, qy)) = nb.f.nearest_boundary(dx, dy) {
                let (sx, sy) = (nb.px + nb.sg * qx, nb.py + nb.sg * qy);
                if (sx - cur.x).abs() <= reach && (sy - cur.y).abs() <= reach {
                    cands.push(lay.clamp_to_bin(prob, rc.os, cur.bin, sx, sy));
                }
            }
        }
        let mut best = cur;
        for &(sx, sy) in &cands {
            if sx == cur.x && sy == cur.y {
                continue;
            }
            let e = rc.eval_local(sx, sy, f64::INFINITY, radius);
            *evals += 1;
            if e < best.e {
                best = Cand { x: sx, y: sy, e, ..cur };
            }
        }
        if best.e < cur.e {
            cur = best;
        } else {
            break;
        }
    }
    cur
}

/// Find a good placement for an item given one context per allowed orientation.
///
/// `current`: the item's present (rotation index, x, y) when it is being repositioned.
/// Returns `None` only when no finite candidate exists (constructive mode).
pub fn search_placement(
    prob: &Problem,
    lay: &Layout,
    rcs: &mut [RotCtx],
    current: Option<(usize, f64, f64)>,
    cfg: &SearchCfg,
    only_bin: Option<usize>,
    rng: &mut Rng,
    evals: &mut u64,
    prof: &mut [f64; 8],
) -> Option<Cand> {
    let tp0 = std::time::Instant::now();
    let nrot = rcs.len();
    let mut best = BestSamples { cap: cfg.n_descents.max(1), v: Vec::with_capacity(cfg.n_descents + 1) };
    let min_dim = rcs.iter().map(|r| r.w.min(r.h)).fold(f64::INFINITY, f64::min);
    let thr = 0.05 * min_dim;

    if let Some((k, x, y)) = current {
        let bin = lay.bin_of(x);
        let e = rcs[k].eval(x, y, f64::INFINITY);
        *evals += 1;
        best.report(Cand { k, bin, x, y, e }, thr);
        // focused samples: bounding-box centre uniformly inside the current bounding box
        let (cw, ch) = (rcs[k].w, rcs[k].h);
        for _ in 0..cfg.n_focus {
            let kk = if nrot > 1 { rng.below(nrot) } else { 0 };
            let rc = &mut rcs[kk];
            if rc.w > lay.bins[bin].len {
                continue;
            }
            let cx = x + rng.f64() * cw;
            let cy = y + rng.f64() * ch;
            let (sx, sy) = lay.clamp_to_bin(prob, rc.os, bin, cx - 0.5 * rc.w, cy - 0.5 * rc.h);
            let e = rc.eval_local(sx, sy, best.bound(), cw.max(ch));
            *evals += 1;
            best.report(Cand { k: kk, bin, x: sx, y: sy, e }, thr);
        }
    }
    let tp1 = std::time::Instant::now();
    prof[1] += (tp1 - tp0).as_secs_f64();
    // container-wide samples
    let nb = lay.bins.len();
    for _ in 0..cfg.n_container {
        let kk = if nrot > 1 { rng.below(nrot) } else { 0 };
        let rc = &rcs[kk];
        let mut b = match only_bin {
            Some(b) => b,
            None if nb > 1 => rng.below(nb),
            None => 0,
        };
        if lay.bins[b].len < rc.w {
            if only_bin.is_some() {
                continue;
            }
            // pick any bin the item fits in
            match (0..nb).find(|&q| lay.bins[q].len >= rc.w) {
                Some(q) => b = q,
                None => continue,
            }
        }
        let bin = &lay.bins[b];
        let sx = rng.range(bin.x0, bin.x0 + bin.len - rc.w);
        let sy = rng.range(0.0, prob.height - rc.h);
        let e = rc.eval(sx, sy, best.bound());
        *evals += 1;
        best.report(Cand { k: kk, bin: b, x: sx, y: sy, e }, thr);
    }
    let tp2 = std::time::Instant::now();
    prof[2] += (tp2 - tp1).as_secs_f64();
    // refine the most promising samples, then polish the winner
    let starts: Vec<Cand> = best.v.clone();
    for s in starts {
        let md = rcs[s.k].w.min(rcs[s.k].h);
        let r = refine(rcs, lay, prob, s, md * cfg.cd_pre.0, md * cfg.cd_pre.1, rng, evals);
        best.report(r, thr);
    }
    let tp3 = std::time::Instant::now();
    prof[3] += (tp3 - tp2).as_secs_f64();
    let b0 = *best.v.first()?;
    if !b0.e.is_finite() {
        return None;
    }
    let md = rcs[b0.k].w.min(rcs[b0.k].h);
    let mut r = refine(rcs, lay, prob, b0, md * cfg.cd_fin.0, md * cfg.cd_fin.1, rng, evals);
    if cfg.snap && (r.e > 0.0 || current.is_none()) {
        r = snap(rcs, lay, prob, r, evals);
    }
    prof[4] += tp3.elapsed().as_secs_f64();
    Some(r)
}
