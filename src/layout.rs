//! Layout state, loss table and collision tracker.

use crate::nfp::sat_depth;
use crate::problem::Problem;

#[derive(Clone, Copy, Debug)]
pub struct Bin {
    pub x0: f64,
    pub len: f64,
}

/// Positions are the bounding-box minimum corners of the oriented shapes.
#[derive(Clone, Debug)]
pub struct Layout {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    /// oriented shape id of every item
    pub os: Vec<u32>,
    pub bins: Vec<Bin>,
}

impl Layout {
    pub fn n(&self) -> usize {
        self.x.len()
    }

    /// Index of the bin an x coordinate belongs to (bins are laid out left to right).
    #[inline]
    pub fn bin_of(&self, x: f64) -> usize {
        if self.bins.len() == 1 {
            return 0;
        }
        let mut b = 0;
        for (k, bin) in self.bins.iter().enumerate() {
            if x >= bin.x0 - 1e-9 {
                b = k;
            }
        }
        b
    }

    /// Right-most extent of the items inside bin `b`, relative to the bin origin.
    pub fn used_len(&self, prob: &Problem, b: usize) -> f64 {
        let mut m: f64 = 0.0;
        for i in 0..self.n() {
            if self.bin_of(self.x[i]) == b {
                m = m.max(self.x[i] + prob.shapes[self.os[i] as usize].w - self.bins[b].x0);
            }
        }
        m
    }

    #[inline]
    pub fn clamp_to_bin(&self, prob: &Problem, os: usize, b: usize, x: f64, y: f64) -> (f64, f64) {
        let s = &prob.shapes[os];
        let bin = &self.bins[b];
        let xh = (bin.x0 + bin.len - s.w).max(bin.x0);
        let yh = (prob.height - s.h).max(0.0);
        (x.clamp(bin.x0, xh), y.clamp(0.0, yh))
    }
}

/// Pair loss `sqrt(k2 * d^2 + k1 * d + k0)` (d = penetration depth) per pair of oriented shapes.
///
/// With `u = d / s` (s = geometric mean of the two shape sizes) this is
/// `s^2 * sqrt(c2 u^2 + c1 u + c0^2)`: a small floor `c0` paid by every collision (so that
/// "resolved" is strictly better than "almost resolved"), a steep start for shallow overlaps
/// and linear growth for deep ones (so that neighbours are pushed apart to make room).  The
/// `s^2` factor makes collisions between big, awkward items weigh more.
pub struct LossTable {
    ns: usize,
    k: Vec<[f64; 3]>,
}

impl LossTable {
    pub fn new(prob: &Problem, c: (f64, f64, f64)) -> Self {
        let ns = prob.n_shapes();
        let mut k = vec![[0.0; 3]; ns * ns];
        for a in 0..ns {
            for b in 0..ns {
                let s2 = prob.shapes[a].size * prob.shapes[b].size;
                let s = s2.sqrt();
                k[a * ns + b] = [c.2 * s2, c.1 * s2 * s, c.0 * c.0 * s2 * s2];
            }
        }
        LossTable { ns, k }
    }

    #[inline(always)]
    pub fn get(&self, a: usize, b: usize) -> [f64; 3] {
        unsafe { *self.k.get_unchecked(a * self.ns + b) }
    }

    #[inline(always)]
    pub fn loss(&self, a: usize, b: usize, pd: f64) -> f64 {
        let k = self.get(a, b);
        ((k[0] * pd + k[1]) * pd + k[2]).sqrt()
    }
}

#[inline]
pub fn pair_loss(prob: &Problem, lt: &LossTable, lay: &Layout, i: usize, j: usize) -> f64 {
    let (si, sj) = (lay.os[i] as usize, lay.os[j] as usize);
    let (f, sg) = prob.pair(sj, si);
    let pd = f.pd(sg * (lay.x[i] - lay.x[j]), sg * (lay.y[i] - lay.y[j]));
    if pd > 0.0 {
        lt.loss(si, sj, pd)
    } else {
        0.0
    }
}

/// Recompute row/column `i` of a symmetric n x n loss matrix.
pub fn recompute_item(loss: &mut [f64], prob: &Problem, lt: &LossTable, lay: &Layout, i: usize) {
    let n = lay.n();
    let (xi, yi) = (lay.x[i], lay.y[i]);
    let si = lay.os[i] as usize;
    let (wi, hi) = (prob.shapes[si].w, prob.shapes[si].h);
    for j in 0..n {
        let mut l = 0.0;
        if j != i {
            let sj = lay.os[j] as usize;
            let s = &prob.shapes[sj];
            // cheap bounding box rejection
            if lay.x[j] < xi + wi && lay.x[j] + s.w > xi && lay.y[j] < yi + hi && lay.y[j] + s.h > yi {
                l = pair_loss(prob, lt, lay, i, j);
            }
        }
        loss[i * n + j] = l;
        loss[j * n + i] = l;
    }
}

#[derive(Clone)]
pub struct Tracker {
    pub n: usize,
    pub loss: Vec<f64>,
    pub wt: Vec<f64>,
}

impl Tracker {
    pub fn new(prob: &Problem, lt: &LossTable, lay: &Layout) -> Self {
        let n = lay.n();
        let mut t = Tracker { n, loss: vec![0.0; n * n], wt: vec![1.0; n * n] };
        t.recompute_all(prob, lt, lay);
        t
    }

    pub fn recompute_all(&mut self, prob: &Problem, lt: &LossTable, lay: &Layout) {
        for i in 0..self.n {
            recompute_item(&mut self.loss, prob, lt, lay, i);
        }
    }

    pub fn reset_weights(&mut self) {
        for w in self.wt.iter_mut() {
            *w = 1.0;
        }
    }

    pub fn total_loss(&self) -> f64 {
        total_loss(&self.loss, self.n)
    }

    pub fn total_wloss(&self) -> f64 {
        total_wloss(&self.loss, &self.wt, self.n)
    }

    /// Guided local search weight update (Voudouris et al.; multipliers as in sparrow).
    pub fn update_weights(&mut self, w_max: f64, w_min: f64, w_decay: f64) {
        let mut mx: f64 = 0.0;
        for &l in &self.loss {
            mx = mx.max(l);
        }
        if mx <= 0.0 {
            return;
        }
        for (w, &l) in self.wt.iter_mut().zip(self.loss.iter()) {
            let m = if l == 0.0 { w_decay } else { w_min + (w_max - w_min) * (l / mx) };
            *w = (*w * m).max(1.0);
        }
    }
}

pub fn total_loss(loss: &[f64], n: usize) -> f64 {
    let mut s = 0.0;
    for i in 0..n {
        for j in i + 1..n {
            s += loss[i * n + j];
        }
    }
    s
}

pub fn total_wloss(loss: &[f64], wt: &[f64], n: usize) -> f64 {
    let mut s = 0.0;
    for i in 0..n {
        for j in i + 1..n {
            s += loss[i * n + j] * wt[i * n + j];
        }
    }
    s
}

#[inline]
pub fn item_loss(loss: &[f64], n: usize, i: usize) -> f64 {
    let mut s = 0.0;
    for l in &loss[i * n..(i + 1) * n] {
        s += *l;
    }
    s
}

/// Result of the independent final check of a layout.
#[derive(Debug, Clone)]
pub struct Verification {
    /// largest penetration depth between convex pieces of two different items
    pub max_overlap: f64,
    /// largest distance by which an item sticks out of its bin
    pub max_outside: f64,
    pub worst_pair: Option<(usize, usize)>,
}

/// Verify a layout without using the NFP fields: separating-axis test on the convex pieces
/// of every pair of items plus container bounds.
pub fn verify(prob: &Problem, lay: &Layout) -> Verification {
    let n = lay.n();
    let mut v = Verification { max_overlap: 0.0, max_outside: 0.0, worst_pair: None };
    for i in 0..n {
        let si = &prob.shapes[lay.os[i] as usize];
        let b = lay.bin_of(lay.x[i]);
        let bin = &lay.bins[b];
        let out = (bin.x0 - lay.x[i])
            .max(lay.x[i] + si.w - (bin.x0 + bin.len))
            .max(-lay.y[i])
            .max(lay.y[i] + si.h - prob.height);
        v.max_outside = v.max_outside.max(out);
        for j in i + 1..n {
            let sj = &prob.shapes[lay.os[j] as usize];
            if lay.x[j] >= lay.x[i] + si.w || lay.x[j] + sj.w <= lay.x[i] || lay.y[j] >= lay.y[i] + si.h || lay.y[j] + sj.h <= lay.y[i] {
                continue;
            }
            let d = crate::geom::p(lay.x[j] - lay.x[i], lay.y[j] - lay.y[i]);
            let dep = sat_depth(&si.cpieces, &sj.cpieces, d);
            if dep > v.max_overlap {
                v.max_overlap = dep;
                v.worst_pair = Some((i, j));
            }
        }
    }
    v
}
