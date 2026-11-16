//! Exact no-fit-polygon (NFP) penetration-depth fields.
//!
//! For two oriented shapes A (fixed) and B (moving) the NFP is the set of displacements
//! `d = pos(B) - pos(A)` for which the interiors of A and B overlap.  With convex
//! decompositions A = U A_k, B = U B_l it is the union of the convex Minkowski sums
//! `A_k (+) (-B_l)`.  We never build that union explicitly.  Instead a `Field` stores
//!
//!  * the convex pieces (as half-plane lists) for an exact "is `d` inside?" test,
//!  * the exact boundary of the union as a segment soup (piece edges clipped against all the
//!    other pieces) for exact penetration depths,
//!  * a uniform grid that classifies every cell as OUT / IN / BAND so a query touches only the
//!    handful of half-planes and segments that matter for its cell.
//!
//! Queries therefore return the *exact* collision decision everywhere and the exact
//! penetration depth (distance to the NFP boundary) wherever the boundary is near, which is
//! what lets the optimiser place items in exact contact.

use crate::geom::*;

#[derive(Clone, Debug)]
pub struct ConvexPoly {
    pub v: Vec<P>,
    /// Half-planes `a*x + b*y + c`: signed distance to each edge line, positive inside.
    pub hp: Vec<[f64; 3]>,
    pub bb: [f64; 4],
}

impl ConvexPoly {
    pub fn new(v: Vec<P>) -> Self {
        let n = v.len();
        let mut hp = Vec::with_capacity(n);
        for i in 0..n {
            let a = v[i];
            let b = v[(i + 1) % n];
            let e = b - a;
            let l = norm(e);
            let (nx, ny) = (-e.y / l, e.x / l);
            hp.push([nx, ny, -(nx * a.x + ny * a.y)]);
        }
        let bb = bbox(&v);
        ConvexPoly { v, hp, bb }
    }

    /// Signed depth of a point: min distance to the edge lines (positive inside).
    #[inline]
    pub fn depth(&self, x: f64, y: f64) -> f64 {
        let mut d = f64::INFINITY;
        for h in &self.hp {
            let s = h[0] * x + h[1] * y + h[2];
            if s < d {
                d = s;
            }
        }
        d
    }
}

fn lowest(v: &[P]) -> usize {
    let mut k = 0;
    for i in 1..v.len() {
        if v[i].y < v[k].y || (v[i].y == v[k].y && v[i].x < v[k].x) {
            k = i;
        }
    }
    k
}

/// Minkowski sum of two convex counter-clockwise polygons.
pub fn minkowski_convex(a: &[P], b: &[P]) -> Vec<P> {
    let (n, m) = (a.len(), b.len());
    let (ia, ib) = (lowest(a), lowest(b));
    let mut out = Vec::with_capacity(n + m);
    let (mut i, mut j) = (0usize, 0usize);
    while i < n || j < m {
        let pa = a[(ia + i) % n];
        let pb = b[(ib + j) % m];
        out.push(pa + pb);
        if i >= n {
            j += 1;
        } else if j >= m {
            i += 1;
        } else {
            let ea = a[(ia + i + 1) % n] - pa;
            let eb = b[(ib + j + 1) % m] - pb;
            let c = cross(ea, eb);
            if c > 0.0 {
                i += 1;
            } else if c < 0.0 {
                j += 1;
            } else {
                i += 1;
                j += 1;
            }
        }
    }
    out
}

#[derive(Clone, Copy, Debug)]
struct Band {
    c0: u32,
    c1: u32,
    s0: u32,
    s1: u32,
}

pub struct Grid {
    pub ox: f64,
    pub oy: f64,
    pub inv_h: f64,
    pub nx: usize,
    pub ny: usize,
    /// cell kinds: 0 = OUT, 1 = IN (bilinear depth), 2 = BAND / near (exact depth)
    pub kind: Vec<u8>,
    /// band record index of the kind-2 cells
    band_of: Vec<u32>,
    /// penetration depth at the grid nodes (only meaningful around IN cells)
    nodes: Vec<f32>,
    /// overlap proxy at the grid nodes (empty unless the field was built with a proxy)
    anodes: Vec<f32>,
    band: Vec<Band>,
    clause_rng: Vec<(u32, u32)>,
    hps: Vec<[f64; 3]>,
    /// x1, y1, ex, ey, 1/len^2
    segs: Vec<[f64; 5]>,
}

pub enum Kind {
    Empty,
    Convex(ConvexPoly),
    Grid(Box<Grid>),
}

pub struct Field {
    /// Bounding box of the NFP: (minx, miny, maxx, maxy).
    pub bb: [f64; 4],
    pub eps: f64,
    pub kind: Kind,
    pub n_pieces: usize,
    pub n_segs: usize,
}

impl Grid {
    #[inline(always)]
    fn bilinear(&self, ix: usize, iy: usize, fx: f64, fy: f64) -> f64 {
        let tx = (fx - ix as f64).clamp(0.0, 1.0);
        let ty = (fy - iy as f64).clamp(0.0, 1.0);
        let w = self.nx + 1;
        let i00 = iy * w + ix;
        let (v00, v10, v01, v11) = unsafe {
            (
                *self.nodes.get_unchecked(i00) as f64,
                *self.nodes.get_unchecked(i00 + 1) as f64,
                *self.nodes.get_unchecked(i00 + w) as f64,
                *self.nodes.get_unchecked(i00 + w + 1) as f64,
            )
        };
        let a = v00 + (v10 - v00) * tx;
        let b = v01 + (v11 - v01) * tx;
        a + (b - a) * ty
    }

    #[inline]
    fn pd(&self, dx: f64, dy: f64, eps: f64) -> f64 {
        let fx = (dx - self.ox) * self.inv_h;
        let fy = (dy - self.oy) * self.inv_h;
        let ix = (fx as usize).min(self.nx - 1);
        let iy = (fy as usize).min(self.ny - 1);
        let c = iy * self.nx + ix;
        let k = unsafe { *self.kind.get_unchecked(c) };
        if k == 0 {
            return 0.0;
        }
        self.pd_cell(c, k, ix, iy, fx, fy, dx, dy, eps)
    }

    fn pd_dir(&self, dx: f64, dy: f64, eps: f64) -> (f64, f64, f64) {
        let fx = (dx - self.ox) * self.inv_h;
        let fy = (dy - self.oy) * self.inv_h;
        let ix = (fx as usize).min(self.nx - 1);
        let iy = (fy as usize).min(self.ny - 1);
        let c = iy * self.nx + ix;
        let k = self.kind[c];
        if k == 0 {
            return (0.0, 0.0, 0.0);
        }
        let bil = |g: &Grid| -> (f64, f64, f64) {
            let tx = (fx - ix as f64).clamp(0.0, 1.0);
            let ty = (fy - iy as f64).clamp(0.0, 1.0);
            let w = g.nx + 1;
            let i00 = iy * w + ix;
            let (v00, v10, v01, v11) = (g.nodes[i00] as f64, g.nodes[i00 + 1] as f64, g.nodes[i00 + w] as f64, g.nodes[i00 + w + 1] as f64);
            let v = (v00 + (v10 - v00) * tx) * (1.0 - ty) + (v01 + (v11 - v01) * tx) * ty;
            let gx = (v10 - v00) * (1.0 - ty) + (v11 - v01) * ty;
            let gy = (v01 - v00) * (1.0 - tx) + (v11 - v10) * tx;
            let l = (gx * gx + gy * gy).sqrt();
            if l > 1e-300 {
                (v.max(eps), -gx / l, -gy / l)
            } else {
                (v.max(eps), 1.0, 0.0)
            }
        };
        if k == 1 {
            return bil(self);
        }
        let b = self.band[self.band_of[c] as usize];
        if b.c0 != u32::MAX {
            let mut inside = false;
            for cl in &self.clause_rng[b.c0 as usize..b.c1 as usize] {
                if self.hps[cl.0 as usize..cl.1 as usize].iter().all(|h| h[0] * dx + h[1] * dy + h[2] > eps) {
                    inside = true;
                    break;
                }
            }
            if !inside {
                return (0.0, 0.0, 0.0);
            }
        }
        let mut best = f64::INFINITY;
        let (mut bx, mut by) = (0.0, 0.0);
        for s in &self.segs[b.s0 as usize..b.s1 as usize] {
            let px = dx - s[0];
            let py = dy - s[1];
            let t = ((px * s[2] + py * s[3]) * s[4]).clamp(0.0, 1.0);
            let qx = px - t * s[2];
            let qy = py - t * s[3];
            let d2 = qx * qx + qy * qy;
            if d2 < best {
                best = d2;
                bx = qx;
                by = qy;
            }
        }
        if !best.is_finite() {
            return bil(self);
        }
        let d = best.sqrt();
        if d > 1e-300 {
            (d.max(eps), -bx / d, -by / d)
        } else {
            (eps, 1.0, 0.0)
        }
    }

    /// Overlap proxy for a query that falls into the non-OUT cell `c` (kind `k`):
    /// the interpolated proxy value (>= 0) when the shapes collide, -1 when they are clear.
    /// The collision decision itself is exact.
    #[inline]
    pub fn proxy_cell(&self, c: usize, k: u8, ix: usize, iy: usize, fx: f64, fy: f64, dx: f64, dy: f64, eps: f64) -> f64 {
        if k == 2 {
            let b = unsafe { *self.band.get_unchecked(*self.band_of.get_unchecked(c) as usize) };
            if b.c0 != u32::MAX {
                let mut inside = false;
                for cl in unsafe { self.clause_rng.get_unchecked(b.c0 as usize..b.c1 as usize) } {
                    let mut ok = true;
                    for h in unsafe { self.hps.get_unchecked(cl.0 as usize..cl.1 as usize) } {
                        if h[0] * dx + h[1] * dy + h[2] <= eps {
                            ok = false;
                            break;
                        }
                    }
                    if ok {
                        inside = true;
                        break;
                    }
                }
                if !inside {
                    return -1.0;
                }
            }
        }
        let tx = (fx - ix as f64).clamp(0.0, 1.0);
        let ty = (fy - iy as f64).clamp(0.0, 1.0);
        let w = self.nx + 1;
        let i00 = iy * w + ix;
        let (v00, v10, v01, v11) = unsafe {
            (
                *self.anodes.get_unchecked(i00) as f64,
                *self.anodes.get_unchecked(i00 + 1) as f64,
                *self.anodes.get_unchecked(i00 + w) as f64,
                *self.anodes.get_unchecked(i00 + w + 1) as f64,
            )
        };
        let a = v00 + (v10 - v00) * tx;
        let b = v01 + (v11 - v01) * tx;
        (a + (b - a) * ty).max(0.0)
    }

    /// Depth for a query that is known to fall into the non-OUT cell `c` (kind `k`).
    #[inline]
    pub fn pd_cell(&self, c: usize, k: u8, ix: usize, iy: usize, fx: f64, fy: f64, dx: f64, dy: f64, eps: f64) -> f64 {
        if k == 1 {
            return self.bilinear(ix, iy, fx, fy).max(eps);
        }
        let b = unsafe { *self.band.get_unchecked(*self.band_of.get_unchecked(c) as usize) };
        if b.c0 != u32::MAX {
            let mut inside = false;
            for cl in unsafe { self.clause_rng.get_unchecked(b.c0 as usize..b.c1 as usize) } {
                let mut ok = true;
                for h in unsafe { self.hps.get_unchecked(cl.0 as usize..cl.1 as usize) } {
                    if h[0] * dx + h[1] * dy + h[2] <= eps {
                        ok = false;
                        break;
                    }
                }
                if ok {
                    inside = true;
                    break;
                }
            }
            if !inside {
                return 0.0;
            }
        }
        let mut best = f64::INFINITY;
        for s in unsafe { self.segs.get_unchecked(b.s0 as usize..b.s1 as usize) } {
            let px = dx - s[0];
            let py = dy - s[1];
            let t = ((px * s[2] + py * s[3]) * s[4]).clamp(0.0, 1.0);
            let qx = px - t * s[2];
            let qy = py - t * s[3];
            let d2 = qx * qx + qy * qy;
            if d2 < best {
                best = d2;
            }
        }
        if best.is_finite() {
            best.sqrt().max(eps)
        } else {
            self.bilinear(ix, iy, fx, fy).max(eps)
        }
    }
}

impl Field {
    /// Penetration depth for displacement (dx, dy): 0 when the shapes do not overlap,
    /// otherwise the distance to the NFP boundary (always > 0).
    #[inline]
    pub fn pd(&self, dx: f64, dy: f64) -> f64 {
        if !(dx > self.bb[0] && dx < self.bb[2] && dy > self.bb[1] && dy < self.bb[3]) {
            return 0.0;
        }
        match &self.kind {
            Kind::Convex(c) => {
                let d = c.depth(dx, dy);
                if d > self.eps {
                    d
                } else {
                    0.0
                }
            }
            Kind::Grid(g) => g.pd(dx, dy, self.eps),
            Kind::Empty => 0.0,
        }
    }

    /// Closest point of the NFP boundary to a displacement, where exact boundary data is at
    /// hand (convex fields, and the band / near cells of gridded fields).  Works on both
    /// sides of the boundary; used to snap placements into exact contact.
    pub fn nearest_boundary(&self, dx: f64, dy: f64) -> Option<(f64, f64)> {
        let best_on = |segs: &mut dyn Iterator<Item = [f64; 5]>| -> Option<(f64, f64)> {
            let mut best = f64::INFINITY;
            let mut q = (0.0, 0.0);
            for s in segs {
                let px = dx - s[0];
                let py = dy - s[1];
                let t = ((px * s[2] + py * s[3]) * s[4]).clamp(0.0, 1.0);
                let (cx, cy) = (s[0] + t * s[2], s[1] + t * s[3]);
                let d2 = (dx - cx) * (dx - cx) + (dy - cy) * (dy - cy);
                if d2 < best {
                    best = d2;
                    q = (cx, cy);
                }
            }
            if best.is_finite() {
                Some(q)
            } else {
                None
            }
        };
        match &self.kind {
            Kind::Convex(c) => {
                let n = c.v.len();
                let mut it = (0..n).map(|i| seg5(c.v[i], c.v[(i + 1) % n]));
                best_on(&mut it)
            }
            Kind::Grid(g) => {
                let fx = (dx - g.ox) * g.inv_h;
                let fy = (dy - g.oy) * g.inv_h;
                if fx < -1.0 || fy < -1.0 || fx > g.nx as f64 + 1.0 || fy > g.ny as f64 + 1.0 {
                    return None;
                }
                let ix = (fx.max(0.0) as usize).min(g.nx - 1);
                let iy = (fy.max(0.0) as usize).min(g.ny - 1);
                let c = iy * g.nx + ix;
                if g.kind[c] != 2 {
                    return None;
                }
                let b = g.band[g.band_of[c] as usize];
                let mut it = g.segs[b.s0 as usize..b.s1 as usize].iter().cloned();
                best_on(&mut it)
            }
            Kind::Empty => None,
        }
    }

    /// Interpolated overlap proxy (>= 0) when colliding, -1 when clear.  Requires a field
    /// built with a proxy.
    #[inline]
    pub fn proxy(&self, dx: f64, dy: f64) -> f64 {
        if !(dx > self.bb[0] && dx < self.bb[2] && dy > self.bb[1] && dy < self.bb[3]) {
            return -1.0;
        }
        match &self.kind {
            Kind::Grid(g) => {
                let fx = (dx - g.ox) * g.inv_h;
                let fy = (dy - g.oy) * g.inv_h;
                let ix = (fx as usize).min(g.nx - 1);
                let iy = (fy as usize).min(g.ny - 1);
                let c = iy * g.nx + ix;
                let k = g.kind[c];
                if k == 0 {
                    return -1.0;
                }
                g.proxy_cell(c, k, ix, iy, fx, fy, dx, dy, self.eps)
            }
            _ => -1.0,
        }
    }

    /// Penetration depth together with the unit direction along which the displacement
    /// leaves the NFP fastest (towards the nearest boundary point).  (0, 0, 0) when clear.
    pub fn pd_dir(&self, dx: f64, dy: f64) -> (f64, f64, f64) {
        if !(dx > self.bb[0] && dx < self.bb[2] && dy > self.bb[1] && dy < self.bb[3]) {
            return (0.0, 0.0, 0.0);
        }
        match &self.kind {
            Kind::Convex(c) => {
                let mut d = f64::INFINITY;
                let mut k = 0;
                for (i, h) in c.hp.iter().enumerate() {
                    let s = h[0] * dx + h[1] * dy + h[2];
                    if s < d {
                        d = s;
                        k = i;
                    }
                }
                if d > self.eps {
                    (d, -c.hp[k][0], -c.hp[k][1])
                } else {
                    (0.0, 0.0, 0.0)
                }
            }
            Kind::Grid(g) => g.pd_dir(dx, dy, self.eps),
            Kind::Empty => (0.0, 0.0, 0.0),
        }
    }

    pub fn memory_bytes(&self) -> usize {
        match &self.kind {
            Kind::Grid(g) => {
                g.kind.len() * 5
                    + g.nodes.len() * 4
                    + g.anodes.len() * 4
                    + g.band.len() * 16
                    + g.clause_rng.len() * 8
                    + g.hps.len() * 24
                    + g.segs.len() * 40
            }
            Kind::Convex(c) => c.hp.len() * 24 + c.v.len() * 16,
            Kind::Empty => 0,
        }
    }
}

/// Portion of segment u->v (as parameter intervals) that lies strictly (by more than `eps`)
/// inside convex polygon `q`.
#[inline]
fn clip_inside(u: P, v: P, q: &ConvexPoly, eps: f64) -> Option<(f64, f64)> {
    let (mut lo, mut hi) = (0.0f64, 1.0f64);
    for h in &q.hp {
        let su = h[0] * u.x + h[1] * u.y + h[2] - eps;
        let sv = h[0] * v.x + h[1] * v.y + h[2] - eps;
        if su <= 0.0 && sv <= 0.0 {
            return None;
        }
        if su > 0.0 && sv > 0.0 {
            continue;
        }
        let t = su / (su - sv);
        if su <= 0.0 {
            if t > lo {
                lo = t;
            }
        } else if t < hi {
            hi = t;
        }
        if lo >= hi {
            return None;
        }
    }
    Some((lo, hi))
}

/// Exact boundary of the union of the (open) convex pieces, as a segment soup.
fn boundary_soup(pieces: &[ConvexPoly], eps: f64) -> Vec<(P, P)> {
    let mut soup = Vec::new();
    let mut iv: Vec<(f64, f64)> = Vec::new();
    for (m, pm) in pieces.iter().enumerate() {
        let n = pm.v.len();
        for e in 0..n {
            let u = pm.v[e];
            let v = pm.v[(e + 1) % n];
            let (ex0, ex1) = (u.x.min(v.x), u.x.max(v.x));
            let (ey0, ey1) = (u.y.min(v.y), u.y.max(v.y));
            iv.clear();
            let mut covered = false;
            for (qi, q) in pieces.iter().enumerate() {
                if qi == m {
                    continue;
                }
                if q.bb[0] >= ex1 || q.bb[2] <= ex0 || q.bb[1] >= ey1 || q.bb[3] <= ey0 {
                    continue;
                }
                if let Some((a, b)) = clip_inside(u, v, q, eps) {
                    if a <= 0.0 && b >= 1.0 {
                        covered = true;
                        break;
                    }
                    iv.push((a, b));
                }
            }
            if covered {
                continue;
            }
            let len = norm(v - u);
            let min_t = (eps * 0.5 / len).min(0.5);
            if iv.is_empty() {
                soup.push((u, v));
                continue;
            }
            iv.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
            let mut cur = 0.0f64;
            for &(a, b) in iv.iter() {
                if a > cur + min_t {
                    soup.push((u + (v - u) * cur, u + (v - u) * a));
                }
                if b > cur {
                    cur = b;
                }
            }
            if cur < 1.0 - min_t {
                soup.push((u + (v - u) * cur, v));
            }
        }
    }
    soup
}

#[inline]
fn seg5(a: P, b: P) -> [f64; 5] {
    let e = b - a;
    let l2 = dot(e, e);
    [a.x, a.y, e.x, e.y, if l2 > 0.0 { 1.0 / l2 } else { 0.0 }]
}

#[inline]
fn dist2_seg5(x: f64, y: f64, s: &[f64; 5]) -> f64 {
    let px = x - s[0];
    let py = y - s[1];
    let t = ((px * s[2] + py * s[3]) * s[4]).clamp(0.0, 1.0);
    let qx = px - t * s[2];
    let qy = py - t * s[3];
    qx * qx + qy * qy
}

/// Build the NFP field of fixed shape A and moving shape B from their convex pieces.
/// `res` is the number of grid cells along the longer side, `eps` the contact tolerance.
/// Inscribed circle sets of the two shapes for the optional overlap-proxy field.
pub struct ProxySpec<'a> {
    pub pa: &'a [[f64; 3]],
    pub pb: &'a [[f64; 3]],
    pub eps: f64,
}

pub fn build_field(a_pieces: &[Vec<P>], b_pieces: &[Vec<P>], res: usize, eps: f64, proxy: Option<ProxySpec>) -> Field {
    // 1. convex pieces A_k (+) (-B_l)
    let mut pieces: Vec<ConvexPoly> = Vec::with_capacity(a_pieces.len() * b_pieces.len());
    for ak in a_pieces {
        for bl in b_pieces {
            let nb: Vec<P> = bl.iter().map(|q| -*q).collect();
            let s = minkowski_convex(ak, &nb);
            if let Some(c) = clean_polygon(&s) {
                pieces.push(ConvexPoly::new(c));
            }
        }
    }
    if pieces.is_empty() {
        return Field { bb: [0.0; 4], eps, kind: Kind::Empty, n_pieces: 0, n_segs: 0 };
    }
    // 2. drop pieces contained in another piece (largest first)
    pieces.sort_by(|p1, p2| signed_area(&p2.v).partial_cmp(&signed_area(&p1.v)).unwrap());
    let np = pieces.len();
    let mut alive = vec![true; np];
    let ctol = eps * 1e-3;
    let contains = |q: &ConvexPoly, m: &ConvexPoly| -> bool {
        if m.bb[0] < q.bb[0] - ctol || m.bb[2] > q.bb[2] + ctol || m.bb[1] < q.bb[1] - ctol || m.bb[3] > q.bb[3] + ctol {
            return false;
        }
        m.v.iter().all(|w| q.depth(w.x, w.y) >= -ctol)
    };
    for m in 0..np {
        for q in 0..np {
            if q == m || !alive[q] {
                continue;
            }
            if contains(&pieces[q], &pieces[m]) && !(q > m && contains(&pieces[m], &pieces[q])) {
                alive[m] = false;
                break;
            }
        }
    }
    let pieces: Vec<ConvexPoly> = pieces.into_iter().zip(alive).filter(|(_, a)| *a).map(|(p, _)| p).collect();
    let mut bb = [f64::INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY];
    for q in &pieces {
        bb[0] = bb[0].min(q.bb[0]);
        bb[1] = bb[1].min(q.bb[1]);
        bb[2] = bb[2].max(q.bb[2]);
        bb[3] = bb[3].max(q.bb[3]);
    }
    if pieces.len() == 1 && proxy.is_none() {
        let c = pieces.into_iter().next().unwrap();
        let ns = c.v.len();
        return Field { bb, eps, kind: Kind::Convex(c), n_pieces: 1, n_segs: ns };
    }

    // 3. exact boundary
    let soup = boundary_soup(&pieces, eps);
    let segs_all: Vec<[f64; 5]> = soup.iter().map(|&(a, b)| seg5(a, b)).collect();

    // 4. grid (coarser for simple NFPs: exactness never depends on the resolution)
    // (a proxy field is interpolated, so there the resolution does matter and is kept)
    let res = if proxy.is_some() { res } else { res.min(((4.0 * (segs_all.len() as f64).sqrt()) as usize).max(12)) };
    let (w, hg) = (bb[2] - bb[0], bb[3] - bb[1]);
    let h = w.max(hg) / res.max(1) as f64;
    let nx = ((w / h).ceil() as usize).max(1);
    let ny = ((hg / h).ceil() as usize).max(1);
    let (ox, oy) = (bb[0], bb[1]);
    let inv_h = 1.0 / h;
    let ncell = nx * ny;
    let pad = eps * 4.0 + 1e-12 * h;
    const F_BAND: u8 = 1;
    const F_IN1: u8 = 2;
    const F_CIN: u8 = 4;
    let mut flag = vec![0u8; ncell];
    let cell_range = |lo: f64, hi: f64, o: f64, n: usize| -> (usize, usize) {
        let a = (((lo - o) * inv_h).floor().max(0.0) as usize).min(n - 1);
        let b = (((hi - o) * inv_h).floor().max(0.0) as usize).min(n - 1);
        (a, b)
    };
    let half_diag = h * std::f64::consts::FRAC_1_SQRT_2;
    // 4a. band cells: crossed by the boundary
    let r_band = half_diag + pad;
    for s in &segs_all {
        let (x0, x1) = (s[0].min(s[0] + s[2]), s[0].max(s[0] + s[2]));
        let (y0, y1) = (s[1].min(s[1] + s[3]), s[1].max(s[1] + s[3]));
        let (ix0, ix1) = cell_range(x0 - pad, x1 + pad, ox, nx);
        let (iy0, iy1) = cell_range(y0 - pad, y1 + pad, oy, ny);
        for iy in iy0..=iy1 {
            let cy = oy + (iy as f64 + 0.5) * h;
            for ix in ix0..=ix1 {
                let cx = ox + (ix as f64 + 0.5) * h;
                if dist2_seg5(cx, cy, s) <= r_band * r_band {
                    flag[iy * nx + ix] |= F_BAND;
                }
            }
        }
    }
    // 4b. clauses
    let mut band_id = vec![u32::MAX; ncell];
    let mut nband = 0u32;
    for c in 0..ncell {
        if flag[c] & F_BAND != 0 {
            band_id[c] = nband;
            nband += 1;
        }
    }
    // per band cell: list of (key, start, len) into tmp_hps
    let mut tmp_hps: Vec<[f64; 3]> = Vec::new();
    let mut tmp_cl: Vec<Vec<(f64, u32, u32)>> = vec![Vec::new(); nband as usize];
    let mut cross_buf: Vec<[f64; 3]> = Vec::new();
    for q in &pieces {
        let (ix0, ix1) = cell_range(q.bb[0], q.bb[2], ox, nx);
        let (iy0, iy1) = cell_range(q.bb[1], q.bb[3], oy, ny);
        for iy in iy0..=iy1 {
            let cy = oy + (iy as f64 + 0.5) * h;
            for ix in ix0..=ix1 {
                let c = iy * nx + ix;
                if flag[c] & F_IN1 != 0 {
                    continue;
                }
                let cx = ox + (ix as f64 + 0.5) * h;
                cross_buf.clear();
                let mut skip = false;
                let mut fmin = f64::INFINITY;
                for hp in &q.hp {
                    let f = hp[0] * cx + hp[1] * cy + hp[2];
                    let r = 0.5 * h * (hp[0].abs() + hp[1].abs());
                    if f + r <= eps {
                        skip = true;
                        break;
                    }
                    if f - r > eps {
                        continue;
                    }
                    cross_buf.push(*hp);
                    if f < fmin {
                        fmin = f;
                    }
                }
                if skip {
                    continue;
                }
                if cross_buf.is_empty() {
                    flag[c] |= F_IN1;
                    continue;
                }
                if flag[c] & F_BAND != 0 {
                    let st = tmp_hps.len() as u32;
                    tmp_hps.extend_from_slice(&cross_buf);
                    tmp_cl[band_id[c] as usize].push((fmin, st, cross_buf.len() as u32));
                } else if fmin > eps {
                    flag[c] |= F_CIN;
                }
            }
        }
    }
    // 4c. final cell kinds + flattened band data
    let mut cells = vec![0u32; ncell];
    let mut band: Vec<Band> = Vec::new();
    let mut clause_rng: Vec<(u32, u32)> = Vec::new();
    let mut hps: Vec<[f64; 3]> = Vec::new();
    let mut final_band_of_cell = vec![u32::MAX; ncell];
    for c in 0..ncell {
        let f = flag[c];
        if f & F_IN1 != 0 {
            cells[c] = 1;
        } else if f & F_BAND != 0 {
            let cl = &mut tmp_cl[band_id[c] as usize];
            if cl.is_empty() {
                cells[c] = 0;
                continue;
            }
            cl.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());
            let c0 = clause_rng.len() as u32;
            for &(_, st, len) in cl.iter() {
                let h0 = hps.len() as u32;
                hps.extend_from_slice(&tmp_hps[st as usize..(st + len) as usize]);
                clause_rng.push((h0, h0 + len));
            }
            let c1 = clause_rng.len() as u32;
            final_band_of_cell[c] = band.len() as u32;
            cells[c] = band.len() as u32 + 2;
            band.push(Band { c0, c1, s0: 0, s1: 0 });
        } else if f & F_CIN != 0 {
            cells[c] = 1;
        }
    }
    drop(tmp_hps);
    drop(tmp_cl);
    // 4d. segment lists: band cells, plus IN cells close to the boundary ("near" cells, which
    //     are always inside but get exact depths as well)
    let r_seg = 3.0 * half_diag + pad + 0.1 * h;
    let r_near = 1.5 * h;
    let r_near_list = r_near + 2.0 * half_diag + pad + 0.1 * h;
    let mut cell_segs: Vec<Vec<(f64, u32)>> = vec![Vec::new(); band.len()];
    let mut near_id = vec![u32::MAX; ncell];
    let mut near_segs: Vec<Vec<(f64, u32)>> = Vec::new();
    for (si, s) in segs_all.iter().enumerate() {
        let (x0, x1) = (s[0].min(s[0] + s[2]), s[0].max(s[0] + s[2]));
        let (y0, y1) = (s[1].min(s[1] + s[3]), s[1].max(s[1] + s[3]));
        let (ix0, ix1) = cell_range(x0 - r_near_list, x1 + r_near_list, ox, nx);
        let (iy0, iy1) = cell_range(y0 - r_near_list, y1 + r_near_list, oy, ny);
        for iy in iy0..=iy1 {
            let cy = oy + (iy as f64 + 0.5) * h;
            for ix in ix0..=ix1 {
                let c = iy * nx + ix;
                let b = final_band_of_cell[c];
                if b == u32::MAX && cells[c] != 1 {
                    continue;
                }
                let cx = ox + (ix as f64 + 0.5) * h;
                let d2 = dist2_seg5(cx, cy, s);
                if b != u32::MAX {
                    if d2 <= r_seg * r_seg {
                        cell_segs[b as usize].push((d2.sqrt(), si as u32));
                    }
                } else if d2 <= r_near_list * r_near_list {
                    if near_id[c] == u32::MAX {
                        near_id[c] = near_segs.len() as u32;
                        near_segs.push(Vec::new());
                    }
                    near_segs[near_id[c] as usize].push((d2.sqrt(), si as u32));
                }
            }
        }
    }
    // A segment can be the nearest one for some point of a cell only if its distance to the
    // cell centre is within one cell diagonal of the smallest such distance.
    let slack = 2.0 * half_diag + pad + 1e-9 * h;
    let mut segs: Vec<[f64; 5]> = Vec::new();
    let push_list = |list: &mut Vec<(f64, u32)>, segs: &mut Vec<[f64; 5]>| -> (u32, u32, f64) {
        let dmin = list.iter().fold(f64::INFINITY, |m, e| m.min(e.0));
        list.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let s0 = segs.len() as u32;
        for &(d, si) in list.iter() {
            if d <= dmin + slack {
                segs.push(segs_all[si as usize]);
            }
        }
        (s0, segs.len() as u32, dmin)
    };
    for (bi, list) in cell_segs.iter_mut().enumerate() {
        let (s0, s1, _) = push_list(list, &mut segs);
        band[bi].s0 = s0;
        band[bi].s1 = s1;
    }
    for c in 0..ncell {
        let k = near_id[c];
        if k == u32::MAX {
            continue;
        }
        let list = &mut near_segs[k as usize];
        let dmin = list.iter().fold(f64::INFINITY, |m, e| m.min(e.0));
        if dmin > r_near {
            continue;
        }
        let (s0, s1, _) = push_list(list, &mut segs);
        cells[c] = band.len() as u32 + 2;
        band.push(Band { c0: u32::MAX, c1: u32::MAX, s0, s1 });
    }
    // 4e. node depths around IN cells (exact distance to the boundary)
    let wn = nx + 1;
    let mut need = vec![false; wn * (ny + 1)];
    for iy in 0..ny {
        for ix in 0..nx {
            if cells[iy * nx + ix] == 1 {
                need[iy * wn + ix] = true;
                need[iy * wn + ix + 1] = true;
                need[(iy + 1) * wn + ix] = true;
                need[(iy + 1) * wn + ix + 1] = true;
            }
        }
    }
    let mut nodes = vec![0f32; wn * (ny + 1)];
    for jy in 0..=ny {
        let y = oy + jy as f64 * h;
        for jx in 0..=nx {
            if !need[jy * wn + jx] {
                continue;
            }
            let x = ox + jx as f64 * h;
            let mut best = f64::INFINITY;
            for s in &segs_all {
                let d2 = dist2_seg5(x, y, s);
                if d2 < best {
                    best = d2;
                }
            }
            nodes[jy * wn + jx] = if best.is_finite() { best.sqrt() as f32 } else { h as f32 };
        }
    }
    // 4f. overlap proxy at the nodes of all non-OUT cells
    let mut anodes: Vec<f32> = Vec::new();
    if let Some(px) = &proxy {
        let mut need = vec![false; wn * (ny + 1)];
        for iy in 0..ny {
            for ix in 0..nx {
                if cells[iy * nx + ix] != 0 {
                    need[iy * wn + ix] = true;
                    need[iy * wn + ix + 1] = true;
                    need[(iy + 1) * wn + ix] = true;
                    need[(iy + 1) * wn + ix + 1] = true;
                }
            }
        }
        let pa = &px.pa[..px.pa.len().min(32)];
        let pb = &px.pb[..px.pb.len().min(32)];
        anodes = vec![0f32; wn * (ny + 1)];
        for jy in 0..=ny {
            let y = oy + jy as f64 * h;
            for jx in 0..=nx {
                if need[jy * wn + jx] {
                    let x = ox + jx as f64 * h;
                    anodes[jy * wn + jx] = pole_proxy(pa, 0.0, 0.0, pb, x, y, px.eps) as f32;
                }
            }
        }
    }
    let n_pieces = pieces.len();
    let n_segs = segs_all.len();
    Field {
        bb,
        eps,
        kind: Kind::Grid(Box::new(Grid {
            ox,
            oy,
            inv_h,
            nx,
            ny,
            kind: cells.iter().map(|&c| c.min(2) as u8).collect(),
            band_of: cells.iter().map(|&c| c.saturating_sub(2)).collect(),
            nodes,
            anodes,
            band,
            clause_rng,
            hps,
            segs,
        })),
        n_pieces,
        n_segs,
    }
}

/// Independent overlap oracle (separating axis theorem on the convex pieces).
/// Returns the largest penetration depth over all pairs of convex pieces of A and B when B is
/// displaced by `d` relative to A (0 when no pair overlaps).
pub fn sat_depth(a_pieces: &[ConvexPoly], b_pieces: &[ConvexPoly], d: P) -> f64 {
    let mut worst = 0.0f64;
    for a in a_pieces {
        for b in b_pieces {
            if b.bb[0] + d.x >= a.bb[2] || b.bb[2] + d.x <= a.bb[0] || b.bb[1] + d.y >= a.bb[3] || b.bb[3] + d.y <= a.bb[1] {
                continue;
            }
            let mut depth = f64::INFINITY;
            // edges of A against vertices of B
            for h in &a.hp {
                let mut mx = f64::NEG_INFINITY;
                for v in &b.v {
                    let s = h[0] * (v.x + d.x) + h[1] * (v.y + d.y) + h[2];
                    if s > mx {
                        mx = s;
                    }
                }
                if mx < depth {
                    depth = mx;
                }
                if depth <= 0.0 {
                    break;
                }
            }
            if depth > 0.0 {
                for h in &b.hp {
                    let mut mx = f64::NEG_INFINITY;
                    for v in &a.v {
                        let s = h[0] * (v.x - d.x) + h[1] * (v.y - d.y) + h[2];
                        if s > mx {
                            mx = s;
                        }
                    }
                    if mx < depth {
                        depth = mx;
                    }
                    if depth <= 0.0 {
                        break;
                    }
                }
            }
            if depth > worst {
                worst = depth;
            }
        }
    }
    worst
}
