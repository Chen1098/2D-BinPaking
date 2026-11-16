//! Basic 2D primitives: points, polygons, hulls, rotations.

use std::ops::{Add, Mul, Neg, Sub};

#[derive(Clone, Copy, Debug, PartialEq, Default)]
pub struct P {
    pub x: f64,
    pub y: f64,
}

#[inline(always)]
pub fn p(x: f64, y: f64) -> P {
    P { x, y }
}

impl Add for P {
    type Output = P;
    #[inline(always)]
    fn add(self, o: P) -> P {
        P { x: self.x + o.x, y: self.y + o.y }
    }
}
impl Sub for P {
    type Output = P;
    #[inline(always)]
    fn sub(self, o: P) -> P {
        P { x: self.x - o.x, y: self.y - o.y }
    }
}
impl Mul<f64> for P {
    type Output = P;
    #[inline(always)]
    fn mul(self, k: f64) -> P {
        P { x: self.x * k, y: self.y * k }
    }
}
impl Neg for P {
    type Output = P;
    #[inline(always)]
    fn neg(self) -> P {
        P { x: -self.x, y: -self.y }
    }
}

#[inline(always)]
pub fn cross(a: P, b: P) -> f64 {
    a.x * b.y - a.y * b.x
}
#[inline(always)]
pub fn dot(a: P, b: P) -> f64 {
    a.x * b.x + a.y * b.y
}
/// Twice the signed area of triangle (a, b, c); > 0 when c is left of a->b.
#[inline(always)]
pub fn orient(a: P, b: P, c: P) -> f64 {
    (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)
}
#[inline(always)]
pub fn norm(a: P) -> f64 {
    (a.x * a.x + a.y * a.y).sqrt()
}

pub fn signed_area(poly: &[P]) -> f64 {
    let n = poly.len();
    if n < 3 {
        return 0.0;
    }
    let o = poly[0];
    let mut s = 0.0;
    for i in 1..n - 1 {
        s += cross(poly[i] - o, poly[i + 1] - o);
    }
    0.5 * s
}

/// (minx, miny, maxx, maxy)
pub fn bbox(poly: &[P]) -> [f64; 4] {
    let mut b = [f64::INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY];
    for q in poly {
        b[0] = b[0].min(q.x);
        b[1] = b[1].min(q.y);
        b[2] = b[2].max(q.x);
        b[3] = b[3].max(q.y);
    }
    b
}

pub fn centroid(poly: &[P]) -> P {
    let n = poly.len();
    let o = poly[0];
    let (mut a, mut cx, mut cy) = (0.0, 0.0, 0.0);
    for i in 1..n - 1 {
        let u = poly[i] - o;
        let v = poly[i + 1] - o;
        let c = cross(u, v);
        a += c;
        cx += c * (u.x + v.x);
        cy += c * (u.y + v.y);
    }
    if a.abs() < 1e-300 {
        return o;
    }
    P { x: o.x + cx / (3.0 * a), y: o.y + cy / (3.0 * a) }
}

/// Remove duplicate and collinear vertices, drop zero-width spikes and return the polygon in
/// counter-clockwise order. Returns `None` when the result is degenerate.
pub fn clean_polygon(pts: &[P]) -> Option<Vec<P>> {
    let b = bbox(pts);
    let scale = (b[2] - b[0]).max(b[3] - b[1]).max(1e-300);
    let dup = 1e-12 * scale;
    let mut v: Vec<P> = Vec::with_capacity(pts.len());
    for &q in pts {
        if let Some(&l) = v.last() {
            if (l.x - q.x).abs() <= dup && (l.y - q.y).abs() <= dup {
                continue;
            }
        }
        v.push(q);
    }
    while v.len() > 1 {
        let (f, l) = (v[0], v[v.len() - 1]);
        if (l.x - f.x).abs() <= dup && (l.y - f.y).abs() <= dup {
            v.pop();
        } else {
            break;
        }
    }
    // iteratively drop collinear vertices / spikes
    loop {
        let n = v.len();
        if n < 3 {
            return None;
        }
        let mut removed = false;
        let mut out: Vec<P> = Vec::with_capacity(n);
        for i in 0..n {
            let a = if let Some(&l) = out.last() { l } else { v[(i + n - 1) % n] };
            let c = v[(i + 1) % n];
            let q = v[i];
            let e1 = q - a;
            let e2 = c - q;
            let cr = cross(e1, e2);
            if cr.abs() <= 1e-12 * norm(e1) * norm(e2) {
                removed = true;
                continue;
            }
            out.push(q);
        }
        v = out;
        if !removed {
            break;
        }
    }
    if v.len() < 3 {
        return None;
    }
    let a = signed_area(&v);
    if a.abs() <= 1e-14 * scale * scale {
        return None;
    }
    if a < 0.0 {
        v.reverse();
    }
    Some(v)
}

/// Convex hull (counter-clockwise, no collinear points).
pub fn convex_hull(pts: &[P]) -> Vec<P> {
    let mut v: Vec<P> = pts.to_vec();
    v.sort_by(|a, b| a.x.partial_cmp(&b.x).unwrap().then(a.y.partial_cmp(&b.y).unwrap()));
    v.dedup();
    let n = v.len();
    if n < 3 {
        return v;
    }
    let mut h: Vec<P> = Vec::with_capacity(2 * n);
    for i in 0..n {
        while h.len() >= 2 && orient(h[h.len() - 2], h[h.len() - 1], v[i]) <= 0.0 {
            h.pop();
        }
        h.push(v[i]);
    }
    let lower = h.len() + 1;
    for i in (0..n - 1).rev() {
        while h.len() >= lower && orient(h[h.len() - 2], h[h.len() - 1], v[i]) <= 0.0 {
            h.pop();
        }
        h.push(v[i]);
    }
    h.pop();
    h
}

pub fn is_convex(poly: &[P]) -> bool {
    let n = poly.len();
    if n < 3 {
        return false;
    }
    for i in 0..n {
        if orient(poly[i], poly[(i + 1) % n], poly[(i + 2) % n]) < 0.0 {
            return false;
        }
    }
    true
}

/// Largest distance between two vertices.
pub fn diameter(pts: &[P]) -> f64 {
    let h = convex_hull(pts);
    let mut d: f64 = 0.0;
    for i in 0..h.len() {
        for j in i + 1..h.len() {
            d = d.max(norm(h[i] - h[j]));
        }
    }
    d
}

/// Rotate a polygon around the origin. Multiples of 90 degrees are exact.
pub fn rotate_poly(poly: &[P], deg: f64) -> Vec<P> {
    let d = deg.rem_euclid(360.0);
    let (c, s) = if d == 0.0 {
        (1.0, 0.0)
    } else if d == 90.0 {
        (0.0, 1.0)
    } else if d == 180.0 {
        (-1.0, 0.0)
    } else if d == 270.0 {
        (0.0, -1.0)
    } else {
        let r = d.to_radians();
        (r.cos(), r.sin())
    };
    poly.iter().map(|q| P { x: c * q.x - s * q.y, y: s * q.x + c * q.y }).collect()
}

/// Squared distance from point `q` to segment `a-b`.
#[inline]
pub fn dist2_point_seg(q: P, a: P, b: P) -> f64 {
    let e = b - a;
    let l2 = dot(e, e);
    let t = if l2 > 0.0 { (dot(q - a, e) / l2).clamp(0.0, 1.0) } else { 0.0 };
    let d = q - (a + e * t);
    dot(d, d)
}

/// Even-odd point in polygon test (boundary points are unspecified).
pub fn point_in_poly(q: P, poly: &[P]) -> bool {
    let n = poly.len();
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let (a, b) = (poly[i], poly[j]);
        if (a.y > q.y) != (b.y > q.y) {
            let x = a.x + (q.y - a.y) * (b.x - a.x) / (b.y - a.y);
            if q.x < x {
                inside = !inside;
            }
        }
        j = i;
    }
    inside
}

/// Greedy set of inscribed circles ("poles of inaccessibility"): each new circle is the largest
/// one that fits inside the polygon without overlapping the previous ones.
/// Stops at `max_n` circles, or earlier once the covered area ratio is high enough
/// (>= 16 circles and 80%, >= 8 circles and 90%), mirroring jagua-rs' surrogate generation.
pub fn poles(poly: &[P], max_n: usize) -> Vec<[f64; 3]> {
    let n = poly.len();
    let bb = bbox(poly);
    let area = signed_area(poly).abs();
    let clearance = |q: P, ps: &[[f64; 3]]| -> f64 {
        if !point_in_poly(q, poly) {
            return -1.0;
        }
        let mut d2 = f64::INFINITY;
        for i in 0..n {
            d2 = d2.min(dist2_point_seg(q, poly[i], poly[(i + 1) % n]));
        }
        let mut c = d2.sqrt();
        for c0 in ps {
            c = c.min(((q.x - c0[0]).powi(2) + (q.y - c0[1]).powi(2)).sqrt() - c0[2]);
        }
        c
    };
    let mut out: Vec<[f64; 3]> = Vec::new();
    let mut covered = 0.0;
    let g = 48usize;
    let (w, h) = (bb[2] - bb[0], bb[3] - bb[1]);
    while out.len() < max_n {
        let mut best = (-1.0f64, p(0.0, 0.0));
        for iy in 0..g {
            for ix in 0..g {
                let q = p(bb[0] + (ix as f64 + 0.5) * w / g as f64, bb[1] + (iy as f64 + 0.5) * h / g as f64);
                let c = clearance(q, &out);
                if c > best.0 {
                    best = (c, q);
                }
            }
        }
        if best.0 <= 0.0 {
            break;
        }
        // pattern-search refinement
        let mut step = w.max(h) / g as f64;
        let (mut c, mut q) = best;
        while step > 1e-6 * w.max(h) {
            let mut improved = false;
            for (dx, dy) in [(1.0, 0.0), (-1.0, 0.0), (0.0, 1.0), (0.0, -1.0), (1.0, 1.0), (-1.0, -1.0), (1.0, -1.0), (-1.0, 1.0)] {
                let q2 = p(q.x + dx * step, q.y + dy * step);
                let c2 = clearance(q2, &out);
                if c2 > c {
                    c = c2;
                    q = q2;
                    improved = true;
                }
            }
            if !improved {
                step *= 0.5;
            }
        }
        if c <= 1e-4 * w.min(h) {
            break;
        }
        out.push([q.x, q.y, c]);
        covered += std::f64::consts::PI * c * c;
        let k = out.len();
        let ratio = covered / area;
        if (k >= 16 && ratio >= 0.8) || (k >= 8 && ratio >= 0.9) {
            break;
        }
    }
    out
}

/// sparrow's overlap proxy between two sets of circles (pole lists translated by a and b).
#[inline]
pub fn pole_proxy(pa: &[[f64; 3]], ax: f64, ay: f64, pb: &[[f64; 3]], bx: f64, by: f64, eps: f64) -> f64 {
    let mut tot = 0.0;
    let e2 = eps * eps;
    for a in pa {
        let (x1, y1, r1) = (a[0] + ax, a[1] + ay, a[2]);
        for b in pb {
            let dx = x1 - (b[0] + bx);
            let dy = y1 - (b[1] + by);
            let pd = r1 + b[2] - (dx * dx + dy * dy).sqrt();
            let dec = if pd >= eps { pd } else { e2 / (2.0 * eps - pd) };
            tot += dec * r1.min(b[2]);
        }
    }
    tot * std::f64::consts::PI
}

/// Offset a counter-clockwise polygon outwards by `d` (mitre joins; very sharp corners get a
/// square cap instead of a long spike).  Meant for small offsets relative to the features.
pub fn inflate(poly: &[P], d: f64) -> Vec<P> {
    let n = poly.len();
    let mut out = Vec::with_capacity(n + 4);
    for i in 0..n {
        let a = poly[(i + n - 1) % n];
        let b = poly[i];
        let c = poly[(i + 1) % n];
        let e1 = (b - a) * (1.0 / norm(b - a));
        let e2 = (c - b) * (1.0 / norm(c - b));
        let n1 = p(e1.y, -e1.x);
        let n2 = p(e2.y, -e2.x);
        let den = 1.0 + dot(n1, n2);
        if den > 0.2 {
            out.push(b + (n1 + n2) * (d / den));
        } else {
            out.push(b + n1 * d + e1 * d);
            out.push(b + n2 * d - e2 * d);
        }
    }
    out
}

/// Conservative simplification of a counter-clockwise simple polygon: vertices are removed
/// only by operations that *grow* the polygon (filling a concave corner, or extending the two
/// neighbours of a short convex edge until they meet), cheapest first, until the added area
/// would exceed `max_ratio` of the original area.  The result always contains the input.
pub fn simplify_outward(poly: &[P], max_ratio: f64) -> Vec<P> {
    let mut v: Vec<P> = poly.to_vec();
    let mut budget = max_ratio * signed_area(poly).abs();
    let empty = |v: &[P], a: P, b: P, c: P, skip: &[usize]| -> bool {
        // no other vertex inside (or on) the counter-clockwise triangle (a, b, c)
        for (j, q) in v.iter().enumerate() {
            if skip.contains(&j) {
                continue;
            }
            if orient(a, b, *q) >= 0.0 && orient(b, c, *q) >= 0.0 && orient(c, a, *q) >= 0.0 {
                return false;
            }
        }
        true
    };
    loop {
        let n = v.len();
        if n <= 4 {
            break;
        }
        // candidates: (added area, kind, index[, new point])
        let mut cands: Vec<(f64, u8, usize, P)> = Vec::with_capacity(2 * n);
        for i in 0..n {
            let (a, b, c) = (v[(i + n - 1) % n], v[i], v[(i + 1) % n]);
            let o = orient(a, b, c);
            if o < 0.0 {
                // concave corner: fill triangle (a, c, b)
                cands.push((-0.5 * o, 0, i, b));
            } else {
                // convex edge b -> c with convex ends: extend a->b and d->c until they meet
                let d = v[(i + 2) % n];
                if orient(b, c, d) > 0.0 {
                    let (e1, e2) = (b - a, c - d);
                    let den = cross(e1, e2);
                    if den.abs() > 1e-300 {
                        let t = cross(d - a, e2) / den;
                        let u = cross(d - a, e1) / den;
                        if t > 1.0 && u > 1.0 {
                            let x = a + e1 * t;
                            let add = 0.5 * orient(b, x, c).abs();
                            cands.push((add, 1, i, x));
                        }
                    }
                }
            }
        }
        cands.sort_by(|p1, p2| p1.0.partial_cmp(&p2.0).unwrap());
        let mut applied = false;
        for &(add, kind, i, x) in &cands {
            if add > budget {
                break;
            }
            let (ia, ic) = ((i + n - 1) % n, (i + 1) % n);
            if kind == 0 {
                if empty(&v, v[ia], v[ic], v[i], &[ia, i, ic]) {
                    v.remove(i);
                    budget -= add;
                    applied = true;
                    break;
                }
            } else if empty(&v, v[i], x, v[ic], &[i, ic]) {
                v[i] = x;
                v.remove(ic);
                budget -= add;
                applied = true;
                break;
            }
        }
        if !applied {
            break;
        }
    }
    v
}
