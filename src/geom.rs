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
