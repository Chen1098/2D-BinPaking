//! Convex decomposition of simple polygons: ear clipping + Hertel-Mehlhorn merging.

use crate::geom::*;
use std::collections::HashMap;

/// Triangulate a simple counter-clockwise polygon by ear clipping.
/// Returns index triples (counter-clockwise).
pub fn triangulate(poly: &[P]) -> Vec<[usize; 3]> {
    let n = poly.len();
    let mut tris: Vec<[usize; 3]> = Vec::with_capacity(n.saturating_sub(2));
    if n < 3 {
        return tris;
    }
    let mut idx: Vec<usize> = (0..n).collect();
    let mut k = 0usize;
    let mut since_clip = 0usize;
    while idx.len() > 3 {
        let m = idx.len();
        k %= m;
        let (ip, ic, inx) = (idx[(k + m - 1) % m], idx[k], idx[(k + 1) % m]);
        let (a, b, c) = (poly[ip], poly[ic], poly[inx]);
        let cr = orient(a, b, c);
        let mut is_ear = cr > 0.0;
        if is_ear {
            for &j in idx.iter() {
                if j == ip || j == ic || j == inx {
                    continue;
                }
                let v = poly[j];
                if v == a || v == b || v == c {
                    continue;
                }
                if orient(a, b, v) >= 0.0 && orient(b, c, v) >= 0.0 && orient(c, a, v) >= 0.0 {
                    is_ear = false;
                    break;
                }
            }
        }
        // a zero-area corner (after earlier clips) can always be dropped
        let degenerate = cr == 0.0;
        if is_ear || degenerate || since_clip > m {
            if is_ear || (since_clip > m && cr > 0.0) {
                tris.push([ip, ic, inx]);
            }
            idx.remove(k);
            since_clip = 0;
            if k > 0 {
                k -= 1;
            }
        } else {
            k += 1;
            since_clip += 1;
        }
    }
    if idx.len() == 3 {
        let (a, b, c) = (idx[0], idx[1], idx[2]);
        if orient(poly[a], poly[b], poly[c]) > 0.0 {
            tris.push([a, b, c]);
        }
    }
    tris
}

/// Decompose a simple counter-clockwise polygon into convex pieces.
pub fn convex_decompose(poly: &[P]) -> Vec<Vec<P>> {
    if is_convex(poly) {
        return vec![poly.to_vec()];
    }
    let tris = triangulate(poly);
    let b = bbox(poly);
    let scale = (b[2] - b[0]).max(b[3] - b[1]);
    let tol = 1e-12 * scale * scale;

    let mut polys: Vec<Option<Vec<usize>>> = tris.iter().map(|t| Some(t.to_vec())).collect();
    let mut emap: HashMap<(usize, usize), usize> = HashMap::new();
    for (id, t) in tris.iter().enumerate() {
        for e in 0..3 {
            emap.insert((t[e], t[(e + 1) % 3]), id);
        }
    }
    // diagonals: undirected edges present in both directions; try the longest first
    let mut diags: Vec<(usize, usize)> = emap
        .keys()
        .filter(|&&(u, v)| u < v && emap.contains_key(&(v, u)))
        .cloned()
        .collect();
    diags.sort_by(|a, b| {
        let la = norm(poly[a.0] - poly[a.1]);
        let lb = norm(poly[b.0] - poly[b.1]);
        lb.partial_cmp(&la).unwrap().then(a.cmp(b))
    });
    for (u, v) in diags {
        let (pa, pb) = match (emap.get(&(u, v)), emap.get(&(v, u))) {
            (Some(&x), Some(&y)) => (x, y),
            _ => continue,
        };
        if pa == pb {
            continue;
        }
        let merged = {
            let a = polys[pa].as_ref().unwrap();
            let bb = polys[pb].as_ref().unwrap();
            let (la, lb) = (a.len(), bb.len());
            let ia = (0..la).find(|&i| a[i] == u && a[(i + 1) % la] == v);
            let ib = (0..lb).find(|&i| bb[i] == v && bb[(i + 1) % lb] == u);
            let (ia, ib) = match (ia, ib) {
                (Some(x), Some(y)) => (x, y),
                _ => continue,
            };
            let a_prev_u = a[(ia + la - 1) % la];
            let a_next_v = a[(ia + 2) % la];
            let b_prev_v = bb[(ib + lb - 1) % lb];
            let b_next_u = bb[(ib + 2) % lb];
            if orient(poly[a_prev_u], poly[u], poly[b_next_u]) < -tol {
                continue;
            }
            if orient(poly[b_prev_v], poly[v], poly[a_next_v]) < -tol {
                continue;
            }
            let mut mg: Vec<usize> = Vec::with_capacity(la + lb - 2);
            for t in 0..la {
                mg.push(a[(ia + 1 + t) % la]);
            }
            for t in 0..lb - 2 {
                mg.push(bb[(ib + 2 + t) % lb]);
            }
            mg
        };
        let id = polys.len();
        for e in 0..merged.len() {
            emap.insert((merged[e], merged[(e + 1) % merged.len()]), id);
        }
        polys[pa] = None;
        polys[pb] = None;
        polys.push(Some(merged));
    }
    let mut out: Vec<Vec<P>> = Vec::new();
    for pl in polys.into_iter().flatten() {
        let pts: Vec<P> = pl.iter().map(|&i| poly[i]).collect();
        if let Some(c) = clean_polygon(&pts) {
            out.push(c);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn check(poly: &[P]) {
        let pieces = convex_decompose(poly);
        let a: f64 = pieces.iter().map(|q| signed_area(q)).sum();
        let a0 = signed_area(poly);
        assert!((a - a0).abs() <= 1e-9 * a0.abs(), "area mismatch {} vs {}", a, a0);
        for q in &pieces {
            assert!(is_convex(q));
        }
    }

    #[test]
    fn l_shape() {
        check(&[p(0., 0.), p(4., 0.), p(4., 1.), p(1., 1.), p(1., 4.), p(0., 4.)]);
    }

    #[test]
    fn comb() {
        check(&[
            p(0., 0.), p(7., 0.), p(7., 5.), p(6., 5.), p(6., 1.), p(5., 1.), p(5., 5.), p(4., 5.),
            p(4., 1.), p(3., 1.), p(3., 5.), p(2., 5.), p(2., 1.), p(1., 1.), p(1., 5.), p(0., 5.),
        ]);
    }

    #[test]
    fn star() {
        let mut v = Vec::new();
        for i in 0..10 {
            let r = if i % 2 == 0 { 5.0 } else { 2.0 };
            let a = std::f64::consts::PI * 2.0 * i as f64 / 10.0;
            v.push(p(r * a.cos(), r * a.sin()));
        }
        check(&v);
    }
}
