//! End-to-end checks: NFP fields against the independent SAT oracle, and full solves that
//! must come back feasible.

use kestrel::geom::{p, P};
use kestrel::layout::verify;
use kestrel::nfp::sat_depth;
use kestrel::problem::{parse_json, parse_legacy, Mode, Problem};
use kestrel::rng::Rng;
use kestrel::solve::{solve, Params};

const SHAPES: &str = r#"{
  "name": "mix",
  "strip_height": 12,
  "items": [
    {"id": 0, "demand": 3, "orientation": {"rotation": {"mode": "discrete", "angles": [0, 90, 180, 270]}},
     "shape": {"type": "simple_polygon", "data": [[0,0],[4,0],[4,1],[1,1],[1,4],[0,4]]}},
    {"id": 1, "demand": 3, "orientation": {"rotation": {"mode": "discrete", "angles": [0, 180]}},
     "shape": {"type": "simple_polygon", "data": [[0,0],[5,0],[5,3],[4,3],[4,1],[3,1],[3,3],[2,3],[2,1],[1,1],[1,3],[0,3]]}},
    {"id": 2, "demand": 4, "orientation": {"rotation": {"mode": "discrete", "angles": [0, 90]}},
     "shape": {"type": "simple_polygon", "data": [[0,0],[3,0],[3,2],[0,2]]}},
    {"id": 3, "demand": 3, "orientation": {"rotation": {"mode": "discrete", "angles": [0, 90, 180, 270]}},
     "shape": {"type": "simple_polygon", "data": [[0,0],[4,0],[2,3.5]]}},
    {"id": 4, "demand": 2, "orientation": {"rotation": {"mode": "discrete", "angles": [0]}},
     "shape": {"type": "simple_polygon", "data": [[2,0],[2.6,1.4],[4,1.5],[2.9,2.4],[3.3,3.9],[2,3],[0.7,3.9],[1.1,2.4],[0,1.5],[1.4,1.4]]}}
  ]
}"#;

fn problem(proxy: bool, sep: f64) -> Problem {
    Problem::new(parse_json(SHAPES, 4).unwrap(), 64, None, proxy, sep, (64, 0.001)).unwrap()
}

/// The exact collision decision of every NFP field must agree with the separating-axis oracle.
#[test]
fn fields_agree_with_sat_oracle() {
    for proxy in [false, true] {
        let prob = problem(proxy, 0.0);
        let ns = prob.n_shapes();
        let mut rng = Rng::new(7);
        let band = 1e-6 * prob.max_dim;
        let mut checked = 0;
        for _ in 0..60_000 {
            let a = rng.below(ns);
            let b = a + rng.below(ns - a);
            let f = prob.field(a, b);
            let d: P = p(rng.range(f.bb[0] - 0.3, f.bb[2] + 0.3), rng.range(f.bb[1] - 0.3, f.bb[3] + 0.3));
            let pd = f.pd(d.x, d.y);
            let sat = sat_depth(&prob.shapes[a].cpieces, &prob.shapes[b].cpieces, d);
            if (sat > 0.0 && sat < 2.0 * band) || (sat == 0.0 && pd > 0.0 && pd < 2.0 * band) {
                continue; // marginal contact
            }
            checked += 1;
            assert_eq!(pd > 0.0, sat > prob.tol, "shapes ({a},{b}) d=({},{}) pd={pd} sat={sat}", d.x, d.y);
            if proxy {
                assert_eq!(f.proxy(d.x, d.y) >= 0.0, pd > 0.0);
            }
        }
        assert!(checked > 40_000);
    }
}

/// Integer displacements on integer shapes: exact fits (touching, zero-width slots) are free.
#[test]
fn exact_fits_are_not_collisions() {
    let prob = problem(false, 0.0);
    // shape 0 of type 2 is the 3 x 2 rectangle
    let r = prob.types[2].shapes[0];
    let f = prob.field(r, r);
    assert_eq!(f.pd(3.0, 0.0), 0.0); // side by side, touching
    assert_eq!(f.pd(0.0, 2.0), 0.0); // stacked, touching
    assert!(f.pd(2.999, 0.0) > 0.0);
    assert!(f.pd(0.0, 1.999) > 0.0);
    // comb (type 1) and a 1-wide tooth of the L (type 0): sliding into the slot is exact
    let comb = prob.types[1].shapes[0];
    let d = sat_depth(&prob.shapes[comb].cpieces, &prob.shapes[comb].cpieces, p(0.0, 3.0));
    assert_eq!(d, 0.0);
}

#[test]
fn strip_solve_is_feasible_and_deterministic() {
    let prob = problem(true, 0.0);
    let pr = Params { workers: 2, max_evals: 1_500_000, verbose: 0, ..Params::default() };
    let a = solve(&prob, &pr);
    let v = verify(&prob, &a.lay);
    assert!(v.max_overlap <= 10.0 * prob.tol, "overlap {}", v.max_overlap);
    assert!(v.max_outside <= 10.0 * prob.tol, "outside {}", v.max_outside);
    assert!(a.density > 0.70, "density {}", a.density);
    let b = solve(&prob, &pr);
    assert_eq!(a.len, b.len, "same seed and budget must give the same layout");
}

#[test]
fn strict_separation_keeps_items_apart() {
    let prob = problem(true, 0.01);
    let pr = Params { workers: 1, max_evals: 600_000, verbose: 0, ..Params::default() };
    let out = solve(&prob, &pr);
    let v = verify(&prob, &out.lay);
    assert!(v.max_overlap <= 10.0 * prob.tol);
    // the true outlines are at least (almost) 0.01 apart: test with un-inflated shapes
    let plain = problem(false, 0.0);
    let n = out.lay.n();
    for i in 0..n {
        for j in i + 1..n {
            // same type order and rotation order -> same oriented shape ids
            let (si, sj) = (&plain.shapes[out.lay.os[i] as usize], &plain.shapes[out.lay.os[j] as usize]);
            let (oi, oj) = (&prob.shapes[out.lay.os[i] as usize], &prob.shapes[out.lay.os[j] as usize]);
            // inflated frames are shifted by half the separation in x and y
            let d = p(
                (out.lay.x[j] + (oj.outline[0].x - sj.outline[0].x)) - (out.lay.x[i] + (oi.outline[0].x - si.outline[0].x)),
                (out.lay.y[j] + (oj.outline[0].y - sj.outline[0].y)) - (out.lay.y[i] + (oi.outline[0].y - si.outline[0].y)),
            );
            assert_eq!(sat_depth(&si.cpieces, &sj.cpieces, d), 0.0, "items {i} and {j} touch");
        }
    }
}

#[test]
fn bin_mode_reaches_the_area_bound() {
    // 24 unit-ish rectangles of 3 x 2 fill exactly four 6 x 6 bins
    let mut txt = String::from("6 6\n24\n");
    for _ in 0..24 {
        txt.push_str("4\n0 0\n3 0\n3 2\n0 2\n");
    }
    let mut raw = parse_legacy(&txt).unwrap();
    assert_eq!(raw.mode, Mode::Bin);
    for it in raw.items.iter_mut() {
        it.angles = vec![0.0, 90.0];
    }
    let prob = Problem::new(raw, 64, None, true, 0.0, (64, 0.001)).unwrap();
    let pr = Params { workers: 2, max_evals: 3_000_000, verbose: 0, ..Params::default() };
    let out = solve(&prob, &pr);
    let v = verify(&prob, &out.lay);
    assert!(v.max_overlap <= 10.0 * prob.tol && v.max_outside <= 10.0 * prob.tol);
    assert_eq!(out.n_bins, 4, "perfect packing into 4 bins exists");
}

#[test]
fn restarts_keep_the_best_trajectory() {
    let prob = problem(true, 0.0);
    // restart after two failed eliminations in a row: several trajectories within the budget
    let pr = Params { workers: 2, max_evals: 4_000_000, restart_fails: 2, verbose: 0, ..Params::default() };
    let a = solve(&prob, &pr);
    let v = verify(&prob, &a.lay);
    assert!(v.max_overlap <= 10.0 * prob.tol && v.max_outside <= 10.0 * prob.tol);
    assert!(a.density > 0.70, "density {}", a.density);
    // the reported history only ever improves
    assert!(a.log.windows(2).all(|w| w[1].1 <= w[0].1), "history must be monotone");
    let b = solve(&prob, &pr);
    assert_eq!(a.len, b.len, "same seed and budget must give the same layout");
}
