//! Consistency checks of the NFP fields against an independent SAT oracle.

use crate::geom::*;
use crate::nfp::{sat_depth, Kind};
use crate::problem::Problem;
use crate::rng::Rng;
use std::time::Instant;

pub fn run(prob: &Problem, threads: usize, samples: usize) -> bool {
    let t0 = Instant::now();
    prob.prewarm(threads);
    let (built, grids, bytes) = prob.field_stats();
    println!(
        "[selftest] {} shapes, {} fields ({} grids) built in {:.2}s, {:.1} MB",
        prob.n_shapes(),
        built,
        grids,
        t0.elapsed().as_secs_f64(),
        bytes as f64 / 1e6
    );
    let ns = prob.n_shapes();
    let mut rng = Rng::new(12345);
    let mut bad = 0usize;
    let mut checked = 0usize;
    let mut colliding = 0usize;
    let mut max_under = 0.0f64;
    let band = 1e-6 * prob.max_dim;
    for _ in 0..samples {
        let a = rng.below(ns);
        let b = a + rng.below(ns - a);
        let f = prob.field(a, b);
        let (w, h) = (f.bb[2] - f.bb[0], f.bb[3] - f.bb[1]);
        let d = p(rng.range(f.bb[0] - 0.05 * w, f.bb[2] + 0.05 * w), rng.range(f.bb[1] - 0.05 * h, f.bb[3] + 0.05 * h));
        let pd = f.pd(d.x, d.y);
        let sat = sat_depth(&prob.shapes[a].cpieces, &prob.shapes[b].cpieces, d);
        if (sat > 0.0 && sat < 2.0 * band) || (sat == 0.0 && pd > 0.0 && pd < 2.0 * band) {
            continue; // marginal contact: either answer is within tolerance
        }
        checked += 1;
        let c1 = pd > 0.0;
        let c2 = sat > prob.tol;
        if c1 != c2 {
            bad += 1;
            if bad <= 10 {
                println!("  MISMATCH shapes ({a},{b}) d=({:.6},{:.6}) pd={pd:.6e} sat={sat:.6e}", d.x, d.y);
            }
        }
        if c1 {
            colliding += 1;
            // the true depth can never be smaller than the deepest convex piece overlap
            if pd < sat {
                max_under = max_under.max((sat - pd) / prob.max_dim);
            }
        }
    }
    println!(
        "[selftest] {checked} random displacements checked, {colliding} colliding, {bad} mismatches, max relative depth underestimate {max_under:.2e}"
    );
    // exact-depth check: a disc of radius pd around a colliding displacement must stay inside the NFP
    let mut ball_bad = 0usize;
    let mut ball_checked = 0usize;
    let mut loose = 0usize;
    for _ in 0..samples / 20 {
        let a = rng.below(ns);
        let b = a + rng.below(ns - a);
        let f = prob.field(a, b);
        if !matches!(f.kind, Kind::Grid(_) | Kind::Convex(_)) {
            continue;
        }
        let d = p(rng.range(f.bb[0], f.bb[2]), rng.range(f.bb[1], f.bb[3]));
        let pd = f.pd(d.x, d.y);
        if pd <= 0.0 || pd > 0.004 * prob.max_dim {
            continue;
        }
        ball_checked += 1;
        let mut any_out = false;
        for k in 0..48 {
            let ang = std::f64::consts::PI * 2.0 * k as f64 / 48.0;
            let u = p(ang.cos(), ang.sin());
            let q_in = d + u * (pd * (1.0 - 1e-3) - 4.0 * band).max(0.0);
            if sat_depth(&prob.shapes[a].cpieces, &prob.shapes[b].cpieces, q_in) <= prob.tol {
                ball_bad += 1;
                if ball_bad <= 5 {
                    println!("  BALL-IN FAIL shapes ({a},{b}) d=({:.6},{:.6}) pd={pd:.6e}", d.x, d.y);
                }
                break;
            }
            let q_out = d + u * (pd * 1.05 + 4.0 * band);
            if sat_depth(&prob.shapes[a].cpieces, &prob.shapes[b].cpieces, q_out) <= prob.tol {
                any_out = true;
            }
        }
        if !any_out {
            loose += 1;
        }
    }
    println!("[selftest] depth check: {ball_checked} samples, {ball_bad} too deep, {loose} looser than a 48-direction probe detects");
    // throughput
    let mut qs = Vec::with_capacity(1 << 16);
    for _ in 0..(1 << 16) {
        let a = rng.below(ns);
        let b = a + rng.below(ns - a);
        let f = prob.field(a, b);
        qs.push((f, rng.range(f.bb[0], f.bb[2]), rng.range(f.bb[1], f.bb[3])));
    }
    let t1 = Instant::now();
    let mut acc = 0.0;
    let reps = 40;
    for _ in 0..reps {
        for (f, x, y) in &qs {
            acc += f.pd(*x, *y);
        }
    }
    let dt = t1.elapsed().as_secs_f64();
    println!(
        "[selftest] field queries: {:.1} M/s single thread (checksum {:.3e})",
        (reps * qs.len()) as f64 / dt / 1e6,
        acc
    );
    bad == 0 && ball_bad == 0
}
