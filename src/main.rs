use kestrel::layout::verify;
use kestrel::problem::{parse_json, parse_legacy, Mode, Problem};
use kestrel::strip::{solve_strip, Params};
use std::path::Path;

fn usage() -> ! {
    eprintln!(
        "kestrel - 2D irregular nesting / bin packing

USAGE:
  kestrel solve -i <instance.json|in.txt> [options]
  kestrel selftest <instance>

OPTIONS:
  -i, --input <file>      instance (JSON, or the legacy plain-text format)
  -t, --time <seconds>    time limit (default 60)
  -s, --seed <n>          random seed (default 1)
  -w, --workers <n>       worker threads (default 3)
  -o, --out <dir>         output directory (default ./out)
      --rotations <list>  override allowed rotations, e.g. 0,90,180,270
      --res <n>           NFP grid resolution (default 64)
  -q, --quiet"
    );
    std::process::exit(2);
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        usage();
    }
    let mut input = String::new();
    let mut out_dir = String::from("out");
    let mut pr = Params::default();
    let mut res = 64usize;
    let mut rots: Option<Vec<f64>> = None;
    let mut kv: Vec<(String, String)> = Vec::new();
    let cmd = args[1].clone();
    let mut i = 2;
    while i < args.len() {
        let a = args[i].clone();
        let mut val = || -> String {
            i += 1;
            match args.get(i) {
                Some(v) => v.clone(),
                None => usage(),
            }
        };
        match a.as_str() {
            "-i" | "--input" => input = val(),
            "-t" | "--time" => pr.time_limit = val().parse().unwrap_or_else(|_| usage()),
            "--evals" => pr.max_evals = (val().parse::<f64>().unwrap_or_else(|_| usage()) * 1e6) as u64,
            "-s" | "--seed" => pr.seed = val().parse().unwrap_or_else(|_| usage()),
            "-w" | "--workers" => pr.workers = val().parse().unwrap_or_else(|_| usage()),
            "-o" | "--out" => out_dir = val(),
            "--res" => res = val().parse().unwrap_or_else(|_| usage()),
            "--rotations" => rots = Some(val().split(',').filter_map(|t| t.trim().parse().ok()).collect()),
            "-q" | "--quiet" => pr.verbose = 0,
            "-v" => pr.verbose = 2,
            "-vv" => pr.verbose = 3,
            "--set" => {
                let v = val();
                let mut it = v.splitn(2, "=");
                kv.push((it.next().unwrap().to_string(), it.next().unwrap_or("").to_string()));
            }
            _ if input.is_empty() && !a.starts_with("-") => input = a.to_string(),
            _ => usage(),
        }
        i += 1;
    }
    for (k, v) in &kv {
        let f: f64 = v.parse().unwrap_or(f64::NAN);
        match k.as_str() {
            "loss_mode" => pr.loss_mode = f as u8,
            "compact" => pr.compact = f != 0.0,
            "relax" => pr.relax = f as usize,
            "finisher" => pr.finisher = f as usize,
            "finish_pd" => pr.finish_pd = f,
            "c0" => pr.loss_c.0 = f,
            "c1" => pr.loss_c.1 = f,
            "c2" => pr.loss_c.2 = f,
            "shrink" => pr.shrink = f,
            "ex_iters" => pr.ex_iters = f as usize,
            "ex_strikes" => pr.ex_strikes = f as usize,
            "cm_iters" => pr.cm_iters = f as usize,
            "cm_strikes" => pr.cm_strikes = f as usize,
            "n_container" => {
                pr.search.n_container = f as usize;
                pr.cmpr_search.n_container = f as usize;
            }
            "n_focus" => {
                pr.search.n_focus = f as usize;
                pr.cmpr_search.n_focus = f as usize;
            }
            "n_descents" => {
                pr.search.n_descents = f as usize;
                pr.cmpr_search.n_descents = f as usize;
            }
            "cd_fin_min" => {
                pr.search.cd_fin.1 = f;
                pr.cmpr_search.cd_fin.1 = f;
            }
            "cm_fin_min" => pr.cmpr_search.cd_fin.1 = f,
            "explore_frac" => pr.explore_frac = f,
            "w_max" => pr.gls.0 = f,
            "w_min" => pr.gls.1 = f,
            "w_decay" => pr.gls.2 = f,
            "pool_std" => pr.pool_std = f,
            "cm_s0" => pr.cm_shrink.0 = f,
            "cm_s1" => pr.cm_shrink.1 = f,
            _ => {
                eprintln!("unknown --set key {k}");
                std::process::exit(2);
            }
        }
    }
    if input.is_empty() {
        usage();
    }
    let text = std::fs::read_to_string(&input).unwrap_or_else(|e| {
        eprintln!("cannot read {input}: {e}");
        std::process::exit(1)
    });
    let raw = if input.to_lowercase().ends_with(".json") { parse_json(&text, 4) } else { parse_legacy(&text) };
    let raw = raw.unwrap_or_else(|e| {
        eprintln!("parse error: {e}");
        std::process::exit(1)
    });
    let prob = Problem::new(raw, res, rots.as_deref(), pr.loss_mode == 2).unwrap_or_else(|e| {
        eprintln!("error: {e}");
        std::process::exit(1)
    });
    match cmd.as_str() {
        "selftest" => {
            let ok = kestrel::selftest::run(&prob, 6, 400_000);
            std::process::exit(if ok { 0 } else { 1 });
        }
        "solve" => {}
        _ => usage(),
    }
    if prob.mode != Mode::Strip {
        eprintln!("bin mode not wired up yet");
        std::process::exit(1);
    }
    let t0 = std::time::Instant::now();
    prob.prewarm(pr.workers.max(1));
    if pr.verbose > 0 {
        let (built, grids, bytes) = prob.field_stats();
        eprintln!(
            "[setup] {} items, {} oriented shapes, {} NFP fields ({} gridded, {:.1} MB) in {:.2}s",
            prob.items.len(),
            prob.n_shapes(),
            built,
            grids,
            bytes as f64 / 1e6,
            t0.elapsed().as_secs_f64()
        );
    }
    let out = solve_strip(&prob, &pr);
    let v = verify(&prob, &out.lay);
    let feasible = v.max_overlap <= 10.0 * prob.tol && v.max_outside <= 10.0 * prob.tol;
    std::fs::create_dir_all(&out_dir).ok();
    let base = Path::new(&out_dir).join(&prob.name);
    let extra = [
        ("strip_length", format!("{:?}", out.len)),
        ("density", format!("{:?}", out.density)),
        ("feasible", format!("{}", feasible)),
        ("max_overlap", format!("{:e}", v.max_overlap)),
        ("seconds", format!("{:.3}", t0.elapsed().as_secs_f64())),
    ];
    std::fs::write(base.with_extension("json"), kestrel::out::to_json(&prob, &out.lay, &extra)).ok();
    let title = format!("{} - length {:.4}, density {:.3}%", prob.name, out.len, 100.0 * out.density);
    std::fs::write(base.with_extension("svg"), kestrel::out::to_svg(&prob, &out.lay, &title)).ok();
    println!(
        "RESULT name={} len={:.6} density={:.4}% feasible={} max_overlap={:.2e} max_outside={:.2e} evals={} moves={} iters={} time={:.1}s",
        prob.name,
        out.len,
        100.0 * out.density,
        feasible,
        v.max_overlap,
        v.max_outside,
        out.evals,
        out.moves,
        out.iters,
        t0.elapsed().as_secs_f64()
    );
}
