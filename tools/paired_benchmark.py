#!/usr/bin/env python3
"""Paired benchmark of kestrel against sparrow (https://github.com/JeroenGar/sparrow).

    python tools/paired_benchmark.py --sparrow /path/to/sparrow --time 1200 --seeds 1

For every instance and seed, one sparrow process and one kestrel process are started at the
same moment, with the same number of threads, the same time limit, the same seed number, the
same (reduced) priority and the same instance file, so that whatever else the machine is
doing affects both alike.  When both have finished, both layouts are checked by
tools/verify.py - exact rational arithmetic, no code shared with either solver - and the
densities in the table are the ones that checker recomputes.

Results are cached per (instance, seed) in the output directory, so an interrupted benchmark
continues where it stopped; `--report` only prints the table of what is there.

No third-party packages.  sparrow has to be built separately
(`cargo build --release --features simd,only_final_svg`).
"""
import argparse
import concurrent.futures as cf
import json
import os
import re
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INSTANCES = ["albano", "dagli", "fu", "jakobs1", "jakobs2", "mao", "marques", "shapes0", "shapes1", "blaz1", "shirts", "swim", "trousers"]
LABEL = {"blaz1": "shapes2 (blaz1)"}
# Gardeyn, Vanden Berghe & Wauters, "An open-source heuristic to reboot 2D nesting research",
# EJOR 2026 (arXiv:2509.13329v3), Tables 2 and 3: expected density and third quartile of 100
# runs of 20 minutes on 3 threads (Ryzen 9 7950X), the best layout the authors ever obtained,
# and the best result of any earlier method.
PAPER = {
    #            E(rho)  Q3     best   previous
    "albano":   (89.47, 89.58, 89.82, 89.58),
    "dagli":    (89.26, 89.50, 90.17, 89.51),
    "fu":       (92.24, 92.40, 92.41, 92.41),
    "jakobs1":  (89.09, 89.09, 89.26, 89.09),
    "jakobs2":  (84.77, 85.25, 87.73, 87.73),
    "mao":      (86.14, 86.45, 86.87, 86.05),
    "marques":  (90.93, 91.02, 92.02, 91.02),
    "shapes0":  (68.60, 68.79, 69.98, 68.79),
    "shapes1":  (75.69, 75.97, 76.73, 76.73),
    "blaz1":    (84.68, 84.88, 86.23, 84.84),
    "shirts":   (89.66, 89.84, 90.92, 88.96),
    "swim":     (78.26, 78.55, 79.83, 75.94),
    "trousers": (91.73, 91.93, 92.62, 91.06),
}
LOW_PRIORITY = 0x00004000 | 0x08000000 if os.name == "nt" else 0  # BELOW_NORMAL | CREATE_NO_WINDOW


def start(cmd, cwd, log):
    kw = dict(cwd=cwd, stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT)
    if os.name == "nt":
        return subprocess.Popen(cmd, creationflags=LOW_PRIORITY, **kw)
    return subprocess.Popen(["nice", "-n", "10"] + cmd, **kw)


def verify(inst_file, sol_file):
    """(density in %, feasible) as recomputed by the independent checker."""
    if not os.path.exists(sol_file):
        return None, False
    r = subprocess.run([sys.executable, os.path.join(ROOT, "tools", "verify.py"), inst_file, sol_file], capture_output=True, text=True)
    m = re.search(r"density\s+: ([0-9.]+)%", r.stdout)
    return (float(m.group(1)) if m else None), "RESULT: FEASIBLE" in r.stdout


def run_pair(a, inst, seed):
    d = os.path.join(a.out, f"{inst}_t{a.time:g}_s{seed}")
    res_file = os.path.join(d, "result.json")
    if os.path.exists(res_file):
        return json.load(open(res_file))
    ds, dk = os.path.join(d, "sparrow"), os.path.join(d, "kestrel")
    os.makedirs(ds, exist_ok=True)
    os.makedirs(dk, exist_ok=True)
    inst_file = os.path.join(a.data, inst + ".json")
    t0 = time.time()
    with open(os.path.join(ds, "stdout.txt"), "w") as ls, open(os.path.join(dk, "stdout.txt"), "w") as lk:
        ps = start([a.sparrow, "-i", inst_file, "-t", f"{a.time:g}", "-s", str(seed), "--workers", str(a.workers)], ds, ls)
        pk = start([a.kestrel, "solve", "-i", inst_file, "-t", f"{a.time:g}", "-s", str(seed), "-w", str(a.workers), "-o", dk, "-q"] + a.kestrel_args.split(), dk, lk)
        ps.wait()
        pk.wait()
    sd, sf = verify(inst_file, os.path.join(ds, "output", f"final_{inst}.json"))
    kd, kf = verify(inst_file, os.path.join(dk, inst + ".json"))
    res = {"instance": inst, "seed": seed, "seconds": a.time, "sparrow": sd, "sparrow_feasible": sf, "kestrel": kd, "kestrel_feasible": kf, "wall": time.time() - t0}
    json.dump(res, open(res_file, "w"))
    return res


def report(a, insts, seeds):
    rows = {}
    for inst in insts:
        for seed in seeds:
            f = os.path.join(a.out, f"{inst}_t{a.time:g}_s{seed}", "result.json")
            if os.path.exists(f):
                r = json.load(open(f))
                if r["kestrel"] is None or r["sparrow"] is None or not r["kestrel_feasible"] or not r["sparrow_feasible"]:
                    print(f"!! {inst} seed {seed}: missing or infeasible result: {r}", file=sys.stderr)
                    continue
                rows.setdefault(inst, []).append(r)
    multi = any(len(v) > 1 for v in rows.values())
    print("| instance | kestrel | sparrow | difference | sparrow paper: expected | Q3 | best ever | best before sparrow |")
    print("|---|---|---|---|---|---|---|---|")
    wins = losses = n = 0
    dsum = 0.0
    above = [0, 0, 0, 0]
    for inst in insts:
        v = rows.get(inst)
        if not v:
            continue
        k = sum(r["kestrel"] for r in v) / len(v)
        s = sum(r["sparrow"] for r in v) / len(v)
        kb, sb = max(r["kestrel"] for r in v), max(r["sparrow"] for r in v)
        n += 1
        # differences below 0.005 points count as ties (the two decimals shown are equal)
        wins += k - s >= 0.005
        losses += s - k >= 0.005
        dsum += k - s
        pe, pq, pb, pp = PAPER[inst]
        for i, ref in enumerate((pe, pq, pb, pp)):
            above[i] += (kb if i >= 2 else k) > ref
        ks = f"**{k:.2f}**" if k - s >= 0.005 else f"{k:.2f}"
        ss = f"**{s:.2f}**" if s - k >= 0.005 else f"{s:.2f}"
        if multi:
            ks += f" (best {kb:.2f})"
            ss += f" (best {sb:.2f})"
        print(f"| {LABEL.get(inst, inst)} | {ks} | {ss} | {k - s:+.2f} | {pe:.2f} | {pq:.2f} | {pb:.2f} | {pp:.2f} |")
    if n:
        print(f"\nkestrel denser on {wins}, sparrow denser on {losses}, tied on {n - wins - losses} of {n} instances; mean difference {dsum / n:+.2f} points.")
        print(f"kestrel above the paper's expected value on {above[0]}, above its third quartile on {above[1]}, "
              f"above the best layout its authors ever obtained on {above[2]}, above the best earlier method on {above[3]} of {n}.")


def main():
    exe = ".exe" if os.name == "nt" else ""
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sparrow", help="path of the sparrow executable")
    ap.add_argument("--kestrel", default=os.path.join(ROOT, "target", "release", "kestrel" + exe))
    ap.add_argument("--kestrel-args", default="--strict", help="extra arguments for kestrel (default: --strict)")
    ap.add_argument("--data", default=os.path.join(ROOT, "data", "strip"))
    ap.add_argument("--out", default=os.path.join(ROOT, "out", "paired"))
    ap.add_argument("--time", type=float, default=1200, help="time limit per run in seconds")
    ap.add_argument("--workers", type=int, default=3, help="threads per solver")
    ap.add_argument("--seeds", default="1")
    ap.add_argument("--inst", default="all", help="comma separated instance names")
    ap.add_argument("--pairs", type=int, default=2, help="pairs running at the same time (each uses 2 x workers threads)")
    ap.add_argument("--report", action="store_true", help="only print the table of the results that exist")
    a = ap.parse_args()
    a.out, a.data, a.kestrel = os.path.abspath(a.out), os.path.abspath(a.data), os.path.abspath(a.kestrel)
    insts = INSTANCES if a.inst == "all" else a.inst.split(",")
    seeds = [int(s) for s in a.seeds.split(",")]
    if not a.report:
        if not a.sparrow:
            ap.error("--sparrow is required (or use --report)")
        a.sparrow = os.path.abspath(a.sparrow)
        with cf.ThreadPoolExecutor(max_workers=a.pairs) as ex:
            futs = [ex.submit(run_pair, a, i, s) for s in seeds for i in insts]
            for f in cf.as_completed(futs):
                r = f.result()
                print(f"{time.strftime('%H:%M:%S')} {r['instance']:9s} seed {r['seed']}: kestrel {r['kestrel']}  sparrow {r['sparrow']}", flush=True)
    report(a, insts, seeds)


if __name__ == "__main__":
    main()
