#!/usr/bin/env python3
"""Run kestrel over a set of instances and seeds and print the densities.

    python tools/bench.py --evals 150 --seeds 1,2,3,4                      # deterministic budget
    python tools/bench.py --time 60 --inst swim,trousers -- --strict       # wall-clock budget
    python tools/bench.py --tag no-proxy --evals 150 -- --set loss_mode=0  # A/B a setting

`--evals N` gives every run a budget of N million evaluated placements; such runs are
bit-for-bit reproducible and independent of machine load, which makes them the right tool
for comparing settings.  Everything after `--` is passed to kestrel unchanged.
Runs are started at below-normal priority on Windows.
"""
import argparse
import concurrent.futures as cf
import os
import re
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CLASSIC = "albano,blaz1,dagli,fu,jakobs1,jakobs2,mao,marques,shapes0,shapes1,shirts,swim,trousers"


def run_one(a, inst, seed, extra, outdir):
    d = os.path.join(outdir, f"{inst}_s{seed}")
    os.makedirs(d, exist_ok=True)
    cmd = [a.bin, "solve", "-i", os.path.join(a.data, inst + ".json"), "-s", str(seed), "-w", str(a.workers), "-o", d, "-q"]
    cmd += ["--evals", str(a.evals), "-t", "1000000"] if a.evals else ["-t", str(a.time)]
    flags = 0x00004000 | 0x08000000 if os.name == "nt" else 0  # BELOW_NORMAL | CREATE_NO_WINDOW
    t0 = time.time()
    p = subprocess.run(cmd + extra, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, creationflags=flags)
    m = re.search(r"density=([0-9.]+)% feasible=(\w+)", p.stdout)
    if not m:
        return inst, seed, None, False, time.time() - t0, (p.stdout + p.stderr)[-400:]
    return inst, seed, float(m.group(1)), m.group(2) == "true", time.time() - t0, ""


def main():
    exe = "kestrel.exe" if os.name == "nt" else "kestrel"
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tag", default="run")
    ap.add_argument("--bin", default=os.path.join(ROOT, "target", "release", exe))
    ap.add_argument("--data", default=os.path.join(ROOT, "data", "strip"))
    ap.add_argument("--out", default=os.path.join(ROOT, "out", "bench"))
    ap.add_argument("--evals", type=float, default=0, help="budget in million evaluations (deterministic)")
    ap.add_argument("--time", type=float, default=60, help="time limit in seconds (if --evals is not given)")
    ap.add_argument("--inst", default="shirts,swim,trousers", help=f"comma separated names, or 'all' ({CLASSIC})")
    ap.add_argument("--seeds", default="1,2,3,4")
    ap.add_argument("--jobs", type=int, default=2, help="runs in parallel")
    ap.add_argument("--workers", type=int, default=3, help="threads per run")
    ap.add_argument("extra", nargs="*")
    a = ap.parse_args()
    insts = (CLASSIC if a.inst == "all" else a.inst).split(",")
    seeds = [int(s) for s in a.seeds.split(",")]
    outdir = os.path.join(a.out, a.tag)
    res = {}
    with cf.ThreadPoolExecutor(max_workers=a.jobs) as ex:
        futs = [ex.submit(run_one, a, i, s, a.extra, outdir) for i in insts for s in seeds]
        for f in cf.as_completed(futs):
            inst, seed, dens, feas, dt, err = f.result()
            if dens is None:
                print(f"  !! {inst} seed {seed} failed: {err}", file=sys.stderr)
                continue
            if not feas:
                print(f"  !! {inst} seed {seed}: INFEASIBLE layout", file=sys.stderr)
            res.setdefault(inst, []).append(dens)
    print(f"{'instance':10s} {'mean':>8s} {'min':>8s} {'max':>8s}   n")
    means = []
    for i in insts:
        v = res.get(i, [])
        if v:
            means.append(sum(v) / len(v))
            print(f"{i:10s} {means[-1]:8.3f} {min(v):8.3f} {max(v):8.3f}   {len(v)}")
    if means:
        print(f"{'MEAN':10s} {sum(means) / len(means):8.3f}")


if __name__ == "__main__":
    main()
