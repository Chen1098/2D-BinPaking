# kestrel - 2D irregular nesting and bin packing

[![CI](https://github.com/Chen1098/2D-BinPaking/actions/workflows/ci.yml/badge.svg)](https://github.com/Chen1098/2D-BinPaking/actions/workflows/ci.yml)

<p align="center">
  <img src="results/strip_20min/swim.svg" height="250" alt="SWIM"/>
  <img src="results/strip_20min/trousers.svg" height="250" alt="TROUSERS"/>
</p>

**kestrel** packs arbitrary polygons - rectangles are just the easy case - into a strip of
minimum length or into as few bins as possible.  It is a from-scratch rewrite (Rust, no
geometry dependencies) of the toy solver this repository used to contain, built to compete
with the best published results.

> **Result in one paragraph.**  The reference is
> [sparrow](https://github.com/JeroenGar/sparrow) (Gardeyn, Vanden Berghe & Wauters,
> *European Journal of Operational Research* 2026), whose best layouts match or beat those
> of every earlier method on all thirteen instances of the classic irregular strip packing
> benchmark.  Both solvers were run **side by side on the same machine**, under that paper's
> protocol (20 minutes, 3 threads), on those thirteen instances.  After 20 minutes kestrel's
> layout was denser on **11** of them, sparrow's on **1**, and one was a tie; the mean
> difference is **+0.48** density points.  kestrel needed a median of **48 seconds** to pass
> the layout sparrow ended with after 20 minutes, and - on a slower machine - it ends above
> the expected value the paper reports for sparrow on 9 instances.  Every layout of both
> solvers was checked by the same independent verifier in exact rational arithmetic.

[中文说明](README.zh-CN.md) · [How it works](docs/ALGORITHM.md) · [Benchmarks](#benchmarks) · [Limitations](#what-this-is-not)

## Quick start

```bash
cargo build --release

# strip packing: minimise the strip length (JSON instance, sparrow / jagua-rs format)
./target/release/kestrel solve -i data/strip/swim.json -t 60

# bin packing: the original plain-text format of this repo (board W H, then polygons)
./target/release/kestrel solve -i data/legacy/in.txt -t 30

# check any result with the independent verifier (plain Python, no packages)
python tools/verify.py data/strip/swim.json out/swim.json
```

Each run writes `out/<name>.json` (rotation and translation of every item) and
`out/<name>.svg` (a drawing of the layout), and prints one summary line:

```
RESULT name=swim bins=1 len=5596.689844 density=79.0411% feasible=true ...
```

Useful options (`kestrel` without arguments prints all of them):

| option | meaning |
|---|---|
| `-t <seconds>` | time limit (default 60) |
| `-w <threads>` | worker threads (default 3) |
| `-s <seed>` | random seed |
| `--evals <millions>` | deterministic budget instead of a time limit: same seed, same result, on any machine |
| `--rotations 0,90,180,270` | override the rotations allowed by the instance |
| `--mode strip\|bin`, `--bin-width <w>` | reinterpret an instance as strip / bin packing |
| `--min-sep <d>`, `--strict` | keep a minimum distance between items (`--strict`: 1e-6 of the item size, which makes the layout overlap free in exact arithmetic) |
| `--simplify <ratio>` | simplify very detailed outlines, only ever growing them |

### Input

* **JSON** (the format of [jagua-rs](https://github.com/JeroenGar/jagua-rs) / sparrow):
  `strip_height` for strip packing or `"bin": {"width": W, "height": H}` for bin packing, and
  a list of items with `demand`, allowed rotations and a simple polygon.  See
  [`data/README.md`](data/README.md).
* **Plain text** - the format of the original version of this repository: `W H`, the number
  of polygons, then for each polygon its vertex count and vertices.  It is solved as a bin
  packing problem (how many `W x H` boards are needed, and how short can the last one be).

## Benchmarks

### The protocol

* sparrow was built from source (commit `67bfaf8`, v0.3.0, `--features simd,only_final_svg`,
  `target-cpu=native`) and run with its default configuration.
* For every instance, sparrow and kestrel were started **at the same moment on the same
  laptop** (Intel Core Ultra 5 125H, Windows 11), each with 3 threads, the same time limit,
  the same seed number, the same priority, reading the same instance file.  Two such pairs
  ran at a time.  Whatever else the machine was doing therefore hit both solvers alike.
* kestrel ran with `--strict`, i.e. its items are strictly separated, like sparrow's.
* Afterwards **both** layouts were checked by [`tools/verify.py`](tools/verify.py), which
  shares no code with either solver and computes all pairwise intersection areas in exact
  rational arithmetic.  All layouts of both solvers passed with an overlap area of exactly
  zero, and the densities below are the ones recomputed by that checker.
* The whole procedure is one script:
  `python tools/paired_benchmark.py --sparrow <path to sparrow> --time 1200 --seeds 1`.

### 20 minutes (the protocol of the sparrow paper)

Density in %, higher is better.

| instance | kestrel | sparrow | difference | kestrel passes sparrow's final result after | sparrow paper: expected | Q3 | best ever |
|---|---|---|---|---|---|---|---|
| albano | **89.59** | 88.89 | +0.69 | 9.0 s | 89.47 | 89.58 | 89.82 |
| dagli | **89.15** | 88.66 | +0.48 | 21 s | 89.26 | 89.50 | 90.17 |
| fu | **92.40** | 91.79 | +0.61 | 1.0 s | 92.24 | 92.40 | 92.41 |
| jakobs1 | 89.08 | 89.08 | +0.00 | - | 89.09 | 89.09 | 89.26 |
| jakobs2 | 84.18 | **84.88** | -0.71 | - | 84.77 | 85.25 | 87.73 |
| mao | **86.54** | 86.32 | +0.22 | 301 s | 86.14 | 86.45 | 86.87 |
| marques | **91.51** | 90.96 | +0.55 | 48 s | 90.93 | 91.02 | 92.02 |
| shapes0 | **68.78** | 67.90 | +0.89 | 19 s | 68.60 | 68.79 | 69.98 |
| shapes1 | **76.72** | 75.22 | +1.49 | 171 s | 75.69 | 75.97 | 76.73 |
| shapes2 (blaz1) | **84.97** | 84.13 | +0.84 | 32 s | 84.68 | 84.88 | 86.23 |
| shirts | **89.79** | 89.74 | +0.05 | 1105 s | 89.66 | 89.84 | 90.92 |
| swim | **78.03** | 77.24 | +0.79 | 325 s | 78.26 | 78.55 | 79.83 |
| trousers | **91.87** | 91.51 | +0.36 | 890 s | 91.73 | 91.93 | 92.62 |

"kestrel passes sparrow's final result after" is the moment in kestrel's run at which its
best layout became denser than the one sparrow ended with after 20 minutes.  The last three
columns are quoted from the sparrow paper: the expected density and the third quartile of 100
runs (on a faster desktop CPU, a Ryzen 9 7950X), and the best layout its authors ever
obtained.

Reading the table:

* **JAKOBS1** is a tie: both solvers end at the same layout length (11.000).
* **JAKOBS2** is sparrow's.  Its pieces have integer coordinates and its best layouts are
  jigsaws of exact fits.  With `--strict`, kestrel keeps a small distance between all items,
  as sparrow always does (the smallest gaps in sparrow's layouts are of the same size), and
  under that rule sparrow found the better layout.  kestrel's default mode lets items touch
  exactly; there it reaches the best known layout length of 22 (87.72 %) in 2 of 6 short
  runs - but a layout of exact contacts is not exactly overlap free in the checker's rational
  arithmetic (it reports an overlap of 2e-8 of an item's area), so those runs are not part of
  this comparison.
* **SHAPES1** (76.72 %) and **FU** (92.40 %) are within 0.01 points of the best layouts the
  sparrow authors report from all their experiments (76.73 % and 92.41 %).
* **SHIRTS** is close to a tie; sparrow was ahead there from the second minute to the
  eighteenth.
* On nine instances kestrel was ahead from the first seconds to the last.

<p align="center">
  <img src="results/convergence_swim.svg" width="640" alt="density over time on SWIM"/>
</p>

### 60 seconds

Same procedure with a time limit of one minute, three pairs per instance (seeds 1, 2, 3).

| instance | kestrel (mean of 3) | sparrow (mean of 3) | difference | pairs in which kestrel was denser |
|---|---|---|---|---|
| albano | **89.20** | 87.80 | +1.41 | 3 of 3 |
| dagli | **88.66** | 86.55 | +2.11 | 2 of 3 |
| fu | **91.99** | 90.57 | +1.42 | 3 of 3 |
| jakobs1 | **89.08** | 89.07 | +0.01 | 3 of 3 |
| jakobs2 | **83.94** | 82.99 | +0.95 | 2 of 3 |
| mao | **85.49** | 84.23 | +1.26 | 3 of 3 |
| marques | **90.93** | 89.75 | +1.18 | 3 of 3 |
| shapes0 | **68.00** | 66.65 | +1.35 | 3 of 3 |
| shapes1 | **75.10** | 73.65 | +1.46 | 3 of 3 |
| shapes2 (blaz1) | **85.00** | 82.80 | +2.20 | 3 of 3 |
| shirts | **89.13** | 88.82 | +0.31 | 3 of 3 |
| swim | **76.87** | 76.20 | +0.66 | 2 of 3 |
| trousers | **91.39** | 91.01 | +0.38 | 2 of 3 |

kestrel was denser on average on all thirteen instances (mean difference +1.13 points) and in
35 of the 39 pairs.

### Where the difference comes from

sparrow tests every candidate placement against the polygons around it.  kestrel precomputes,
for every pair of shapes, the exact region of relative positions in which they collide (the
no-fit polygon) together with a map of how bad the collision is, and then answers the same
question with one table lookup.  In the runs above it evaluated about six times as many
candidate placements per second as sparrow (between 3.7 and 7.9 times, depending on the
instance) - on the same search scheme, which is sparrow's.
On top of that come moves that need the exact geometry: sliding all items into exact
contact, and squeezing a layout until it jams.  [docs/ALGORITHM.md](docs/ALGORITHM.md)
describes all of it, including a list of ideas that were tried and did not help.

The layouts and solution files of both solvers are in
[`results/strip_20min/`](results/strip_20min);
`python tools/verify.py data/strip/<name>.json results/strip_20min/kestrel/<name>.json`
(or `.../sparrow/<name>.json`) repeats the check.

<p align="center">
  <img src="results/strip_20min/shirts.svg" height="170" alt="SHIRTS"/>
  <img src="results/strip_20min/albano.svg" height="170" alt="ALBANO"/>
</p>
<p align="center">
  <img src="results/strip_20min/mao.svg" height="200" alt="MAO"/>
  <img src="results/strip_20min/marques.svg" height="200" alt="MARQUES"/>
  <img src="results/strip_20min/dagli.svg" height="200" alt="DAGLI"/>
</p>

## What this is not

* **Not a proof of optimality.**  It is a heuristic; the lower bounds for these instances are
  far below any known layout.
* **Statistics are thin.**  The paper reports 100 runs per instance; the tables above rest on
  one pair of runs per instance for the 20-minute table and three for the 60-second table, because every pair of 20-minute runs costs 20 minutes of a laptop.  Single
  runs of this kind of search scatter by several tenths of a point (the paper reports
  interquartile ranges of 0.2-0.6), so differences of that size in one row mean little; the
  picture across the thirteen rows is what counts.
* **Tuned on the benchmark.**  The few parameters that are not sparrow's were chosen on these
  same instances.
* **The gap closes with time.**  kestrel is far ahead in the first minute and less far after twenty.  Its
  speed lets it reach in seconds what sparrow reaches in minutes, but the last tenths of a
  point on these instances are a matter of finding a different arrangement, and there both
  solvers search in the same way.  On SHIRTS and JAKOBS2 that was not enough to be ahead.
* **Discrete rotations only.**  Free rotation is approximated by a fixed number of angles;
  sparrow supports continuous rotation natively.
* **Memory grows with the square of the number of distinct oriented shapes** (one field per
  pair).  Instances with a few hundred distinct shapes are fine; thousands are not.
* **Simple polygons, rectangular containers.**  Holes inside items are ignored (treated as
  solid), containers with irregular outlines or defects are not supported.
* **Bin packing has no external comparison here.**  It is tested (perfect packings are
  found, the result is verified exactly), but it was not benchmarked against the bin packing
  literature.

## Repository layout

```
src/geom.rs       points, polygons, hulls, offsets, simplification
src/decomp.rs     convex decomposition (ear clipping + Hertel-Mehlhorn)
src/nfp.rs        exact no-fit-polygon fields
src/problem.rs    instances, oriented shapes, the field table
src/search.rs     placement search: sampling, coordinate descent, snap
src/separator.rs  guided local search, push relaxation
src/solve.rs      strip and bin packing drivers, compaction, squeeze
src/pool.rs       spin-waiting fork-join pool
tools/verify.py            independent exact checker (kestrel and sparrow solution files)
tools/paired_benchmark.py  the side-by-side benchmark against sparrow
tools/bench.py             deterministic A/B runs of kestrel settings
data/             benchmark instances
results/          the layouts behind the tables above
```

`cargo test --release` runs the unit and end-to-end tests; `kestrel selftest <instance>`
cross-checks the fields of an instance against an independent oracle.

## History

The first version of this repository (2023) wrapped every polygon in its bounding box and
scanned the board one unit at a time for a free spot; it is still in the git history.  On
its own test file (`data/legacy/in.txt`: 57 random convex polygons on a 700 x 700 board) its
boxes spread over the whole board (699 x 658).  kestrel packs the same polygons into the
first 432 units of the board's length in 30 seconds - 62 % of the board.

## References and credits

* J. Gardeyn, G. Vanden Berghe, T. Wauters, *An open-source heuristic to reboot 2D nesting
  research*, EJOR 2026 ([arXiv:2509.13329](https://arxiv.org/abs/2509.13329)).  kestrel's
  search loop (sampling moves, guided local search weights, exploration and compression
  phases) follows this paper, and its collision measure is the paper's overlap proxy,
  tabulated.  That sparrow is open source is what made a same-machine comparison possible at
  all.
* The benchmark instances are from [ESICUP](https://www.euro-online.org/websites/esicup/data-sets/),
  in the JSON conversion shipped with sparrow (MIT License).
* Further references are listed in [docs/ALGORITHM.md](docs/ALGORITHM.md).

MIT License.
