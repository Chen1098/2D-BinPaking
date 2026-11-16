# How kestrel works

kestrel solves 2D irregular packing ("nesting") problems: given polygons and the rotations
each of them may use, place all of them without overlap

* in a strip of fixed height and minimum length (**strip packing**), or
* in as few identical rectangular bins as possible (**bin packing**; the least filled bin is
  additionally made as short as possible).

The search follows the scheme that currently defines the state of the art for this problem -
*allow overlaps, then eliminate them* with a guided local search - as in sparrow (Gardeyn,
Vanden Berghe & Wauters, EJOR 2026) and its predecessors (Egeblad et al. 2007, Umetani et al.
2009, Imamichi et al. 2009, Elkeran 2013, Sato et al. 2019).  What is different is the
machinery underneath and around that loop:

| | sparrow | kestrel |
|---|---|---|
| collision test | polygon against polygon through a quadtree, per query | **exact** point lookup in a precomputed no-fit-polygon field |
| collision severity | sum over pairs of inscribed circles, per query | the same quantity **tabulated as a field**, interpolated |
| contacts | whatever precision the sampling reaches | **exact contact** by projection onto the NFP boundary; exact fits are feasible |
| moves of the whole layout | none (single items are re-sampled) | **gravity compaction**, **push relaxation**, **jam squeeze** |
| problem types | strip | strip **and bins** |
| candidate placements evaluated per second | 1 | about 6 (3.7 to 7.9 across the benchmark; same machine, same thread count) |

The rest of this document goes through the pieces in the order a query flows through them,
and ends with the list of things that were tried and did *not* help.

## 1. Exact no-fit-polygon fields

For a fixed shape `A` and a moving shape `B` the *no-fit polygon* NFP(A, B) is the set of
displacements `d = pos(B) - pos(A)` for which the interiors of the two shapes overlap.
Deciding "do these two placed items collide?" is then a point-in-region test, and the
distance from `d` to the boundary of the NFP is the *penetration depth*: the shortest
translation that separates the two items.

NFPs of non-convex polygons are notoriously awkward to compute robustly, so kestrel never
builds one as a polygon.  Every oriented shape is decomposed into convex pieces (ear clipping
+ Hertel-Mehlhorn merging, `src/decomp.rs`), and

```
NFP(A, B) = union over k, l of   A_k (+) (-B_l)        (convex Minkowski sums)
```

A `Field` (`src/nfp.rs`) stores this union *implicitly*:

1. **The convex pieces** as lists of half-planes.  `d` is inside the NFP iff it is strictly
   inside at least one piece - an exact test that needs no polygon clipping library.
2. **The exact boundary of the union as a segment soup.**  Every edge of every piece is
   clipped against all other pieces; what survives is exactly the boundary.  The pieces are
   treated as *open* sets, so zero-width channels (a tab sliding into a slot of exactly the
   same width) and single touching points survive as boundary: exact fits stay feasible.
3. **A uniform grid** over the bounding box of the NFP.  Every cell is classified once:
   * `OUT`  - no piece reaches the cell: the answer is "no collision" after one byte lookup;
   * `IN`   - some piece contains the whole cell: collision, depth by bilinear interpolation
     of exact node values;
   * `BAND` - the boundary crosses the cell: the cell carries the few half-planes that
     matter inside it and the few boundary segments that can be nearest to a point in it;
   * `near` - an `IN` cell close to the boundary: always colliding, but with the exact
     segment list too, so depths are exact wherever they are small.

The collision *decision* is therefore exact everywhere, and the depth is exact wherever the
boundary is near.  Neither depends on the grid resolution; the grid only decides how many
half-planes and segments one query looks at.  One thread answers roughly ten million
queries per second.

Fields are built in parallel, one per unordered pair of oriented shapes; the 210 fields of
SWIM take about 0.3 s and 32 MB.  `kestrel selftest <instance>` compares them against an
independent separating-axis oracle (400 000 random displacements; no disagreement on any of
the benchmark instances).

Tolerance: two items count as colliding when one convex piece penetrates another by more than
`1e-9 x` the problem scale.  With `--strict` (or `--min-sep d`) the collision shapes are
inflated instead, and the resulting layouts are overlap free in exact rational arithmetic.
Outlines with more than 64 vertices are first simplified, growing by at most 0.1 % of their
area (`--simplify`); the solution is still drawn and checked with the true outlines.

## 2. What a collision costs: the proxy field

The guided local search needs a number for "how bad is this collision".  The obvious choice
with NFPs at hand - the exact penetration depth - turns out to be a mediocre guide.  At equal
numbers of evaluated placements, sparrow's measure (the summed, softened overlap of the two
shapes' inscribed circles, which behaves like the square root of the overlap *area* scaled
by the shape sizes) eliminates collisions clearly better:

| pair loss (150 M evaluations per run; SHIRTS, SWIM, TROUSERS; 4 seeds each) | mean density |
|---|---|
| penetration depth with floor | 84.88 % |
| circle-overlap proxy, computed per query (7x slower per evaluation) | 85.42 % |
| **circle-overlap proxy as a field (kestrel)** | **85.52 %** |

That proxy is a function of the relative displacement only - so it can be tabulated.
kestrel evaluates it once at the nodes of each NFP grid and interpolates.  A query is then
one exact collision decision (section 1) and, if the pair collides, one bilinear lookup.
This is where the speed comes from: the guide is as good as the per-query one and costs a
fraction of it.

`--set loss_mode=0` switches to the depth-based loss, `--set loss_mode=1` evaluates the proxy
directly (both kept for ablation).

## 3. Eliminating collisions

`src/separator.rs`, `src/search.rs`.  For a fixed container size:

* Every colliding item, in random order, is moved to the best placement found by
  **sampling** (the current position, 25 samples around it, 50 anywhere in the container, in
  any allowed rotation), **coordinate descent** on the three most promising samples, and a
  final finer descent.  A placement is scored by the weighted sum of its pair losses and
  an evaluation stops as soon as the running sum exceeds the best known score.
* **Guided local search**: after every pass the weight of each pair that still collides is
  multiplied (x1.2 ... x2, by severity), all others decay by 5 %.  Items that keep fighting
  over the same spot eventually move elsewhere.
* Several workers run the pass from the same state with different random orders; the best
  result is kept.  The workers live in a small spin-waiting pool (`src/pool.rs`): one pass
  takes a few hundred microseconds, which is less than the wake-up latency of a general
  purpose thread pool.
* The best state seen is restored after 200 passes without progress; three such rounds
  without progress count as failure.

These rules and their parameters are sparrow's.  Two additions use the exact geometry:

* **Snap.**  Sampling and descent only get close to a contact.  The final placement is
  projected onto the exact NFP boundary of the neighbours it collides with (alternating
  between them for corners).  This is what makes exact fits reachable at all.
* **Finisher.**  Whenever a new best state has only shallow overlaps left, the push relaxation
  (below) tries to dissolve them into whatever slack the layout has.

## 4. Moving the whole layout

Re-sampling single items is bad at one thing: using slack that is spread thinly over many
gaps.  kestrel has deterministic, sampling-free operations for that.

* **Gravity compaction** (`compact`).  Items slide left, then down, into exact contact, in
  order of their position.
* **Push relaxation** (`relax`).  Every colliding pair is pushed apart along its exact
  penetration direction, half the depth each, Gauss-Seidel style, with positions clamped to
  the container.  Pressure travels along contact chains, items slide along slanted edges, and
  small overlaps are absorbed by slack anywhere in the layout - or the relaxation detects that
  the total overlap stopped shrinking and gives up (the layout is *jammed*).
* **Jam squeeze** (`squeeze`).  Shrink the container by ever smaller steps, relaxing each
  time, until nothing gives any more.  The result is a local minimum of the container length
  under continuous deformation of the layout.

After every successful elimination the layout is compacted and squeezed.  The sampling search
is thus mostly asked the hard question - *find a different arrangement* - instead of
re-discovering slack that a deterministic procedure can collect.  This is what makes the
first minute so much faster than sparrow's (see the README); in very long runs the advantage
shrinks, because both solvers then spend their time on the same hard question.

## 5. The driver

`src/solve.rs`.

**Strip packing.**  A bottom-left-fill construction gives a first layout.  *Exploration*
(80 % of the budget): shrink the strip by 0.1 %, eliminate the collisions, tighten, repeat.
When the elimination fails, the least infeasible layouts are kept in a pool; one of the
better ones is restored, two large items are swapped (items that end up covered by a swapped
item move along with it), and the elimination is tried again.  *Compression* (20 %): tiny
shrink steps at a random cut position, always from the best layout.  Finally the layout is
squeezed at full precision.

**Bin packing.**  Construction fills bins first-fit.  *Phase A* dissolves the emptiest bin:
its items are dropped into the other bins and the collisions are eliminated; repeated until
it fails or the area bound is reached.  *Phase B* treats the emptiest bin as a strip and
shortens it with the strip exploration - items migrate into the other bins as it shrinks, and
if it runs empty it is dropped.  All bins live side by side in one layout, so every move
(sampling, swaps, relaxation) works across bins without special cases.

## 6. Reproducibility and verification

* `--evals N` replaces the time limit by a budget of N million evaluated placements.  The run
  is then bit-for-bit deterministic and independent of machine load; everything in section 7
  was measured this way.
* Every result is checked twice: by a separating-axis test on the convex pieces inside the
  solver (independent of the fields), and by `tools/verify.py`, which shares no code with the
  solver and recomputes all pairwise intersection areas in exact rational arithmetic from
  the instance file and the solution file alone.  The checker also reads sparrow's solution
  files, so both sides of the comparison in the README went through it.

## 7. What was tried and did not help

Single runs of this kind of search are noisy: with a fixed budget of 500 M evaluations and
nothing changed but the seed, MAO ends anywhere between 85.4 % and 86.6 % and DAGLI between
87.6 % and 90.0 % (8 seeds each).  Most of a long run is spent *stuck* - a typical trajectory
on a small instance makes its last improvement within the first few hundred million
evaluations.  The obvious question is how to use the stuck time.  All of the following were
implemented and measured against the default on the same seeds; none made a difference that
survived the noise, and none is on by default.  (The restarts are still available as
`--set restart_fails=30`; the code of the others is in the commit tagged `exp-snapshot`.)

| idea | measured (mean density, default -> variant) |
|---|---|
| restart from a fresh construction after 30 failed eliminations in a row | better on five instances, worse on four: ALBANO 88.76 -> 89.24, JAKOBS2 87.17 -> 87.72, SHAPES1 76.00 -> 76.36, MAO 85.90 -> 86.13, SHIRTS 89.38 -> 89.60, but SHAPES2 84.69 -> 84.36, DAGLI 89.04 -> 88.78, TROUSERS 91.93 -> 91.70, MARQUES 90.96 -> 90.89 (1500 M or 3000 M evaluations, 2 to 4 seeds) |
| halve the shrink step after a failure, double it after a success | MAO 85.91 -> 85.73, DAGLI 88.41 -> 88.18 (8 seeds) |
| cut the layout at a random position instead of the middle when shrinking | MAO 86.01 -> 85.74, DAGLI 88.30 -> 88.31 (6 seeds) |
| give up an elimination early when it looks worse than earlier failed ones | MAO 86.01 -> 85.40, DAGLI 88.30 -> 87.92 (6 seeds) |
| basin hopping: loosen by 3 %, swap, tighten, accept if within 1 % of the best | MAO 86.01 -> 86.05, DAGLI 88.30 -> 88.33 (6 seeds) |
| twice as many samples and descents per move | SHIRTS 88.87 -> 88.63 (600 M, 6 seeds) |
| half of the budget for the compression phase instead of 20 % | SHIRTS 89.38 -> 89.38, TROUSERS 91.93 -> 91.98 (3000 M, 2 seeds) |
| explore without tightening (tighten only the end result) | SHIRTS 88.87 -> 89.19 (600 M, 6 seeds), but DAGLI 89.24 -> 88.33, TROUSERS 91.70 -> 91.20, ALBANO 88.49 -> 88.95 (2 seeds); 3000 M: SHIRTS 89.38 -> 89.54, TROUSERS 91.93 -> 91.91 |

(These runs let items touch, the default; the comparison in the README uses `--strict`.)
The one decision that did matter is the one in section 2: which quantity the search
minimises.

One more thing was learned the hard way.  On JAKOBS2, whose pieces have integer coordinates,
the default mode reaches the best known length 22 in a third of all short runs, because its
layouts may contain exact fits; with any positive separation between the items - the rule
sparrow plays by - the same search gets stuck at length 23 almost every time, and sparrow
does better there.

## References

* J. Gardeyn, G. Vanden Berghe, T. Wauters. *An open-source heuristic to reboot 2D nesting
  research.* European Journal of Operational Research, 2026. arXiv:2509.13329.
* A. K. Sato, T. C. Martins, A. M. Gomes, M. S. G. Tsuzuki. *Raster penetration map applied
  to the irregular packing problem.* EJOR 279(2), 2019.
* A. Elkeran. *A new approach for sheet nesting problem using guided cuckoo search and
  pairwise clustering.* EJOR 231(3), 2013.
* T. Imamichi, M. Yagiura, H. Nagamochi. *An iterated local search algorithm based on
  nonlinear programming for the irregular strip packing problem.* Discrete Optimization
  6(4), 2009.
* S. Umetani et al. *Solving the irregular strip packing problem via guided local search for
  overlap minimization.* ITOR 16(6), 2009.
* J. Egeblad, B. K. Nielsen, A. Odgaard. *Fast neighborhood search for two- and
  three-dimensional nesting problems.* EJOR 183(3), 2007.
* J. A. Bennell, J. F. Oliveira. *The geometry of nesting problems: a tutorial.* EJOR
  184(2), 2008.
* C. Voudouris, E. Tsang, A. Alsheddy. *Guided local search.* Handbook of Metaheuristics,
  2010.
