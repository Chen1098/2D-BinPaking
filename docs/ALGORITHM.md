# How kestrel works

kestrel solves 2D irregular packing ("nesting") problems: given polygons and the rotations
each of them may use, place all of them without overlap

* in a strip of fixed height and minimum length (**strip packing**), or
* in as few identical rectangular bins as possible (**bin packing**; the least filled bin is
  additionally made as short as possible).

The search follows the scheme that currently defines the state of the art for this problem -
*allow overlaps, then eliminate them* with a guided local search, as in sparrow (Gardeyn,
Vanden Berghe & Wauters, EJOR 2026) and its predecessors (Egeblad et al. 2007, Umetani et al.
2009, Imamichi et al. 2009, Elkeran 2013, Sato et al. 2019).  What is different is everything
underneath and around that loop:

| | sparrow | kestrel |
|---|---|---|
| collision test | polygon/polygon via quadtree, conservative | **exact**, O(1) lookup in a precomputed no-fit-polygon field |
| collision severity | sum over pairs of inscribed circles, evaluated per query | the same quantity **precomputed as a field** and interpolated |
| contact precision | ~0.1 % of the item size (sampling) | **exact contact** (projection onto the NFP boundary) |
| continuous moves | none (only sampled single-item moves) | **gravity compaction** and **push relaxation** of the whole layout |
| problem types | strip | strip **and bins** |
| throughput (3 threads, SWIM) | ~1.4 M placements/s on the test machine | ~8 M placements/s |

The rest of this document goes through the pieces in the order a query flows through them.

## 1. Exact no-fit-polygon fields

For a fixed shape `A` and a moving shape `B` the *no-fit polygon* NFP(A, B) is the set of
displacements `d = pos(B) - pos(A)` for which the interiors of the two shapes overlap.
Deciding "do these two placed items collide?" is then a point-in-region test, and the
distance from `d` to the boundary of the NFP is the *penetration depth*: the shortest
translation that separates the two items.

NFPs of non-convex polygons are notoriously awkward to compute robustly, so kestrel never
builds one.  Every oriented shape is decomposed into convex pieces (ear clipping +
Hertel-Mehlhorn merging, `src/decomp.rs`), and

```
NFP(A, B) = union over k, l of   A_k (+) (-B_l)        (convex Minkowski sums)
```

A `Field` (`src/nfp.rs`) stores this union *implicitly*:

1. **The convex pieces** as lists of half-planes.  `d` is inside the NFP iff it is strictly
   inside at least one piece - an exact test that needs no polygon clipping library.
2. **The exact boundary of the union as a segment soup.**  Every edge of every piece is
   clipped against all other pieces; what survives is exactly the boundary.  The pieces are
   treated as *open* sets, so zero-width channels (a tab sliding into a slot of exactly the
   same width) and single touching points survive as boundary - exact fits stay feasible.
3. **A uniform grid** over the bounding box of the NFP.  Every cell is classified once:
   * `OUT`  - no piece reaches the cell: the answer is "no collision", after one byte lookup;
   * `IN`   - some piece contains the whole cell: collision, depth by bilinear interpolation
     of exact node values;
   * `BAND` - the boundary crosses the cell: the cell carries the few half-planes that
     matter inside it and the few boundary segments that can be nearest to a point in it;
   * `near` - an `IN` cell close to the boundary: always colliding, but with the exact
     segment list too, so depths are exact wherever they are small.

The collision *decision* is therefore exact everywhere, and the depth is exact wherever the
boundary is near.  Neither depends on the grid resolution; the grid only decides how many
half-planes and segments one query looks at.  A query costs 20-70 ns.

All fields are built lazily and in parallel; the 210 fields of SWIM take about 0.2 s.
`kestrel selftest <instance>` compares the fields against an independent separating-axis
oracle (400 000 random displacements; zero disagreements on all benchmark instances).

Tolerance: two items count as colliding when one convex piece penetrates another by more than
`1e-9 x` the problem scale.  With `--strict` (or `--min-sep d`) the collision shapes are
inflated instead, and the resulting layouts are overlap free in exact rational arithmetic.

## 2. What a collision costs: the proxy field

The guided local search needs a number for "how bad is this collision".  The obvious choice
with NFPs at hand - the exact penetration depth - turns out to be a mediocre guide.  At equal
numbers of evaluated placements, sparrow's measure (the summed, softened overlap of the two
shapes' inscribed circles, which behaves like the square root of the overlap *area* scaled
by the shape sizes) eliminates collisions clearly better:

| pair loss (150 M evaluations, SHIRTS / SWIM / TROUSERS, 4 seeds each) | mean density |
|---|---|
| penetration depth with floor | 84.88 % |
| circle-overlap proxy, computed per query (7x slower) | 85.42 % |
| **circle-overlap proxy as a field (kestrel)** | **85.52 %** |

But that proxy is a function of the relative displacement only - so it can be tabulated.
kestrel evaluates it once at the nodes of each NFP grid and interpolates.  A query is then:
one exact collision decision (section 1), and if it collides, one bilinear lookup.  This is
both the better guide and the faster one, and it is the single most important reason kestrel
gets more out of every second than a per-query evaluation can.

`--set loss_mode=0` switches to the depth-based loss, `--set loss_mode=1` evaluates the proxy
directly (both kept for ablation).

## 3. Eliminating collisions

`src/separator.rs`, `src/search.rs`.  For a fixed container size:

* Every colliding item, in random order, is moved to the best placement found by
  **sampling** (the current position, 25 samples around it, 50 anywhere in the container, in
  any allowed rotation), **coordinate descent** on the three most promising samples, and a
  final finer descent.  A placement is scored by the weighted sum of its pair losses and
  evaluation stops as soon as the running sum exceeds the best known score.
* **Guided local search**: after every pass the weight of each pair that still collides is
  multiplied (x1.2 ... x2, by severity), all others decay by 5 %.  Items that keep fighting
  over the same spot eventually jump elsewhere.
* Several workers run the pass from the same state with different random orders; the best
  result is kept.  The workers live in a small spin-waiting pool (`src/pool.rs`): one pass is
  only a few hundred microseconds, which is less than the wake-up latency of a general
  purpose thread pool.
* The best state seen is restored after 200 passes without progress; three such attempts
  without progress count as failure.

These rules and their parameters are sparrow's.  Two kestrel additions make use of the exact
geometry:

* **Snap.**  Sampling and descent only get close to a contact.  The final placement is
  projected onto the exact NFP boundary of the neighbours it collides with (alternating
  between them for corners).  This is what makes exact fits reachable at all.
* **Finisher.**  Whenever a new best state has only shallow overlaps left, the push relaxation
  (below) tries to dissolve them into whatever slack the layout has.

## 4. Moving the whole layout

Single-item moves are bad at one thing: using slack that is spread thinly over many gaps.
kestrel has two deterministic, sampling-free operations for that.

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

After every successful separation the layout is compacted and squeezed.  The stochastic
search is thus only ever asked the hard question - *find a different arrangement* - and never
spends its time re-discovering slack that a deterministic procedure can collect.

## 5. The driver

`src/solve.rs`.

**Strip packing.**  A bottom-left-fill construction gives a first layout.  *Exploration*
(80 % of the budget): shrink the strip by 0.1 %, eliminate the collisions, tighten, repeat.
When the elimination fails, the least infeasible layouts are kept in a pool; one of the
better ones is restored, two large items are swapped (items that end up covered by a swapped
item move along with it), and the elimination is tried again.  *Compression* (20 %): tiny
shrink steps, always restarting from the best layout.  Finally the layout is squeezed at full
precision.

**Bin packing.**  Construction fills bins first-fit.  *Phase A* dissolves the emptiest bin:
its items are dropped into the other bins and the collisions are eliminated; repeated until
it fails or the area bound is reached.  *Phase B* treats the emptiest bin as a strip and
shortens it with the strip exploration - items migrate into the other bins as it shrinks, and
if it runs empty it is dropped.  All bins live side by side in one layout, so every move
(sampling, swaps, relaxation) works across bins without special cases.

## 6. Reproducibility and verification

* `--evals N` replaces the time limit by a budget of N million evaluated placements.  The run
  is then bit-for-bit deterministic and independent of machine load; all tuning in this
  repository was done this way.
* Every result is checked twice: by a separating-axis test on the convex pieces inside the
  solver (independent of the fields), and by `tools/verify.py`, which shares no code with the
  solver and recomputes all pairwise intersection areas in exact rational arithmetic from
  the instance file and the solution file alone.

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
