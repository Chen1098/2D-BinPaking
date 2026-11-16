# Results

The layouts behind the tables in the main README.

```
strip_20min/<instance>.svg            kestrel's layout (drawing)
strip_20min/kestrel/<instance>.json   kestrel's layout (rotation and translation of every item)
strip_20min/sparrow/<instance>.json   sparrow's layout from the same paired run
convergence_swim.svg                  density over time of both solvers on SWIM
```

All of them come from one run of

```
python tools/paired_benchmark.py --sparrow <path to sparrow> --time 1200 --seeds 1
```

on an Intel Core Ultra 5 125H laptop (Windows 11): for every instance one sparrow process and
one kestrel process, started at the same moment, 3 threads and 20 minutes each, two such pairs
at a time.  sparrow: commit `67bfaf8` (v0.3.0), default configuration.  kestrel: commit
`19cb992` of this repository, `--strict`.

Any file can be re-checked with the independent verifier, which reads both formats:

```
python tools/verify.py data/strip/swim.json results/strip_20min/kestrel/swim.json
python tools/verify.py data/strip/swim.json results/strip_20min/sparrow/swim.json
```

It prints the exact overlap area of the worst pair (0 for every file here), the strip length
and the density.
