# Instances

## `strip/` - the classic irregular strip packing benchmark

The thirteen instances used in the literature to compare nesting heuristics
(ALBANO, DAGLI, FU, JAKOBS1/2, MAO, MARQUES, SHAPES0/1, SHAPES2 = `blaz1`, SHIRTS, SWIM,
TROUSERS).  The shapes originate from the [ESICUP](https://www.euro-online.org/websites/esicup/data-sets/)
data sets; the JSON files are the conversions shipped with
[sparrow](https://github.com/JeroenGar/sparrow) (MIT License, Copyright (c) 2025 Jeroen
Gardeyn, KU Leuven) and are included unchanged so that both solvers read byte-identical
input.

Note that, exactly as in the sparrow paper, the strip height of six instances with integer
coordinates (`blaz1`, `fu`, `jakobs1`, `jakobs2`, `shapes0`, `shapes1`) is inflated by 0.01 %
(e.g. 40 -> 40.004).  The files are kept that way to compare on equal terms; all densities
in this repository are computed with the height that is written in the file.

Format (jagua-rs / sparrow style):

```json
{
  "name": "albano",
  "strip_height": 4900.0,
  "items": [
    {"id": 0, "demand": 2,
     "orientation": {"rotation": {"mode": "discrete", "angles": [0, 180]}},
     "shape": {"type": "simple_polygon", "data": [[0.0, 86.0], [966.0, 142.0], "..."]}}
  ]
}
```

## `bin/` - bin packing example

`trousers_bins.json` is TROUSERS with the strip replaced by bins:
`"bin": {"width": 100, "height": 79}`.

## `legacy/` - the original input of this repository

`in.txt` is the test file of the 2023 version of this project: a 700 x 700 board and 57
random convex polygons, in the plain-text format

```
W H
n
k          <- vertices of polygon 1
x y        <- k lines
...
```

`tools/gen_random.py` generates more instances of this kind.
