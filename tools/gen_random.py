#!/usr/bin/env python3
"""Generate a random instance in the plain-text format of the original project.

    python tools/gen_random.py [--n 57] [--board 700 700] [--size 100] [--points 10] [--seed 1] > in.txt

Every polygon is the convex hull of `--points` random points in a `--size` x `--size` box
(the same recipe as the original `数据随机生成.cpp`).  Output format:

    W H
    n
    k            <- number of vertices of the first polygon
    x y          <- k lines
    ...
"""
import argparse
import random


def hull(pts):
    pts = sorted(set(pts))
    if len(pts) < 3:
        return pts

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower, upper = [], []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=57)
    ap.add_argument("--board", type=float, nargs=2, default=[700, 700])
    ap.add_argument("--size", type=float, default=100.0)
    ap.add_argument("--points", type=int, default=10)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    print(f"{a.board[0]:g} {a.board[1]:g}")
    print(a.n)
    for _ in range(a.n):
        h = []
        while len(h) < 3:
            h = hull([(round(rng.uniform(0, a.size), 2), round(rng.uniform(0, a.size), 2)) for _ in range(a.points)])
        print(len(h))
        for x, y in h:
            print(f"{x:.2f} {y:.2f}")
        print()


if __name__ == "__main__":
    main()
