#!/usr/bin/env python3
"""Independent solution checker (no third-party packages, exact rational arithmetic).

    python tools/verify.py <instance.json|in.txt> <solution.json>

Reads the solution files written by kestrel and those written by sparrow / jagua-rs.

The checker shares no code with the solver.  It rebuilds every placed polygon from the
*instance* file and the (rotation, translation) pairs of the *solution* file, then

  * checks that every item of the instance is placed exactly `demand` times,
  * checks that every polygon lies inside its bin / the strip,
  * computes the exact intersection area of every pair of polygons (triangulation by ear
    clipping + exact convex clipping, all in `fractions.Fraction`),
  * recomputes the strip length and the packing density.

Exit code 0 means the layout is feasible within the tolerance printed at the end.
"""
import json
import math
import sys
from fractions import Fraction as F


def load_instance(path):
    txt = open(path, encoding="utf-8").read()
    if path.lower().endswith(".json"):
        d = json.loads(txt, parse_float=F, parse_int=F)
        items = {}
        for k, it in enumerate(d["items"]):
            data = it["shape"]["data"]
            if isinstance(data, dict):
                if "outer" in data:
                    data = data["outer"]
                else:
                    x0, y0 = data.get("x_min", F(0)), data.get("y_min", F(0))
                    w, h = data["width"], data["height"]
                    data = [[x0, y0], [x0 + w, y0], [x0 + w, y0 + h], [x0, y0 + h]]
            items[int(it.get("id", k))] = (int(it.get("demand", 1)), [(F(p[0]), F(p[1])) for p in data])
        if "bin" in d:
            return items, F(d["bin"]["height"]), F(d["bin"]["width"])
        return items, F(d["strip_height"]), None
    tok = txt.split()
    w, h, n = F(tok[0]), F(tok[1]), int(tok[2])
    pos = 3
    items = {}
    for i in range(n):
        k = int(tok[pos]); pos += 1
        pts = []
        for _ in range(k):
            pts.append((F(tok[pos]), F(tok[pos + 1]))); pos += 2
        items[i] = (1, pts)
    return items, h, w


def rot(p, deg):
    d = deg % 360
    x, y = p
    if d == 0:
        return (x, y)
    if d == 90:
        return (-y, x)
    if d == 180:
        return (-x, -y)
    if d == 270:
        return (y, -x)
    c, s = F(math.cos(math.radians(float(d)))), F(math.sin(math.radians(float(d))))
    return (c * x - s * y, s * x + c * y)


def area2(poly):
    return sum(poly[i][0] * poly[(i + 1) % len(poly)][1] - poly[(i + 1) % len(poly)][0] * poly[i][1] for i in range(len(poly)))


def orient(a, b, c):
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def triangulate(poly):
    """Ear clipping of a simple polygon (made counter-clockwise first)."""
    pts = list(poly)
    if area2(pts) < 0:
        pts.reverse()
    idx = list(range(len(pts)))
    tris = []
    guard = 0
    while len(idx) > 3 and guard < 100000:
        guard += 1
        m = len(idx)
        clipped = False
        for k in range(m):
            ia, ib, ic = idx[k - 1], idx[k], idx[(k + 1) % m]
            a, b, c = pts[ia], pts[ib], pts[ic]
            cr = orient(a, b, c)
            if cr == 0:
                idx.pop(k)
                clipped = True
                break
            if cr < 0:
                continue
            ok = True
            for j in idx:
                if j in (ia, ib, ic):
                    continue
                v = pts[j]
                if v in (a, b, c):
                    continue
                if orient(a, b, v) >= 0 and orient(b, c, v) >= 0 and orient(c, a, v) >= 0:
                    ok = False
                    break
            if ok:
                tris.append((a, b, c))
                idx.pop(k)
                clipped = True
                break
        if not clipped:
            raise ValueError("polygon is not simple")
    if len(idx) == 3:
        a, b, c = pts[idx[0]], pts[idx[1]], pts[idx[2]]
        if orient(a, b, c) > 0:
            tris.append((a, b, c))
    return tris


def clip_area2(tri_a, tri_b):
    """Twice the area of the intersection of two counter-clockwise triangles (exact)."""
    out = list(tri_a)
    for i in range(3):
        p, q = tri_b[i], tri_b[(i + 1) % 3]
        if not out:
            return F(0)
        inp, out = out, []
        for j in range(len(inp)):
            cur, nxt = inp[j], inp[(j + 1) % len(inp)]
            sc, sn = orient(p, q, cur), orient(p, q, nxt)
            if sc >= 0:
                out.append(cur)
            if (sc > 0 and sn < 0) or (sc < 0 and sn > 0):
                t = sc / (sc - sn)
                out.append((cur[0] + t * (nxt[0] - cur[0]), cur[1] + t * (nxt[1] - cur[1])))
    if len(out) < 3:
        return F(0)
    return abs(area2(out))


def bbox(pts):
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    return min(xs), min(ys), max(xs), max(ys)


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    items, height, bin_w = load_instance(sys.argv[1])
    sol = json.loads(open(sys.argv[2], encoding="utf-8").read(), parse_float=F, parse_int=F)
    if "solution" in sol:
        # a jagua-rs / sparrow solution file
        claimed = sol["solution"].get("strip_width", 0)
        placements = [
            {"item_id": p["item_id"], "rotation": p["transformation"]["rotation"], "translation": p["transformation"]["translation"]}
            for p in sol["solution"]["layout"]["placed_items"]
        ]
    else:
        claimed = sol.get("strip_length", 0)
        placements = sol["placements"]
    placed = []
    count = {}
    for pl in placements:
        iid = int(pl["item_id"])
        count[iid] = count.get(iid, 0) + 1
        tx, ty = F(pl["translation"][0]), F(pl["translation"][1])
        poly = [rot(p, pl["rotation"]) for p in items[iid][1]]
        poly = [(x + tx, y + ty) for x, y in poly]
        placed.append((iid, int(pl.get("bin", 0)), poly))
    ok = True
    for iid, (demand, _) in items.items():
        if count.get(iid, 0) != demand:
            print(f"FAIL: item {iid} placed {count.get(iid, 0)} times, demand {demand}")
            ok = False
    scale = max(height, max(max(bbox(p)[2] - bbox(p)[0], bbox(p)[3] - bbox(p)[1]) for _, _, p in placed))
    tol_len = scale * F(1, 10**7)
    nbins = max(b for _, b, _ in placed) + 1
    used = [F(0)] * nbins
    max_out = F(0)
    for iid, b, poly in placed:
        x0, y0, x1, y1 = bbox(poly)
        used[b] = max(used[b], x1)
        out = max(-x0, -y0, y1 - height, (x1 - bin_w) if bin_w is not None else F(0))
        max_out = max(max_out, out)
    if max_out > tol_len:
        print(f"FAIL: an item sticks out of the container by {float(max_out):.3e}")
        ok = False
    tris = [triangulate(p) for _, _, p in placed]
    areas = [abs(area2(p)) / 2 for _, _, p in placed]
    for a, t in zip(areas, tris):
        ta = sum(abs(area2(x)) for x in t) / 2
        if ta != a:
            print("FAIL: triangulation does not reproduce the polygon area (non-simple polygon?)")
            ok = False
    boxes = [bbox(p) for _, _, p in placed]
    tboxes = [[bbox(t) for t in ts] for ts in tris]
    worst = (F(0), None)
    worst_rel = F(0)
    n = len(placed)
    pairs = 0
    for i in range(n):
        for j in range(i + 1, n):
            if placed[i][1] != placed[j][1]:
                continue
            a, b = boxes[i], boxes[j]
            if a[2] <= b[0] or b[2] <= a[0] or a[3] <= b[1] or b[3] <= a[1]:
                continue
            pairs += 1
            s = F(0)
            for ta, ba in zip(tris[i], tboxes[i]):
                for tb, bb in zip(tris[j], tboxes[j]):
                    if ba[2] <= bb[0] or bb[2] <= ba[0] or ba[3] <= bb[1] or bb[3] <= ba[1]:
                        continue
                    s += clip_area2(ta, tb)
            s /= 2
            if s > worst[0]:
                worst = (s, (i, j))
            worst_rel = max(worst_rel, s / min(areas[i], areas[j]))
    total = sum(areas)
    print(f"items placed        : {n} in {nbins} bin(s); {pairs} close pairs checked exactly")
    print(f"max outside         : {float(max_out):.3e}")
    print(f"max pair overlap    : area {float(worst[0]):.3e}  ({float(worst_rel):.3e} of the smaller item)  pair {worst[1]}")
    if worst_rel > F(1, 10**6):
        print("FAIL: overlapping items")
        ok = False
    if bin_w is None:
        length = used[0]
        dens = total / (height * length)
        print(f"strip length        : {float(length):.9f}   (solution file says {float(claimed):.9f})")
        print(f"density             : {float(dens) * 100:.4f}%")
    else:
        frac = (nbins - 1) + used[-1] / bin_w
        print(f"bins used           : {nbins}   (fractional: {float(frac):.4f})")
        print(f"density per bin     : " + ", ".join(
            f"{float(sum(a for a, (_, b, _) in zip(areas, placed) if b == k) / (height * bin_w)) * 100:.2f}%" for k in range(nbins)))
    print("RESULT:", "FEASIBLE" if ok else "INFEASIBLE", "(tolerances: outside 1e-7 x scale, overlap area 1e-6 x item area)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
