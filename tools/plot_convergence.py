#!/usr/bin/env python3
"""Draw density-over-time curves as an SVG (no third-party packages).

    python tools/plot_convergence.py out.svg "title" label1=solution1.json [label2=sparrow_log.txt ...]

A kestrel solution file contributes its `history` (seconds, strip length); a sparrow log
contributes its "[EXPL] feasible solution found" / "[CMPR] success" lines.
"""
import json
import re
import sys

COLORS = ["#c2410c", "#1d4ed8", "#15803d", "#7e22ce"]


def load(path):
    txt = open(path, encoding="utf-8", errors="replace").read()
    pts = []
    if path.lower().endswith(".json"):
        d = json.loads(txt)
        total = d["density"] * d["strip_length"]  # density * length is constant
        for t, length in d.get("history", []):
            pts.append((t, 100.0 * total / length))
        pts.append((d["seconds"], 100.0 * d["density"]))
    else:
        for l in txt.splitlines():
            m = re.match(r"\[INFO\] \[(\d+):(\d+):(\d+)\].*(?:feasible solution found! \(width: [\d.]+, dens: ([\d.]+)%\)|\[CMPR\] success at [\d.]+% \([\d.]+ \| ([\d.]+)%\))", l)
            if m:
                t = int(m[1]) * 3600 + int(m[2]) * 60 + int(m[3])
                pts.append((t, float(m[4] or m[5])))
    best = 0.0
    out = []
    for t, d in sorted(pts):
        if d > best:
            best = d
            out.append((t, d))
    return out


def main():
    if len(sys.argv) < 4:
        print(__doc__)
        return 2
    out, title = sys.argv[1], sys.argv[2]
    series = []
    for a in sys.argv[3:]:
        label, path = a.split("=", 1)
        series.append((label, load(path)))
    tmax = max(p[-1][0] for _, p in series if p)
    dmax = max(p[-1][1] for _, p in series if p)
    # start the y axis where the slower solver is after 1% of the time
    dmin = min(max((d for t, d in p if t <= 0.02 * tmax), default=p[0][1]) for _, p in series if p)
    dmin = dmin - 0.1 * (dmax - dmin)
    W, H, L, R, T, B = 760, 380, 62, 20, 36, 44
    X = lambda t: L + (W - L - R) * t / tmax
    Y = lambda d: T + (H - T - B) * (1 - (d - dmin) / (dmax + 0.05 * (dmax - dmin) - dmin))
    o = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" font-family="sans-serif" font-size="12">',
         f'<rect width="{W}" height="{H}" fill="white"/>',
         f'<text x="{W/2}" y="20" text-anchor="middle" font-size="14" font-weight="bold">{title}</text>']
    step = 0.5 if dmax - dmin < 6 else 1.0
    d = (int(dmin / step) + 1) * step
    while d <= dmax + 1e-9:
        o.append(f'<line x1="{L}" y1="{Y(d):.1f}" x2="{W-R}" y2="{Y(d):.1f}" stroke="#e5e7eb"/>')
        o.append(f'<text x="{L-6}" y="{Y(d)+4:.1f}" text-anchor="end">{d:.1f}</text>')
        d += step
    for k in range(0, 5):
        t = tmax * k / 4
        o.append(f'<line x1="{X(t):.1f}" y1="{T}" x2="{X(t):.1f}" y2="{H-B}" stroke="#f3f4f6"/>')
        o.append(f'<text x="{X(t):.1f}" y="{H-B+16}" text-anchor="middle">{t/60:.0f} min</text>')
    o.append(f'<text x="14" y="{(T+H-B)/2}" text-anchor="middle" transform="rotate(-90 14 {(T+H-B)/2})">density (%)</text>')
    o.append(f'<rect x="{L}" y="{T}" width="{W-L-R}" height="{H-T-B}" fill="none" stroke="#9ca3af"/>')
    for i, (label, pts) in enumerate(series):
        if not pts:
            continue
        c = COLORS[i % len(COLORS)]
        path = []
        prev = None
        for t, dd in pts:
            if dd < dmin:
                continue
            if prev is not None:
                path.append(f"L{X(t):.1f},{Y(prev):.1f}")
            path.append(("M" if not path else "L") + f"{X(t):.1f},{Y(dd):.1f}")
            prev = dd
        path.append(f"L{X(tmax):.1f},{Y(prev):.1f}")
        o.append(f'<path d="{" ".join(path)}" fill="none" stroke="{c}" stroke-width="2"/>')
        # legend in the lower right corner, where no curve of a converging search runs
        ly = H - B - 14 - 18 * (len(series) - 1 - i)
        o.append(f'<line x1="{W-R-170}" y1="{ly-4}" x2="{W-R-142}" y2="{ly-4}" stroke="{c}" stroke-width="2"/>')
        o.append(f'<text x="{W-R-136}" y="{ly}">{label}: {pts[-1][1]:.2f}%</text>')
    o.append("</svg>")
    open(out, "w", encoding="utf-8").write("\n".join(o) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
