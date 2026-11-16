//! Problem definition, instance parsing and the lazily built NFP field table.

use crate::decomp::convex_decompose;
use crate::geom::*;
use crate::nfp::{build_field, ConvexPoly, Field, ProxySpec};
use serde_json::Value;
use std::sync::OnceLock;

/// One item type in one allowed orientation, translated so its bounding box starts at (0, 0).
pub struct Shape {
    pub type_id: usize,
    pub angle: f64,
    /// collision shape (the item outline, inflated by half the minimum separation)
    pub poly: Vec<P>,
    /// the true item outline in the same frame
    pub outline: Vec<P>,
    pub pieces: Vec<Vec<P>>,
    pub cpieces: Vec<ConvexPoly>,
    pub w: f64,
    pub h: f64,
    /// bounding-box minimum of the rotated original polygon: final = R(angle) * orig + (pos - off)
    pub off: P,
    pub area: f64,
    pub centroid: P,
    /// sqrt of the convex hull area
    pub size: f64,
    pub hull_area: f64,
    pub diam: f64,
    /// inscribed circles (x, y, r), used by the optional pole-proxy loss
    pub poles: Vec<[f64; 3]>,
}

pub struct ItemType {
    pub ext_id: i64,
    pub demand: usize,
    pub orig: Vec<P>,
    /// oriented shape ids, one per allowed (deduplicated) rotation
    pub shapes: Vec<usize>,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Mode {
    Strip,
    Bin,
}

pub struct Problem {
    pub name: String,
    pub mode: Mode,
    /// fixed extent along y
    pub height: f64,
    /// extent of one bin along x (bin mode only)
    pub bin_width: f64,
    pub types: Vec<ItemType>,
    pub shapes: Vec<Shape>,
    /// type index of every item instance
    pub items: Vec<usize>,
    pub tol: f64,
    pub res: usize,
    /// build overlap-proxy fields (needed by the proxy-field loss)
    pub proxy: bool,
    pub total_area: f64,
    pub max_dim: f64,
    fields: Vec<OnceLock<Field>>,
}

pub struct RawItem {
    pub ext_id: i64,
    pub demand: usize,
    pub poly: Vec<P>,
    pub angles: Vec<f64>,
}

pub struct RawInstance {
    pub name: String,
    /// minimum distance between items requested by the instance (0 if absent)
    pub min_sep: f64,
    pub mode: Mode,
    pub height: f64,
    pub bin_width: f64,
    pub items: Vec<RawItem>,
}

fn parse_points(v: &Value) -> Result<Vec<P>, String> {
    let arr = v.as_array().ok_or("polygon data must be an array of [x, y]")?;
    let mut out = Vec::with_capacity(arr.len());
    for q in arr {
        let xy = q.as_array().ok_or("point must be [x, y]")?;
        if xy.len() < 2 {
            return Err("point must be [x, y]".into());
        }
        out.push(p(xy[0].as_f64().ok_or("bad x")?, xy[1].as_f64().ok_or("bad y")?));
    }
    Ok(out)
}

/// Parse the jagua-rs / sparrow style JSON instance format (plus a small `bin` extension).
pub fn parse_json(text: &str, continuous_steps: usize) -> Result<RawInstance, String> {
    let v: Value = serde_json::from_str(text).map_err(|e| format!("invalid JSON: {e}"))?;
    let name = v.get("name").and_then(|x| x.as_str()).unwrap_or("instance").to_string();
    let mut mode = Mode::Strip;
    let mut height = v.get("strip_height").and_then(|x| x.as_f64()).unwrap_or(0.0);
    let mut bin_width = 0.0;
    if let Some(b) = v.get("bin") {
        mode = Mode::Bin;
        bin_width = b.get("width").and_then(|x| x.as_f64()).ok_or("bin.width missing")?;
        height = b.get("height").and_then(|x| x.as_f64()).ok_or("bin.height missing")?;
    } else if let Some(bins) = v.get("bins").and_then(|x| x.as_array()) {
        if let Some(b) = bins.first() {
            let d = b.get("shape").and_then(|s| s.get("data")).ok_or("bins[0].shape.data missing")?;
            if let (Some(w), Some(h)) = (d.get("width").and_then(|x| x.as_f64()), d.get("height").and_then(|x| x.as_f64())) {
                mode = Mode::Bin;
                bin_width = w;
                height = h;
            } else if let Ok(pts) = parse_points(d) {
                let bb = bbox(&pts);
                mode = Mode::Bin;
                bin_width = bb[2] - bb[0];
                height = bb[3] - bb[1];
            }
        }
    }
    if !(height > 0.0) {
        return Err("instance needs a positive `strip_height` or a `bin` {width, height}".into());
    }
    let items_v = v.get("items").and_then(|x| x.as_array()).ok_or("`items` missing")?;
    let mut items = Vec::new();
    for (k, it) in items_v.iter().enumerate() {
        let ext_id = it.get("id").and_then(|x| x.as_i64()).unwrap_or(k as i64);
        let demand = it.get("demand").and_then(|x| x.as_u64()).unwrap_or(1) as usize;
        let shape = it.get("shape").ok_or("item.shape missing")?;
        let data = shape.get("data").ok_or("item.shape.data missing")?;
        let poly = if data.is_array() {
            parse_points(data)?
        } else if let Some(o) = data.get("outer") {
            parse_points(o)?
        } else if let (Some(w), Some(h)) = (data.get("width").and_then(|x| x.as_f64()), data.get("height").and_then(|x| x.as_f64())) {
            let x0 = data.get("x_min").and_then(|x| x.as_f64()).unwrap_or(0.0);
            let y0 = data.get("y_min").and_then(|x| x.as_f64()).unwrap_or(0.0);
            vec![p(x0, y0), p(x0 + w, y0), p(x0 + w, y0 + h), p(x0, y0 + h)]
        } else {
            return Err("unsupported item shape".into());
        };
        let mut angles: Vec<f64> = Vec::new();
        if let Some(rot) = it.get("orientation").and_then(|o| o.get("rotation")) {
            match rot.get("mode").and_then(|m| m.as_str()).unwrap_or("discrete") {
                "discrete" => {
                    if let Some(a) = rot.get("angles").and_then(|a| a.as_array()) {
                        angles = a.iter().filter_map(|x| x.as_f64()).collect();
                    }
                }
                "stepped" => {
                    let step = rot.get("step").and_then(|x| x.as_f64()).unwrap_or(90.0);
                    let mut a = 0.0;
                    while a < 360.0 - 1e-9 {
                        angles.push(a);
                        a += step;
                    }
                }
                "continuous" => {
                    let n = continuous_steps.max(1);
                    angles = (0..n).map(|i| 360.0 * i as f64 / n as f64).collect();
                }
                _ => {}
            }
        } else if let Some(a) = it.get("allowed_orientations").and_then(|a| a.as_array()) {
            angles = a.iter().filter_map(|x| x.as_f64()).collect();
        }
        if angles.is_empty() {
            angles.push(0.0);
        }
        items.push(RawItem { ext_id, demand, poly, angles });
    }
    let min_sep = v.get("min_item_separation").and_then(|x| x.as_f64()).unwrap_or(0.0);
    Ok(RawInstance { name, min_sep, mode, height, bin_width, items })
}

/// Parse the original plain-text format of this repository:
/// `W H`, `n`, then for every polygon `k` followed by `k` lines `x y`.
pub fn parse_legacy(text: &str) -> Result<RawInstance, String> {
    let mut it = text.split_whitespace().map(|t| t.parse::<f64>());
    let mut next = |what: &str| -> Result<f64, String> {
        match it.next() {
            Some(Ok(v)) => Ok(v),
            _ => Err(format!("legacy input: expected {what}")),
        }
    };
    let w = next("board width")?;
    let h = next("board height")?;
    let n = next("polygon count")? as usize;
    let mut items = Vec::with_capacity(n);
    for i in 0..n {
        let k = next("vertex count")? as usize;
        let mut poly = Vec::with_capacity(k);
        for _ in 0..k {
            let x = next("x")?;
            let y = next("y")?;
            poly.push(p(x, y));
        }
        items.push(RawItem { ext_id: i as i64, demand: 1, poly, angles: vec![0.0] });
    }
    Ok(RawInstance { name: "legacy".into(), min_sep: 0.0, mode: Mode::Bin, height: h, bin_width: w, items })
}

fn same_polygon(a: &[P], b: &[P], tol: f64) -> bool {
    if a.len() != b.len() {
        return false;
    }
    let n = a.len();
    'shift: for s in 0..n {
        for i in 0..n {
            let q = b[(i + s) % n];
            if (a[i].x - q.x).abs() > tol || (a[i].y - q.y).abs() > tol {
                continue 'shift;
            }
        }
        return true;
    }
    false
}

impl Problem {
    /// `min_sep`: minimum distance to keep between any two items (0 = items may touch).
    pub fn new(raw: RawInstance, res: usize, rot_override: Option<&[f64]>, proxy: bool, min_sep: f64) -> Result<Problem, String> {
        let mut types = Vec::new();
        let mut shapes: Vec<Shape> = Vec::new();
        let mut items = Vec::new();
        let mut total_area = 0.0;
        let mut max_dim: f64 = 0.0;
        for r in &raw.items {
            let orig = clean_polygon(&r.poly).ok_or_else(|| format!("item {} has a degenerate polygon", r.ext_id))?;
            let bb0 = bbox(&orig);
            let dtol = 1e-9 * (bb0[2] - bb0[0]).max(bb0[3] - bb0[1]);
            let angles: &[f64] = rot_override.unwrap_or(&r.angles);
            let mut ids = Vec::new();
            for &ang in angles {
                let rp0 = rotate_poly(&orig, ang);
                let rp = if min_sep > 0.0 {
                    clean_polygon(&inflate(&rp0, 0.5 * min_sep)).ok_or_else(|| format!("item {}: cannot apply the minimum separation", r.ext_id))?
                } else {
                    rp0.clone()
                };
                let bb = bbox(&rp);
                let off = p(bb[0], bb[1]);
                let poly: Vec<P> = rp.iter().map(|q| *q - off).collect();
                let outline: Vec<P> = rp0.iter().map(|q| *q - off).collect();
                let (w, h) = (bb[2] - bb[0], bb[3] - bb[1]);
                if h > raw.height * (1.0 + 1e-12) {
                    continue;
                }
                if raw.mode == Mode::Bin && w > raw.bin_width * (1.0 + 1e-12) {
                    continue;
                }
                if ids.iter().any(|&s: &usize| same_polygon(&shapes[s].poly, &poly, dtol)) {
                    continue;
                }
                let pieces = convex_decompose(&poly);
                let cpieces: Vec<ConvexPoly> = pieces.iter().map(|q| ConvexPoly::new(q.clone())).collect();
                let hull = convex_hull(&poly);
                let hull_area = signed_area(&hull).abs();
                let area_c = signed_area(&poly);
                let area = signed_area(&outline);
                let pa: f64 = pieces.iter().map(|q| signed_area(q)).sum();
                if (pa - area_c).abs() > 1e-7 * area_c {
                    return Err(format!(
                        "item {}: convex decomposition failed (area {} vs {}); is the polygon simple?",
                        r.ext_id, pa, area_c
                    ));
                }
                ids.push(shapes.len());
                shapes.push(Shape {
                    type_id: types.len(),
                    angle: ang,
                    centroid: centroid(&poly),
                    diam: diameter(&poly),
                    poles: poles(&poly, 64),
                    poly,
                    outline,
                    pieces,
                    cpieces,
                    w,
                    h,
                    off,
                    area,
                    size: hull_area.sqrt(),
                    hull_area,
                });
                max_dim = max_dim.max(w).max(h);
            }
            if ids.is_empty() {
                return Err(format!("item {} does not fit the container in any allowed orientation", r.ext_id));
            }
            let a = shapes[ids[0]].area;
            total_area += a * r.demand as f64;
            for _ in 0..r.demand {
                items.push(types.len());
            }
            types.push(ItemType { ext_id: r.ext_id, demand: r.demand, orig, shapes: ids });
        }
        if items.is_empty() {
            return Err("no items".into());
        }
        let ns = shapes.len();
        let mut fields = Vec::with_capacity(ns * (ns + 1) / 2);
        fields.resize_with(ns * (ns + 1) / 2, OnceLock::new);
        let scale = raw.height.max(max_dim).max(if raw.mode == Mode::Bin { raw.bin_width } else { total_area / raw.height });
        Ok(Problem {
            name: raw.name,
            mode: raw.mode,
            height: raw.height,
            bin_width: raw.bin_width,
            types,
            shapes,
            items,
            tol: 1e-9 * scale,
            res,
            proxy,
            total_area,
            max_dim,
            fields,
        })
    }

    #[inline(always)]
    pub fn n_shapes(&self) -> usize {
        self.shapes.len()
    }

    #[inline(always)]
    fn fidx(&self, a: usize, b: usize) -> usize {
        debug_assert!(a <= b);
        let ns = self.shapes.len();
        a * ns - (a * a.wrapping_sub(1)) / 2 + (b - a)
    }

    /// Field of the unordered shape pair (a <= b): displacement is pos(b) - pos(a).
    #[inline]
    pub fn field(&self, a: usize, b: usize) -> &Field {
        let i = self.fidx(a, b);
        self.fields[i].get_or_init(|| {
            let (sa, sb) = (&self.shapes[a], &self.shapes[b]);
            let px = if self.proxy { Some(ProxySpec { pa: &sa.poles, pb: &sb.poles, eps: 0.01 * sa.diam.max(sb.diam) }) } else { None };
            build_field(&sa.pieces, &sb.pieces, self.res, self.tol, px)
        })
    }

    /// Field and sign for a fixed shape `sf` and a moving shape `sm`:
    /// query with `sign * (pos_m - pos_f)`.
    #[inline(always)]
    pub fn pair(&self, sf: usize, sm: usize) -> (&Field, f64) {
        if sf <= sm {
            (self.field(sf, sm), 1.0)
        } else {
            (self.field(sm, sf), -1.0)
        }
    }

    /// Build every field up front on `threads` threads.
    pub fn prewarm(&self, threads: usize) {
        let ns = self.shapes.len();
        let mut pairs = Vec::with_capacity(ns * (ns + 1) / 2);
        for a in 0..ns {
            for b in a..ns {
                pairs.push((a, b));
            }
        }
        let next = std::sync::atomic::AtomicUsize::new(0);
        std::thread::scope(|sc| {
            for _ in 0..threads.max(1) {
                sc.spawn(|| loop {
                    let k = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    if k >= pairs.len() {
                        break;
                    }
                    let (a, b) = pairs[k];
                    self.field(a, b);
                });
            }
        });
    }

    pub fn field_stats(&self) -> (usize, usize, usize) {
        let mut built = 0;
        let mut bytes = 0;
        let mut grids = 0;
        for f in &self.fields {
            if let Some(f) = f.get() {
                built += 1;
                bytes += f.memory_bytes();
                if matches!(f.kind, crate::nfp::Kind::Grid(_)) {
                    grids += 1;
                }
            }
        }
        (built, grids, bytes)
    }
}
