//! Solution export: JSON placements and an SVG drawing.

use crate::geom::*;
use crate::layout::Layout;
use crate::problem::Problem;
use std::fmt::Write;

/// Final polygon of item `i` in container coordinates.
pub fn placed_polygon(prob: &Problem, lay: &Layout, i: usize) -> Vec<P> {
    let s = &prob.shapes[lay.os[i] as usize];
    s.outline.iter().map(|q| p(q.x + lay.x[i], q.y + lay.y[i])).collect()
}

pub fn to_json(prob: &Problem, lay: &Layout, extra: &[(&str, String)]) -> String {
    let mut o = String::new();
    let _ = write!(o, "{{\n  \"name\": {:?},\n", prob.name);
    for (k, v) in extra {
        let _ = write!(o, "  {:?}: {},\n", k, v);
    }
    let _ = write!(o, "  \"height\": {:?},\n  \"bins\": [", prob.height);
    for (k, b) in lay.bins.iter().enumerate() {
        let _ = write!(o, "{}{{\"x0\": {:?}, \"length\": {:?}}}", if k > 0 { ", " } else { "" }, b.x0, b.len);
    }
    o.push_str("],\n  \"placements\": [\n");
    for i in 0..lay.n() {
        let s = &prob.shapes[lay.os[i] as usize];
        let t = &prob.types[prob.items[i]];
        let b = lay.bin_of(lay.x[i]);
        // final = R(rotation) * original + translation   (translation relative to the bin origin)
        let tx = lay.x[i] - s.off.x - lay.bins[b].x0;
        let ty = lay.y[i] - s.off.y;
        let _ = write!(
            o,
            "    {{\"item_id\": {}, \"bin\": {}, \"rotation\": {:?}, \"translation\": [{:?}, {:?}]}}{}\n",
            t.ext_id,
            b,
            s.angle,
            tx,
            ty,
            if i + 1 < lay.n() { "," } else { "" }
        );
    }
    o.push_str("  ]\n}\n");
    o
}

fn color(k: usize, n: usize) -> String {
    let h = 360.0 * (k as f64 * 0.618_033_988_75).fract();
    let l = if n > 12 && k % 2 == 1 { 62 } else { 72 };
    format!("hsl({:.0},55%,{}%)", h, l)
}

pub fn to_svg(prob: &Problem, lay: &Layout, title: &str) -> String {
    let h = prob.height;
    let gap = if lay.bins.len() > 1 { 0.04 * h } else { 0.0 };
    // bins are stored far apart; draw them next to each other
    let mut draw_x0 = Vec::with_capacity(lay.bins.len());
    let mut cur = 0.0;
    for b in &lay.bins {
        draw_x0.push(cur);
        cur += b.len + gap;
    }
    let total_w = (cur - gap).max(1e-9);
    let px = 1200.0 / total_w.max(h * 0.5);
    let sw = 0.0015 * total_w.max(h);
    let mut o = String::new();
    let _ = write!(
        o,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"{} {} {} {}\" width=\"{:.0}\" height=\"{:.0}\">\n",
        -2.0 * sw,
        -2.0 * sw,
        total_w + 4.0 * sw,
        h + 4.0 * sw,
        (total_w + 4.0 * sw) * px,
        (h + 4.0 * sw) * px
    );
    let _ = write!(o, "<title>{}</title>\n", title);
    for (k, b) in lay.bins.iter().enumerate() {
        let _ = write!(
            o,
            "<rect x=\"{}\" y=\"0\" width=\"{}\" height=\"{}\" fill=\"#f4f4f1\" stroke=\"#333\" stroke-width=\"{}\"/>\n",
            draw_x0[k], b.len, h, sw
        );
    }
    let nt = prob.types.len();
    for i in 0..lay.n() {
        let b = lay.bin_of(lay.x[i]);
        let dx = draw_x0[b] - lay.bins[b].x0;
        let poly = placed_polygon(prob, lay, i);
        let mut pts = String::new();
        for q in &poly {
            let _ = write!(pts, "{:.6},{:.6} ", q.x + dx, h - q.y);
        }
        let _ = write!(
            o,
            "<polygon points=\"{}\" fill=\"{}\" stroke=\"#222\" stroke-width=\"{}\" stroke-linejoin=\"round\"/>\n",
            pts.trim_end(),
            color(prob.items[i], nt),
            sw * 0.6
        );
    }
    o.push_str("</svg>\n");
    o
}
