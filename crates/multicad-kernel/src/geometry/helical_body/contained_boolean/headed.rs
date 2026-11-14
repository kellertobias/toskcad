//! Whole-box membership against the retained headed-screw/nut components.
//! Hex walls use the same source-owned finite vertices as their boundary.
use super::*;
use BoxClassification::{Inside, Outside, UnresolvedBoundary};
use multicad_rust_kernel::headed_screw::{ScrewHead, hex_vertices};

fn intersection(a: BoxClassification, b: BoxClassification) -> BoxClassification {
    match (a, b) {
        (Outside, _) | (_, Outside) => Outside,
        (Inside, Inside) => Inside,
        _ => UnresolvedBoundary,
    }
}
fn subtract(a: BoxClassification, b: BoxClassification) -> BoxClassification {
    match (a, b) {
        (Outside, _) | (_, Inside) => Outside,
        (Inside, Outside) => Inside,
        _ => UnresolvedBoundary,
    }
}
fn slab(z: B, low: f64, high: f64) -> BoxClassification {
    if !z.valid() {
        UnresolvedBoundary
    } else if z.hi < low || z.lo > high || low == high {
        Outside
    } else if z.lo > low && z.hi < high {
        Inside
    } else {
        UnresolvedBoundary
    }
}
fn circle([x, y]: [B; 2], radius: f64) -> BoxClassification {
    let rho = x.square() + y.square();
    let limit = B::point(radius).square();
    if !rho.valid() {
        UnresolvedBoundary
    } else if rho.lo > limit.hi {
        Outside
    } else if rho.hi < limit.lo {
        Inside
    } else {
        UnresolvedBoundary
    }
}
fn hexagon([x, y]: [B; 2], width: f64) -> BoxClassification {
    let vertices = hex_vertices(width);
    let mut inside = true;
    for (a, b) in vertices.iter().zip(vertices.iter().cycle().skip(1)).take(6) {
        let [ax, ay] = a.map(B::point);
        let [bx, by] = b.map(B::point);
        let cross = (bx - ax) * (y - ay) - (by - ay) * (x - ax);
        if !cross.valid() {
            return UnresolvedBoundary;
        }
        if cross.hi < 0. {
            return Outside;
        }
        inside &= cross.lo > 0.;
    }
    if inside { Inside } else { UnresolvedBoundary }
}
pub(super) fn screw_relation(source: &Screw, bounds: [[f64; 2]; 3]) -> BoxClassification {
    let [x, y, z] = bounds.map(|[lo, hi]| B::new(lo, hi));
    if [x, y, z].iter().any(|v| !v.valid()) {
        return UnresolvedBoundary;
    }
    let spec = source.spec();
    let head_outer = match spec.head {
        ScrewHead::SocketCap { radius_mm, .. } => circle([x, y], radius_mm),
        ScrewHead::ExternalHex {
            across_flats_mm, ..
        } => hexagon([x, y], across_flats_mm),
    };
    let mut head = intersection(
        head_outer,
        slab(z, source.head_top_z(), source.head_seat_z()),
    );
    if let Some((width, depth)) = spec.head.socket() {
        let cavity = intersection(
            hexagon([x, y], width),
            slab(z, source.head_top_z(), source.head_top_z() + depth),
        );
        head = subtract(head, cavity);
    }
    let shank = intersection(
        circle([x, y], spec.shank_radius_mm),
        slab(z, source.head_seat_z(), 0.),
    );
    let thread = source
        .thread()
        .classify_box(bounds)
        .unwrap_or(UnresolvedBoundary);
    let regions = [head, shank, thread];
    if regions.contains(&Inside) {
        Inside
    } else if regions.iter().all(|&region| region == Outside) {
        Outside
    } else {
        UnresolvedBoundary
    }
}
pub(super) fn nut_relation(source: &Nut, bounds: [[f64; 2]; 3]) -> BoxClassification {
    let [x, y, z] = bounds.map(|[lo, hi]| B::new(lo, hi));
    if [x, y, z].iter().any(|v| !v.valid()) {
        return UnresolvedBoundary;
    }
    let stock = intersection(
        hexagon([x, y], source.across_flats_mm()),
        slab(z, 0., source.thread().spec().length_mm),
    );
    subtract(
        stock,
        source
            .thread()
            .classify_box(bounds)
            .unwrap_or(UnresolvedBoundary),
    )
}

#[cfg(test)]
mod tests;
