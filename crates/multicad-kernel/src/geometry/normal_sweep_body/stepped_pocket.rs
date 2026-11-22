//! Cut a box pocket into one cap of an unpierced normal sweep.
//!
//! A box with one axis along the sweep direction projects onto a rectangle
//! of the profile plane. When it reaches past one cap, ends strictly between
//! the caps, and the rectangle widened by a margin meets no profile ring
//! except holes that lie inside the rectangle shrunk by the same margin, the
//! rectangle's boundary lies wholly in material or wholly outside it; one of
//! its corners decides which. Inside, the difference above the floor is the
//! sweep of the profile minus the rectangle and below it the unchanged
//! sweep. Its exact boundary is the profile's own walls and caps, the
//! pocket's four planes and the floor plane, and the enclosed holes' walls
//! cut off at the floor. That boundary is written as an exchange solid and
//! read back by the retained-boundary reader, which validates it, and its
//! integrated volume must equal the analytic one.
use super::opening_clearance::Rectangle;
use super::*;

fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| a[i] - b[i])
}

fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

/// Whether `edge` lies in the closed rectangle; arcs are judged by their
/// whole circle.
fn edge_within(rectangle: &Rectangle, edge: &Edge) -> bool {
    let closed = |p: [f64; 2], reach: f64| {
        let [s, t] = rectangle.local(p);
        s - reach >= 0. && s + reach <= rectangle.length && t.abs() + reach <= rectangle.half_width
    };
    match *edge {
        Edge::Line { a, b } => closed(a, 0.) && closed(b, 0.),
        Edge::Arc { center, a, .. } => {
            closed(center, (a[0] - center[0]).portable_hypot(a[1] - center[1]))
        }
        Edge::Circle {
            center, radius_mm, ..
        } => closed(center, radius_mm),
    }
}

/// The rectangle with every side moved outward by `offset`.
fn offset(rectangle: Rectangle, offset: f64) -> Rectangle {
    Rectangle {
        base: std::array::from_fn(|i| rectangle.base[i] - offset * rectangle.axis[i]),
        axis: rectangle.axis,
        length: rectangle.length + 2. * offset,
        half_width: rectangle.half_width + offset,
    }
}

impl Solid {
    pub(in crate::geometry) fn stepped_pocket_cut(
        &self,
        tool: &Self,
        operation: BooleanOperation,
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        if operation != BooleanOperation::Subtract
            || policy.boolean_fuzzy_mm != 0.
            || self.policy != policy
            || tool.policy != policy
        {
            return None;
        }
        let (body, placement) = self.normal_sweep_body()?;
        let (prism, frame) = tool.rust_prism_frame()?;
        if body.openings.is_some()
            || prism.rounded_end().is_some()
            || prism.inner_loops().next().is_some()
        {
            return None;
        }
        let size = prism.dimensions().map(Rational::as_f64);
        let normal = placement.normal();
        let axes = [frame.u, frame.v, frame.normal()];
        let height = body.sweep.height_mm();
        let magnitude = frame
            .origin
            .into_iter()
            .chain(placement.origin)
            .chain(size)
            .map(f64::abs)
            .fold(height.max(1.), f64::max);
        let rounding = 64. * f64::EPSILON * magnitude;
        let margin = policy.minimum_feature_mm.max(rounding);
        let along = (0..3).find(|&k| norm(cross(axes[k], normal)) <= 64. * f64::EPSILON)?;
        let [i, j] = match along {
            0 => [1, 2],
            1 => [0, 2],
            _ => [0, 1],
        };
        let start = dot(sub3(frame.origin, placement.origin), normal);
        let end = start + dot(axes[along], normal) * size[along];
        let [low, high] = [start.min(end), start.max(end)];
        let levels = if body.reverse {
            [-height, 0.]
        } else {
            [0., height]
        };
        let inner = |level: f64| level > levels[0] + margin && level < levels[1] - margin;
        let (floor, top) = if high >= levels[1] - rounding && inner(low) {
            (low, true)
        } else if low <= levels[0] + rounding && inner(high) {
            (high, false)
        } else {
            return None;
        };
        let project = |offset: [f64; 3]| {
            let point: [f64; 3] = std::array::from_fn(|n| frame.origin[n] + offset[n]);
            let d = sub3(point, placement.origin);
            [dot(d, placement.u), dot(d, placement.v)]
        };
        let side = |k: usize| axes[k].map(|x| x * size[k]);
        let a = project([0.; 3]);
        let b = project(side(i));
        let c = project(std::array::from_fn(|n| side(i)[n] + side(j)[n]));
        let d = project(side(j));
        let length = (b[0] - a[0]).portable_hypot(b[1] - a[1]);
        let width = (d[0] - a[0]).portable_hypot(d[1] - a[1]);
        if !(length > 2. * margin && width > 2. * margin) {
            return None;
        }
        let rectangle = Rectangle {
            base: std::array::from_fn(|n| (a[n] + d[n]) / 2.),
            axis: [(b[0] - a[0]) / length, (b[1] - a[1]) / length],
            length,
            half_width: width / 2.,
        };
        let widened = offset(rectangle, margin);
        let shrunk = offset(rectangle, -margin);
        let rings = body.sweep.rings();
        if rings[0].iter().any(|edge| widened.meets(edge)) {
            return None;
        }
        let mut enclosed = Vec::with_capacity(rings.len() - 1);
        for ring in &rings[1..] {
            if ring.iter().all(|edge| edge_within(&shrunk, edge)) {
                enclosed.push(true);
            } else if ring.iter().all(|edge| !widened.meets(edge)) {
                enclosed.push(false);
            } else {
                return None;
            }
        }
        // No ring crosses the rectangle's boundary, so one corner in the
        // material puts the whole boundary there.
        if !body.sweep.contains([a[0], a[1], height / 2.], 0.) {
            return None;
        }
        let mut pocket: Vec<Edge> = [a, b, c, d]
            .iter()
            .zip([b, c, d, a].iter())
            .map(|(&a, &b)| Edge::Line { a, b })
            .collect();
        if normal_sweep::signed_area(&pocket) > 0. {
            pocket = pocket.iter().rev().map(|edge| edge.reversed()).collect();
        }
        Some(self.write_stepped_pocket(body, placement, &pocket, &enclosed, floor, top, policy))
    }

    #[allow(clippy::too_many_arguments)]
    fn write_stepped_pocket(
        &self,
        body: &NormalSweepBody,
        placement: PlaneFrame,
        pocket: &[Edge],
        enclosed: &[bool],
        floor: f64,
        top: bool,
        policy: GeometryPolicy,
    ) -> Result<Option<Self>, GeometryError> {
        const OPERATION: &str = "stepped pocket cut";
        let fail = |kind, message: String| error(OPERATION, kind, message);
        let mut writer = multicad_rust_kernel::step::Writer::new();
        let text = writer
            .append_stepped_pocket_sweep(
                &body.sweep,
                pocket,
                enclosed,
                floor,
                top,
                placement.origin,
                placement.u,
                placement.v,
                body.reverse,
            )
            .and_then(|()| writer.finish())
            .map_err(|cause| fail(GeometryErrorKind::Unsupported, format!("{cause:?}")))?;
        let result = Self::import_rust_brep(&text, body.feature, policy)?;
        // The pocket removes its outline less the enclosed holes, from the
        // floor to the cap it enters.
        let height = body.sweep.height_mm();
        let levels = if body.reverse {
            [-height, 0.]
        } else {
            [0., height]
        };
        let depth = if top {
            levels[1] - floor
        } else {
            floor - levels[0]
        };
        let rings = body.sweep.rings();
        let holes: f64 = rings[1..]
            .iter()
            .zip(enclosed)
            .filter(|(_, inside)| **inside)
            .map(|(ring, _)| normal_sweep::signed_area(ring).abs())
            .sum();
        let removed = (normal_sweep::signed_area(pocket).abs() - holes) * depth;
        let expected = body.sweep.volume_mm3() - removed;
        let actual = result.volume_mm3()?;
        if !(removed > 0. && (actual - expected).abs() <= 1e-9 * expected.abs().max(1.)) {
            return Err(fail(
                GeometryErrorKind::InvalidResult,
                format!("The pocketed boundary encloses {actual} mm³, not {expected} mm³"),
            ));
        }
        Ok(Some(result))
    }
}
