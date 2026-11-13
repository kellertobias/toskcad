//! Exact plane sections of rectangular stock minus finite circular cylinders.
//!
//! A plane meets the stock box in a convex polygon and each cylinder in an
//! ellipse clipped by its two cap planes. The section region is the polygon
//! minus those convex regions. Its boundary is published from one set of
//! junction points: every point shared by two boundary pieces is computed
//! once, from the same root, so the loops close exactly. Ellipse pieces keep
//! their analytic form; the hatch boundary replaces them by tangent polygons
//! that lie in the void, so no stroke enters a bore.
use super::*;
use multicad_rust_kernel::curved_projection::{Curve, Edge as CurvedEdge, Ellipse};
use std::f64::consts::TAU;

type P2 = [f64; 2];
type P3 = [f64; 3];

/// A rectangular box in world coordinates: `origin + Σ size[i] * axes[i]`.
pub(crate) struct Stock {
    pub(crate) origin: P3,
    pub(crate) axes: [P3; 3],
    pub(crate) size: P3,
}

/// A finite circular cylinder: base disk center, unit axis, radius, height.
pub(crate) struct Bore {
    pub(crate) base: P3,
    pub(crate) axis: P3,
    pub(crate) radius: f64,
    pub(crate) height: f64,
}

/// Exact section boundary and its hatch segments in drawing coordinates.
pub(crate) struct Section {
    pub(crate) curves: Vec<Curve>,
    pub(crate) hatch: Vec<[P2; 2]>,
    /// The section polygon of the stock box before bores are removed.
    pub(crate) polygon: Vec<P2>,
}

fn dot3(a: P3, b: P3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn sub3(a: P3, b: P3) -> P3 {
    std::array::from_fn(|i| a[i] - b[i])
}

fn unresolved() -> GeometryError {
    err(
        GeometryErrorKind::Unsupported,
        "Rust oblique section cannot resolve a near-tangent stock or bore boundary",
    )
}

/// One bore's region in the section plane: an ellipse `center + major
/// cos θ u + minor sin θ v`, trimmed to `low <= axial(θ) <= high`, where
/// `axial(θ) = station + reach cos θ` measures along the bore axis.
struct Region {
    ellipse: Ellipse,
    station: f64,
    reach: f64,
    low: f64,
    high: f64,
}

impl Region {
    fn point(&self, angle: f64) -> P2 {
        self.ellipse.point(angle)
    }

    fn axial(&self, angle: f64) -> f64 {
        self.station + self.reach * angle.cos()
    }

    /// Material outside the bore at `point`: outside the ellipse or the slab.
    fn outside(&self, point: P2, tolerance: f64) -> bool {
        let e = &self.ellipse;
        let d = [point[0] - e.center[0], point[1] - e.center[1]];
        let x = (d[0] * e.u[0] + d[1] * e.u[1]) / e.major;
        let y = (d[0] * e.v[0] + d[1] * e.v[1]) / e.minor;
        let axial = self.station + self.reach * x;
        x.hypot(y) > 1. + tolerance || axial < self.low - tolerance || axial > self.high + tolerance
    }
}

/// A bore's void in the section plane: a trimmed ellipse where the plane
/// crosses the axis, or the rectangle between the two generator lines and
/// the cap lines where the plane runs along it.
enum Void {
    Ellipse(Region),
    Strip([P2; 4]),
}

impl Void {
    fn outside(&self, point: P2, tolerance: f64) -> bool {
        match self {
            Self::Ellipse(region) => region.outside(point, tolerance),
            Self::Strip(corners) => (0..4).any(|i| {
                let [a, b] = [corners[i], corners[(i + 1) % 4]];
                (b[0] - a[0]) * (point[1] - a[1]) - (b[1] - a[1]) * (point[0] - a[0])
                    < -tolerance * (b[0] - a[0]).hypot(b[1] - a[1])
            }),
        }
    }
}

/// The first recorded junction within `tolerance` of `point`, or `point`
/// itself as a new junction. Two root computations of one physical
/// junction, such as a cap chord ending on a stock face, then share a
/// single coordinate so the boundary loops close exactly.
fn snap(anchors: &mut Vec<P2>, point: P2, tolerance: f64) -> P2 {
    if let Some(anchor) = anchors
        .iter()
        .find(|a| (a[0] - point[0]).hypot(a[1] - point[1]) <= tolerance)
    {
        return *anchor;
    }
    anchors.push(point);
    point
}

/// Angles where `a cos θ + b sin θ = c`, or nothing. A double root is a
/// tangency that cannot be classified and refuses.
fn roots(a: f64, b: f64, c: f64, tolerance: f64) -> Result<Vec<f64>, GeometryError> {
    let length = a.hypot(b);
    if length <= tolerance {
        return Ok(Vec::new());
    }
    let ratio = c / length;
    if (ratio.abs() - 1.).abs() <= tolerance / length {
        return Err(unresolved());
    }
    if ratio.abs() > 1. {
        return Ok(Vec::new());
    }
    let phase = b.atan2(a);
    let spread = ratio.acos();
    Ok(vec![
        (phase + spread).rem_euclid(TAU),
        (phase - spread).rem_euclid(TAU),
    ])
}

/// The convex section polygon of the stock box in drawing coordinates.
fn polygon(
    stock: &Stock,
    axes: [P3; 3],
    origin: P3,
    depth: f64,
    tolerance: f64,
) -> Result<Vec<P2>, GeometryError> {
    let corner = |bits: usize| -> P3 {
        std::array::from_fn(|i| {
            stock.origin[i]
                + (0..3)
                    .map(|k| {
                        f64::from(u8::from(bits & (1 << k) != 0)) * stock.size[k] * stock.axes[k][i]
                    })
                    .sum::<f64>()
        })
    };
    let height = |p: P3| dot3(sub3(p, origin), axes[2]) - depth;
    // A plane containing a whole stock face has no transversal section.
    for axis in 0..3 {
        for high in [0, 1 << axis] {
            let face: Vec<_> = (0..8).filter(|bits| bits & (1 << axis) == high).collect();
            if face
                .iter()
                .all(|&bits| height(corner(bits)).abs() <= tolerance)
            {
                return Err(err(
                    GeometryErrorKind::Unsupported,
                    "Rust oblique sections do not cut along a stock face",
                ));
            }
        }
    }
    let mut points: Vec<P3> = Vec::new();
    for bits in 0..8 {
        let a = corner(bits);
        let ha = height(a);
        if ha.abs() <= tolerance {
            points.push(a);
        }
        for k in 0..3 {
            if bits & (1 << k) != 0 {
                continue;
            }
            let b = corner(bits | (1 << k));
            let hb = height(b);
            if (ha < -tolerance && hb > tolerance) || (ha > tolerance && hb < -tolerance) {
                let t = ha / (ha - hb);
                points.push(std::array::from_fn(|i| a[i] + (b[i] - a[i]) * t));
            }
        }
    }
    let mut flat: Vec<P2> = points
        .into_iter()
        .map(|p| {
            [
                dot3(sub3(p, origin), axes[0]),
                dot3(sub3(p, origin), axes[1]),
            ]
        })
        .collect();
    flat.dedup_by(|a, b| (a[0] - b[0]).hypot(a[1] - b[1]) <= tolerance);
    let mut unique: Vec<P2> = Vec::new();
    for p in flat {
        if unique
            .iter()
            .all(|q| (p[0] - q[0]).hypot(p[1] - q[1]) > tolerance)
        {
            unique.push(p);
        }
    }
    if unique.len() < 3 {
        return Ok(Vec::new());
    }
    let n = unique.len() as f64;
    let center = [0, 1].map(|i| unique.iter().map(|p| p[i]).sum::<f64>() / n);
    unique.sort_by(|a, b| {
        (a[1] - center[1])
            .atan2(a[0] - center[0])
            .total_cmp(&(b[1] - center[1]).atan2(b[0] - center[0]))
    });
    Ok(unique)
}

/// Clip a counterclockwise convex polygon to `[lo, hi]` along one axis.
fn clip(polygon: Vec<P2>, axis: usize, [lo, hi]: [f64; 2]) -> Vec<P2> {
    let mut result = polygon;
    for (level, sign) in [(lo, 1.), (hi, -1.)] {
        let inside = |p: P2| sign * (p[axis] - level) >= 0.;
        let input = std::mem::take(&mut result);
        for (i, &a) in input.iter().enumerate() {
            let b = input[(i + 1) % input.len()];
            if inside(a) {
                result.push(a);
            }
            if inside(a) != inside(b) {
                let t = (level - a[axis]) / (b[axis] - a[axis]);
                result.push([a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t]);
            }
        }
        if result.is_empty() {
            break;
        }
    }
    result
}

fn region(
    bore: &Bore,
    axes: [P3; 3],
    plane: P3,
    origin: P3,
    tolerance: f64,
) -> Result<Option<Void>, GeometryError> {
    let w = bore.axis;
    let along = dot3(w, axes[2]);
    if along.abs() <= 1e-9 {
        // The plane runs along the axis at distance `offset`: the void is the
        // strip between two generators, trimmed by the cap lines.
        let offset = dot3(sub3(bore.base, plane), axes[2]);
        let gap = bore.radius - offset.abs();
        if gap.abs() <= tolerance {
            return Err(unresolved());
        }
        if gap < 0. {
            return Ok(None);
        }
        let half = (bore.radius * bore.radius - offset * offset).sqrt();
        let foot: P3 = std::array::from_fn(|i| bore.base[i] - offset * axes[2][i]);
        let start = [
            dot3(sub3(foot, origin), axes[0]),
            dot3(sub3(foot, origin), axes[1]),
        ];
        let d = [dot3(w, axes[0]), dot3(w, axes[1])];
        let length = d[0].hypot(d[1]);
        let d = [d[0] / length, d[1] / length];
        let side = [-d[1] * half, d[0] * half];
        let end = [start[0] + d[0] * bore.height, start[1] + d[1] * bore.height];
        // Counterclockwise: start - side, end - side, end + side, start + side.
        return Ok(Some(Void::Strip([
            [start[0] - side[0], start[1] - side[1]],
            [end[0] - side[0], end[1] - side[1]],
            [end[0] + side[0], end[1] + side[1]],
            [start[0] + side[0], start[1] + side[1]],
        ])));
    }
    let station = dot3(sub3(plane, bore.base), axes[2]) / along;
    let center3: P3 = std::array::from_fn(|i| bore.base[i] + station * w[i]);
    let center = [
        dot3(sub3(center3, origin), axes[0]),
        dot3(sub3(center3, origin), axes[1]),
    ];
    let p = [dot3(w, axes[0]), dot3(w, axes[1])];
    let length = p[0].hypot(p[1]);
    let (u, v) = if length <= tolerance {
        ([1., 0.], [0., 1.])
    } else {
        let u = [p[0] / length, p[1] / length];
        (u, [-u[1], u[0]])
    };
    let major = bore.radius / along.abs();
    Ok(Some(Void::Ellipse(Region {
        ellipse: Ellipse {
            center,
            u,
            v,
            major,
            minor: bore.radius,
            trim: [0., TAU],
        },
        station,
        reach: length * major,
        low: 0.,
        high: bore.height,
    })))
}

/// Section the stock minus `bores` in the plane at `depth` along the view
/// normal. `Ok(None)` means the plane misses all material.
pub(crate) fn section(
    stock: &Stock,
    bores: &[Bore],
    options: DrawingProjection,
    depth: f64,
    tolerance: f64,
) -> Result<Option<Section>, GeometryError> {
    let axes = planar::view_axes(options.frame);
    let origin = options.frame.origin;
    let plane: P3 = std::array::from_fn(|i| origin[i] + axes[2][i] * depth);
    let mut outline = polygon(stock, axes, origin, depth, tolerance)?;
    for (axis, window) in [options.section_u_mm, options.section_v_mm]
        .into_iter()
        .enumerate()
    {
        if let Some(window) = window
            && !outline.is_empty()
        {
            outline = clip(outline, axis, window);
        }
    }
    if outline.len() < 3 {
        return Ok(None);
    }
    let voids: Vec<Void> = bores
        .iter()
        .map(|bore| region(bore, axes, plane, origin, tolerance))
        .collect::<Result<Vec<_>, _>>()?
        .into_iter()
        .flatten()
        .collect();
    let count = outline.len();
    // Inward half-planes of the polygon: `normal · point >= offset`.
    let sides: Vec<(P2, f64, P2, P2)> = (0..count)
        .map(|i| {
            let [a, b] = [outline[i], outline[(i + 1) % count]];
            let d = [b[0] - a[0], b[1] - a[1]];
            let length = d[0].hypot(d[1]);
            let normal = [-d[1] / length, d[0] / length];
            (normal, normal[0] * a[0] + normal[1] * a[1], a, b)
        })
        .collect();
    let inside_polygon = |point: P2| {
        sides.iter().all(|(normal, offset, ..)| {
            normal[0] * point[0] + normal[1] * point[1] >= offset - tolerance
        })
    };
    // A cap chord on the polygon boundary separates the bore from empty
    // space, not from material, so only chords strictly inside count.
    let strictly_inside = |point: P2| {
        sides.iter().all(|(normal, offset, ..)| {
            normal[0] * point[0] + normal[1] * point[1] > offset + tolerance
        })
    };
    let outside_bores = |point: P2, skip: Option<usize>| {
        voids
            .iter()
            .enumerate()
            .all(|(i, void)| Some(i) == skip || void.outside(point, tolerance))
    };

    let mut curves = Vec::new();
    let mut hatch = Vec::new();
    let mut anchors: Vec<P2> = Vec::new();
    // Junctions on each polygon side, as (parameter along the side, point).
    let mut junctions: Vec<Vec<(f64, P2)>> = vec![Vec::new(); count];
    let along_side = |side: usize, point: P2| {
        let (_, _, a, b) = sides[side];
        let d = [b[0] - a[0], b[1] - a[1]];
        ((point[0] - a[0]) * d[0] + (point[1] - a[1]) * d[1]) / (d[0] * d[0] + d[1] * d[1])
    };
    for (index, void) in voids.iter().enumerate() {
        let region = match void {
            Void::Ellipse(region) => region,
            Void::Strip(corners) => {
                // Strip sides are exact lines; keep the parts strictly inside
                // the polygon and outside every other bore.
                for i in 0..4 {
                    let [a, b] = [corners[i], corners[(i + 1) % 4]];
                    let [a, b] = [a, b].map(|p| snap(&mut anchors, p, tolerance));
                    let mut stops = vec![(0., a), (1., b)];
                    for (side, &(normal, offset, ..)) in sides.iter().enumerate() {
                        let [pa, pb] = [a, b].map(|p| normal[0] * p[0] + normal[1] * p[1] - offset);
                        if (pa < -tolerance && pb > tolerance)
                            || (pa > tolerance && pb < -tolerance)
                        {
                            let t = pa / (pa - pb);
                            let point = snap(
                                &mut anchors,
                                [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t],
                                tolerance,
                            );
                            let s = along_side(side, point);
                            if (-tolerance..=1. + tolerance).contains(&s) {
                                junctions[side].push((s, point));
                            }
                            stops.push((t, point));
                        }
                    }
                    stops.sort_by(|x, y| x.0.total_cmp(&y.0));
                    for pair in stops.windows(2) {
                        let [p, q] = [pair[0].1, pair[1].1];
                        let middle = [(p[0] + q[0]) / 2., (p[1] + q[1]) / 2.];
                        if (q[0] - p[0]).hypot(q[1] - p[1]) > tolerance
                            && strictly_inside(middle)
                            && outside_bores(middle, Some(index))
                        {
                            curves.push(Curve::Line([p, q]));
                            hatch.push([p, q]);
                        }
                    }
                }
                continue;
            }
        };
        let e = &region.ellipse;
        // Arc breakpoints: the slab limits and every polygon side.
        let mut breaks = Vec::new();
        for level in [region.low, region.high] {
            breaks.extend(roots(region.reach, 0., level - region.station, tolerance)?);
        }
        for (side, &(normal, offset, ..)) in sides.iter().enumerate() {
            let center = normal[0] * e.center[0] + normal[1] * e.center[1];
            let a = e.major * (normal[0] * e.u[0] + normal[1] * e.u[1]);
            let b = e.minor * (normal[0] * e.v[0] + normal[1] * e.v[1]);
            for angle in roots(a, b, offset - center, tolerance)? {
                let point = snap(&mut anchors, region.point(angle), tolerance);
                let t = along_side(side, point);
                let axial = region.axial(angle);
                if (-tolerance..=1. + tolerance).contains(&t)
                    && axial >= region.low - tolerance
                    && axial <= region.high + tolerance
                {
                    junctions[side].push((t, point));
                }
                breaks.push(angle);
            }
        }
        breaks.sort_by(f64::total_cmp);
        breaks.dedup_by(|a, b| (*a - *b).abs() <= 1e-12);
        let spans: Vec<[f64; 2]> = if breaks.is_empty() {
            vec![[0., TAU]]
        } else {
            (0..breaks.len())
                .map(|i| {
                    let start = breaks[i];
                    let end = if i + 1 < breaks.len() {
                        breaks[i + 1]
                    } else {
                        breaks[0] + TAU
                    };
                    [start, end]
                })
                .collect()
        };
        for [start, end] in spans {
            let middle = region.point((start + end) / 2.);
            let axial = region.axial((start + end) / 2.);
            if end - start <= 1e-12
                || axial < region.low
                || axial > region.high
                || !inside_polygon(middle)
                || !outside_bores(middle, Some(index))
            {
                continue;
            }
            // Published trims stay within one turn starting at zero.
            for trim in [[start, end.min(TAU)], [0., end - TAU]] {
                if trim[1] - trim[0] > 1e-12 {
                    curves.push(Curve::Ellipse(Ellipse { trim, ..*e }));
                }
            }
            // The void lies inside the ellipse: chords plus tangent tips
            // outside it keep every hatch stroke in material.
            let step_limit = 4. * (tolerance / (4. * (e.major + tolerance))).sqrt().asin();
            let pieces = ((end - start) / step_limit.clamp(1e-3, 0.5)).ceil().max(1.) as usize;
            if pieces > 100_000 {
                return Err(err(
                    GeometryErrorKind::Limit,
                    "Rust oblique section exceeded its hatch budget",
                ));
            }
            let step = (end - start) / pieces as f64;
            let closed = end - start >= TAU - 1e-12;
            let first = if closed {
                region.point(start)
            } else {
                snap(&mut anchors, region.point(start), tolerance)
            };
            let mut last = first;
            for i in 0..pieces {
                let middle = start + step * (i as f64 + 0.5);
                let scale = 1. / (step / 2.).cos();
                let (s, c) = middle.sin_cos();
                let tip = std::array::from_fn(|k| {
                    e.center[k] + scale * (e.major * c * e.u[k] + e.minor * s * e.v[k])
                });
                // A closed ellipse ends exactly where it starts.
                let next = if i + 1 == pieces && closed {
                    first
                } else if i + 1 == pieces {
                    snap(&mut anchors, region.point(end), tolerance)
                } else {
                    region.point(start + step * (i + 1) as f64)
                };
                hatch.push([last, tip]);
                hatch.push([tip, next]);
                last = next;
            }
        }
        // Cap chords run along `v` inside the ellipse and the polygon.
        for level in [region.low, region.high] {
            let roots = roots(region.reach, 0., level - region.station, tolerance)?;
            let [Some(&first), Some(&second)] = [roots.first(), roots.get(1)] else {
                continue;
            };
            let [a, b] =
                [first, second].map(|angle| snap(&mut anchors, region.point(angle), tolerance));
            let mut stops = vec![(0., a), (1., b)];
            for (side, &(normal, offset, ..)) in sides.iter().enumerate() {
                let [pa, pb] = [a, b].map(|p| normal[0] * p[0] + normal[1] * p[1] - offset);
                if (pa < -tolerance && pb > tolerance) || (pa > tolerance && pb < -tolerance) {
                    let t = pa / (pa - pb);
                    let point = snap(
                        &mut anchors,
                        [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t],
                        tolerance,
                    );
                    let s = along_side(side, point);
                    if (-tolerance..=1. + tolerance).contains(&s) {
                        junctions[side].push((s, point));
                    }
                    stops.push((t, point));
                }
            }
            stops.sort_by(|x, y| x.0.total_cmp(&y.0));
            for pair in stops.windows(2) {
                let [p, q] = [pair[0].1, pair[1].1];
                let middle = [(p[0] + q[0]) / 2., (p[1] + q[1]) / 2.];
                if (q[0] - p[0]).hypot(q[1] - p[1]) > tolerance
                    && strictly_inside(middle)
                    && outside_bores(middle, Some(index))
                {
                    curves.push(Curve::Line([p, q]));
                    hatch.push([p, q]);
                }
            }
        }
    }
    for (side, mut stops) in junctions.into_iter().enumerate() {
        let (_, _, a, b) = sides[side];
        stops.push((0., snap(&mut anchors, a, tolerance)));
        stops.push((1., snap(&mut anchors, b, tolerance)));
        stops.sort_by(|x, y| x.0.total_cmp(&y.0));
        for pair in stops.windows(2) {
            let [p, q] = [pair[0].1, pair[1].1];
            let middle = [(p[0] + q[0]) / 2., (p[1] + q[1]) / 2.];
            if (q[0] - p[0]).hypot(q[1] - p[1]) > tolerance && outside_bores(middle, None) {
                curves.push(Curve::Line([p, q]));
                hatch.push([p, q]);
            }
        }
    }
    if curves.is_empty() {
        return Ok(None);
    }
    Ok(Some(Section {
        curves,
        hatch,
        polygon: outline,
    }))
}

/// Publish an exact section as Section vectors with its hatch.
pub(crate) fn publish(
    section: Section,
    hatch: Option<DrawingHatch>,
    options: DrawingProjection,
    source: usize,
    tolerance: f64,
) -> Result<Vec<ProjectedVector>, GeometryError> {
    let edges = section
        .curves
        .into_iter()
        .map(|curve| CurvedEdge {
            hidden: false,
            curve,
        })
        .collect();
    let mut vectors = curved::publish(edges, options, source)?;
    for vector in &mut vectors {
        vector.role = VectorRole::Section;
    }
    if let Some(hatch) = hatch.filter(|h| h.pattern != HatchPattern::None) {
        let mut curves = (options.maximum_curves as usize).saturating_sub(vectors.len());
        let mut points = (options.maximum_points as usize)
            .saturating_sub(vectors.iter().map(|v| v.preview_mm.len()).sum());
        vectors.extend(publish_hatch(
            source,
            hatch,
            &section.hatch,
            tolerance,
            &mut curves,
            &mut points,
            &mut section_hatch::Budget::default(),
        )?);
    }
    Ok(vectors)
}

/// Sections of the oblique-bore solid: a stock box minus one finite oblique
/// cylinder. Section-only cuts take any transversal plane. Retained and
/// finite-depth cuts are admitted when the view runs along the bore axis and
/// every retained stock corner projects inside the section polygon: then all
/// retained stock lies behind the cut face except the bore, which is seen
/// along its axis, so every visible edge behind the plane coincides with the
/// section boundary and is published there.
pub(crate) fn project_oblique_bore(
    input: &DrawingSolid<'_>,
    options: DrawingProjection,
    depth: f64,
) -> Option<Result<Vec<ProjectedVector>, GeometryError>> {
    let (body, frame) = input.solid.oblique_bore_body()?;
    let tolerance = input.solid.geometry_policy().linear_tolerance_mm;
    let geometry = body.drawing_box_and_tool(frame);
    let (origin, axes, size) = (geometry.origin, geometry.axes, geometry.size);
    let (base, axis, radius, height) = (
        geometry.base,
        geometry.axis,
        geometry.radius,
        geometry.height,
    );
    let stock = Stock { origin, axes, size };
    let bore = Bore {
        base,
        axis,
        radius,
        height,
    };
    Some((|| {
        if options.detail.is_some() {
            return Err(err(
                GeometryErrorKind::Unsupported,
                "Rust oblique bore sections do not support a detail crop",
            ));
        }
        let Some(section) = section(
            &stock,
            std::slice::from_ref(&bore),
            options,
            depth,
            tolerance,
        )?
        else {
            return Ok(Vec::new());
        };
        let view = planar::view_axes(options.frame);
        let retained_boundary = if options.section_only {
            false
        } else {
            let along = dot3(axis, view[2]);
            if (1. - along.abs()) > 1e-12 {
                return Err(err(
                    GeometryErrorKind::Unsupported,
                    "Rust retained oblique bore sections require a view along the bore axis",
                ));
            }
            // Corners behind the plane, within any finite depth, must project
            // inside the convex section polygon.
            let count = section.polygon.len();
            let far = options.section_depth_mm.map(|range| depth - range);
            for bits in 0..8 {
                let corner: P3 = std::array::from_fn(|i| {
                    origin[i]
                        + (0..3)
                            .map(|k| {
                                f64::from(u8::from(bits & (1 << k) != 0)) * size[k] * axes[k][i]
                            })
                            .sum::<f64>()
                });
                let offset = sub3(corner, options.frame.origin);
                let level = dot3(offset, view[2]);
                if level > depth + tolerance || far.is_some_and(|far| level < far - tolerance) {
                    continue;
                }
                let point = [dot3(offset, view[0]), dot3(offset, view[1])];
                let inside = (0..count).all(|i| {
                    let [a, b] = [section.polygon[i], section.polygon[(i + 1) % count]];
                    (b[0] - a[0]) * (point[1] - a[1]) - (b[1] - a[1]) * (point[0] - a[0])
                        >= -tolerance
                });
                if !inside {
                    return Err(err(
                        GeometryErrorKind::Unsupported,
                        "Rust retained oblique bore sections require the retained stock behind the cut face",
                    ));
                }
            }
            true
        };
        let curves = section.curves.clone();
        let mut vectors = publish(section, input.hatch, options, 0, tolerance)?;
        if retained_boundary {
            let mut remaining = options;
            remaining.maximum_curves = remaining
                .maximum_curves
                .saturating_sub(vectors.len() as u32);
            let edges = curves
                .into_iter()
                .map(|curve| CurvedEdge {
                    hidden: false,
                    curve,
                })
                .collect();
            vectors.extend(curved::publish(edges, remaining, 0)?);
        }
        Ok(vectors)
    })())
}
