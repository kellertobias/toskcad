//! Exact drawing members of a stepped shaft with one perpendicular
//! cross-bore. The shaft draws as its coaxial bands; the bored band carries
//! the exact tool void for occlusion. The bored-wall rims are published as
//! chords within the drawing chord tolerance of the exact quartic, the
//! bore's inner outlines as their exact generator segments and blind
//! bottoms as their exact circles.
use super::*;
use multicad_rust_kernel::{
    cross_bore::End,
    curved_projection::Ellipse,
    display::Tessellation,
    planar_projection::{LineCandidate, SpatialEllipseCandidate},
    round_occlusion::CrossVoid,
};
use std::f64::consts::TAU;

type Point = [f64; 3];

pub(super) struct Bored {
    pub rounds: Vec<RoundOccluder>,
    pub lines: Vec<LineCandidate>,
    pub bottoms: Vec<SpatialEllipseCandidate>,
}

/// A spatial circle `center + ρ cos t a + ρ sin t b` as an ellipse in its
/// principal paper axes with harmonic depth.
fn spatial_circle(
    source: usize,
    center: Point,
    [a, b]: [Point; 2],
    radius: f64,
) -> SpatialEllipseCandidate {
    let [a, b] = [a.map(|c| c * radius), b.map(|c| c * radius)];
    let [pa, pb] = [[a[0], a[1]], [b[0], b[1]]];
    let turn = 0.5
        * (2. * (pa[0] * pb[0] + pa[1] * pb[1]))
            .atan2(pa[0] * pa[0] + pa[1] * pa[1] - pb[0] * pb[0] - pb[1] * pb[1]);
    let (sin, cos) = turn.sin_cos();
    let major: [f64; 2] = std::array::from_fn(|i| cos * pa[i] + sin * pb[i]);
    let minor: [f64; 2] = std::array::from_fn(|i| -sin * pa[i] + cos * pb[i]);
    let [lengths, depths] = [
        [major[0].hypot(major[1]), minor[0].hypot(minor[1])],
        [cos * a[2] + sin * b[2], -sin * a[2] + cos * b[2]],
    ];
    let u = major.map(|c| c / lengths[0]);
    let v = if lengths[1] > 1e-12 * lengths[0] {
        minor.map(|c| c / lengths[1])
    } else {
        [-u[1], u[0]]
    };
    SpatialEllipseCandidate {
        source,
        ellipse: Ellipse {
            center: [center[0], center[1]],
            u,
            v,
            major: lengths[0],
            minor: lengths[1],
            trim: [0., TAU],
        },
        depth: center[2],
        depth_cosine: depths[0],
        depth_sine: depths[1],
    }
}

pub(super) fn scene(
    solid: &Solid,
    owner: usize,
    source: usize,
    options: DrawingProjection,
) -> Result<Option<Bored>, GeometryError> {
    let Some((body, placement)) = solid.cross_bore_body() else {
        return Ok(None);
    };
    let geometry = body.geometry();
    let unsupported = |message| err(GeometryErrorKind::Unsupported, message);
    let [ring] = geometry.shaft().rings() else {
        return Err(unsupported(
            "Rust cross-bore drawings need one stepped meridian",
        ));
    };
    let bands = stepped::axial_bands(ring, options.tangent_edges)
        .ok_or_else(|| unsupported("Rust cross-bore drawings need a stepped shaft"))?;
    let screen = planar::view_axes(options.frame);
    let vector = |p: Point| screen.map(|axis| dot(axis, p));
    let point = |p: Point| vector(std::array::from_fn(|i| p[i] - options.frame.origin[i]));
    let world = |local: Point| crate::rust_prism_topology::position(placement, local);
    let direction = |local: Point| {
        let axes = [placement.u, placement.v, placement.normal()];
        vector(std::array::from_fn(|i| {
            (0..3).map(|k| axes[k][i] * local[k]).sum()
        }))
    };
    let (_, wall, span) = geometry.wall();
    let ([travel, across], offset, level, stations) = geometry.bore_chart();
    let ends = geometry.ends();
    let axes = [
        direction(travel),
        direction(across),
        direction([0., 0., 1.]),
    ];
    let void = CrossVoid::new(
        point(world([0.; 3])),
        axes,
        wall,
        offset,
        level,
        geometry.tool().2,
        stations,
        ends.map(|end| end == End::Through),
    )
    .map_err(planar::projection_error)?;
    let axis = vector(unit(placement.normal()));
    let mut rounds = Vec::with_capacity(bands.len());
    for [radius, low, high] in bands {
        let mut round =
            RoundOccluder::new(point(world([0., 0., low])), axis, radius, high - low, None)
                .map_err(planar::projection_error)?
                .with_contour_owner(owner);
        if radius == wall && low <= span[0] && span[1] <= high {
            round = round.with_cross_void(void);
        }
        rounds.push(round);
    }
    // Bored-wall rims as chords of the exact curve.
    let chords = Tessellation {
        deflection: options.chord_tolerance_mm,
        angle: 0.25,
        triangle_limit: options.maximum_curves as usize,
    };
    let mut lines = Vec::new();
    for (side, end) in ends.into_iter().enumerate() {
        if end != End::Through {
            continue;
        }
        let samples = geometry.sample_rim(side, chords).map_err(|_| {
            err(
                GeometryErrorKind::Limit,
                "Rust cross-bore rim exceeded the drawing curve budget",
            )
        })?;
        let samples: Vec<Point> = samples.into_iter().map(|p| point(world(p))).collect();
        lines.extend(samples.windows(2).map(|pair| LineCandidate {
            source,
            points: [pair[0], pair[1]],
        }));
    }
    // Inner outlines: generators where the line of sight is tangent to the
    // tool wall, between the retained ends of the bore.
    let toward = screen[2];
    let local_view = [placement.u, placement.v, placement.normal()].map(|axis| dot(axis, toward));
    let view = [travel, across, [0., 0., 1.]].map(|axis| dot(axis, local_view));
    if view[1].hypot(view[2]) > 1e-9 {
        let angle = (-view[1]).atan2(view[2]);
        for angle in [angle, angle + std::f64::consts::PI] {
            let [from, to] = [0, 1].map(|side| point(world(geometry.rim_point(side, angle))));
            lines.push(LineCandidate {
                source,
                points: [from, to],
            });
        }
    }
    let mut bottoms = Vec::new();
    for side in 0..2 {
        if let Some((center, _, radius)) = geometry.cap(side) {
            bottoms.push(spatial_circle(
                source,
                point(world(center)),
                [direction(across), direction([0., 0., 1.])],
                radius,
            ));
        }
    }
    Ok(Some(Bored {
        rounds,
        lines,
        bottoms,
    }))
}
