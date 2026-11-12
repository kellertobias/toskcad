//! Compare complete finite imported carriers with an owned screw source.
use super::*;
use crate::{head_profile::Curve, profiled_screw::Screw};

pub fn profiled_screw_import_matches(
    imported: &ImportSolid,
    source: &Screw,
    frame: [Point; 3],
    tolerance: f64,
) -> bool {
    let Some((other, head_volume_error_mm3)) = imported.profiled_screw_geometry() else {
        return false;
    };
    let reserve =
        16. * (source.profile().rounding_error_mm() + other.profile().rounding_error_mm());
    if !tolerance.is_finite()
        || tolerance <= 0.
        || !head_volume_error_mm3.is_finite()
        || head_volume_error_mm3 <= 0.
        || reserve > tolerance / 4.
    {
        return false;
    }
    let other_frame = [imported.origin_mm, imported.u, imported.v];
    if !thread_nets_match(
        source.thread(),
        other.thread(),
        [None, None],
        [frame, other_frame],
        tolerance,
    ) {
        return false;
    }
    let world = |[origin, u, v]: [Point; 3], r: Point2, z: f64| {
        let n = cross(u, v);
        std::array::from_fn(|i| origin[i] + r[0] * u[i] + r[1] * v[i] + z * n[i])
    };
    let shifts = [source.spec().shank_length_mm, other.spec().shank_length_mm];
    let near = |a: Point, b: Point| length(sub(a, b)) + reserve <= tolerance / 2.;
    let station = |a: Point2, b: Point2| {
        (a[0] - b[0]).abs()
            + length(sub(
                world(frame, [0.; 2], a[1] - shifts[0]),
                world(other_frame, [0.; 2], b[1] - shifts[1]),
            ))
            + reserve
            <= tolerance / 2.
    };
    let a: Vec<_> = source.profile().segments().collect();
    let b: Vec<_> = other.profile().segments().collect();
    if a.len() != b.len() {
        return false;
    }
    let meridian = a.iter().zip(&b).all(|(a, b)| {
        let curve = match (a.curve, b.curve) {
            (Curve::Line, Curve::Line) => true,
            (
                Curve::Arc {
                    center: ca,
                    radius_mm: ra,
                    sweep: sa,
                },
                Curve::Arc {
                    center: cb,
                    radius_mm: rb,
                    sweep: sb,
                },
            ) => ra == rb && station(ca, cb) && ra * (sa - sb).abs() + reserve <= tolerance / 2.,
            _ => false,
        };
        a.role == b.role && station(a.start, b.start) && station(a.end, b.end) && curve
    });
    let points = [source, other].map(|s| s.drive().corners_mm());
    let drive = (0..6).any(|start| {
        [0., 1.].into_iter().all(|floor| {
            points[0].into_iter().enumerate().all(|(i, a)| {
                let b = points[1][(i + start) % 6];
                near(
                    world(
                        frame,
                        a,
                        source.head_top_z() + floor * source.spec().socket_depth_mm,
                    ),
                    world(
                        other_frame,
                        b,
                        other.head_top_z() + floor * other.spec().socket_depth_mm,
                    ),
                )
            })
        })
    });
    meridian && drive
}
