//! Generated worm-wheel face to a finite B-rep point, by certified branch
//! and bound over the exact generating patch.
use super::*;
use multicad_core::references::{AnalyticGeometry, TopologyKind};
use multicad_precision::PortableMath as _;
use multicad_rust_kernel::worm_gear::wheel_distance::Error;

type Point = [f64; 3];
/// Largest declared gap between the returned distance and its proven lower
/// bound, whatever the model tolerance.
const MAX_GAP_MM: f64 = 1e-5;

fn dot(a: Point, b: Point) -> f64 {
    (0..3).map(|i| a[i] * b[i]).sum()
}
fn norm(a: Point) -> f64 {
    a[0].portable_hypot(a[1]).portable_hypot(a[2])
}
fn sub(a: Point, b: Point) -> Point {
    std::array::from_fn(|i| a[i] - b[i])
}

pub(super) fn measure(
    a: &Solid,
    sa: BRepSelection,
    b: &Solid,
    sb: BRepSelection,
    tolerance: f64,
) -> Option<Result<MinimumDistance, GeometryError>> {
    let (solid, face, probe, index, reverse) = match (sa, sb) {
        (BRepSelection::Face(f), BRepSelection::Vertex(i)) => (a, f, b, i, false),
        (BRepSelection::Vertex(i), BRepSelection::Face(f)) => (b, f, a, i, true),
        _ => return None,
    };
    let (wheel, frame, sector, piece) = solid.rust_wheel_patch(face)?;
    Some((|| {
        let unsupported = |message: &str| error(GeometryErrorKind::Unsupported, message);
        let topology = probe
            .rust_topology()
            .ok_or_else(|| unsupported("Point selection has no exact Rust topology"))??;
        let point = topology
            .elements
            .iter()
            .find(|e| e.key.kind == TopologyKind::Vertex && e.key.ordinal == index)
            .and_then(|e| match e.geometry {
                Some(AnalyticGeometry::Point { position_mm }) => Some(position_mm),
                _ => None,
            })
            .ok_or_else(|| {
                error(
                    GeometryErrorKind::Input,
                    "Selected Rust B-rep vertex is missing",
                )
            })?;
        let unit = |v: Point| {
            let n = norm(v);
            v.map(|x| x / n)
        };
        let u = unit(frame.u);
        let n = unit(frame.normal());
        let v = [
            n[1] * u[2] - n[2] * u[1],
            n[2] * u[0] - n[0] * u[2],
            n[0] * u[1] - n[1] * u[0],
        ];
        let delta = sub(point, frame.origin);
        let reach = norm(delta);
        let magnitude = point
            .into_iter()
            .chain(frame.origin)
            .map(f64::abs)
            .fold(reach.max(1.), f64::max);
        // Chart drift moves both the query and the placed witness.
        let placement =
            (norm(sub(frame.u, u)) + norm(sub(frame.v, v)) + norm(sub(frame.normal(), n))) * reach
                + 64. * f64::EPSILON * magnitude;
        let budget = tolerance.clamp(1e-7, MAX_GAP_MM);
        if !placement.is_finite() || placement > budget / 8. {
            return Err(unsupported(
                "Worm-wheel distance placement cannot be certified within model tolerance",
            ));
        }
        let m = wheel.module_mm();
        let (s, c) = (-wheel.sector_angle(sector)).portable_sin_cos();
        let local = [dot(delta, u), dot(delta, v), dot(delta, n)].map(|x| x / m);
        let query = [
            c * local[0] - s * local[1],
            s * local[0] + c * local[1],
            local[2],
        ];
        let found = wheel
            .patch_distance(piece, query, budget / 2. / m)
            .map_err(|cause| match cause {
                Error::Limit => error(
                    GeometryErrorKind::Limit,
                    "Worm-wheel face distance exceeded its cell budget",
                ),
                Error::Invalid | Error::Unresolved => {
                    unsupported("Worm-wheel face distance bounds could not be certified")
                }
            })?;
        let (s, c) = wheel.sector_angle(sector).portable_sin_cos();
        let [x, y, z] = found.point.map(|x| x * m);
        let on_face =
            crate::rust_prism_topology::position(frame, [c * x - s * y, s * x + c * y, z]);
        let distance_mm = norm(sub(point, on_face));
        // The witness distance is an upper bound; the proven lower bound
        // lies within the declared gap below it.
        if !distance_mm.is_finite()
            || (distance_mm - found.distance * m).abs() > budget / 4.
            || distance_mm - found.lower * m > budget
        {
            return Err(unsupported(
                "Worm-wheel face distance witness exceeds the numerical budget",
            ));
        }
        Ok(MinimumDistance {
            distance_mm,
            witnesses: vec![if reverse {
                DistanceWitness {
                    a_mm: point,
                    b_mm: on_face,
                }
            } else {
                DistanceWitness {
                    a_mm: on_face,
                    b_mm: point,
                }
            }],
            witnesses_complete: false,
            total_solutions: None,
            inner_solution: false,
        })
    })())
}
