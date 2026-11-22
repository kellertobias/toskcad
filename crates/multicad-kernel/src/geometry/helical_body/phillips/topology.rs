//! Actual head/thread junction and retained finite cutter incidence.
use super::super::{boundary::Boundary, *};
use crate::geometry::internal_thread_treatment_body::interval::Interval as B;
use multicad_core::references::Origin;
use multicad_rust_kernel::phillips_drive::{
    BoundaryTopology, CurveKey, FaceKey, PhillipsScrew, ScrewFace as PhiFace, VertexKey,
};
const PREFIX: &str = "screw.drive/phillips/";
pub(in crate::geometry::helical_body) fn face_name(
    source: &PhillipsScrew,
    role: PhiFace,
) -> String {
    match role {
        PhiFace::ThreadShoulder => "screw.thread/start".into(),
        PhiFace::ThreadTip => "screw.thread/end".into(),
        PhiFace::ThreadWall(index) => format!(
            "thread/start/{}/{}",
            index / 4,
            BAND_NAMES[index as usize % 4]
        ),
        PhiFace::HeadTop => "screw.head/top.land".into(),
        PhiFace::ProfileWall(index) => format!(
            "screw.head/{}",
            source
                .profile()
                .segments()
                .nth(index as usize)
                .expect("retained profile face")
                .role
                .source_role()
        ),
        PhiFace::Cavity(key) => format!("{PREFIX}face/{key:?}"),
    }
}
fn curve_name(topology: &BoundaryTopology<'_>, key: CurveKey) -> String {
    if let CurveKey::Mouth { edge } = key {
        let (quarter, role) = topology.events().mouth().edges()[edge].source_identity();
        format!("{PREFIX}curve/mouth/{quarter}/{role:?}")
    } else {
        format!("{PREFIX}curve/{key:?}")
    }
}
pub(super) fn curve_origin(
    body: &HelicalBody,
    topology: &BoundaryTopology<'_>,
    key: CurveKey,
) -> Origin {
    body.origin(curve_name(topology, key))
}
fn vertex_name(key: VertexKey) -> String {
    format!("{PREFIX}vertex/{key:?}")
}
pub(in crate::geometry::helical_body) fn is_cavity_edge(element: &TopologyElement) -> bool {
    element.key.kind == TopologyKind::Edge
        && element
            .origins
            .iter()
            .any(|origin| origin.role.starts_with(PREFIX))
}
/// Literal stored placement columns, including its raw third column. The
/// midpoint is a display/reference representative of this original source,
/// accepted only with its full directed coordinate error inside policy.
pub(super) fn placed_point(
    frame: PlaneFrame,
    point: [[f64; 2]; 3],
    head_top: f64,
    policy: GeometryPolicy,
) -> Result<Point, GeometryError> {
    placed_point_with_error(frame, point, head_top, policy).map(|(point, _)| point)
}
pub(super) fn placed_point_with_error(
    frame: PlaneFrame,
    point: [[f64; 2]; 3],
    head_top: f64,
    policy: GeometryPolicy,
) -> Result<(Point, f64), GeometryError> {
    let mut point = point.map(|[lo, hi]| B::new(lo, hi));
    point[2] = point[2] + B::point(head_top);
    let axes = [frame.u, frame.v, frame.normal()];
    let world: [B; 3] = std::array::from_fn(|i| {
        (0..3).fold(B::point(frame.origin[i]), |sum, j| {
            sum + B::point(axes[j][i]) * point[j]
        })
    });
    if world.into_iter().any(|v| !v.valid()) {
        return Err(super::failure(
            multicad_rust_kernel::phillips_drive::Error::Unresolved,
        ));
    }
    let midpoint = world.map(|v| v.lo.midpoint(v.hi));
    let error = world.into_iter().zip(midpoint).fold(0_f64, |sum, (v, p)| {
        (sum + (p - v.lo).next_up().max((v.hi - p).next_up())).next_up()
    });
    if error > policy.linear_tolerance_mm / 8. {
        return Err(error_message(
            "Source vertex placement reaches the precision floor",
        ));
    }
    policy.point(midpoint)?;
    Ok((midpoint, error))
}
fn error_message(message: &str) -> GeometryError {
    error(
        "Rust Phillips topology",
        GeometryErrorKind::Unsupported,
        message,
    )
}
fn next_ordinal(elements: &[TopologyElement], kind: TopologyKind) -> u32 {
    elements
        .iter()
        .filter(|e| e.key.kind == kind)
        .map(|e| e.key.ordinal)
        .max()
        .map_or(0, |v| v + 1)
}
fn face_index(source: &PhillipsScrew, key: FaceKey) -> usize {
    let role = if key == FaceKey::Mouth {
        PhiFace::HeadTop
    } else {
        PhiFace::Cavity(key)
    };
    source
        .faces()
        .iter()
        .position(|candidate| *candidate == role)
        .expect("same immutable cutter side owners")
}
/// Extend the original helical graph. Mouth is replaced by the actual head's
/// top land, retaining every original cavity edge and vertex without a cap.
pub(in crate::geometry::helical_body) fn extend(
    body: &HelicalBody,
    source: &PhillipsScrew,
    frame: PlaneFrame,
    policy: GeometryPolicy,
    elements: &mut Vec<TopologyElement>,
) -> Result<(), GeometryError> {
    policy.validate()?;
    frame.checked(policy)?;
    let roles = source.faces();
    let count = body.rails.len();
    let segments: Vec<_> = source.profile().segments().collect();
    let base_faces = elements
        .iter()
        .filter(|e| e.key.kind == TopologyKind::Face)
        .collect::<Vec<_>>();
    if base_faces.len() != 2 + count
        || base_faces
            .iter()
            .enumerate()
            .any(|(i, e)| e.key.ordinal != i as u32)
    {
        return Err(error_message(
            "Original helical face handoff is incomplete or stale",
        ));
    }

    for (ordinal, &role) in roles.iter().enumerate().skip(2 + count) {
        let origin = body.origin(face_name(source, role));
        let plane = super::plane(source, role, frame, policy)?;
        elements.push(TopologyElement {
            key: TopologyKey {
                kind: TopologyKind::Face,
                ordinal: ordinal as u32,
            },
            origins: BTreeSet::from([origin.clone()]),
            boundary: BTreeSet::new(),
            frames: plane
                .into_iter()
                .map(|frame| AnchoredFrame {
                    origin: origin.clone(),
                    frame,
                })
                .collect(),
            geometry: plane.map(|frame| AnalyticGeometry::Plane { frame }),
        });
    }
    // The common helper consumes the actual profile and the original thread
    // rails. At an equal-radius junction crest arcs belong to the final head
    // wall; the remaining finite pockets retain their shoulder ownership.
    {
        let mut add = Boundary::new(body, frame, &roles, elements, false);
        let last = PhiFace::ProfileWall(segments.len() as u32 - 1);
        let shoulder = PhiFace::ThreadShoulder;
        let rate = source.thread().helical_rate();
        let start = source.visible_thread_interval_mm()[0];
        if source.profile().end()[0] > source.thread().spec().major_radius_mm {
            add.arc(
                "screw.head/thread.junction/outer".into(),
                start,
                source.profile().end()[0],
                [0., TAU],
                &[shoulder, last],
            );
        } else {
            for i in 0..count {
                if (body.band(i) - 2) % 4 == 1 {
                    let edge = body.rim_origin(0, i);
                    add.face(shoulder).boundary.remove(&edge);
                    add.face(last).boundary.insert(edge);
                } else {
                    let phase = body.rails[i];
                    add.arc(
                        format!("screw.head/thread.junction/{}", body.rail_name(i)),
                        start,
                        source.thread().spec().major_radius_mm,
                        [(start - phase) / rate, (phase - body.next_phase(i)) / rate],
                        &[shoulder, last],
                    );
                }
            }
        }
        for (i, segment) in segments.iter().enumerate() {
            let before = if i == 0 {
                PhiFace::HeadTop
            } else {
                PhiFace::ProfileWall(i as u32 - 1)
            };
            add.arc(
                format!("screw.head/{}/start", segment.role.source_role()),
                segment.start[1] - source.spec().shank_length_mm,
                segment.start[0],
                [0., TAU],
                &[before, PhiFace::ProfileWall(i as u32)],
            );
        }
    }
    let cutter = source
        .cavity()
        .source()
        .boundary_topology()
        .map_err(super::failure)?;
    let face_slots: Vec<_> = elements
        .iter()
        .enumerate()
        .filter_map(|(i, e)| (e.key.kind == TopologyKind::Face).then_some(i))
        .collect();
    // The validated base plus the complete immutable source-role suffix above
    // supplies exactly the face table consumed by Boundary and this mapping.
    let first_edge = next_ordinal(elements, TopologyKind::Edge);
    let first_vertex = next_ordinal(elements, TopologyKind::Vertex);
    let mut edge_names = Vec::with_capacity(cutter.edges().len());
    let mut face_uses = vec![Vec::<usize>::new(); cutter.edges().len()];
    for face in cutter.faces() {
        let actual = face_slots[face_index(source, face.source_identity())];
        for use_ in face.edge_uses() {
            face_uses[use_.edge].push(actual);
        }
    }
    for (i, edge) in cutter.edges().iter().enumerate() {
        let origin = body.origin(curve_name(&cutter, edge.source_identity()));
        let mut adjacent = BTreeSet::new();
        for &face in &face_uses[i] {
            elements[face].boundary.insert(origin.clone());
            adjacent.extend(elements[face].origins.iter().cloned());
        }
        edge_names.push(origin.clone());
        elements.push(TopologyElement {
            key: TopologyKey {
                kind: TopologyKind::Edge,
                ordinal: first_edge + i as u32,
            },
            origins: BTreeSet::from([origin]),
            boundary: adjacent,
            frames: vec![],
            geometry: None,
        });
    }
    let mut corners = vec![BTreeSet::<Origin>::new(); cutter.vertices().len()];
    for (index, edge) in cutter.edges().iter().enumerate() {
        for vertex in edge.vertices() {
            corners[vertex].insert(edge_names[index].clone());
        }
    }
    for (i, vertex) in cutter.vertices().iter().enumerate() {
        let point = placed_point(frame, vertex.point_mm_bounds(), source.head_top_z(), policy)?;
        elements.push(TopologyElement {
            key: TopologyKey {
                kind: TopologyKind::Vertex,
                ordinal: first_vertex + i as u32,
            },
            origins: BTreeSet::from([body.origin(vertex_name(vertex.source_identity()))]),
            boundary: corners[i].clone(),
            frames: vec![],
            geometry: Some(AnalyticGeometry::Point { position_mm: point }),
        });
    }
    Ok(())
}

#[cfg(test)]
mod tests;
