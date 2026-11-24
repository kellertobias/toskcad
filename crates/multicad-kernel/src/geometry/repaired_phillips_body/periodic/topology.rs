//! Finite topology of the retained periodic material, independent of meshing.
//! The four periodic side bands share their original phase nodes. The crest
//! meets the final head wall directly: no zero-width shoulder is introduced.
use super::*;
use multicad_core::references::{
    AnalyticGeometry, AnchoredFrame, Origin, TopologyElement, TopologyKey, TopologyKind,
    TopologySnapshot,
};
use multicad_rust_kernel::{
    phillips_drive::{FaceKey, Inner, WallKind},
    step_read::PhillipsMaterialFaceSource as FaceSource,
};
use std::collections::BTreeSet;

struct Graph {
    feature: EntityId,
    elements: Vec<TopologyElement>,
    edges: u32,
    vertices: u32,
}
impl Graph {
    fn origin(&self, role: impl Into<String>) -> Origin {
        Origin {
            source: self.feature,
            role: role.into(),
        }
    }
    fn vertex(&mut self, role: String, point: Point) -> usize {
        let index = self.elements.len();
        self.elements.push(TopologyElement {
            key: TopologyKey {
                kind: TopologyKind::Vertex,
                ordinal: self.vertices,
            },
            origins: BTreeSet::from([self.origin(role)]),
            boundary: BTreeSet::new(),
            frames: vec![],
            geometry: Some(AnalyticGeometry::Point { position_mm: point }),
        });
        self.vertices += 1;
        index
    }
    fn edge(&mut self, role: String, faces: [usize; 2], vertices: &[usize]) {
        let origin = self.origin(role);
        let mut boundary = BTreeSet::new();
        for face in faces {
            self.elements[face].boundary.insert(origin.clone());
            boundary.extend(self.elements[face].origins.iter().cloned());
        }
        for &vertex in vertices {
            self.elements[vertex].boundary.insert(origin.clone());
        }
        self.elements.push(TopologyElement {
            key: TopologyKey {
                kind: TopologyKind::Edge,
                ordinal: self.edges,
            },
            origins: BTreeSet::from([origin]),
            boundary,
            frames: vec![],
            geometry: None,
        });
        self.edges += 1;
    }
}
pub(super) fn cavity_curve_name(
    topology: &multicad_rust_kernel::phillips_drive::BoundaryTopology<'_>,
    key: multicad_rust_kernel::phillips_drive::CurveKey,
) -> String {
    if let multicad_rust_kernel::phillips_drive::CurveKey::Mouth { edge } = key {
        let (quarter, role) = topology.events().mouth().edges()[edge].source_identity();
        format!("screw.drive/phillips/curve/mouth/{quarter}/{role:?}")
    } else {
        format!("screw.drive/phillips/curve/{key:?}")
    }
}
fn face_name(source: &OwnedPeriodicMaterial, face: FaceSource) -> String {
    match face {
        FaceSource::HeadMeridian { segment } => format!(
            "screw.head/{}",
            source
                .head()
                .view()
                .profile_proposal()
                .segments()
                .nth(segment)
                .expect("retained head face")
                .role
                .source_role()
        ),
        FaceSource::ThreadSide { face } => format!(
            "thread/start/0/{}",
            ["upper", "root", "lower", "crest"][face]
        ),
        FaceSource::ThreadTip => "screw.thread/end".into(),
        FaceSource::JoinAnnulus => "screw.thread/start".into(),
        FaceSource::TopAnnulus => "screw.head/top.land".into(),
        FaceSource::ReversedCavity { face } => format!("screw.drive/phillips/face/{face:?}"),
    }
}
// A retained inner wall is the physical C-C curve swept through its corner.
// Cutter's validated Inner law defines its tilted cylinder (Rounded) or the
// signed pair of planes (ShallowV); envelopes, outer walls and head faces do
// not inherit this construction. Multiple true patches may share ancestry.
fn inner_curve_ancestry(source: &OwnedPeriodicMaterial, face: FaceSource) -> Option<String> {
    let FaceSource::ReversedCavity {
        face:
            FaceKey::Wall {
                kind: WallKind::Inner,
                quarter,
                lower,
            },
    } = face
    else {
        return None;
    };
    let curve = match source.cavity().spec().inner {
        Inner::Rounded { .. } => "inner.round",
        Inner::ShallowV { .. } => {
            if lower {
                "inner.lower"
            } else {
                "inner.upper"
            }
        }
    };
    Some(format!("screw.phillips.corner.{quarter}/{curve}"))
}
impl RepairedPhillipsBody {
    pub(crate) fn periodic_topology(
        &self,
        frame: PlaneFrame,
        owner: EntityId,
        revision: u64,
    ) -> Result<TopologySnapshot, GeometryError> {
        let source = self.source.periodic().ok_or_else(|| {
            unsupported("Retained polynomial source topology has not been admitted")
        })?;
        self.validate_frame(frame)?;
        let placement = self.placement(frame)?;
        let reserve = self.policy.linear_tolerance_mm / 8.;
        let point = |p| placement.bounded_point(p, reserve);
        let mut graph = Graph {
            feature: self.feature,
            elements: vec![],
            edges: 0,
            vertices: 0,
        };
        for index in 0..source.face_count() {
            let face = source
                .face_source(index)
                .expect("immutable complete face inventory");
            let plane = super::topology_planes::plane(source, face, placement, self.policy)?;
            let origin = graph.origin(face_name(source, face));
            graph.elements.push(TopologyElement {
                key: TopologyKey {
                    kind: TopologyKind::Face,
                    ordinal: index as u32,
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
        let head = source.head().view();
        let h = head.segment_count();
        let tip = h + 4;
        let join = h + 5;
        let top = h + 6;
        let caps = source
            .cap_expressions_mm()
            .map(|p| I::point(p[0]) + I::point(p[1]));
        let mut cap_vertices = [[0; 4]; 2];
        let mut outer = [0; 4];
        for cap in 0..2 {
            for band in 0..4 {
                let bounds = source
                    .thread()
                    .bounds(
                        band,
                        [0., 0.],
                        [
                            caps[cap].lo.max(0.),
                            caps[cap].hi.min(source.thread().length_mm()),
                        ],
                        16384,
                        None,
                    )
                    .map_err(source_failure)?;
                let p = bounds.position_mm.map(|p| I::new(p[0], p[1]));
                cap_vertices[cap][band] =
                    graph.vertex(format!("thread/cap/{cap}/vertex/{band}"), point(p)?);
                if cap == 0 {
                    outer[band] = if band == 0 || band == 3 {
                        cap_vertices[cap][band]
                    } else {
                        let ratio = I::point(source.thread().radius_nodes_mm()[0])
                            / I::point(source.thread().radius_nodes_mm()[band]);
                        graph.vertex(
                            format!("screw.head/thread.junction/vertex/{band}"),
                            point([p[0] * ratio, p[1] * ratio, p[2]])?,
                        )
                    };
                }
            }
        }
        for band in 0..4 {
            let next = (band + 1) % 4;
            graph.edge(
                format!("thread/rail/{band}"),
                [h + band, h + (band + 3) % 4],
                &[cap_vertices[0][band], cap_vertices[1][band]],
            );
            graph.edge(
                format!("thread/cap/1/rim/{band}"),
                [h + band, tip],
                &[cap_vertices[1][band], cap_vertices[1][next]],
            );
            graph.edge(
                format!("thread/cap/0/rim/{band}"),
                [h + band, if band == 3 { h - 1 } else { join }],
                &[cap_vertices[0][band], cap_vertices[0][next]],
            );
            if band != 3 {
                graph.edge(
                    format!("screw.head/thread.junction/{band}"),
                    [join, h - 1],
                    &[outer[band], outer[next]],
                );
            }
        }
        // Each immutable repaired meridian endpoint defines one exact common
        // circle. Closed circles have no artificial parameter-seam vertex.
        for segment in 0..h {
            graph.edge(
                format!(
                    "screw.head/{}/start",
                    head.profile_proposal()
                        .segments()
                        .nth(segment)
                        .expect("retained segment")
                        .role
                        .source_role()
                ),
                [if segment == 0 { top } else { segment - 1 }, segment],
                &[],
            );
        }
        let cavity = source
            .cavity()
            .boundary_topology()
            .map_err(|e| unsupported(format!("Actual cavity topology: {e:?}")))?;
        let face_index = |key| {
            if key == FaceKey::Mouth {
                top
            } else {
                (0..source.face_count()).find(|&i| matches!(source.face_source(i), Some(FaceSource::ReversedCavity { face }) if face == key)).expect("same retained cavity inventory")
            }
        };
        let mouth = source.cavity_mouth_expression_mm();
        let mouth = I::point(mouth[0]) + I::point(mouth[1]);
        let mut vertices = Vec::with_capacity(cavity.vertices().len());
        for vertex in cavity.vertices() {
            let mut p = vertex.point_mm_bounds().map(|p| I::new(p[0], p[1]));
            p[2] = p[2] + mouth;
            vertices.push(graph.vertex(
                format!("screw.drive/phillips/vertex/{:?}", vertex.source_identity()),
                point(p)?,
            ));
        }
        let mut uses = vec![Vec::new(); cavity.edges().len()];
        for face in cavity.faces() {
            for edge in face.edge_uses() {
                uses[edge.edge].push(face_index(face.source_identity()));
            }
        }
        for (i, edge) in cavity.edges().iter().enumerate() {
            let adjacent: [usize; 2] = uses[i]
                .as_slice()
                .try_into()
                .map_err(|_| unsupported("Actual cavity edge is not two-sided"))?;
            let v = edge.vertices().map(|v| vertices[v]);
            graph.edge(
                cavity_curve_name(&cavity, edge.source_identity()),
                adjacent,
                &v,
            );
        }
        // Add semantic construction ancestry after the finite graph has used
        // each adjacent face's unchanged primary source identity. This keeps
        // exact edge/face incidence while retaining both physical meanings.
        for index in 0..source.face_count() {
            if let Some(role) = inner_curve_ancestry(
                source,
                source
                    .face_source(index)
                    .expect("same immutable face inventory"),
            ) {
                let origin = graph.origin(role);
                graph.elements[index].origins.insert(origin);
            }
        }
        let snapshot = TopologySnapshot {
            owner,
            revision,
            elements: graph.elements.into(),
        };
        snapshot.validate().map_err(unsupported)?;
        Ok(snapshot)
    }
}
