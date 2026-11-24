//! Finite topology of the actual owned degree-eleven material.
//! Original cyclic side nets retain their cap endpoints and shared rails;
//! the real join annulus meets one complete head-bottom circle.
use super::*;
use multicad_core::references::{
    AnalyticGeometry, AnchoredFrame, Origin, TopologyElement, TopologyKey, TopologyKind,
    TopologySnapshot,
};
use multicad_rust_kernel::{
    phillips_drive::{FaceKey, Inner, WallKind},
    step_read::PhillipsDegree11MaterialFace as FaceSource,
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
fn face_name(source: &Owned, face: FaceSource) -> String {
    match face {
        FaceSource::HeadMeridian { segment, .. } => format!(
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
        FaceSource::ThreadSide { index, .. } => format!("thread/source/{index}"),
        FaceSource::ThreadTip { .. } => "screw.thread/end".into(),
        FaceSource::JoinAnnulus { .. } => "screw.thread/start".into(),
        FaceSource::TopAnnulus { .. } => "screw.head/top.land".into(),
        FaceSource::ReversedCavity { face, .. } => format!("screw.drive/phillips/face/{face:?}"),
    }
}
// A retained inner wall is the physical C-C curve swept through its corner.
// Cutter's validated Inner law defines its tilted cylinder (Rounded) or the
// signed pair of planes (ShallowV); envelopes, outer walls and head faces do
// not inherit this construction. Multiple true patches may share ancestry.
fn inner_curve_ancestry(source: &Owned, face: FaceSource) -> Option<String> {
    let FaceSource::ReversedCavity {
        face:
            FaceKey::Wall {
                kind: WallKind::Inner,
                quarter,
                lower,
            },
        ..
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
    pub(crate) fn degree11_topology(
        &self,
        frame: PlaneFrame,
        owner: EntityId,
        revision: u64,
    ) -> Result<TopologySnapshot, GeometryError> {
        let source = self
            .source
            .degree11()
            .ok_or_else(|| unsupported("Solid has no actual degree-eleven material topology"))?;
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
            let plane = super::planes::plane(source, index, placement, self.policy)?;
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
        // Resolve every face against the actual immutable inventory. No
        // four-band layout or authored section recipe defines these ordinals.
        let find = |predicate: &dyn Fn(FaceSource) -> bool| {
            (0..source.face_count())
                .find(|&i| source.face_source(i).is_some_and(predicate))
                .ok_or_else(|| unsupported("Incomplete actual degree-eleven face inventory"))
        };
        let tip = find(&|f| matches!(f, FaceSource::ThreadTip { .. }))?;
        let join = find(&|f| matches!(f, FaceSource::JoinAnnulus { .. }))?;
        let top = find(&|f| matches!(f, FaceSource::TopAnnulus { .. }))?;
        let head_faces = (0..h)
            .map(|segment| {
                find(&|f| matches!(f, FaceSource::HeadMeridian { segment: i, .. } if i == segment))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let count = source.thread_faces().len();
        if count == 0 || h == 0 {
            return Err(unsupported("Empty actual degree-eleven boundary"));
        }
        let sides = (0..count)
            .map(|index| {
                find(&|f| matches!(f, FaceSource::ThreadSide { index: i, .. } if i == index))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let mut cap_vertices = [Vec::with_capacity(count), Vec::with_capacity(count)];
        for (cap, vertices) in cap_vertices.iter_mut().enumerate() {
            for (index, face) in source.thread_faces().iter().enumerate() {
                let net = face.source();
                let height = net.dimensions()[1];
                let p = net.poles()[if cap == 0 { 0 } else { height - 1 }];
                vertices.push(graph.vertex(
                    format!("thread/cap/{cap}/vertex/{index}"),
                    point(p.map(I::point))?,
                ));
            }
        }
        for index in 0..count {
            let next = (index + 1) % count;
            graph.edge(
                format!("thread/rail/{index}"),
                [sides[index], sides[(index + count - 1) % count]],
                &[cap_vertices[0][index], cap_vertices[1][index]],
            );
            graph.edge(
                format!("thread/cap/1/rim/{index}"),
                [sides[index], tip],
                &[cap_vertices[1][index], cap_vertices[1][next]],
            );
            graph.edge(
                format!("thread/cap/0/rim/{index}"),
                [sides[index], join],
                &[cap_vertices[0][index], cap_vertices[0][next]],
            );
        }
        // Closed physical circles have no artificial parameter seam vertex.
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
                [
                    if segment == 0 {
                        top
                    } else {
                        head_faces[segment - 1]
                    },
                    head_faces[segment],
                ],
                &[],
            );
        }
        graph.edge(
            "screw.head/thread.junction".into(),
            [head_faces[h - 1], join],
            &[],
        );
        let cavity = source
            .cavity()
            .boundary_topology()
            .map_err(|e| unsupported(format!("Actual cavity topology: {e:?}")))?;
        let face_index = |key| {
            if key == FaceKey::Mouth {
                top
            } else {
                (0..source.face_count()).find(|&i| matches!(source.face_source(i), Some(FaceSource::ReversedCavity { face, .. }) if face == key)).expect("same retained cavity inventory")
            }
        };
        let mouth = I::point(source.chart().outer_top_z);
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
