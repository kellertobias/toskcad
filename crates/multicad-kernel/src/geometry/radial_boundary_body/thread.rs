//! Original authored cone-spine law and semantic ownership, never inferred.
use super::*;
use crate::HelicalThread;
use multicad_core::references::{TopologyElement, TopologyKind};
use multicad_rust_kernel::{
    helical_thread::Spec,
    radial_tensor::{AuthoredBoundary, BoundaryInput},
    tapered_thread::Thread,
};

#[cfg(test)]
mod tests;

const CAPS: [&str; 2] = ["cap.start", "cap.end"];
const BANDS: [&str; 4] = ["lower", "crest", "upper", "root"];

/// Original legacy Custom cone-spine transport, distinct from axial-flat NPT
/// and rounded BSPT laws. Construction and occurrence frames stay separate.
#[derive(Clone, Debug, PartialEq)]
pub struct RetainedHelicalThreadSource {
    pub profile: HelicalThread,
    pub transport: Thread,
    pub construction_frame: PlaneFrame,
    pub cap_faces: [usize; 2],
    /// Original start order; lower flank, crest, upper flank.
    pub band_faces: Vec<[usize; 3]>,
    pub gap_faces: Vec<usize>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum FaceRole {
    Cap(usize),
    Band { start: usize, band: usize },
}
impl FaceRole {
    fn name(self) -> String {
        match self {
            Self::Cap(end) => CAPS[end].into(),
            Self::Band { start, band } => format!("thread/start/{start}/{}", BANDS[band]),
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum EdgeRole {
    Rail {
        start: usize,
        band: usize,
    },
    Rim {
        end: usize,
        start: usize,
        band: usize,
        segment: Option<usize>,
    },
}
impl EdgeRole {
    fn name(self) -> String {
        match self {
            Self::Rail { start, band } => {
                format!("thread/start/{start}/{}/leading_rail", BANDS[band])
            }
            Self::Rim {
                end,
                start,
                band,
                segment,
            } => format!(
                "{}/rim/start/{start}/{}{}",
                CAPS[end],
                BANDS[band],
                segment.map_or_else(String::new, |i| format!("/segment/{i}"))
            ),
        }
    }
}

pub(super) struct AuthoredThread {
    source: RetainedHelicalThreadSource,
    faces: Vec<FaceRole>,
    edges: Vec<EdgeRole>,
    vertices: Vec<String>,
}
fn invalid() -> GeometryError {
    error(
        "retained thread",
        GeometryErrorKind::InvalidResult,
        "Authored thread ownership does not match its finite boundary",
    )
}

impl Solid {
    /// Direct owned legacy Custom transport. Axial pipe and rounded-root laws
    /// require their own source providers; they are never approximated here.
    pub fn rust_retained_helical_thread_with_interrupt(
        profile: HelicalThread,
        feature: EntityId,
        frame: PlaneFrame,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, GeometryError> {
        policy.validate()?;
        checkpoint(interrupt)?;
        frame.checked(policy)?;
        if feature.is_nil() || policy.boolean_fuzzy_mm != 0. {
            return Err(error(
                "retained thread",
                GeometryErrorKind::Input,
                "Authored retained threads need a source identity and exact geometry policy",
            ));
        }
        let transport = Thread::new_with_interrupt(
            Spec {
                major_radius_mm: profile.major_radius_mm,
                pitch_mm: profile.pitch_mm,
                length_mm: profile.length_mm,
                depth_mm: profile.depth_mm,
                crest_width_mm: profile.crest_width_mm,
                lower_flank_deg: profile.lower_flank_deg,
                upper_flank_deg: profile.upper_flank_deg,
                starts: profile.starts,
                left_hand: profile.left_hand,
                root_radius_mm: None,
            },
            profile.radial_taper,
            policy.linear_tolerance_mm,
            policy.minimum_feature_mm,
            interrupt,
        )
        .map_err(failure)?;
        let authored = transport
            .retained_boundary_with_interrupt(
                [frame.origin, frame.u, frame.v],
                policy.linear_tolerance_mm,
                policy.minimum_feature_mm,
                interrupt,
            )
            .map_err(failure)?;
        let metadata = AuthoredThread::new(profile, transport, frame, &authored)?;
        checkpoint(interrupt)?;
        Self::from_certified_radial_boundary(
            authored.boundary,
            feature,
            policy,
            interrupt,
            Some(Rc::new(metadata)),
        )
    }

    /// This metadata exists only for direct authored construction. Blind STEP
    /// imports retain source geometry without acquiring a guessed thread law.
    pub fn rust_retained_helical_thread_source(
        &self,
    ) -> Option<(&RetainedHelicalThreadSource, PlaneFrame)> {
        let (body, occurrence) = self.radial_boundary_body()?;
        Some((&body.authored_thread.as_ref()?.source, occurrence))
    }
}

impl AuthoredThread {
    fn new(
        profile: HelicalThread,
        transport: Thread,
        construction_frame: PlaneFrame,
        authored: &AuthoredBoundary,
    ) -> Result<Self, GeometryError> {
        let input = authored.boundary.input();
        let starts = profile.starts as usize;
        if authored.band_faces.len() != starts || authored.gap_faces.len() != starts {
            return Err(invalid());
        }
        let mut roles: Vec<Option<FaceRole>> = vec![None; input.faces.len()];
        let mut assign = |index: usize, role: FaceRole| -> Result<(), GeometryError> {
            if roles
                .get_mut(index)
                .ok_or_else(invalid)?
                .replace(role)
                .is_some()
            {
                return Err(invalid());
            }
            Ok(())
        };
        for (end, index) in authored.cap_faces.into_iter().enumerate() {
            assign(index, FaceRole::Cap(end))?;
        }
        for (start, bands) in authored.band_faces.iter().enumerate() {
            for (band, &index) in bands.iter().enumerate() {
                assign(index, FaceRole::Band { start, band })?;
            }
            assign(authored.gap_faces[start], FaceRole::Band { start, band: 3 })?;
        }
        let faces = roles
            .into_iter()
            .collect::<Option<Vec<_>>>()
            .ok_or_else(invalid)?;
        for (index, role) in faces.iter().enumerate() {
            let support = &input.surfaces[input.faces[index].surface];
            let compatible = matches!(
                (role, support),
                (FaceRole::Cap(_), Surface::Plane(_))
                    | (FaceRole::Band { band: 0..=2, .. }, Surface::Tensor(_))
                    | (FaceRole::Band { band: 3, .. }, Surface::Cone { .. })
            );
            if !compatible {
                return Err(invalid());
            }
        }
        let mut edge_faces = vec![BTreeSet::new(); input.edges.len()];
        let mut vertex_faces = vec![BTreeSet::new(); input.vertices.len()];
        let mut vertex_edges = vec![Vec::new(); input.vertices.len()];
        for (face, data) in input.faces.iter().enumerate() {
            for &(edge, _) in &data.outer {
                edge_faces[edge].insert(face);
                for vertex in input.edges[edge].vertices {
                    vertex_faces[vertex].insert(face);
                }
            }
        }
        for (edge, data) in input.edges.iter().enumerate() {
            for vertex in data.vertices {
                vertex_edges[vertex].push(edge);
            }
        }
        let edges = (0..input.edges.len())
            .map(|edge| Self::edge_role(input, edge, &edge_faces, &vertex_faces, &faces, starts))
            .collect::<Result<Vec<_>, _>>()?;
        let vertices = (0..input.vertices.len())
            .map(|vertex| {
                Self::vertex_role(&vertex_edges[vertex], &vertex_faces[vertex], &faces, &edges)
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self {
            source: RetainedHelicalThreadSource {
                profile,
                transport,
                construction_frame,
                cap_faces: authored.cap_faces,
                band_faces: authored.band_faces.clone(),
                gap_faces: authored.gap_faces.clone(),
            },
            faces,
            edges,
            vertices,
        })
    }

    fn edge_role(
        input: &BoundaryInput,
        edge: usize,
        edge_faces: &[BTreeSet<usize>],
        vertex_faces: &[BTreeSet<usize>],
        roles: &[FaceRole],
        starts: usize,
    ) -> Result<EdgeRole, GeometryError> {
        let adjacent = edge_faces[edge]
            .iter()
            .map(|&i| roles[i])
            .collect::<Vec<_>>();
        if adjacent.len() != 2 {
            return Err(invalid());
        }
        if let Some(end) = adjacent.iter().find_map(|r| {
            if let FaceRole::Cap(end) = r {
                Some(*end)
            } else {
                None
            }
        }) {
            let (start, band) = adjacent
                .iter()
                .find_map(|r| {
                    if let FaceRole::Band { start, band } = r {
                        Some((*start, *band))
                    } else {
                        None
                    }
                })
                .ok_or_else(invalid)?;
            let segment = if band == 3 {
                let linked = input.edges[edge]
                    .vertices
                    .into_iter()
                    .flat_map(|v| vertex_faces[v].iter().copied())
                    .map(|i| roles[i])
                    .collect::<Vec<_>>();
                match (
                    linked.contains(&FaceRole::Band { start, band: 2 }),
                    linked.contains(&FaceRole::Band {
                        start: (start + 1) % starts,
                        band: 0,
                    }),
                ) {
                    (true, false) => Some(0),
                    (false, true) => Some(1),
                    _ => return Err(invalid()),
                }
            } else {
                None
            };
            return Ok(EdgeRole::Rim {
                end,
                start,
                band,
                segment,
            });
        }
        let band = adjacent.iter().find_map(|r| {
            if let FaceRole::Band { start, band: 0 } = r {
                Some((*start, 0))
            } else {
                None
            }
        });
        if let Some((start, _)) = band {
            if adjacent.contains(&FaceRole::Band { start, band: 1 }) {
                return Ok(EdgeRole::Rail { start, band: 1 });
            }
            if adjacent.contains(&FaceRole::Band {
                start: (start + starts - 1) % starts,
                band: 3,
            }) {
                return Ok(EdgeRole::Rail { start, band: 0 });
            }
        }
        if let Some(start) = adjacent.iter().find_map(|r| {
            if let FaceRole::Band { start, band: 2 } = r {
                Some(*start)
            } else {
                None
            }
        }) {
            if adjacent.contains(&FaceRole::Band { start, band: 1 }) {
                return Ok(EdgeRole::Rail { start, band: 2 });
            }
            if adjacent.contains(&FaceRole::Band { start, band: 3 }) {
                return Ok(EdgeRole::Rail { start, band: 3 });
            }
        }
        Err(invalid())
    }

    fn vertex_role(
        vertex_edges: &[usize],
        vertex_faces: &BTreeSet<usize>,
        faces: &[FaceRole],
        edges: &[EdgeRole],
    ) -> Result<String, GeometryError> {
        let end = vertex_faces
            .iter()
            .find_map(|&i| {
                if let FaceRole::Cap(end) = faces[i] {
                    Some(end)
                } else {
                    None
                }
            })
            .ok_or_else(invalid)?;
        if let Some((start, band)) = vertex_edges.iter().find_map(|&i| {
            if let EdgeRole::Rail { start, band } = edges[i] {
                Some((start, band))
            } else {
                None
            }
        }) {
            return Ok(format!(
                "{}/corner/start/{start}/{}",
                CAPS[end], BANDS[band]
            ));
        }
        let start = vertex_edges
            .iter()
            .find_map(|&i| {
                if let EdgeRole::Rim { start, band: 3, .. } = edges[i] {
                    Some(start)
                } else {
                    None
                }
            })
            .ok_or_else(invalid)?;
        Ok(format!("{}/root-seam/start/{start}", CAPS[end]))
    }

    pub(super) fn face_origins(&self, feature: EntityId, index: usize) -> Vec<Origin> {
        let mut result = vec![Origin {
            source: feature,
            role: self.faces[index].name(),
        }];
        if matches!(self.faces[index], FaceRole::Band { band: 3, .. }) {
            result.push(Origin {
                source: feature,
                role: "thread.root".into(),
            });
        }
        result
    }
    pub(super) fn apply_origin(&self, feature: EntityId, element: &mut TopologyElement) {
        let role = match element.key.kind {
            TopologyKind::Face => return,
            TopologyKind::Edge => self.edges[element.key.ordinal as usize].name(),
            TopologyKind::Vertex => self.vertices[element.key.ordinal as usize].clone(),
        };
        element.origins = BTreeSet::from([Origin {
            source: feature,
            role,
        }]);
    }
}
