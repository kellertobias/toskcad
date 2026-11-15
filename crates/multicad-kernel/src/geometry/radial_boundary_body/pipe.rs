//! Direct authored axial pipe laws, distinct from retained/imported supports.
use super::*;
use multicad_core::references::{TopologyElement, TopologyKind};
use multicad_rust_kernel::axial_pipe_thread::authored::{
    AuthoredSection, AuthoredSpec, AuthoredThread,
};
#[cfg(test)]
mod tests;
/// Construction provenance exists only for an explicit authored factory.
/// Independent STEP imports never acquire a guessed pipe recipe or phase.
#[derive(Clone, Debug, PartialEq)]
pub struct RetainedAxialPipeThreadSource {
    pub transport: AuthoredThread,
    pub construction_frame: PlaneFrame,
    pub cap_faces: [usize; 2],
    pub band_faces: Vec<usize>,
    pub root_faces: Vec<usize>,
    /// Authored law to retained net Euclidean bound in the local source chart.
    pub local_source_error_mm: f64,
    /// Complete retained-body radial bound in the identity construction chart.
    /// Present only when an unchanged stock limit was independently proved inactive.
    pub retained_radius_upper_bound_mm: Option<f64>,
}
pub(super) struct AuthoredPipe {
    source: RetainedAxialPipeThreadSource,
    faces: Vec<String>,
    edges: Vec<String>,
    vertices: Vec<String>,
}
impl Solid {
    pub fn rust_authored_axial_pipe_thread_with_interrupt(
        spec: AuthoredSpec,
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
                "authored pipe thread",
                GeometryErrorKind::Input,
                "Authored pipe threads require a source identity and exact geometry policy",
            ));
        }
        let transport =
            AuthoredThread::new(spec, policy.linear_tolerance_mm, policy.minimum_feature_mm)
                .map_err(failure)?;
        if spec.crest_radius_limit_mm.is_some() && frame != PlaneFrame::default() {
            return Err(error(
                "authored pipe thread",
                GeometryErrorKind::Unsupported,
                "Stock containment for a posed pipe construction is not yet certified",
            ));
        }
        match transport.stock_clip_plan(Default::default(), interrupt) {
            Ok(source) => return Self::from_clipped_npt_source(source, feature, policy, interrupt),
            Err(multicad_rust_kernel::axial_pipe_thread::Error::Unsupported) => {}
            Err(cause) => return Err(failure(cause)),
        }
        // The raw builder explicitly refuses active moving OD clips. A stock
        // limit is forwarded unchanged; no positive material is silently lost.
        let raw = transport
            .raw_boundary_ownership(policy.linear_tolerance_mm, interrupt)
            .map_err(|cause| {
                if cause == multicad_rust_kernel::radial_tensor::Error::Unsupported
                    && spec.crest_radius_limit_mm.is_some()
                {
                    error(
                        "authored pipe thread",
                        GeometryErrorKind::Unsupported,
                        "Moving or unresolved pipe stock clipping is not yet certified",
                    )
                } else {
                    failure(cause)
                }
            })?;
        let local = raw
            .ownership
            .bind_surfaces(policy.linear_tolerance_mm, interrupt)
            .and_then(|b| {
                b.certify_axial(
                    policy.linear_tolerance_mm,
                    policy.minimum_feature_mm,
                    interrupt,
                )
            })
            .and_then(|b| b.certify_caps(interrupt))
            .and_then(|b| b.certify_ordering(interrupt))
            .and_then(|b| b.certify_material(policy.linear_tolerance_mm, interrupt))
            .map_err(failure)?;
        let retained_radius_upper_bound_mm = if let Some(limit) = spec.crest_radius_limit_mm {
            let uncapped = AuthoredThread::new(
                AuthoredSpec {
                    crest_radius_limit_mm: None,
                    ..spec
                },
                policy.linear_tolerance_mm,
                policy.minimum_feature_mm,
            )
            .map_err(failure)?;
            let author = uncapped
                .radial_bounds([0., spec.pitch_mm], [0., spec.length_mm])
                .map_err(failure)?[1];
            // Positive additions round outward. The first receipt covers the
            // entire authored-law/net transfer and cap projection; the second
            // covers any closure admitted by the occupied-material certificate.
            let upper =
                ((author + raw.source_error_mm).next_up() + local.boundary_error_mm()).next_up();
            if !upper.is_finite() || upper > limit {
                return Err(error(
                    "authored pipe thread",
                    GeometryErrorKind::Unsupported,
                    "The complete retained pipe body is not proved inside the stock OD",
                ));
            }
            Some(upper)
        } else {
            None
        };
        let boundary = if frame == PlaneFrame::default() {
            local
        } else {
            local
                .transformed_with_interrupt(
                    [frame.origin, frame.u, frame.v],
                    policy.linear_tolerance_mm,
                    policy.minimum_feature_mm,
                    interrupt,
                )
                .map_err(failure)?
        };
        let source = RetainedAxialPipeThreadSource {
            transport,
            construction_frame: frame,
            cap_faces: raw.cap_faces,
            band_faces: raw.band_faces,
            root_faces: raw.root_faces,
            local_source_error_mm: raw.source_error_mm,
            retained_radius_upper_bound_mm,
        };
        let metadata = AuthoredPipe::new(source, boundary.input())?;
        Self::from_certified_radial_boundary(
            boundary,
            feature,
            policy,
            interrupt,
            None,
            Some(Rc::new(metadata)),
        )
    }
    pub fn rust_authored_axial_pipe_thread_source(
        &self,
    ) -> Option<(&RetainedAxialPipeThreadSource, PlaneFrame)> {
        let (body, occurrence) = self.radial_boundary_body()?;
        Some((&body.authored_pipe.as_ref()?.source, occurrence))
    }
}
fn invalid() -> GeometryError {
    error(
        "authored pipe thread",
        GeometryErrorKind::InvalidResult,
        "Authored pipe semantic incidence does not match its finite boundary",
    )
}
impl AuthoredPipe {
    fn new(
        source: RetainedAxialPipeThreadSource,
        input: &multicad_rust_kernel::radial_tensor::BoundaryInput,
    ) -> Result<Self, GeometryError> {
        let mut faces = vec![None; input.faces.len()];
        for (end, face) in source.cap_faces.into_iter().enumerate() {
            if !matches!(input.surfaces[input.faces[face].surface], Surface::Plane(_)) {
                return Err(invalid());
            }
            faces[face] = Some(if end == 0 { "cap.start" } else { "cap.end" }.to_string());
        }
        let names: &[&str] = match source.transport.spec().section {
            AuthoredSection::NptFlat { .. } => &["lower", "crest", "upper"],
            AuthoredSection::BsptRounded { .. } => {
                &["crest.leading", "upper", "root", "lower", "crest.trailing"]
            }
        };
        if source.band_faces.len() != names.len() {
            return Err(invalid());
        }
        for (&face, name) in source.band_faces.iter().zip(names) {
            if !matches!(
                input.surfaces[input.faces[face].surface],
                Surface::Tensor(_)
            ) || faces[face]
                .replace(format!("thread/start/0/{name}"))
                .is_some()
            {
                return Err(invalid());
            }
        }
        if source.root_faces.len() != 1
            || matches!(
                source.transport.spec().section,
                AuthoredSection::BsptRounded { .. }
            ) && source.root_faces[0] != source.band_faces[2]
        {
            return Err(invalid());
        }
        for &face in &source.root_faces {
            if faces[face].is_none() {
                if !matches!(
                    input.surfaces[input.faces[face].surface],
                    Surface::Cone { .. }
                ) {
                    return Err(invalid());
                }
                faces[face] = Some("thread/start/0/root".into());
            }
        }
        let faces = faces
            .into_iter()
            .collect::<Option<Vec<_>>>()
            .ok_or_else(invalid)?;
        let mut edge_faces = vec![BTreeSet::new(); input.edges.len()];
        let mut vertex_faces = vec![BTreeSet::new(); input.vertices.len()];
        for (i, face) in input.faces.iter().enumerate() {
            for &(edge, _) in &face.outer {
                edge_faces[edge].insert(faces[i].clone());
                for vertex in input.edges[edge].vertices {
                    vertex_faces[vertex].insert(faces[i].clone());
                }
            }
        }
        let join = |set: &BTreeSet<String>| set.iter().cloned().collect::<Vec<_>>().join("|");
        let vertices = vertex_faces
            .iter()
            .map(|set| format!("pipe/corner/{}", join(set)))
            .collect::<Vec<_>>();
        let edges = input
            .edges
            .iter()
            .enumerate()
            .map(|(i, edge)| {
                let mut endpoints = edge.vertices.map(|v| vertices[v].clone());
                endpoints.sort();
                format!(
                    "pipe/edge/{}/{}/{}",
                    join(&edge_faces[i]),
                    endpoints[0],
                    endpoints[1]
                )
            })
            .collect::<Vec<_>>();
        if edges.iter().collect::<BTreeSet<_>>().len() != edges.len()
            || vertices.iter().collect::<BTreeSet<_>>().len() != vertices.len()
        {
            return Err(invalid());
        }
        Ok(Self {
            source,
            faces,
            edges,
            vertices,
        })
    }
    pub(super) fn face_origins(&self, feature: EntityId, index: usize) -> Vec<Origin> {
        let mut roles = vec![Origin {
            source: feature,
            role: self.faces[index].clone(),
        }];
        if self.source.root_faces.contains(&index) {
            roles.push(Origin {
                source: feature,
                role: "thread.root".into(),
            });
        }
        roles
    }
    pub(super) fn apply_origin(&self, feature: EntityId, element: &mut TopologyElement) {
        let role = match element.key.kind {
            TopologyKind::Face => return,
            TopologyKind::Edge => &self.edges[element.key.ordinal as usize],
            TopologyKind::Vertex => &self.vertices[element.key.ordinal as usize],
        };
        element.origins = BTreeSet::from([Origin {
            source: feature,
            role: role.clone(),
        }]);
    }
}
