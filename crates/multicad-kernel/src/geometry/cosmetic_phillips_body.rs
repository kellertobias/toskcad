//! An actual Phillips head and recess on an explicitly cosmetic smooth shaft.
//!
//! The physical source remains the retained Phillips screw: its head and
//! recess faces, edges and vertices are the physical ones. Only the visible
//! thread interval is replaced by its nominal cylindrical envelope. No
//! manufacturing mesh or STEP admission follows from this conversion.
use super::*;
use crate::rust_prism_topology::position;
use multicad_core::references::{
    AnalyticGeometry, AnchoredFrame, Origin, TopologyElement, TopologyKey, TopologyKind,
};
use multicad_rust_kernel::{
    cosmetic_profiled_screw::PhillipsEnvelope,
    phillips_drive::{PhillipsScrew, ScrewFace as Role},
};
use std::collections::BTreeSet;
use std::f64::consts::{PI, TAU};

const OPERATION: &str = "cosmetic Phillips screw";

pub(crate) struct CosmeticPhillipsBody {
    pub(super) feature: EntityId,
    geometry: PhillipsEnvelope,
    /// The retained physical source placed at `frame`; it owns the head and
    /// recess boundary, never the shaft of this envelope.
    physical: Solid,
    frame: PlaneFrame,
}

fn is_thread(role: Role) -> bool {
    matches!(
        role,
        Role::ThreadShoulder | Role::ThreadTip | Role::ThreadWall(_)
    )
}

impl Solid {
    /// Display the actual Phillips head and recess over the nominal smooth
    /// envelope of the visible thread interval.
    pub fn rust_cosmetic_phillips_screw(
        source: PhillipsScrew,
        feature: EntityId,
        frame: PlaneFrame,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        let physical =
            Self::rust_phillips_screw_from_source(source.clone(), feature, frame, policy)?;
        let body = CosmeticPhillipsBody {
            feature,
            geometry: PhillipsEnvelope::new(source),
            physical,
            frame,
        };
        frame.checked_metric_bounds(body.local_bounds(), policy)?;
        let bounds = body.bounds();
        policy.point(bounds.min)?;
        policy.point(bounds.max)?;
        Ok(Self {
            body: SolidBody::CosmeticPhillips(Rc::new(body)),
            policy,
            provenance: vec![],
            rust_origins: vec![],
            retained_placement: None,
        })
    }
    pub(crate) fn cosmetic_phillips_body(&self) -> Option<&CosmeticPhillipsBody> {
        match &self.body {
            SolidBody::CosmeticPhillips(body) => Some(body),
            _ => None,
        }
    }
}

impl CosmeticPhillipsBody {
    pub(super) fn transformed(&self, transform: RigidTransform) -> Result<Self, GeometryError> {
        Ok(Self {
            feature: self.feature,
            geometry: self.geometry.clone(),
            physical: self.physical.transformed(transform)?,
            frame: transform.frame(self.frame),
        })
    }
    fn source(&self) -> &PhillipsScrew {
        self.geometry.head_source()
    }
    pub(super) fn volume_mm3(&self) -> f64 {
        self.geometry.volume_mm3()
    }
    fn local_bounds(&self) -> Bounds {
        self.placed_bounds(PlaneFrame::default())
    }
    fn placed_bounds(&self, frame: PlaneFrame) -> Bounds {
        let axes = [frame.u, frame.v, frame.normal()];
        let support = |i, sign: f64| self.geometry.support(axes.map(|a| sign * a[i]));
        Bounds {
            min: std::array::from_fn(|i| frame.origin[i] - support(i, -1.)),
            max: std::array::from_fn(|i| frame.origin[i] + support(i, 1.)),
        }
    }
    pub(super) fn bounds(&self) -> Bounds {
        self.placed_bounds(self.frame)
    }
    /// Index of a retained role in the physical source face order.
    fn physical_index(&self, role: Role) -> usize {
        self.source()
            .faces()
            .iter()
            .position(|r| *r == role)
            .expect("retained physical Phillips face")
    }
    fn index(&self, role: Role) -> usize {
        self.geometry
            .faces()
            .iter()
            .position(|r| *r == role)
            .expect("cosmetic Phillips face")
    }
    fn cap_plane(&self, z: f64) -> ReferenceFrame {
        ReferenceFrame {
            origin_mm: position(self.frame, [0., 0., z]),
            u: self.frame.u,
            v: self.frame.v,
        }
    }
    fn last_profile_wall(&self) -> Role {
        Role::ProfileWall(self.source().profile().segments().count() as u32 - 1)
    }
    pub(super) fn faces(&self) -> Result<Vec<Face>, GeometryError> {
        let physical = self.physical.faces()?;
        let [a, b] = self.geometry.shaft_interval_mm();
        let r = self.geometry.shaft_radius_mm();
        let end = self.source().profile().end()[0];
        Ok(self
            .geometry
            .faces()
            .iter()
            .enumerate()
            .map(|(i, &role)| {
                let index = i as u32;
                let cap = |z: f64, area| {
                    let plane = self.cap_plane(z);
                    Face {
                        index,
                        plane: Some((plane.origin_mm, plane.normal())),
                        center_mm: plane.origin_mm,
                        area_mm2: area,
                    }
                };
                match role {
                    Role::ThreadShoulder => cap(a, PI * (end - r) * (end + r)),
                    Role::ThreadTip => cap(b, PI * r * r),
                    Role::ThreadWall(_) => Face {
                        index,
                        plane: None,
                        center_mm: position(self.frame, [0., 0., a.midpoint(b)]),
                        area_mm2: TAU * r * (b - a),
                    },
                    role => Face {
                        index,
                        ..physical[self.physical_index(role)].clone()
                    },
                }
            })
            .collect())
    }
    fn origin(&self, role: impl Into<String>) -> Origin {
        Origin {
            source: self.feature,
            role: role.into(),
        }
    }
    fn circle(&self, z: f64, radius: f64) -> AnalyticGeometry {
        AnalyticGeometry::CircularArc {
            frame: self.cap_plane(z),
            radius_mm: radius,
            start_rad: 0.,
            sweep_rad: TAU,
        }
    }
    /// The physical head and recess topology, with the helical thread faces
    /// and every edge bounding them replaced by the smooth envelope.
    pub(crate) fn topology(
        &self,
        owner: EntityId,
        revision: u64,
    ) -> Result<TopologySnapshot, GeometryError> {
        let physical = self.physical.topology_snapshot(owner, revision)?;
        let roles = self.source().faces();
        let face_role = |element: &TopologyElement| {
            (element.key.kind == TopologyKind::Face)
                .then(|| roles.get(element.key.ordinal as usize).copied())
                .flatten()
        };
        let replaced: BTreeSet<Origin> = physical
            .elements
            .iter()
            .filter(|e| face_role(e).is_some_and(is_thread))
            .flat_map(|e| e.boundary.iter().cloned())
            .collect();
        let mut faces: Vec<Option<TopologyElement>> = vec![None; self.geometry.faces().len()];
        let mut edges = Vec::new();
        let mut vertices = Vec::new();
        for mut element in physical.elements {
            match element.key.kind {
                TopologyKind::Face => {
                    let Some(role) = face_role(&element) else {
                        return Err(error(
                            OPERATION,
                            GeometryErrorKind::InvalidResult,
                            "Physical Phillips topology has an unowned face",
                        ));
                    };
                    if is_thread(role) {
                        continue;
                    }
                    element.boundary.retain(|origin| !replaced.contains(origin));
                    let index = self.index(role);
                    element.key.ordinal = index as u32;
                    faces[index] = Some(element);
                }
                TopologyKind::Edge => {
                    if !element.origins.iter().any(|o| replaced.contains(o)) {
                        edges.push(element);
                    }
                }
                _ => {
                    // Helical rail endpoints belong to the replaced thread.
                    if element
                        .origins
                        .iter()
                        .all(|o| !o.role.starts_with("thread/"))
                    {
                        vertices.push(element);
                    }
                }
            }
        }
        let [a, b] = self.geometry.shaft_interval_mm();
        let r = self.geometry.shaft_radius_mm();
        for (role, name, plane) in [
            (Role::ThreadShoulder, "cap.start", Some(self.cap_plane(a))),
            (Role::ThreadTip, "cap.end", Some(self.cap_plane(b))),
            (Role::ThreadWall(0), "thread.cosmetic.envelope", None),
        ] {
            if !self.geometry.faces().contains(&role) {
                continue;
            }
            let origin = self.origin(name);
            let index = self.index(role);
            faces[index] = Some(TopologyElement {
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
        let shoulder = self.geometry.faces().contains(&Role::ThreadShoulder);
        let mut added = vec![
            (
                "cosmetic.shaft/tip",
                self.circle(b, r),
                [Role::ThreadTip, Role::ThreadWall(0)],
            ),
            (
                "cosmetic.shaft/start",
                self.circle(a, r),
                [
                    Role::ThreadWall(0),
                    if shoulder {
                        Role::ThreadShoulder
                    } else {
                        self.last_profile_wall()
                    },
                ],
            ),
        ];
        if shoulder {
            added.push((
                "screw.head/thread.junction/outer",
                self.circle(a, self.source().profile().end()[0]),
                [Role::ThreadShoulder, self.last_profile_wall()],
            ));
        }
        let mut faces: Vec<TopologyElement> = faces
            .into_iter()
            .map(|face| {
                face.ok_or_else(|| {
                    error(
                        OPERATION,
                        GeometryErrorKind::InvalidResult,
                        "Physical Phillips topology lacks a retained head face",
                    )
                })
            })
            .collect::<Result<_, _>>()?;
        for (name, geometry, owners) in added {
            let origin = self.origin(name);
            for role in owners {
                faces[self.index(role)].boundary.insert(origin.clone());
            }
            edges.push(TopologyElement {
                key: TopologyKey {
                    kind: TopologyKind::Edge,
                    ordinal: 0,
                },
                origins: BTreeSet::from([origin]),
                boundary: BTreeSet::new(),
                frames: vec![],
                geometry: Some(geometry),
            });
        }
        for (ordinal, edge) in edges.iter_mut().enumerate() {
            edge.key.ordinal = ordinal as u32;
        }
        for (ordinal, vertex) in vertices.iter_mut().enumerate() {
            vertex.key.ordinal = ordinal as u32;
        }
        let mut elements = faces;
        elements.extend(edges);
        elements.extend(vertices);
        Ok(TopologySnapshot {
            owner,
            revision,
            elements,
        })
    }
    pub(super) fn mesh(
        &self,
        options: multicad_rust_kernel::display::Tessellation,
    ) -> Result<Mesh, GeometryError> {
        use multicad_rust_kernel::helical_thread_mesh::PhillipsMeshError;
        let local = self
            .geometry
            .mesh(options, &|| false)
            .map_err(|cause| match cause {
                PhillipsMeshError::Geometry(cause) => cause,
                PhillipsMeshError::Cancelled => {
                    multicad_rust_kernel::display::MeshError::UnresolvedBoundary
                }
            });
        crate::rust_display::placed_mesh(local, self.frame, false)
    }
}

#[cfg(test)]
mod tests;
