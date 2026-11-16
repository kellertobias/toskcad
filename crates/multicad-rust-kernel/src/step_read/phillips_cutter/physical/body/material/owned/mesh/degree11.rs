//! Actual D11 owned material candidate; never enters the D9 source providers.
use super::*;
use crate::step_read::phillips_cutter::physical::body::{
    Degree11MaterialFace, OwnedDegree11Material,
};
mod thread;
#[derive(Clone)]
pub struct JoinedDegree11MaterialCandidate<'s> {
    source: &'s OwnedDegree11Material,
    mesh: DisplayMesh,
    connectivity: crate::manufacturing::MeshProof,
    distance: f64,
    quality: Tessellation,
    work: usize,
    remaining: usize,
    peak: usize,
    delegated: DelegatedWorkspaceBounds,
    combined_peak: usize,
}
impl<'s> JoinedDegree11MaterialCandidate<'s> {
    pub fn source(&self) -> &'s OwnedDegree11Material {
        self.source
    }
    pub fn mesh(&self) -> &DisplayMesh {
        &self.mesh
    }
    pub fn connectivity(&self) -> crate::manufacturing::MeshProof {
        self.connectivity
    }
    pub fn certified_tessellation(&self) -> Tessellation {
        self.quality
    }
    pub fn new_surface_error_mm(&self) -> f64 {
        self.distance
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn remaining_work(&self) -> usize {
        self.remaining
    }
    pub fn accounted_peak_bytes(&self) -> usize {
        self.peak
    }
    pub fn delegated_workspace_bounds(&self) -> DelegatedWorkspaceBounds {
        self.delegated
    }
    pub fn combined_live_upper_bound_bytes(&self) -> usize {
        self.combined_peak
    }
}
fn head(
    source: &OwnedDegree11Material,
    angular: &mut Vec<AngularCurve<'_, 12>>,
    quality: Tessellation,
    work: &mut Work<'_>,
) -> R<Head> {
    let n = source.head().view().segment_count();
    work.storage(n * std::mem::size_of::<usize>())?;
    let mut counts = vec![1usize; n];
    work.held += counts.capacity() * std::mem::size_of::<usize>();
    for _ in 0..64 {
        match head_attempt_curve(source.head(), angular, &counts, quality, work)? {
            Ok(head) => return Ok(head),
            Err(Split::T(part)) => {
                counts[part] = counts[part]
                    .checked_mul(2)
                    .ok_or(PhillipsReadError::Limit)?
            }
            Err(Split::U) => {
                let old = angular.capacity() * std::mem::size_of::<AngularCurve<'_, 12>>();
                let count = angular
                    .len()
                    .checked_mul(2)
                    .ok_or(PhillipsReadError::Limit)?;
                work.storage(count * std::mem::size_of::<AngularCurve<'_, 12>>())?;
                let mut children = Vec::with_capacity(count);
                work.storage(children.capacity() * std::mem::size_of::<AngularCurve<'_, 12>>())?;
                for &a in angular.iter() {
                    work.charge(1)?;
                    let m = a.range[0] + (a.range[1] - a.range[0]) / 2.;
                    if m <= a.range[0] || m >= a.range[1] {
                        return Err(PhillipsReadError::Unsupported);
                    }
                    children.extend([
                        AngularCurve {
                            range: [a.range[0], m],
                            ..a
                        },
                        AngularCurve {
                            range: [m, a.range[1]],
                            ..a
                        },
                    ]);
                }
                work.held = work.held - old
                    + children.capacity() * std::mem::size_of::<AngularCurve<'_, 12>>();
                work.storage(old)?;
                *angular = children;
            }
        }
    }
    Err(PhillipsReadError::Limit)
}
fn cavity_face(source: &OwnedDegree11Material, key: CavityFace) -> Option<u32> {
    (0..source.face_count()).find_map(|i| match source.face_source(i)? {
        Degree11MaterialFace::ReversedCavity { face, .. } if face == key => Some(i as u32),
        _ => None,
    })
}
impl OwnedDegree11Material {
    /// Facets of this SAME detached degree11 material. The whole source C0,
    /// actual normal and closed indexed connectivity gates are mandatory.
    /// A later consumer must still earn global literal embedding and physical
    /// manufacturing/publication; no nominal source is reconstructed here.
    pub fn joined_mesh_candidate(
        &self,
        query: JoinedMeshOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> R<JoinedDegree11MaterialCandidate<'_>> {
        let options = query.tessellation;
        options.check().map_err(|_| PhillipsReadError::Invalid)?;
        if query.maximum_work > MESH_WORK || query.maximum_bytes > MESH_BYTES {
            return Err(PhillipsReadError::Invalid);
        }
        let quality = source_quality(options)?;
        let mut work = Work {
            remaining: query.maximum_work,
            held: self.retained_budget_bytes(),
            peak: 0,
            maximum_bytes: query.maximum_bytes,
            maximum_triangles: options.triangle_limit,
            delegated: Default::default(),
            interrupt,
        };
        work.storage(thread::SCRATCH)?;
        let thread = thread::build(self, options, &mut work)?;
        work.held += mesh_bytes(&thread.mesh)
            + (thread.entry.capacity() + thread.tip.capacity()) * 4
            + thread.curves.capacity() * std::mem::size_of::<[P; 12]>()
            + thread.divisions.capacity() * std::mem::size_of::<usize>();
        let angular_count: usize = thread.divisions.iter().sum();
        work.storage(angular_count * std::mem::size_of::<AngularCurve<'_, 12>>())?;
        let mut angular = Vec::with_capacity(angular_count);
        for (curve, &divisions) in thread.curves.iter().zip(&thread.divisions) {
            for k in 0..divisions {
                angular.push(AngularCurve {
                    poles: curve,
                    spans: 1,
                    reverse: false,
                    range: [
                        k as f64 / divisions as f64,
                        (k + 1) as f64 / divisions as f64,
                    ],
                });
            }
        }
        work.held += angular.capacity() * std::mem::size_of::<AngularCurve<'_, 12>>();
        work.storage(0)?;
        let head = head(self, &mut angular, options, &mut work)?;
        work.held += mesh_bytes(&head.mesh) + (head.top.capacity() + head.bottom.capacity()) * 4;
        work.storage(1024 * 1024)?;
        work.charge(4096)?;
        let atlas = self.cavity().surface_patches().map_err(cavity_error)?;
        work.held += atlas.retained_budget_bytes();
        work.storage(0)?;
        work.delegated.cutter_bytes = delegated::cutter(options.triangle_limit)?;
        let cavity = atlas
            .boundary_mesh_with_mouth_trace(
                Tessellation {
                    deflection: options.deflection * 0.99,
                    angle: quality.source_angle,
                    ..options
                },
                crate::phillips_drive::MouthTraceOptions {
                    maximum_work: work.remaining,
                    ..Default::default()
                },
                &|| interrupt.is_some_and(|f| f()),
            )
            .map_err(cavity_error)?;
        work.charge(cavity.work_used())?;
        work.held += mesh_bytes(cavity.boundary_mesh().mesh()) + cavity.retained_metadata_bytes();
        work.storage(0)?;
        let h = self.head().view().segment_count();
        let t = self.thread_faces().len();
        let mut mesh = DisplayMesh {
            vertices: Vec::new(),
            indices: Vec::new(),
            face_indices: Vec::new(),
        };
        let thread_offset = append_mesh(
            &mut mesh,
            &thread.mesh,
            &mut |f| Some(h as u32 + f),
            0.,
            &mut work,
        )?;
        let head_offset = append_mesh(&mut mesh, &head.mesh, &mut Some, 0., &mut work)?;
        let cavity_faces = atlas.topology().faces();
        work.storage(cavity_faces.len() * std::mem::size_of::<Option<u32>>())?;
        work.charge(cavity_faces.len() * (self.face_count() + 1))?;
        let mut cavity_map = Vec::with_capacity(cavity_faces.len());
        for face in cavity_faces {
            let key = face.source_identity();
            cavity_map.push(if key == CavityFace::Mouth {
                None
            } else {
                Some(cavity_face(self, key).ok_or(PhillipsReadError::Unsupported)?)
            });
        }
        work.held += cavity_map.capacity() * std::mem::size_of::<Option<u32>>();
        let cavity_offset = append_mesh(
            &mut mesh,
            cavity.boundary_mesh().mesh(),
            &mut |f| cavity_map[f as usize],
            self.chart().outer_top_z,
            &mut work,
        )?;
        let translated = translation_transfer(
            cavity.boundary_mesh().mesh(),
            &mesh,
            cavity_offset,
            self.chart().outer_top_z,
            options.deflection * 0.01,
            quality.angular_reserve,
            &mut work,
        )?;
        for (tri, &face) in mesh
            .indices
            .as_chunks_mut::<3>()
            .0
            .iter_mut()
            .zip(&mesh.face_indices)
        {
            work.charge(1)?;
            if face as usize >= h + t + 3 {
                tri.swap(1, 2);
            }
        }
        let ring_bytes = (head.bottom.len()
            + head.top.len()
            + thread.entry.len()
            + thread.tip.len()
            + cavity.mouth().len())
            * 4;
        work.storage(mesh_bytes(&mesh) + ring_bytes)?;
        let bottom: Vec<_> = head.bottom.iter().map(|&i| i + head_offset).collect();
        let top: Vec<_> = head.top.iter().map(|&i| i + head_offset).collect();
        let entry: Vec<_> = thread.entry.iter().map(|&i| i + thread_offset).collect();
        let tip: Vec<_> = thread.tip.iter().map(|&i| i + thread_offset).collect();
        let mouth: Vec<_> = cavity
            .mouth()
            .iter()
            .map(|p| p.vertex + cavity_offset)
            .collect();
        work.held += (bottom.capacity()
            + top.capacity()
            + entry.capacity()
            + tip.capacity()
            + mouth.capacity())
            * 4;
        work.storage(mesh_bytes(&mesh))?;
        let tip_error = plane(
            &mut mesh,
            &tip,
            &[],
            false,
            PlaneProof {
                face: (h + t) as u32,
                positive: true,
                source_z: self.cap_stations_mm()[1],
                boundary_error: thread.error,
                options,
            },
            &mut work,
        )?;
        let join_error = plane(
            &mut mesh,
            &bottom,
            &entry,
            false,
            PlaneProof {
                face: (h + t + 1) as u32,
                positive: true,
                source_z: self.cap_stations_mm()[0],
                boundary_error: thread.error.max(head.distance),
                options,
            },
            &mut work,
        )?;
        let cavity_error =
            (I::point(cavity.boundary_mesh().source_distance_mm()) + I::point(translated)).hi;
        let top_error = plane(
            &mut mesh,
            &top,
            &mouth,
            true,
            PlaneProof {
                face: (h + t + 2) as u32,
                positive: false,
                source_z: self.chart().outer_top_z,
                boundary_error: cavity_error.max(head.distance),
                options,
            },
            &mut work,
        )?;
        let mesh = compact(&mesh, &mut work)?;
        // Same exact connectivity checks, with every operation, cancellation
        // point and actual temporary capacity charged to THIS remaining query.
        // These bytes are in producer peak, not an opaque delegated allowance.
        let mesh_storage = mesh_bytes(&mesh);
        let connectivity = crate::manufacturing::certify_counted(&mesh, &mut |steps, scratch| {
            work.charge(steps)?;
            if let Some(bytes) = scratch {
                work.storage(mesh_storage + bytes)?;
            }
            Ok(())
        })
        .map_err(|error| match error {
            crate::manufacturing::CountedError::Mesh(error) => unresolved_mesh(error),
            crate::manufacturing::CountedError::Interrupted(error) => error,
        })?;
        work.storage(mesh_bytes(&mesh))?;
        if mesh.face_indices.len() > options.triangle_limit {
            return Err(PhillipsReadError::Limit);
        }
        work.charge(1)?;
        Ok(JoinedDegree11MaterialCandidate {
            source: self,
            mesh,
            connectivity,
            distance: thread
                .error
                .max(head.distance)
                .max(tip_error)
                .max(join_error)
                .max(top_error)
                .max(cavity_error),
            quality: options,
            work: query.maximum_work - work.remaining,
            remaining: work.remaining,
            peak: work.peak,
            delegated: work.delegated,
            combined_peak: work
                .peak
                .checked_add(work.delegated.maximum())
                .ok_or(PhillipsReadError::Limit)?,
        })
    }
}
#[cfg(test)]
mod tests;
