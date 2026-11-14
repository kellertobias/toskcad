use super::*;
use crate::geometry::entry_chamfer_body::placement::measured::{self, SourceMesh};
use multicad_rust_kernel::{
    display::{DisplayMesh, Tessellation},
    manufacturing,
};
impl HangerMaterial {
    pub fn mesh(
        &self,
        options: MeshOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Mesh, GeometryError> {
        checkpoint(interrupt)?;
        let quality = Tessellation {
            deflection: options.linear_deflection_mm,
            angle: options.angular_deflection_rad,
            triangle_limit: options.maximum_triangles as usize,
        };
        let mesh = measured::qualify(self.frame, quality, &mut |quality| {
            let candidate = self
                .source
                .tessellate_candidate(quality, 16_000_000, interrupt)
                .map_err(failure)?;
            let chord_mm = candidate.source_chord_bound_mm();
            let normal_rad = candidate.source_normal_angle_bound_rad();
            let mesh = candidate.into_mesh();
            Ok(SourceMesh {
                chord_mm,
                normal_rad,
                mesh: Mesh {
                    vertices: mesh.vertices,
                    indices: mesh.indices,
                    face_indices: mesh.face_indices,
                },
            })
        })?;
        checkpoint(interrupt)?;
        Ok(mesh)
    }
    pub fn manufacturing_mesh(
        &self,
        options: MeshOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<crate::ManufacturingMesh, GeometryError> {
        let Mesh {
            vertices,
            indices,
            face_indices,
        } = self.mesh(options, interrupt)?;
        let mesh = DisplayMesh {
            vertices,
            indices,
            face_indices,
        };
        let proof = certify(&mesh)?;
        checkpoint(interrupt)?;
        Ok(crate::ManufacturingMesh {
            mesh: Mesh {
                vertices: mesh.vertices,
                indices: mesh.indices,
                face_indices: mesh.face_indices,
            },
            maximum_displacement_mm: 0.,
            edge_pairs: proof.edge_pairs,
            vertex_links: proof.vertex_links,
            removed_native_pole_triangles: 0,
            removed_exact_zero_triangles: 0,
            exact_vertex_pairs: 0,
            internal_constraint_edges: 0,
            removed_meshing_internal_edges: 0,
            removed_meshing_constraint_wire_containers: 0,
            removed_meshing_internal_vertices: 0,
            remeshed_without_internal_constraints: false,
        })
    }
}

pub(super) fn certify(mesh: &DisplayMesh) -> Result<manufacturing::MeshProof, GeometryError> {
    manufacturing::certify(mesh).map_err(|cause| {
        error(
            "Rust authored hanger manufacturing mesh",
            GeometryErrorKind::InvalidResult,
            format!("{cause:?}"),
        )
    })
}
