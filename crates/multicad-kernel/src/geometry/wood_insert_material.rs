//! Native paired-thread material. Solid topology/Model admission is separate.
mod faces;
use super::entry_chamfer_body::placement::measured::{self, SourceMesh};
use super::internal_thread_treatment_body::interval::Interval as I;
use super::*;
pub use faces::WoodInsertFaceProperties;
use multicad_rust_kernel::{
    display::{DisplayMesh, Tessellation},
    manufacturing,
    wood_insert::{self as source, Insert, QueryOptions},
};
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WoodInsertMaterial {
    source: Insert,
    frame: PlaneFrame,
}
fn failure(cause: source::Error) -> GeometryError {
    error(
        "Rust wood insert",
        match cause {
            source::Error::Invalid => GeometryErrorKind::Input,
            source::Error::Unsupported | source::Error::Unresolved => {
                GeometryErrorKind::Unsupported
            }
            source::Error::Limit => GeometryErrorKind::Limit,
            source::Error::Cancelled => GeometryErrorKind::Cancelled,
        },
        format!("Paired wood-insert source proof failed: {cause:?}"),
    )
}
fn checkpoint(interrupt: Option<&dyn Fn() -> bool>) -> Result<(), GeometryError> {
    if interrupt.is_some_and(|f| f()) {
        Err(failure(source::Error::Cancelled))
    } else {
        Ok(())
    }
}
impl WoodInsertMaterial {
    pub fn new(
        source: Insert,
        frame: PlaneFrame,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, GeometryError> {
        checkpoint(interrupt)?;
        policy.validate()?;
        let r = source.outer().maximum_radius_mm();
        policy.length(source.minimum_wall_mm())?;
        frame.checked_metric_extent([r, r, source.full_length_mm()], policy)?;
        let body = Self { source, frame };
        let bounds = body.bounds()?;
        policy.point(bounds.min)?;
        policy.point(bounds.max)?;
        checkpoint(interrupt)?;
        Ok(body)
    }
    pub fn source(&self) -> &Insert {
        &self.source
    }
    pub fn frame(&self) -> PlaneFrame {
        self.frame
    }
    pub fn bounds(&self) -> Result<Bounds, GeometryError> {
        placed_bounds(&self.source, self.frame)
    }
    /// Actual paired source volume under the literal affine frame. The mesh
    /// signed volume is a separate value and is never substituted here.
    pub fn volume_bounds(
        &self,
        options: QueryOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<[f64; 2], GeometryError> {
        let v = self
            .source
            .volume(options, interrupt)
            .map_err(failure)?
            .material_mm3;
        let [u, vaxis, n] =
            [self.frame.u, self.frame.v, self.frame.normal()].map(|a| a.map(I::point));
        let c = [(1, 2), (2, 0), (0, 1)].map(|(j, k)| u[j] * vaxis[k] - u[k] * vaxis[j]);
        let det = (0..3).fold(I::ZERO, |s, i| s + c[i] * n[i]);
        let actual = I::new(v[0], v[1]) * det;
        if !actual.valid() || actual.lo <= 0. || actual.width() > options.volume_error_mm3 {
            return Err(failure(source::Error::Unresolved));
        }
        checkpoint(interrupt)?;
        Ok([actual.lo, actual.hi])
    }
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
            let c = self
                .source
                .tessellate_candidate(quality, interrupt)
                .map_err(failure)?;
            let chord_mm = c.source_chord_bound_mm();
            let normal_rad = c.source_normal_angle_bound_rad();
            let m = c.into_mesh();
            Ok(SourceMesh {
                chord_mm,
                normal_rad,
                mesh: Mesh {
                    vertices: m.vertices,
                    indices: m.indices,
                    face_indices: m.face_indices,
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
        let m = DisplayMesh {
            vertices,
            indices,
            face_indices,
        };
        // Source residual and actual literal embedding are already proved.
        let p = certify(&m)?;
        checkpoint(interrupt)?;
        Ok(crate::ManufacturingMesh {
            mesh: Mesh {
                vertices: m.vertices,
                indices: m.indices,
                face_indices: m.face_indices,
            },
            maximum_displacement_mm: 0.,
            edge_pairs: p.edge_pairs,
            vertex_links: p.vertex_links,
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

pub(super) fn certify(m: &DisplayMesh) -> Result<manufacturing::MeshProof, GeometryError> {
    manufacturing::certify(m).map_err(|cause| {
        error(
            "Rust wood insert manufacturing mesh",
            GeometryErrorKind::InvalidResult,
            format!("{cause:?}"),
        )
    })
}

pub(super) fn placed_bounds(source: &Insert, frame: PlaneFrame) -> Result<Bounds, GeometryError> {
    let b = source.bounds_mm().map(|[lo, hi]| I::new(lo, hi));
    let axes = [frame.u, frame.v, frame.normal()];
    let placed: [I; 3] = std::array::from_fn(|i| {
        (0..3).fold(I::point(frame.origin[i]), |s, j| {
            s + I::point(axes[j][i]) * b[j]
        })
    });
    if placed.into_iter().any(|b| !b.valid()) {
        return Err(failure(source::Error::Unresolved));
    }
    Ok(Bounds {
        min: placed.map(|b| b.lo),
        max: placed.map(|b| b.hi),
    })
}
