//! Print-quality meshes of authored tapered pipe laws. The dedicated law
//! tessellation proves each facet's distance to the authored law and its
//! normal angle; callers add their own law-to-boundary bound before
//! admitting the mesh, and map law roles to their own B-rep faces.
use super::*;
use multicad_rust_kernel::axial_pipe_thread::Error as Cause;

/// B-rep face of every law role; `extra` is the stock wall of a sleeve.
/// Law role indices are `band < bands`, Root, Clip, StartCap, EndCap, then
/// the extra wall, as in `LawMesh::role_index`.
#[derive(Clone, Copy, Debug)]
pub(crate) struct LawFaces {
    pub(crate) band_count: u32,
    pub(crate) bands: [u32; 5],
    pub(crate) root: u32,
    pub(crate) clip: u32,
    pub(crate) caps: [u32; 2],
    pub(crate) extra: u32,
}
impl LawFaces {
    pub(crate) fn extra_role(&self) -> u32 {
        self.band_count + 4
    }
    pub(crate) fn remap(&self, faces: &mut [u32]) -> Result<(), GeometryError> {
        let n = self.band_count;
        for face in faces {
            *face = match *face {
                b if b < n && b < 5 => self.bands[b as usize],
                r if r == n => self.root,
                c if c == n + 1 => self.clip,
                s if s == n + 2 => self.caps[0],
                e if e == n + 3 => self.caps[1],
                x if x == n + 4 => self.extra,
                _ => {
                    return Err(error(
                        "authored pipe print mesh",
                        GeometryErrorKind::InvalidResult,
                        "A law facet has no B-rep face",
                    ));
                }
            };
        }
        Ok(())
    }
}
pub(crate) fn failure(cause: Cause) -> GeometryError {
    error(
        "authored pipe print mesh",
        match cause {
            Cause::Invalid => GeometryErrorKind::Input,
            Cause::Limit => GeometryErrorKind::Limit,
            Cause::Cancelled => GeometryErrorKind::Cancelled,
            Cause::Unsupported | Cause::Unresolved => GeometryErrorKind::Unsupported,
        },
        format!("Authored law tessellation failed: {cause:?}"),
    )
}
/// A law facet bound plus the actual boundary's distance to its law must
/// still meet the requested deflection.
pub(crate) fn admit(
    deviation: f64,
    boundary_error: f64,
    deflection: f64,
) -> Result<(), GeometryError> {
    let total = (deviation + boundary_error).next_up();
    if !total.is_finite() || total > deflection {
        return Err(error(
            "authored pipe print mesh",
            GeometryErrorKind::Unsupported,
            format!(
                "Law facets are {total} mm from the actual boundary, above the requested deflection"
            ),
        ));
    }
    Ok(())
}
/// Closed oriented manufacturing receipt of an already placed mesh, with
/// its volume compared to the solid's exact volume within the facet bound.
pub(crate) fn receipt(
    mesh: Mesh,
    exact_volume: [f64; 2],
    distance: f64,
) -> Result<crate::ManufacturingMesh, GeometryError> {
    let display = multicad_rust_kernel::display::DisplayMesh {
        vertices: mesh.vertices,
        indices: mesh.indices,
        face_indices: mesh.face_indices,
    };
    let proof = multicad_rust_kernel::manufacturing::certify(&display).map_err(|cause| {
        error(
            "authored pipe print mesh",
            GeometryErrorKind::InvalidResult,
            format!("Law print mesh is not one closed oriented shell: {cause:?}"),
        )
    })?;
    let mut area = 0f64;
    for t in display.indices.as_chunks::<3>().0 {
        let p = t.map(|i| {
            let i = 3 * i as usize;
            [
                display.vertices[i],
                display.vertices[i + 1],
                display.vertices[i + 2],
            ]
        });
        let a: [f64; 3] = std::array::from_fn(|k| p[1][k] - p[0][k]);
        let b: [f64; 3] = std::array::from_fn(|k| p[2][k] - p[0][k]);
        let c = [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ];
        area += 0.5 * (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt();
    }
    // Every facet lies within `distance` of the actual boundary.
    let reserve = distance * area * 1.01 + (exact_volume[1] - exact_volume[0]) + 1e-9;
    let middle = exact_volume[0] + (exact_volume[1] - exact_volume[0]) * 0.5;
    let gap = (proof.signed_volume_mm3 - middle).abs();
    if gap.is_nan() || gap > reserve {
        return Err(error(
            "authored pipe print mesh",
            GeometryErrorKind::InvalidResult,
            "Law print mesh volume differs from the exact solid beyond its facet bound",
        ));
    }
    Ok(crate::ManufacturingMesh {
        mesh: Mesh {
            vertices: display.vertices,
            indices: display.indices,
            face_indices: display.face_indices,
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
