//! Conformal source-chart approximation. This stage certifies shared indices,
//! source distance and facet normals; global embedding is a separate gate.
use super::*;
use crate::display::{DisplayMesh, MeshError, Tessellation};
use crate::manufacturing::{self, MeshProof};
use std::collections::HashMap;
mod schema;
use schema::{Complex, Side};

/// A closed indexed source approximation. This is not a manufacturing receipt:
/// global geometric embedding and the head/cavity union remain separate gates.
#[derive(Clone, Debug)]
pub struct BoundaryMesh {
    mesh: DisplayMesh,
    connectivity: MeshProof,
    source_distance_mm: f64,
    normal_chord: f64,
    refinement: TopologyCertificate,
}
impl BoundaryMesh {
    pub fn mesh(&self) -> &DisplayMesh {
        &self.mesh
    }
    pub fn connectivity(&self) -> MeshProof {
        self.connectivity
    }
    pub fn source_distance_mm(&self) -> f64 {
        self.source_distance_mm
    }
    pub fn normal_chord(&self) -> f64 {
        self.normal_chord
    }
    pub fn refinement(&self) -> TopologyCertificate {
        self.refinement
    }
}
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
enum PointKey {
    Anchor(usize),
    Edge(usize, u64),
    Interior(usize, u64, u64),
}
struct Builder<'a, 's> {
    atlas: &'a PatchAtlas<'s>,
    complex: Complex,
    points: HashMap<PointKey, usize>,
    mesh: DisplayMesh,
}
fn groups(law: Law) -> [usize; 2] {
    match law {
        Law::Envelope(i, _) => [i, 4],
        Law::Wall(i, _, _) => [i, 5],
        Law::Trace(i, _, _) => [i, 6],
        Law::Gap(_) => [4, 6],
        Law::Floor(_) => [7, 4],
    }
}
impl<'a, 's> Builder<'a, 's> {
    fn new(atlas: &'a PatchAtlas<'s>, cancel: &dyn Fn() -> bool) -> Result<Self, Error> {
        Ok(Self {
            atlas,
            complex: Complex::build(atlas, cancel)?,
            points: HashMap::new(),
            mesh: DisplayMesh {
                vertices: Vec::new(),
                indices: Vec::new(),
                face_indices: Vec::new(),
            },
        })
    }
    fn point(&mut self, patch: usize, uv: [f64; 2]) -> Result<usize, Error> {
        let side = if uv[0] == 0. {
            Some((0, uv[1]))
        } else if uv[0] == 1. {
            Some((1, uv[1]))
        } else if uv[1] == 0. {
            Some((2, uv[0]))
        } else if uv[1] == 1. {
            Some((3, uv[0]))
        } else {
            None
        };
        let key = match side {
            Some((side, t)) => match self.complex.sides[patch][side] {
                Side::Collapsed(node) => PointKey::Anchor(node),
                Side::Curve(edge) if t == 0. || t == 1. => {
                    PointKey::Anchor(self.complex.edges[edge].ends[t as usize])
                }
                Side::Curve(edge) => PointKey::Edge(edge, t.to_bits()),
            },
            None => PointKey::Interior(patch, uv[0].to_bits(), uv[1].to_bits()),
        };
        if let Some(&id) = self.points.get(&key) {
            return Ok(id);
        }
        let point = match key {
            PointKey::Anchor(node) => self.complex.points[node],
            PointKey::Edge(edge, bits) => {
                self.complex
                    .edge_point(self.atlas, edge, f64::from_bits(bits))?
            }
            PointKey::Interior(..) => self.atlas.state(patch, uv.map(I::point))?.0.map(|p| p.v),
        };
        let id = self.mesh.vertices.len() / 3;
        self.mesh.vertices.extend(point.map(I::midpoint));
        self.points.insert(key, id);
        Ok(id)
    }
    fn literal(&self, indices: [usize; 4]) -> [[f64; 3]; 4] {
        indices.map(|id| std::array::from_fn(|i| self.mesh.vertices[3 * id + i]))
    }
}
fn cross(a: [I; 3], b: [I; 3]) -> [I; 3] {
    std::array::from_fn(|i| a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3])
}
fn facet_chord(
    points: [[f64; 3]; 3],
    center: [I; 3],
    source_spread: f64,
    law: Law,
) -> Result<f64, Error> {
    let [a, b, c] = points.map(|p| p.map(I::point));
    let n = cross(
        std::array::from_fn(|i| b[i] - a[i]),
        std::array::from_fn(|i| c[i] - a[i]),
    );
    let length = norm(n);
    if !length.valid() || length.lo <= 0. {
        return Err(Error::Unresolved);
    }
    let projected = match law {
        Law::Wall(_, q, _) => rotate([n[0], n[1]], (4 - q) % 4)[0],
        Law::Trace(..) | Law::Gap(_) => -n[2],
        _ => n[2],
    };
    if projected.lo <= 0. {
        return Err(Error::Unresolved);
    }
    let chord = norm(std::array::from_fn(|i| n[i] / length - center[i])) + I::point(source_spread);
    if !chord.valid() {
        return Err(Error::Unresolved);
    }
    Ok(chord.hi)
}
fn mesh_error(error: MeshError) -> Error {
    match error {
        MeshError::InvalidOptions => Error::Invalid,
        MeshError::BudgetExceeded => Error::Limit,
        MeshError::UnresolvedBoundary => Error::Unresolved,
    }
}
fn resolve_facet(
    result: Result<f64, Error>,
    law: Law,
    refine: &mut [bool; 8],
    groups: [usize; 2],
) -> Option<f64> {
    match result {
        Ok(chord) => Some(chord),
        Err(_) => {
            let [u, v] = groups;
            // Tangential contact needs angular refinement independently of the
            // radial parameter; refining both preserves a folded aspect ratio.
            if !matches!(law, Law::Envelope(..)) {
                refine[u] = true;
            }
            refine[v] = true;
            None
        }
    }
}
fn quad_normal_score(
    builder: &mut Builder<'_, '_>,
    index: usize,
    uv: [[f64; 2]; 2],
) -> Result<f64, Error> {
    let pairs = [
        [uv[0][0], uv[1][0]],
        [uv[0][1], uv[1][0]],
        [uv[0][1], uv[1][1]],
        [uv[0][0], uv[1][1]],
    ];
    let [a, b, c, d] = pairs.map(|p| builder.point(index, p));
    let ids = [a?, b?, c?, d?];
    let literal = builder.literal(ids);
    let bound = builder.atlas.certify_quad(index, uv, literal)?;
    let center = builder
        .atlas
        .state(index, uv.map(|[a, b]| I::point(a + (b - a) / 2.)))?
        .1;
    let mut score = 0_f64;
    for mut tri in [[0, 1, 2], [0, 2, 3]] {
        if builder.atlas.patch_winding(index)? < 0 {
            tri.swap(1, 2);
        }
        let keys = tri.map(|i| ids[i]);
        if keys[0] == keys[1] || keys[1] == keys[2] || keys[2] == keys[0] {
            continue;
        }
        score = score.max(facet_chord(
            tri.map(|i| literal[i]),
            center,
            bound.normal_radius_chord,
            builder.atlas.patches[index].law,
        )?);
    }
    Ok(score)
}
fn normal_direction(builder: &Builder<'_, '_>, index: usize, uv: [[f64; 2]; 2]) -> usize {
    // These disposable prospective cells use the identical semantic corner
    // registry and proof. Their vertices never enter the published mesh.
    let mut preview = Builder {
        atlas: builder.atlas,
        complex: builder.complex.clone(),
        points: HashMap::new(),
        mesh: DisplayMesh {
            vertices: Vec::new(),
            indices: Vec::new(),
            face_indices: Vec::new(),
        },
    };
    let mut scores = [0_f64; 2];
    for (axis, score) in scores.iter_mut().enumerate() {
        let [a, b] = uv[axis];
        let middle = a + (b - a) / 2.;
        for interval in [[a, middle], [middle, b]] {
            let mut child = uv;
            child[axis] = interval;
            *score =
                score.max(quad_normal_score(&mut preview, index, child).unwrap_or(f64::INFINITY));
        }
    }
    usize::from(scores[1] < scores[0])
}
fn budget(counts: [usize; 8], patches: &[SurfacePatch], limit: usize) -> Result<(), Error> {
    let cells = patches
        .iter()
        .try_fold(0usize, |sum, p| {
            let [a, b] = groups(p.law);
            counts[a]
                .checked_mul(counts[b])
                .and_then(|n| n.checked_mul(2))
                .and_then(|n| sum.checked_add(n))
        })
        .ok_or(Error::Limit)?;
    if cells > limit {
        return Err(Error::Limit);
    }
    Ok(())
}
impl PatchAtlas<'_> {
    /// Produce a conformal approximation without welding coordinates or repairing
    /// winding. The result does not assert a globally embedded manufacturing solid.
    pub fn boundary_mesh(&self, options: Tessellation) -> Result<BoundaryMesh, Error> {
        self.boundary_mesh_with_cancel(options, &|| false)
    }
    pub fn boundary_mesh_with_cancel(
        &self,
        options: Tessellation,
        cancel: &dyn Fn() -> bool,
    ) -> Result<BoundaryMesh, Error> {
        self.refine_from_grid(options, cancel, [1; 8])
    }
    fn refine_from_grid(
        &self,
        options: Tessellation,
        cancel: &dyn Fn() -> bool,
        mut counts: [usize; 8],
    ) -> Result<BoundaryMesh, Error> {
        options.check().map_err(mesh_error)?;
        if counts.iter().any(|count| !count.is_power_of_two()) {
            return Err(Error::Invalid);
        }
        let allowed = (I::point(options.angle) / I::point(2.)).sin() * I::point(2.);
        loop {
            checkpoint(cancel)?;
            budget(counts, &self.patches, options.triangle_limit)?;
            let mut builder = Builder::new(self, cancel)?;
            let mut refine = [false; 8];

            let mut distance = 0_f64;
            let mut normals = 0_f64;
            for (index, patch) in self.patches.iter().enumerate() {
                checkpoint(cancel)?;
                let [gu, gv] = groups(patch.law);
                let [nu, nv] = [counts[gu], counts[gv]];
                for u in 0..nu {
                    for v in 0..nv {
                        checkpoint(cancel)?;
                        if refine[gu] || refine[gv] {
                            continue;
                        }
                        let uv = [
                            [u as f64 / nu as f64, (u + 1) as f64 / nu as f64],
                            [v as f64 / nv as f64, (v + 1) as f64 / nv as f64],
                        ];
                        let pairs = [
                            [uv[0][0], uv[1][0]],
                            [uv[0][1], uv[1][0]],
                            [uv[0][1], uv[1][1]],
                            [uv[0][0], uv[1][1]],
                        ];
                        let ids = pairs.map(|p| builder.point(index, p));
                        let [a, b, c, d] = ids;
                        let ids = [a?, b?, c?, d?];
                        let literal = builder.literal(ids);
                        let bound = self.certify_quad(index, uv, literal)?;
                        let center = self
                            .state(index, uv.map(|[a, b]| I::point(a + (b - a) / 2.)))?
                            .1;
                        let mut valid = bound.positional_error_mm() <= options.deflection;
                        for mut tri in [[0, 1, 2], [0, 2, 3]] {
                            if self.patch_winding(index)? < 0 {
                                tri.swap(1, 2);
                            }
                            let indices = tri.map(|i| ids[i]);
                            if indices[0] == indices[1]
                                || indices[1] == indices[2]
                                || indices[2] == indices[0]
                            {
                                continue;
                            }
                            let chord = facet_chord(
                                tri.map(|i| literal[i]),
                                center,
                                bound.normal_radius_chord,
                                patch.law,
                            );
                            let Some(chord) =
                                resolve_facet(chord, patch.law, &mut refine, [gu, gv])
                            else {
                                continue;
                            };
                            normals = normals.max(chord);
                            valid &= chord <= allowed.lo;
                            builder.mesh.triangle(indices, patch.face as u32);
                        }
                        distance = distance.max(bound.positional_error_mm());
                        if !valid {
                            // Shared source grids refine globally; no hanging edge vertices.
                            // Prefer the coarser parameter direction, retaining full mixed error.
                            if bound.positional_error_mm() <= options.deflection {
                                let axis = normal_direction(&builder, index, uv);
                                refine[[gu, gv][axis]] = true;
                            } else if nu <= nv {
                                refine[gu] = true;
                            } else {
                                refine[gv] = true;
                            }
                        }
                    }
                }
            }

            if !refine.into_iter().any(|v| v) {
                let connectivity = manufacturing::certify(&builder.mesh).map_err(mesh_error)?;
                return Ok(BoundaryMesh {
                    mesh: builder.mesh,
                    connectivity,
                    source_distance_mm: distance,
                    normal_chord: normals,
                    refinement: builder.complex.certificate,
                });
            }
            for (count, refine) in counts.iter_mut().zip(refine) {
                if refine {
                    *count = count.checked_mul(2).ok_or(Error::Limit)?;
                }
            }
        }
    }
}
#[cfg(test)]
mod tests;
