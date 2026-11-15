//! Immutable actual literal geometry. No parser borrow or authored recipe.
use super::*;
mod measure;
pub use measure::{FaceMeasure, FaceMeasureOptions};
mod mesh;
mod query;
pub use mesh::{
    LiteralEmbeddingProof, LiteralGeometryOptions, LiteralGeometryPlan, LiteralMesh,
    LiteralMeshCandidate, LiteralMeshOptions,
};
pub use query::{Location, MaterialQuery, QueryOptions};
use std::{ops::Range, sync::Arc};
#[derive(Clone)]
pub struct OwnedMaterial(Arc<Source>);
struct StoredScalar {
    value: f64,
    bounds: [f64; 2],
    spelling: Range<usize>,
}
struct StoredPoint {
    id: u32,
    coordinates: [StoredScalar; 3],
}
struct StoredRecord {
    id: u32,
    components: Range<usize>,
}
struct Source {
    literal: String,
    records: Vec<StoredRecord>,
    components: Vec<Range<usize>>,
    points: Vec<StoredPoint>,
    vertices: Vec<Vertex>,
    curves: Vec<Curve>,
    tensors: Vec<Tensor>,
    edges: Vec<Edge>,
    faces: Vec<Face>,
    member: Member,
    linear: [[I; 3]; 3],
    inverse: [[I; 3]; 3],
    origin: [I; 3],
    inverse_origin: [I; 3],
    determinant: I,
    local_bounds: [[f64; 2]; 3],
    world_bounds: [[f64; 2]; 3],
    local_volume: [f64; 2],
    world_volume: [f64; 2],
    local_areas: [[f64; 2]; 2],
    cap_heights: [[f64; 2]; 2],
    radial_squared: [f64; 2],
    work: usize,
    retained: usize,
    peak: usize,
    document_bytes: usize,
}
fn bounds(i: I) -> [f64; 2] {
    [i.lo, i.hi]
}
fn interval([lo, hi]: [f64; 2]) -> I {
    I::new(lo, hi)
}
fn hold(budget: &mut Budget<'_>, bytes: usize) -> Result<()> {
    budget.step(1)?;
    if bytes > budget.limit - budget.bytes {
        return Err(Error::Limit);
    }
    budget.bytes += bytes;
    budget.peak = budget.peak.max(budget.bytes);
    Ok(())
}
fn copy<T: Copy>(values: &[T], budget: &mut Budget<'_>) -> Result<Vec<T>> {
    let mut out = budget.allocate(values.len())?;
    out.extend_from_slice(values);
    Ok(out)
}
impl OwnedMaterial {
    pub fn member(&self) -> Member {
        self.0.member
    }
    pub fn point(&self, index: usize) -> Option<Point<'_>> {
        let point = self.0.points.get(index)?;
        Some(Point {
            id: point.id,
            coordinates: std::array::from_fn(|axis| {
                let scalar = &point.coordinates[axis];
                numeric::retained_scalar(
                    scalar.value,
                    scalar.bounds,
                    &self.0.literal[scalar.spelling.clone()],
                )
            }),
        })
    }
    pub fn point_count(&self) -> usize {
        self.0.points.len()
    }
    pub fn vertices(&self) -> &[Vertex] {
        &self.0.vertices
    }
    pub fn curves(&self) -> &[Curve] {
        &self.0.curves
    }
    pub fn tensors(&self) -> &[Tensor] {
        &self.0.tensors
    }
    pub fn edges(&self) -> &[Edge] {
        &self.0.edges
    }
    pub fn faces(&self) -> &[Face] {
        &self.0.faces
    }
    pub fn source_entity_count(&self) -> usize {
        self.0.records.len()
    }
    pub fn original_component(&self, entity: u32, component: usize) -> Option<&str> {
        let record = &self.0.records[self
            .0
            .records
            .binary_search_by_key(&entity, |r| r.id)
            .ok()?];
        if component >= record.components.len() {
            return None;
        }
        Some(&self.0.literal[self.0.components[record.components.start + component].clone()])
    }
    pub fn local_enclosure_mm(&self) -> [[f64; 2]; 3] {
        self.0.local_bounds
    }
    pub fn world_enclosure_mm(&self) -> [[f64; 2]; 3] {
        self.0.world_bounds
    }
    pub fn local_volume_mm3(&self) -> [f64; 2] {
        self.0.local_volume
    }
    pub fn world_volume_mm3(&self) -> [f64; 2] {
        self.0.world_volume
    }
    pub fn local_cap_areas_mm2(&self) -> [[f64; 2]; 2] {
        self.0.local_areas
    }
    pub fn determinant_bounds(&self) -> [f64; 2] {
        bounds(self.0.determinant)
    }
    pub fn linear_columns_bounds(&self) -> [[[f64; 2]; 3]; 3] {
        self.0.linear.map(|v| v.map(bounds))
    }
    pub fn inverse_columns_bounds(&self) -> [[[f64; 2]; 3]; 3] {
        self.0.inverse.map(|v| v.map(bounds))
    }
    pub fn origin_mm(&self) -> [[f64; 2]; 3] {
        self.0.origin.map(bounds)
    }
    pub fn inverse_origin_mm(&self) -> [[f64; 2]; 3] {
        self.0.inverse_origin.map(bounds)
    }
    pub fn map_box_mm(&self, input: [[f64; 2]; 3]) -> Result<[[f64; 2]; 3]> {
        placement::transport(self.0.linear, self.0.origin, input)
    }
    pub fn inverse_box_mm(&self, input: [[f64; 2]; 3]) -> Result<[[f64; 2]; 3]> {
        placement::transport(self.0.inverse, self.0.inverse_origin, input)
    }
    pub fn work_used(&self) -> usize {
        self.0.work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.0.retained
    }
    /// Borrowed neutral arena plus construction query peak were simultaneously
    /// live. This is separate from the final immutable retained allocation.
    pub fn construction_query_peak_bytes(&self) -> usize {
        self.0.peak
    }
    pub fn original_document_budget_bytes(&self) -> usize {
        self.0.document_bytes
    }
    pub fn cell(
        &self,
        tensor: usize,
        index: [usize; 2],
    ) -> Result<crate::tapered_thread::TensorCell11> {
        let tensor = self.0.tensors.get(tensor).ok_or(Error::Invalid)?;
        if index.into_iter().zip(tensor.spans).any(|(i, n)| i >= n) {
            return Err(Error::Invalid);
        }
        let height = 11 * tensor.spans[1] + 1;
        let controls = std::array::from_fn(|i| {
            let point = &self.0.points
                [tensor.poles[(11 * index[0] + i / 12) * height + 11 * index[1] + i % 12]];
            std::array::from_fn(|axis| interval(point.coordinates[axis].bounds))
        });
        Ok(crate::tapered_thread::calculation_cell11(controls)
            .expect("retained validated literal cell"))
    }
}
impl Selected<'_, '_> {
    pub fn retain(
        &self,
        options: MaterialOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<OwnedMaterial> {
        if options.maximum_work == 0
            || options.maximum_work > 64_000_000
            || options.maximum_bytes == 0
            || options.maximum_bytes > 64 * 1024 * 1024
            || options.maximum_nodes > 65_536
            || options.maximum_depth > 32
        {
            return Err(Error::Invalid);
        }
        let pose = self.placement(
            Options {
                maximum_work: options.maximum_work.min(16_000_000),
                maximum_bytes: options.maximum_bytes,
            },
            interrupt,
        )?;
        let pose_work = pose.work_used() - self.parsed.work;
        let maximum_work = options
            .maximum_work
            .checked_sub(pose_work)
            .filter(|n| *n > 0)
            .ok_or(Error::Limit)?;
        let maximum_bytes = options
            .maximum_bytes
            .checked_sub(size_of::<Placement<'_, '_, '_>>())
            .filter(|n| *n > 0)
            .ok_or(Error::Limit)?;
        let proof_options = MaterialOptions {
            maximum_work,
            maximum_bytes,
            ..options
        };
        let occupied = self
            .parsed
            .projection(proof_options, interrupt)?
            .occupied(proof_options, interrupt)?;
        let used = occupied.work_used() + pose_work;
        let bytes = occupied.query_budget_bytes() + size_of::<Placement<'_, '_, '_>>();
        let mut budget = Budget {
            used,
            remaining: options.maximum_work.checked_sub(used).ok_or(Error::Limit)?,
            bytes,
            peak: bytes.max(pose.query_budget_bytes()),
            limit: options.maximum_bytes,
            interrupt,
        };
        // Graph work is included exactly once in occupied; pose contributes
        // only its additional query work. Neutral document is a separate scope.
        let fixed = size_of::<OwnedMaterial>()
            + size_of::<Source>()
            + 2 * size_of::<usize>()
            + 2 * std::mem::align_of::<usize>();
        hold(&mut budget, fixed)?;
        let arena = self.parsed.document.records();
        let mut total = 0usize;
        let mut count = 0usize;
        for id in self.entity_ids() {
            budget.step(32)?;
            for entity in &arena[id] {
                let source = entity.source.expect("neutral original component");
                budget.step(8)?;
                total = total.checked_add(source.len()).ok_or(Error::Limit)?;
                count += 1;
            }
        }
        let mut literal = budget.allocate::<u8>(total)?;
        let mut components = budget.allocate(count)?;
        let mut records = budget.allocate(self.entity_ids().len())?;
        for id in self.entity_ids() {
            budget.step(32)?;
            let first = components.len();
            for entity in &arena[id] {
                let source = entity.source.expect("neutral original component");
                budget.step(source.len() + 16)?;
                let start = literal.len();
                literal.extend_from_slice(source.as_bytes());
                components.push(start..literal.len());
            }
            records.push(StoredRecord {
                id: *id,
                components: first..components.len(),
            });
        }
        budget.step(literal.len() + 1)?;
        let literal = String::from_utf8(literal).expect("concatenated original UTF8 components");
        let mut points = budget.allocate(self.parsed.points.len())?;
        for point in &self.parsed.points {
            // Original-arena and owned-record lookups plus three checked
            // scalar spelling offsets, cache copies and publication.
            budget.step(128)?;
            let record = &records[records
                .binary_search_by_key(&point.id, |r| r.id)
                .expect("selected point owner")];
            let span = &components[record.components.start];
            let original = arena[&point.id][0]
                .source
                .expect("validated point component");
            let coordinates = point.coordinates.map(|scalar| {
                let offset = scalar.spelling().as_ptr() as usize - original.as_ptr() as usize;
                assert!(
                    offset + scalar.spelling().len() <= original.len(),
                    "scalar belongs to its original component"
                );
                StoredScalar {
                    value: scalar.value(),
                    bounds: scalar.bounds(),
                    spelling: span.start + offset..span.start + offset + scalar.spelling().len(),
                }
            });
            points.push(StoredPoint {
                id: point.id,
                coordinates,
            });
        }
        let vertices = copy(&self.parsed.vertices, &mut budget)?;
        let edges = copy(&self.parsed.edges, &mut budget)?;
        let mut faces = budget.allocate(self.parsed.faces.len())?;
        for face in &self.parsed.faces {
            budget.step(16)?;
            faces.push(Face {
                id: face.id,
                support: face.support,
                coedges: face.coedges,
                tensor_sides: face.tensor_sides,
                forward: face.forward,
            });
        }
        let mut curves = budget.allocate(self.parsed.curves.len())?;
        for curve in &self.parsed.curves {
            budget.step(16)?;
            curves.push(Curve {
                id: curve.id,
                spans: curve.spans,
                poles: copy(&curve.poles, &mut budget)?,
            });
        }
        let mut tensors = budget.allocate(self.parsed.tensors.len())?;
        for tensor in &self.parsed.tensors {
            budget.step(16)?;
            tensors.push(Tensor {
                id: tensor.id,
                spans: tensor.spans,
                poles: copy(&tensor.poles, &mut budget)?,
            });
        }
        let mut local = [I::new(f64::INFINITY, f64::NEG_INFINITY); 3];
        for tensor in &tensors {
            for point in &tensor.poles {
                budget.step(6)?;
                for (axis, range) in local.iter_mut().enumerate() {
                    range.lo = range.lo.min(points[*point].coordinates[axis].bounds[0]);
                    range.hi = range.hi.max(points[*point].coordinates[axis].bounds[1]);
                }
            }
        }
        // Fixed matrix/cache transport, literal mass scaling and the bounded
        // descriptor/capacity publication ledger below (eight nets,20curves).
        budget.step(2048)?;
        let linear = pose.linear_columns_bounds().map(|v| v.map(interval));
        let inverse = pose.inverse_columns_bounds().map(|v| v.map(interval));
        let origin = pose.origin_mm().map(interval);
        let inverse_origin = pose.inverse_origin_mm().map(interval);
        let determinant = interval(pose.determinant_bounds());
        let world_bounds = placement::transport(linear, origin, local.map(bounds))?;
        let world_volume = interval(occupied.volume_mm3()) * determinant;
        if !world_volume.valid() || world_volume.lo <= 0. {
            return Err(Error::Unsupported);
        }
        budget.step(1)?;
        let retained = fixed
            + literal.capacity()
            + records.capacity() * size_of::<StoredRecord>()
            + components.capacity() * size_of::<Range<usize>>()
            + points.capacity() * size_of::<StoredPoint>()
            + vertices.capacity() * size_of::<Vertex>()
            + edges.capacity() * size_of::<Edge>()
            + faces.capacity() * size_of::<Face>()
            + curves.capacity() * size_of::<Curve>()
            + tensors.capacity() * size_of::<Tensor>()
            + curves
                .iter()
                .map(|c| c.poles.capacity() * size_of::<usize>())
                .sum::<usize>()
            + tensors
                .iter()
                .map(|t| t.poles.capacity() * size_of::<usize>())
                .sum::<usize>();
        Ok(OwnedMaterial(Arc::new(Source {
            literal,
            records,
            components,
            points,
            vertices,
            curves,
            tensors,
            edges,
            faces,
            member: self.member,
            linear,
            inverse,
            origin,
            inverse_origin,
            determinant,
            local_bounds: local.map(bounds),
            world_bounds,
            local_volume: occupied.volume_mm3(),
            world_volume: bounds(world_volume),
            local_areas: occupied.cap_areas_mm2(),
            cap_heights: occupied.projection().cap_height_bounds_mm(),
            radial_squared: occupied.projection().lateral_radius_squared_bounds(),
            work: budget.used,
            retained,
            peak: budget.peak,
            document_bytes: self.parsed.document.retained_budget_bytes(),
        })))
    }
}
#[cfg(test)]
mod tests;
