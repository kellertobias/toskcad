//! Accounted retained source bytes, not RSS or allocator fragmentation.
//! No geometry is changed, sampled or reconstructed by this metadata query.
use super::*;
fn payload<T>(values: &Vec<T>) -> usize {
    values.capacity().saturating_mul(std::mem::size_of::<T>())
}
struct PoleOwners {
    keys: [(*const Point, usize); 102],
    count: usize,
}
impl PoleOwners {
    fn new() -> Self {
        Self {
            keys: [(std::ptr::null(), 0); 102],
            count: 0,
        }
    }
    fn allocation(&mut self, tensor: &RetainedTensor) -> usize {
        let poles = tensor.poles();
        let key = (poles.as_ptr(), poles.len());
        if self.keys[..self.count].contains(&key) {
            return 0;
        }
        // BoundaryOwnership admits <=34 surfaces/faces. Only input surfaces,
        // BoundTensor and ChartTensor retain Arc pole copies: <=3*34 owners.
        // All three families are traversed exactly once below, without cloning.
        self.keys[self.count] = key;
        self.count += 1;
        // Arc<[Point]> has its exact slice length (no hidden Vec capacity).
        // Two counters plus alignment slack conservatively cover its header.
        poles
            .len()
            .saturating_mul(std::mem::size_of::<Point>())
            .saturating_add(2 * std::mem::size_of::<usize>() + 2 * std::mem::align_of::<Point>())
    }
}
impl RadialBoundary {
    /// Conservative accounted retained bytes of this immutable certified source:
    /// all actual Vec capacities/nested headers/projected controls, plus each
    /// shared pole allocation once. This excludes query scratch, unrelated
    /// owners, allocator fragmentation and process RSS. Work is bounded metadata
    /// only (<=34 faces/surfaces,112 edges/curves,80 vertices and102 Arc refs);
    /// no control-coordinate scan, heap allocation or source cache is involved.
    pub fn retained_budget_bytes(&self) -> usize {
        let caps = &self.ordered.caps;
        let axial = &caps.axial;
        let trimmed = &axial.boundary;
        let owned = &trimmed.ownership;
        let input = &owned.input;
        let mut owners = PoleOwners::new();
        // The nested fixed objects and every Vec/Arc descriptor header are
        // already contained in Self or in their enclosing dynamic elements.
        let mut bytes = std::mem::size_of::<Self>();
        let mut add = |n: usize| bytes = bytes.saturating_add(n);
        add(payload(&input.vertices));
        add(payload(&input.curves));
        for curve in &input.curves {
            if let Curve::Spline { poles, .. } = curve {
                add(payload(poles));
            }
        }
        add(payload(&input.edges));
        add(payload(&input.surfaces));
        for surface in &input.surfaces {
            if let Surface::Tensor(t) = surface {
                add(owners.allocation(t));
            }
        }
        add(payload(&input.faces));
        for face in &input.faces {
            add(payload(&face.outer));
        }
        add(payload(&owned.topology.edge_faces));
        add(payload(&owned.topology.vertex_edges));
        for edges in &owned.topology.vertex_edges {
            add(payload(edges));
        }
        add(payload(&trimmed.tensor_edges));
        add(payload(&trimmed.tensors));
        for t in &trimmed.tensors {
            add(owners.allocation(&t.surface));
        }
        add(payload(&axial.tensors));
        for t in &axial.tensors {
            add(payload(&t.controls));
            add(owners.allocation(&t.surface));
        }
        add(payload(&axial.cones));
        add(payload(&caps.sides));
        add(payload(&caps.kinds));
        add(payload(&caps.rails));
        for (_, rail) in &caps.rails {
            add(payload(rail));
        }
        for rims in &caps.rims {
            add(payload(rims));
            for rim in rims {
                add(payload(rim));
            }
        }
        add(payload(&self.ordered.rails));
        for rail in &self.ordered.rails {
            add(payload(rail));
        }
        bytes
    }
}
#[cfg(test)]
mod tests;
