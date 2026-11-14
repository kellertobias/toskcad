//! Shared source-owned INCIDENCE only; no coordinate or geometry adapter.
use super::*;
pub(in crate::step::insert) trait Input {
    fn domain(&self) -> &DomainPlan<'_>;
    fn row_count(&self) -> usize;
    fn row_key(&self, index: usize) -> RowKey;
    fn row_ends(&self, index: usize) -> [usize; 2];
    fn patch_count(&self) -> usize;
    fn patch(&self, index: usize) -> SurfacePatch;
    fn vertex_count(&self) -> usize;
    fn vertex_owner(&self, index: usize) -> VertexOwner;
}
impl Input for SurfacePlan<'_, '_> {
    fn domain(&self) -> &DomainPlan<'_> {
        SurfacePlan::domain(self)
    }
    fn row_count(&self) -> usize {
        self.rows().len()
    }
    fn row_key(&self, index: usize) -> RowKey {
        self.rows()[index].key()
    }
    fn row_ends(&self, index: usize) -> [usize; 2] {
        self.rows()[index].endpoints()
    }
    fn patch_count(&self) -> usize {
        self.patches().len()
    }
    fn patch(&self, index: usize) -> SurfacePatch {
        self.patches()[index]
    }
    fn vertex_count(&self) -> usize {
        self.vertices().len()
    }
    fn vertex_owner(&self, index: usize) -> VertexOwner {
        self.vertices()[index].owner()
    }
}
pub(in crate::step::insert) struct Built {
    pub edges: Vec<ShellEdge>,
    pub faces: Vec<ShellFace>,
    pub used_vertices: usize,
}
pub(in crate::step::insert) fn prepare(
    surface: &dyn Input,
    memory: &mut Memory,
    w: &mut Work<'_>,
) -> Result<Built, Error> {
    let mut edges = memory.allocate(surface.row_count() + 2 * surface.patch_count(), w)?;
    for i in 0..surface.row_count() {
        w.charge(1)?;
        edges.push(ShellEdge {
            ends: surface.row_ends(i),
            geometry: EdgeGeometry::Row(i),
        });
    }
    let mut ruling_keys = memory.allocate(surface.patch_count() * 2, w)?;
    let mut faces = memory.allocate(surface.patch_count() + 2, w)?;
    for patch_index in 0..surface.patch_count() {
        let patch = surface.patch(patch_index);
        let slab = &surface.domain().slabs()[patch.slab];
        let [a, b] = patch.rows.map(|i| surface.row_ends(i));
        let mut wire = memory.allocate(4, w)?;
        if a[0] != b[0] {
            wire.push(ruling(
                surface,
                slab,
                0,
                [a[0], b[0]],
                &mut ruling_keys,
                &mut edges,
                w,
            )?);
        }
        wire.push(ShellCoedge {
            edge: patch.rows[1],
            forward: true,
        });
        if a[1] != b[1] {
            wire.push(ruling(
                surface,
                slab,
                1,
                [b[1], a[1]],
                &mut ruling_keys,
                &mut edges,
                w,
            )?);
        }
        wire.push(ShellCoedge {
            edge: patch.rows[0],
            forward: false,
        });
        if !patch.same_sense {
            wire.reverse();
            for c in &mut wire {
                c.forward = !c.forward;
            }
        }
        let mut loops = memory.allocate(1, w)?;
        loops.push(wire);
        faces.push(ShellFace {
            geometry: FaceGeometry::Tensor(patch_index),
            same_sense: patch.same_sense,
            loops,
        });
    }
    for (end, key) in [(0, FaceKey::StartCap), (1, FaceKey::EndCap)] {
        let source_loops = surface
            .domain()
            .boundary()
            .face_loops(key)
            .expect("actual paired annular cap");
        let mut loops = memory.allocate(source_loops.len(), w)?;
        for source_loop in source_loops {
            // Each source curve contributes only its retained shared rows;
            // no disk scaffold or circular-envelope substitute is emitted.
            w.charge(source_loop.len() * surface.row_count())?;
            let count = source_loop
                .iter()
                .map(|c| {
                    (0..surface.row_count())
                        .filter(|&i| surface.row_key(i).curve == c.curve)
                        .count()
                })
                .sum();
            let mut wire = memory.allocate(count, w)?;
            for c in source_loop {
                let begin = bound(surface, c.curve, false);
                let finish = bound(surface, c.curve, true);
                w.charge(2 * (usize::BITS - surface.row_count().leading_zeros() + 1) as usize)?;
                let original = surface
                    .domain()
                    .boundary()
                    .curve_endpoints(c.curve)
                    .expect("actual cap curve endpoints");
                let first = surface.vertex_owner(surface.row_ends(begin)[0]);
                let last = surface.vertex_owner(surface.row_ends(finish - 1)[1]);
                let forward = cap_forward([first, last], original, c.forward)?;
                for j in 0..(finish - begin) {
                    w.charge(1)?;
                    wire.push(ShellCoedge {
                        edge: begin + if forward { j } else { (finish - begin) - 1 - j },
                        forward,
                    });
                }
            }
            loops.push(wire);
        }
        faces.push(ShellFace {
            geometry: FaceGeometry::Cap(end),
            same_sense: end == 1,
            loops,
        });
    }
    memory.release(&ruling_keys);
    drop(ruling_keys);
    let used_vertices = verify(surface.vertex_count(), &edges, &faces, memory, w)?;
    w.charge(1)?;
    Ok(Built {
        edges,
        faces,
        used_vertices,
    })
}
fn bound(input: &dyn Input, key: BoundaryCurveKey, upper: bool) -> usize {
    let mut lo = 0;
    let mut hi = input.row_count();
    while lo < hi {
        let mid = lo + (hi - lo) / 2;
        let value = input.row_key(mid).curve;
        if value < key || (upper && value == key) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    lo
}
