//! Closed source-owned incidence over the ACTUAL coupled coefficient rows.
//! No float-coordinate adapter; literal embedding/material remain separate.
use super::*;
use crate::step::insert::{
    SurfacePatch, VertexOwner,
    topology::{
        ShellEdge, ShellFace,
        incidence::{self, Input},
    },
};
#[derive(Clone, Copy, Debug)]
pub struct VertexBinding {
    owner: VertexOwner,
    row: usize,
    endpoint: usize,
}
impl VertexBinding {
    pub fn owner(&self) -> VertexOwner {
        self.owner
    }
    pub fn coefficient_reference(&self) -> [usize; 2] {
        [self.row, self.endpoint * 11]
    }
}
struct Metadata<'p, 's> {
    plan: &'p CoupledPlan<'s>,
    vertices: Vec<VertexBinding>,
    ends: Vec<[usize; 2]>,
    patches: Vec<SurfacePatch>,
}
impl Input for Metadata<'_, '_> {
    fn domain(&self) -> &DomainPlan<'_> {
        self.plan.domain()
    }
    fn row_count(&self) -> usize {
        self.plan.rows().len()
    }
    fn row_key(&self, index: usize) -> RowKey {
        self.plan.rows()[index].key()
    }
    fn row_ends(&self, index: usize) -> [usize; 2] {
        self.ends[index]
    }
    fn patch_count(&self) -> usize {
        self.patches.len()
    }
    fn patch(&self, index: usize) -> SurfacePatch {
        self.patches[index]
    }
    fn vertex_count(&self) -> usize {
        self.vertices.len()
    }
    fn vertex_owner(&self, index: usize) -> VertexOwner {
        self.vertices[index].owner
    }
}
pub struct CoupledTopology<'p, 's> {
    metadata: Metadata<'p, 's>,
    built: incidence::Built,
    work: usize,
    bytes: usize,
}
impl<'p, 's> CoupledTopology<'p, 's> {
    pub fn plan(&self) -> &'p CoupledPlan<'s> {
        self.metadata.plan
    }
    pub fn vertices(&self) -> &[VertexBinding] {
        &self.metadata.vertices
    }
    pub fn row_endpoints(&self, index: usize) -> [usize; 2] {
        self.metadata.ends[index]
    }
    pub fn patches(&self) -> &[SurfacePatch] {
        &self.metadata.patches
    }
    pub fn edges(&self) -> &[ShellEdge] {
        &self.built.edges
    }
    pub fn faces(&self) -> &[ShellFace] {
        &self.built.faces
    }
    pub fn used_vertices(&self) -> usize {
        self.built.used_vertices
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_bytes(&self) -> usize {
        self.bytes
    }
    pub fn prepare(
        plan: &'p CoupledPlan<'s>,
        maximum_work: usize,
        maximum_bytes: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, Error> {
        let mut w = Work::new(maximum_work, interrupt)?;
        w.charge(plan.work_used() + 1)?;
        if maximum_bytes == 0 || maximum_bytes > 64 * 1024 * 1024 {
            return Err(Error::Invalid);
        }
        let mut memory = Memory::new(plan.retained_bytes(), maximum_bytes)?;
        memory.add(size_of::<Self>())?;
        let source = plan.source();
        let count = [Side::Outer, Side::Inner]
            .into_iter()
            .map(|side| source.source_boundary(side).points().len())
            .sum::<usize>();
        let mut owners = memory.allocate(count + 2 * plan.rows().len(), &mut w)?;
        for side in [Side::Outer, Side::Inner] {
            for i in 0..source.source_boundary(side).points().len() {
                w.charge(1)?;
                owners.push(VertexOwner::Original(vertex_key(side, i as u32)));
            }
        }
        let mut representatives = memory.allocate(count + 2 * plan.rows().len(), &mut w)?;
        representatives.resize(count, None::<(usize, usize)>);
        let mut ends: Vec<[usize; 2]> = memory.allocate(plan.rows().len(), &mut w)?;
        for (index, row) in plan.rows().iter().enumerate() {
            let mut pair = [0; 2];
            for (i, owner) in row.endpoint_owners().into_iter().enumerate() {
                w.charge(1)?;
                let vertex = match owner {
                    VertexOwner::Original(key) => vertex_index(source, key),
                    VertexOwner::Station { .. } => {
                        if i == 0
                            && index > 0
                            && plan.rows()[index - 1].endpoint_owners()[1] == owner
                        {
                            ends[index - 1][1]
                        } else {
                            owners.push(owner);
                            representatives.push(None);
                            owners.len() - 1
                        }
                    }
                };
                representatives[vertex].get_or_insert((index, i));
                pair[i] = vertex;
            }
            ends.push(pair);
        }
        let mut vertices = memory.allocate(owners.len(), &mut w)?;
        for (owner, representative) in owners.iter().zip(&representatives) {
            w.charge(1)?;
            let (row, endpoint) = representative
                .expect("CoupledPlan certifies complete original and shared station ownership");
            vertices.push(VertexBinding {
                owner: *owner,
                row,
                endpoint,
            });
        }
        let mut patches = memory.allocate(plan.domain().slabs().len(), &mut w)?;
        for (slab_index, slab) in plan.domain().slabs().iter().enumerate() {
            w.charge(2 * (usize::BITS - plan.rows().len().leading_zeros() + 1) as usize)?;
            let rows = slab.rows.map(|curve| {
                plan.rows()
                    .binary_search_by_key(
                        &RowKey {
                            curve,
                            interval: slab.stations[0],
                        },
                        CoupledRow::key,
                    )
                    .expect("registered exact source slab row")
            });
            let left = source.source_boundary(slab.side).source().spec().left_hand;
            patches.push(SurfacePatch {
                slab: slab_index,
                rows,
                same_sense: left != (slab.side == Side::Inner),
            });
        }
        memory.release(&owners);
        memory.release(&representatives);
        drop(owners);
        drop(representatives);
        let metadata = Metadata {
            plan,
            vertices,
            ends,
            patches,
        };
        let built = incidence::prepare(&metadata, &mut memory, &mut w)?;
        Ok(Self {
            metadata,
            built,
            work: w.used,
            bytes: memory.used,
        })
    }
}
#[cfg(test)]
mod tests;
