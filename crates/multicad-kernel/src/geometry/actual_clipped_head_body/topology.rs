//! Original physical IDs and finite incidence; no derived verification face.
use super::*;
use multicad_core::references::{
    AnalyticGeometry, AnchoredFrame, TopologyElement, TopologyKey, TopologyKind,
};
use std::collections::BTreeSet;
impl ActualClippedHeadBody {
    fn origin(&self, kind: &str, id: u32) -> Origin {
        Origin {
            source: self.feature,
            role: format!("step.actual-clipped-head/{kind}/#{id}"),
        }
    }
    fn vertex_ids(&self) -> Vec<u32> {
        let local = self.source.local_geometry();
        local
            .auxiliary_thread()
            .vertices()
            .iter()
            .map(|v| v.entity())
            .chain(local.physical_roles().head_vertex_entities)
            .collect()
    }
    fn edge_ids(&self) -> Vec<u32> {
        let local = self.source.local_geometry();
        local
            .auxiliary_thread()
            .edges()
            .iter()
            .map(|v| v.entity())
            .chain(local.physical_roles().head_edge_entities)
            .collect()
    }
    fn edge_vertices(&self, index: usize) -> Option<[usize; 2]> {
        let local = self.source.local_geometry();
        let aux = local.auxiliary_thread();
        if index < aux.edges().len() {
            Some(aux.edges()[index].vertex_indices())
        } else {
            local
                .head()
                .edges
                .get(index - aux.edges().len())
                .map(|v| v.map(|i| i + aux.vertices().len()))
        }
    }
    fn face_loops(&self, index: usize) -> Result<Vec<Vec<(usize, bool)>>, GeometryError> {
        let local = self.source.local_geometry();
        let aux = local.auxiliary_thread();
        if index >= 17 {
            return Err(failure(read::Error::Invalid));
        }
        if index < 9 {
            let id = local.physical_roles().thread_face_entities[index];
            let f = aux
                .faces()
                .iter()
                .find(|f| f.entity() == id)
                .ok_or_else(|| failure(read::Error::Invalid))?;
            return Ok(vec![
                f.coedges()
                    .iter()
                    .map(|c| (c.edge_index(), c.forward()))
                    .collect(),
            ]);
        }
        let f = &local.head().faces[index - 9];
        let mut loops = vec![
            f.outer
                .iter()
                .map(|&(e, d)| (aux.edges().len() + e, d))
                .collect(),
        ];
        if !f.inner.is_empty() {
            let inner = f
                .inner
                .iter()
                .map(|&(id, d)| {
                    aux.edges()
                        .iter()
                        .position(|e| e.entity() == id)
                        .map(|e| (e, d))
                        .ok_or_else(|| failure(read::Error::Invalid))
                })
                .collect::<Result<Vec<_>, _>>()?;
            loops.push(inner);
        }
        Ok(loops)
    }
    pub(crate) fn topology(
        &self,
        owner: EntityId,
        revision: u64,
    ) -> Result<TopologySnapshot, GeometryError> {
        let local = self.source.local_geometry();
        let roles = local.physical_roles();
        let vertices = self.vertex_ids();
        let edges = self.edge_ids();
        let mut elements = Vec::with_capacity(17 + edges.len() + vertices.len());
        for index in 0..17 {
            let id = if index < 9 {
                roles.thread_face_entities[index]
            } else {
                roles.head_face_entities[index - 9]
            };
            let origin = self.origin("face", id);
            let plane = self.plane(index as u32)?;
            let mut boundary = BTreeSet::new();
            for ring in self.face_loops(index)? {
                for (edge, _) in ring {
                    boundary.insert(self.origin("edge", edges[edge]));
                    for v in self
                        .edge_vertices(edge)
                        .ok_or_else(|| failure(read::Error::Invalid))?
                    {
                        boundary.insert(self.origin("vertex", vertices[v]));
                    }
                }
            }
            elements.push(TopologyElement {
                key: TopologyKey {
                    kind: TopologyKind::Face,
                    ordinal: index as u32,
                },
                origins: BTreeSet::from([origin.clone()]),
                boundary,
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
        for (index, &id) in edges.iter().enumerate() {
            let v = self
                .edge_vertices(index)
                .ok_or_else(|| failure(read::Error::Invalid))?;
            let geometry =
                (index >= local.auxiliary_thread().edges().len()).then(|| AnalyticGeometry::Line {
                    a_mm: self.vertices[v[0]],
                    b_mm: self.vertices[v[1]],
                });
            elements.push(TopologyElement {
                key: TopologyKey {
                    kind: TopologyKind::Edge,
                    ordinal: index as u32,
                },
                origins: BTreeSet::from([self.origin("edge", id)]),
                boundary: v
                    .map(|i| self.origin("vertex", vertices[i]))
                    .into_iter()
                    .collect(),
                frames: vec![],
                geometry,
            });
        }
        for (index, &id) in vertices.iter().enumerate() {
            elements.push(TopologyElement {
                key: TopologyKey {
                    kind: TopologyKind::Vertex,
                    ordinal: index as u32,
                },
                origins: BTreeSet::from([self.origin("vertex", id)]),
                boundary: BTreeSet::new(),
                frames: vec![],
                geometry: Some(AnalyticGeometry::Point {
                    position_mm: self.vertices[index],
                }),
            });
        }
        Ok(TopologySnapshot {
            owner,
            revision,
            elements,
        })
    }
    pub(crate) fn planar_loops(&self, index: u32) -> Vec<Vec<(u32, bool)>> {
        if !matches!(self.plane(index), Ok(Some(_))) {
            return vec![];
        }
        self.face_loops(index as usize)
            .unwrap_or_default()
            .into_iter()
            .map(|v| v.into_iter().map(|(i, d)| (i as u32, !d)).collect())
            .collect()
    }
    pub(in crate::geometry) fn sample_edge(
        &self,
        index: u32,
        q: Tessellation,
    ) -> Option<Result<Vec<[f64; 3]>, MeshError>> {
        let index = index as usize;
        self.edge_vertices(index)?;
        Some((|| {
            if !q.deflection.is_finite()
                || q.deflection <= 0.
                || !(0.01..=1.).contains(&q.angle)
                || q.triangle_limit == 0
            {
                return Err(MeshError::InvalidOptions);
            }
            let local = self.source.local_geometry();
            let aux = local.auxiliary_thread();
            if index >= aux.edges().len() {
                let v = self.edge_vertices(index).ok_or(MeshError::InvalidOptions)?;
                return Ok(vec![self.vertices[v[0]], self.vertices[v[1]]]);
            }
            let count = q.triangle_limit.min(131_072);
            let max_output = count.checked_mul(24).ok_or(MeshError::BudgetExceeded)?;
            let external = self
                .source
                .receipt()
                .retained_bytes
                .checked_sub(aux.retained_budget_bytes())
                .and_then(|v| v.checked_add(self.external_bytes()))
                .and_then(|v| {
                    v.checked_add(
                        self.occurrence
                            .as_ref()
                            .map_or(0, |a| a.retained_budget_bytes()),
                    )
                })
                .and_then(|v| v.checked_add(max_output))
                .ok_or(MeshError::BudgetExceeded)?;
            let available = BYTES
                .checked_sub(external)
                .ok_or(MeshError::BudgetExceeded)?;
            let scale = local.axial_sewing_transfer().positive_axial_scale_bounds()[1].max(1.);
            let sigma = (I::point(scale)
                * I::point(self.source.actual_frame().scale_bounds()[1])
                * I::point(
                    self.occurrence
                        .as_ref()
                        .map_or(1., |a| a.linear_l2_norm_upper()),
                ))
            .hi;
            let deflection = (I::point(q.deflection) / I::point(sigma)).lo;
            let sampled = aux
                .sample_edge(
                    index,
                    deflection,
                    count,
                    read::Options {
                        maximum_bytes: available,
                        ..Default::default()
                    },
                    None,
                )
                .map_err(|e| match e {
                    read::Error::Invalid => MeshError::InvalidOptions,
                    read::Error::Limit => MeshError::BudgetExceeded,
                    _ => MeshError::UnresolvedBoundary,
                })?;
            if sampled
                .points_mm()
                .len()
                .checked_mul(256)
                .is_none_or(|v| v > sampled.remaining_work())
            {
                return Err(MeshError::BudgetExceeded);
            }
            let mut output = Vec::new();
            output
                .try_reserve_exact(sampled.points_mm().len())
                .map_err(|_| MeshError::BudgetExceeded)?;
            if output.capacity() * 24 > max_output {
                return Err(MeshError::BudgetExceeded);
            }
            let mut rounding = 0f64;
            for &p in sampled.points_mm() {
                let bound = local
                    .axial_sewing_transfer()
                    .map_raw_box_mm(p.map(|v| [v, v]))
                    .map_err(|_| MeshError::UnresolvedBoundary)?;
                let bound = self
                    .map_local_box(bound)
                    .map_err(|_| MeshError::UnresolvedBoundary)?;
                let literal = bound.map(|v| v[0].midpoint(v[1]));
                let error = ranges(bound)
                    .into_iter()
                    .zip(literal)
                    .fold(I::ZERO, |s, (v, p)| {
                        s + I::point((v - I::point(p)).magnitude())
                    });
                rounding = rounding.max(error.hi);
                output.push(literal);
            }
            let error =
                I::point(sigma) * I::point(sampled.position_error_mm()) + I::point(rounding);
            if !error.valid() || error.hi > q.deflection {
                return Err(MeshError::UnresolvedBoundary);
            }
            Ok(output)
        })())
    }
}
