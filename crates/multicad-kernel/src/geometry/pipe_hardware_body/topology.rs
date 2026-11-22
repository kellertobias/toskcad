//! Original thread identities are retained; only real hex-stock owners are
//! added. No mesh vertices or facets are substituted for finite topology.
use super::*;
fn origin(source: &Source, role: String) -> Origin {
    Origin {
        source: source.feature,
        role: format!("pipe.hardware/{role}"),
    }
}
fn reverse(ring: &[(u32, bool)]) -> Vec<(u32, bool)> {
    ring.iter().rev().map(|(e, f)| (*e, !*f)).collect()
}
pub(super) fn extend(
    source: &mut Source,
    base: u32,
    cap_loops: [Vec<Vec<(u32, bool)>>; 2],
    head: bool,
) -> Result<(), GeometryError> {
    // Both actual authored source builders earn one finite cap cycle. This
    // assertion depends on their typed construction, never on file headers.
    assert!(
        cap_loops.iter().all(|r| r.len() == 1),
        "sealed authored pipe has one cap cycle"
    );
    for i in 0..12 {
        let own = origin(source, format!("hex.vertex.{i}"));
        source.topology.elements.push(TopologyElement {
            key: TopologyKey {
                kind: TopologyKind::Vertex,
                ordinal: source.vertex_offset + i as u32,
            },
            origins: BTreeSet::from([own]),
            boundary: BTreeSet::new(),
            frames: vec![],
            geometry: Some(AnalyticGeometry::Point {
                position_mm: midpoint(source.hex_vertices[i]),
            }),
        });
    }
    for i in 0..18 {
        let (a, b) = if i < 12 {
            let start = if i < 6 { 0 } else { 6 };
            (i, start + (i + 1) % 6)
        } else {
            (i - 12, i - 6)
        };
        let own = origin(source, format!("hex.edge.{i}"));
        let ends = [
            origin(source, format!("hex.vertex.{a}")),
            origin(source, format!("hex.vertex.{b}")),
        ];
        source.topology.elements.push(TopologyElement {
            key: TopologyKey {
                kind: TopologyKind::Edge,
                ordinal: source.edge_offset + i as u32,
            },
            origins: BTreeSet::from([own]),
            boundary: ends.into_iter().collect(),
            frames: vec![],
            geometry: Some(AnalyticGeometry::Line {
                a_mm: midpoint(source.hex_vertices[a]),
                b_mm: midpoint(source.hex_vertices[b]),
            }),
        });
    }
    let e = source.edge_offset;
    for i in 0..6u32 {
        let j = (i + 1) % 6;
        source.loops.insert(
            base + i,
            vec![vec![
                (e + i, true),
                (e + 12 + j, true),
                (e + 6 + i, false),
                (e + 12 + i, false),
            ]],
        );
    }
    let bottom = (0..6).rev().map(|i| (e + i, false)).collect::<Vec<_>>();
    let top = (0..6).map(|i| (e + 6 + i, true)).collect::<Vec<_>>();
    if head {
        source
            .loops
            .insert(base + 6, vec![bottom, cap_loops[1][0].clone()]);
        source.loops.insert(base + 7, vec![top]);
    } else {
        source
            .loops
            .insert(base + 6, vec![bottom, reverse(&cap_loops[0][0])]);
        source
            .loops
            .insert(base + 7, vec![top, reverse(&cap_loops[1][0])]);
    }
    for i in 0..8 {
        let index = base + i;
        let own = origin(source, format!("hex.face.{i}"));
        let mut boundary = BTreeSet::new();
        for &(edge, _) in source.loops[&index].iter().flatten() {
            let element = source
                .topology
                .elements
                .iter()
                .find(|v| v.key.kind == TopologyKind::Edge && v.key.ordinal == edge)
                .expect("same finite source edge owner");
            boundary.extend(element.origins.iter().cloned());
            boundary.extend(element.boundary.iter().cloned());
        }
        let frame = source.faces[index as usize]
            .frame
            .expect("actual added stock plane");
        source.topology.elements.push(TopologyElement {
            key: TopologyKey {
                kind: TopologyKind::Face,
                ordinal: index,
            },
            origins: BTreeSet::from([own.clone()]),
            boundary,
            frames: vec![AnchoredFrame { origin: own, frame }],
            geometry: Some(AnalyticGeometry::Plane { frame }),
        });
    }
    source.topology.elements.sort_by_key(|e| e.key);
    Ok(())
}
impl PipeHardwareBody {
    pub(crate) fn sample_edge(
        &self,
        index: u32,
        options: Tessellation,
    ) -> Option<Result<Vec<[f64; 3]>, multicad_rust_kernel::display::MeshError>> {
        use multicad_rust_kernel::display::MeshError;
        if index >= self.source.edge_offset {
            return None;
        }
        let local = Tessellation {
            deflection: (I::point(options.deflection) * I::point(0.5)
                / I::point(self.placement.sigma[1]))
            .lo,
            triangle_limit: options.triangle_limit.min(131072),
            ..options
        };
        let points = self.source.thread.sample_non_analytic_edge(index, local)?;
        Some(points.and_then(|points| {
            points
                .into_iter()
                .map(|p| {
                    let (point, error) = self
                        .placement
                        .point(p)
                        .map_err(|_| MeshError::UnresolvedBoundary)?;
                    if error > options.deflection / 4. {
                        return Err(MeshError::UnresolvedBoundary);
                    }
                    Ok(point)
                })
                .collect()
        }))
    }
    fn frame(
        &self,
        frame: ReferenceFrame,
        policy: GeometryPolicy,
    ) -> Result<ReferenceFrame, GeometryError> {
        let origin = self.qualified_point(frame.origin_mm, policy)?;
        let u = midpoint(self.placement.linear_bounds(frame.u.map(I::point)));
        let v = midpoint(self.placement.linear_bounds(frame.v.map(I::point)));
        Ok(ReferenceFrame {
            origin_mm: origin,
            u,
            v,
        })
    }
    fn qualified_point(
        &self,
        p: [f64; 3],
        policy: GeometryPolicy,
    ) -> Result<[f64; 3], GeometryError> {
        let (point, error) = self.placement.point(p)?;
        if error > policy.linear_tolerance_mm / 8. {
            return Err(failure(Cause::Unresolved));
        }
        policy.point(point)?;
        Ok(point)
    }
    pub(crate) fn topology(
        &self,
        owner: EntityId,
        revision: u64,
        policy: GeometryPolicy,
    ) -> Result<TopologySnapshot, GeometryError> {
        let mut result = self.source.topology.clone();
        result.owner = owner;
        result.revision = revision;
        for element in &mut result.elements {
            for anchored in &mut element.frames {
                anchored.frame = self.frame(anchored.frame, policy)?;
            }
            element.geometry = match element.geometry {
                Some(AnalyticGeometry::Point { position_mm }) => Some(AnalyticGeometry::Point {
                    position_mm: self.qualified_point(position_mm, policy)?,
                }),
                Some(AnalyticGeometry::Line { a_mm, b_mm }) => Some(AnalyticGeometry::Line {
                    a_mm: self.qualified_point(a_mm, policy)?,
                    b_mm: self.qualified_point(b_mm, policy)?,
                }),
                Some(AnalyticGeometry::Plane { frame }) => Some(AnalyticGeometry::Plane {
                    frame: self.frame(frame, policy)?,
                }),
                other if self.placement.is_identity() => other,
                // A near-rigid literal affine map need not preserve circles.
                // Finite owner/sampling remains authoritative; no analytic
                // circle metadata is invented from a normalized proxy.
                _ => None,
            };
        }
        Ok(result)
    }
    pub(crate) fn planar_loops(&self, index: u32) -> Vec<Vec<(u32, bool)>> {
        self.source.loops.get(&index).cloned().unwrap_or_default()
    }
    pub(crate) fn display_edges(
        &self,
        deflection: f64,
        policy: GeometryPolicy,
    ) -> Result<Vec<DisplayEdge>, GeometryError> {
        let local = (I::point(deflection) * I::point(0.5) / I::point(self.placement.sigma[1])).lo;
        let original = self.source.thread.topology_snapshot(self.feature, 0)?;
        let edges = crate::rust_display::edges(Some(&self.source.thread), &original, local)?;
        let mut result = Vec::with_capacity(edges.len() + 18);
        for edge in edges {
            let mut points = Vec::with_capacity(edge.points_mm.len());
            for p in edge.points_mm {
                let (mapped, error) = self.placement.point(p)?;
                if error > deflection / 4. {
                    return Err(failure(Cause::Unresolved));
                }
                policy.point(mapped)?;
                points.push(mapped);
            }
            result.push(DisplayEdge {
                index: edge.index,
                points_mm: points,
            });
        }
        for i in 0..18 {
            let (a, b) = if i < 12 {
                let start = if i < 6 { 0 } else { 6 };
                (i, start + (i + 1) % 6)
            } else {
                (i - 12, i - 6)
            };
            let mut points = Vec::with_capacity(2);
            for index in [a, b] {
                let exact = self
                    .placement
                    .point_bounds(self.source.hex_vertices[index])?;
                let point = midpoint(exact);
                let error = norm(std::array::from_fn(|k| exact[k] - I::point(point[k])));
                if !error.valid() || error.hi > deflection {
                    return Err(failure(Cause::Unresolved));
                }
                policy.point(point)?;
                points.push(point);
            }
            result.push(DisplayEdge {
                index: self.source.edge_offset + i as u32,
                points_mm: points,
            });
        }
        Ok(result)
    }
}
