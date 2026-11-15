use super::*;
use crate::geometry::circular_pocket_body::face_origin;
impl OrthogonalSeatBody {
    pub(crate) fn topology(
        &self,
        frame: PlaneFrame,
        owner: EntityId,
        revision: u64,
    ) -> Result<TopologySnapshot, GeometryError> {
        let mut snapshot = self.original.topology.as_ref().clone();
        snapshot.owner = owner;
        snapshot.revision = revision;
        let stock_faces = self.source.boundary().faces.len();
        let wall_origins: Vec<_> = self
            .original
            .through_origins
            .iter()
            .map(|o| BTreeSet::from([o.clone(), face_origin(o, "wall")]))
            .collect();
        let cone_origins: Vec<_> = self
            .origins
            .iter()
            .map(|o| BTreeSet::from([o.clone(), face_origin(o, "cone")]))
            .collect();
        let mut walls = Vec::new();
        let mut cones = Vec::new();
        for (i, bore) in self.source.bores().iter().enumerate() {
            let mut adjacent = BTreeSet::new();
            for end in 0..2 {
                let cap = bore.faces[end];
                let normal = self.source.boundary().faces[cap].normal;
                let center = bore.centers_mm[end];
                let seat = self
                    .source
                    .seats()
                    .iter()
                    .enumerate()
                    .find(|(_, s)| s.bore == i && usize::from(s.end) == end);
                let exterior_owner = seat.map_or(&wall_origins[i], |(j, _)| &cone_origins[j]);
                snapshot.elements[cap]
                    .boundary
                    .extend(exterior_owner.iter().cloned());
                let cap_origins = snapshot.elements[cap].origins.clone();
                let circle =
                    |ordinal: usize, center: Point, radius: f64, origins: BTreeSet<Origin>| {
                        let a = (0..3)
                            .find(|&j| normal[j] != 0.)
                            .expect("qualified cap axis");
                        let u = std::array::from_fn(|j| f64::from(j == (a + 1) % 3));
                        let mut edge = element(
                            TopologyKind::Edge,
                            ordinal,
                            origins,
                            AnalyticGeometry::CircularArc {
                                frame: ReferenceFrame {
                                    origin_mm: center,
                                    u,
                                    v: cross(normal, u),
                                },
                                radius_mm: radius,
                                start_rad: 0.,
                                sweep_rad: TAU,
                            },
                        );
                        edge.boundary = edge.origins.clone();
                        edge
                    };
                snapshot.elements.push(circle(
                    self.original.edge_count as usize + 2 * i + end,
                    center,
                    seat.map_or(bore.radius_mm, |(_, s)| s.radius_mm),
                    cap_origins.union(exterior_owner).cloned().collect(),
                ));
                if let Some((j, s)) = seat {
                    let inner = std::array::from_fn(|k| center[k] - normal[k] * s.depth_mm);
                    snapshot.elements.push(circle(
                        self.original.edge_count as usize + 2 * self.source.bores().len() + j,
                        inner,
                        bore.radius_mm,
                        wall_origins[i].union(&cone_origins[j]).cloned().collect(),
                    ));
                    adjacent.extend(cone_origins[j].iter().cloned());
                    cones.push(TopologyElement {
                        key: TopologyKey {
                            kind: TopologyKind::Face,
                            ordinal: (stock_faces + self.source.bores().len() + j) as u32,
                        },
                        origins: cone_origins[j].clone(),
                        boundary: cap_origins.union(&wall_origins[i]).cloned().collect(),
                        frames: vec![],
                        geometry: None,
                    });
                } else {
                    adjacent.extend(cap_origins);
                }
            }
            walls.push(TopologyElement {
                key: TopologyKey {
                    kind: TopologyKind::Face,
                    ordinal: (stock_faces + i) as u32,
                },
                origins: wall_origins[i].clone(),
                boundary: adjacent,
                frames: vec![],
                geometry: None,
            });
        }
        snapshot.elements.extend(walls);
        snapshot.elements.extend(cones);
        for element in &mut snapshot.elements {
            let point = |p| transform_point(frame, p);
            let vector = |v| transform_vector(frame, v);
            element.geometry = element.geometry.map(|g| g.transformed(point, vector));
            for anchor in &mut element.frames {
                anchor.frame = ReferenceFrame {
                    origin_mm: point(anchor.frame.origin_mm),
                    u: vector(anchor.frame.u),
                    v: vector(anchor.frame.v),
                };
            }
        }
        snapshot
            .validate()
            .map_err(|e| invalid(format!("Countersunk original stock topology: {e}")))?;
        Ok(snapshot)
    }
}
