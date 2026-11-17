use super::*;
use multicad_core::references::{AnchoredFrame, ReferenceFrame};
impl RoundedBlindBoreBody {
    fn origin(&self, role: impl Into<String>) -> Origin {
        Origin {
            source: self.feature,
            role: role.into(),
        }
    }
    pub(crate) fn topology(
        &self,
        frame: PlaneFrame,
        owner: EntityId,
        revision: u64,
    ) -> Result<TopologySnapshot, GeometryError> {
        let mut result = self.inner.topology(frame, owner, revision)?;
        let inward = |origin: &Origin| Origin {
            source: origin.source,
            role: if origin.role == "cap.end" {
                "receiving.floor".into()
            } else if origin.role.starts_with("thread/") || origin.role == "thread.root" {
                format!("receiving.blind/{}", origin.role)
            } else {
                origin.role.clone()
            },
        };
        for element in &mut result.elements {
            element.origins = element.origins.iter().map(inward).collect();
            element.boundary = element.boundary.iter().map(inward).collect();
            for anchor in &mut element.frames {
                anchor.origin = inward(&anchor.origin);
            }
            if element.key.kind == TopologyKind::Face && element.key.ordinal == 1 {
                if let Some(AnalyticGeometry::Plane { frame }) = &mut element.geometry {
                    frame.v = frame.v.map(|v| -v);
                }
                for anchor in &mut element.frames {
                    anchor.frame.v = anchor.frame.v.map(|v| -v);
                }
            }
        }
        let wall = self.origin("receiving.outer_wall");
        let top = self.origin("receiving.stock.end");
        let rims = [
            self.origin("receiving.outer_rim.start"),
            self.origin("receiving.outer_rim.end"),
        ];
        result.elements[0].boundary.insert(rims[0].clone());
        for (cap, rim) in rims.iter().enumerate() {
            let z = if cap == 0 {
                0.
            } else {
                self.geometry.stock_height_mm().as_f64()
            };
            result.elements.push(TopologyElement {
                key: TopologyKey {
                    kind: TopologyKind::Edge,
                    ordinal: 15 + cap as u32,
                },
                origins: BTreeSet::from([rim.clone()]),
                boundary: BTreeSet::from([
                    if cap == 0 {
                        self.origin("cap.start")
                    } else {
                        top.clone()
                    },
                    wall.clone(),
                ]),
                frames: vec![],
                geometry: Some(AnalyticGeometry::CircularArc {
                    frame: ReferenceFrame {
                        origin_mm: position(frame, [0., 0., z]),
                        u: frame.u,
                        v: frame.v,
                    },
                    radius_mm: self.geometry.stock_radius_mm().as_f64(),
                    start_rad: 0.,
                    sweep_rad: std::f64::consts::TAU,
                }),
            });
        }
        result.elements.push(TopologyElement {
            key: TopologyKey {
                kind: TopologyKind::Face,
                ordinal: 7,
            },
            origins: BTreeSet::from([wall]),
            boundary: rims.iter().cloned().collect(),
            frames: vec![],
            geometry: None,
        });
        let top_frame = ReferenceFrame {
            origin_mm: position(frame, [0., 0., self.geometry.stock_height_mm().as_f64()]),
            u: frame.u,
            v: frame.v,
        };
        result.elements.push(TopologyElement {
            key: TopologyKey {
                kind: TopologyKind::Face,
                ordinal: 8,
            },
            origins: BTreeSet::from([top.clone()]),
            boundary: BTreeSet::from([rims[1].clone()]),
            frames: vec![AnchoredFrame {
                origin: top,
                frame: top_frame,
            }],
            geometry: Some(AnalyticGeometry::Plane { frame: top_frame }),
        });
        if let Some((source, authored, reverse)) = &self.inherited {
            let chart = helical_body::stock_frame::Chart::new(*authored);
            for element in &mut result.elements {
                let old = match (element.key.kind, element.key.ordinal) {
                    (TopologyKind::Face, 0) => {
                        Some((TopologyKind::Face, usize::from(*reverse) as u32))
                    }
                    (TopologyKind::Face, 8) => {
                        Some((TopologyKind::Face, usize::from(!*reverse) as u32))
                    }
                    (TopologyKind::Face, 7) => Some((TopologyKind::Face, 2)),
                    (TopologyKind::Edge, rim @ 15..=16) => Some((
                        TopologyKind::Edge,
                        if *reverse { 16 - rim } else { rim - 15 },
                    )),
                    _ => None,
                };
                if let Some((kind, ordinal)) = old
                    && let Some(old) = source
                        .elements
                        .iter()
                        .find(|e| e.key.kind == kind && e.key.ordinal == ordinal)
                {
                    element.origins.extend(old.origins.iter().cloned());
                    for anchor in &old.frames {
                        let origin = chart.point(anchor.frame.origin_mm);
                        let u = chart.vector(anchor.frame.u);
                        let v = chart.vector(anchor.frame.v);
                        let direction = |p: Point| {
                            std::array::from_fn(|i| {
                                frame.u[i] * p[0] + frame.v[i] * p[1] + frame.normal()[i] * p[2]
                            })
                        };
                        element.frames.push(AnchoredFrame {
                            origin: anchor.origin.clone(),
                            frame: ReferenceFrame {
                                origin_mm: position(frame, origin),
                                u: direction(u),
                                v: direction(v),
                            },
                        });
                    }
                }
            }
        }
        result
            .validate()
            .map_err(|e| failure(GeometryErrorKind::InvalidResult, e))?;
        Ok(result)
    }
    pub(in crate::geometry) fn display_edges(
        &self,
        frame: PlaneFrame,
        deflection: f64,
    ) -> Result<Vec<DisplayEdge>, GeometryError> {
        let mut edges = self.inner.display_edges(frame, deflection)?;
        let mut topology = self.topology(frame, self.feature, 0)?;
        topology
            .elements
            .retain(|e| e.key.kind == TopologyKind::Edge && e.key.ordinal >= 15);
        edges.extend(crate::rust_display::edges(None, &topology, deflection)?);
        Ok(edges)
    }
}
