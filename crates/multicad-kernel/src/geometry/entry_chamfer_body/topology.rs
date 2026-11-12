use super::*;
use source::{Face as Role, SourceCurve};
impl EntryChamferBody {
    fn face_role(&self, role: Role) -> String {
        match role {
            Role::EntryCap => "cap.start".into(),
            Role::EndCap if self.geometry.terminal_relief().is_some() => {
                "thread.end-relief/shoulder".into()
            }
            Role::EndCap | Role::FarCap => "cap.end".into(),
            Role::ThreadBand { start, band } => {
                format!("thread/{}", self.references.thread_faces[&(start, band)])
            }
            Role::ConeHalf { start, half } => format!("thread.entry-chamfer/{start}/{half:?}"),
            Role::NeckWall => "thread.end-relief/neck.wall".into(),
        }
    }
    fn reference_edge_role(&self, edge: &EntryChamferEdge) -> String {
        if let Some(role) = &edge.rail {
            return role.clone();
        }
        if let Some((start, band)) = edge.original {
            format!(
                "{}/rim/{}",
                self.face_role(Role::EndCap),
                helical_body::phase_role(self.geometry.source_thread(), start, band)
            )
        } else {
            self.edge_role(self.topology.edges()[edge.pieces[0]].source())
        }
    }
    fn reference_edge_geometry(
        &self,
        edge: &EntryChamferEdge,
        frame: PlaneFrame,
    ) -> Option<AnalyticGeometry> {
        if let Some((start, band)) = edge.original {
            let tool = self.geometry.source_thread();
            let original = (tool.face_index(i64::from(start), band) - 2) % 4;
            if original != 1 && (original != 3 || tool.spec().root_radius_mm.is_some()) {
                return None;
            }
            let knots = tool.profile_breaks();
            let pitch = tool.spec().pitch_mm;
            let a = f64::from(start) * pitch + knots[band];
            let b = if band == 4 {
                f64::from(start + 1) * pitch + knots[1]
            } else {
                f64::from(start) * pitch + knots[band + 1]
            };
            let radius_mm = if original == 1 {
                tool.spec().major_radius_mm
            } else {
                tool.root_radius_mm()
            };
            let z = self.geometry.material_length_mm();
            Some(AnalyticGeometry::CircularArc {
                frame: ReferenceFrame {
                    origin_mm: position(frame, [0., 0., z]),
                    u: frame.u,
                    v: frame.v,
                },
                radius_mm,
                start_rad: (z - a) / tool.helical_rate(),
                sweep_rad: (a - b) / tool.helical_rate(),
            })
        } else if let SourceCurve::Circle { radius_mm, z_mm } =
            self.topology.edges()[edge.pieces[0]].source()
        {
            Some(AnalyticGeometry::CircularArc {
                frame: ReferenceFrame {
                    origin_mm: position(frame, [0., 0., z_mm]),
                    u: frame.u,
                    v: frame.v,
                },
                radius_mm,
                start_rad: 0.,
                sweep_rad: std::f64::consts::TAU,
            })
        } else {
            None
        }
    }
    fn edge_role(&self, curve: SourceCurve) -> String {
        let endpoint = |p: source::PhaseEndpoint| format!("{:?}", p.constraint());
        match curve {
            SourceCurve::Intersection {
                start,
                band,
                phase,
                active,
            } => format!(
                "entry.intersection/{start}/{band}/{active}/{}/{}",
                endpoint(phase[0]),
                endpoint(phase[1])
            ),
            SourceCurve::EntryArc { start, phase } => format!(
                "entry.cap-rim/{start}/{}/{}",
                endpoint(phase[0]),
                endpoint(phase[1])
            ),
            SourceCurve::EndProfile { start, band, phase } => format!(
                "entry.end-rim/{start}/{band}/{}/{}",
                endpoint(phase[0]),
                endpoint(phase[1])
            ),
            SourceCurve::WallRail {
                start,
                band,
                phase,
                active,
            } => format!(
                "entry.wall-rail/{start}/{band}/{active}/{}",
                endpoint(phase)
            ),
            SourceCurve::ConeRail { start, band, phase } => {
                format!("entry.cone-rail/{start}/{band}/{}", endpoint(phase))
            }
            SourceCurve::Circle { z_mm, .. } => {
                if z_mm == self.geometry.full_length_mm() {
                    "cap.end/neck-rim".into()
                } else {
                    "thread.end-relief/neck.entry-rim".into()
                }
            }
            SourceCurve::NeckSeam { .. } => "thread.end-relief/neck.seam".into(),
        }
    }
    pub(crate) fn topology(
        &self,
        frame: PlaneFrame,
        owner: EntityId,
        revision: u64,
    ) -> Result<TopologySnapshot, GeometryError> {
        let graph = &self.topology;
        let faces = graph
            .faces()
            .iter()
            .map(|f| self.origin(self.face_role(f.role())))
            .collect::<Vec<_>>();
        let edges = self
            .references
            .edges
            .iter()
            .map(|e| self.origin(self.reference_edge_role(e)))
            .collect::<Vec<_>>();
        let mut elements = Vec::new();
        let mut adjacent = vec![BTreeSet::new(); edges.len()];
        let mut corners = vec![BTreeSet::new(); graph.vertices().len()];
        for (i, f) in graph.faces().iter().enumerate() {
            let mut origins = BTreeSet::from([faces[i].clone()]);
            if matches!(f.role(), Role::ThreadBand { band: 1, .. }) {
                origins.insert(self.origin("thread.root"));
            }
            let boundary = f
                .loops()
                .iter()
                .flatten()
                .map(|u| {
                    let edge = self.references.raw_to_edge[u.edge];
                    adjacent[edge].insert(faces[i].clone());
                    edges[edge].clone()
                })
                .collect();
            let plane = self.cap_frame(frame, f.role());
            elements.push(TopologyElement {
                key: TopologyKey {
                    kind: TopologyKind::Face,
                    ordinal: i as u32,
                },
                origins,
                boundary,
                frames: plane
                    .into_iter()
                    .map(|f| AnchoredFrame {
                        origin: faces[i].clone(),
                        frame: f,
                    })
                    .collect(),
                geometry: plane.map(|frame| AnalyticGeometry::Plane { frame }),
            });
        }
        if let Some(stock) = &self.receiving {
            stock.extend_topology(&self.geometry, frame, &mut elements);
        }
        let mut physical_vertices = vec![false; graph.vertices().len()];
        for (i, e) in graph.edges().iter().enumerate() {
            let edge = self.references.raw_to_edge[i];
            for v in e.vertices() {
                corners[v].insert(edges[edge].clone());
                if !matches!(
                    e.source(),
                    SourceCurve::EndProfile { .. }
                        | SourceCurve::Circle { .. }
                        | SourceCurve::NeckSeam { .. }
                ) {
                    physical_vertices[v] = true;
                }
            }
        }
        for (i, e) in self.references.edges.iter().enumerate() {
            let geometry = self.reference_edge_geometry(e, frame);
            elements.push(TopologyElement {
                key: TopologyKey {
                    kind: TopologyKind::Edge,
                    ordinal: i as u32,
                },
                origins: BTreeSet::from([edges[i].clone()]),
                boundary: adjacent[i].clone(),
                frames: vec![],
                geometry,
            });
        }
        let mut vertex_index = 0;
        for (i, point) in graph.vertices().iter().enumerate() {
            let inherited = self.references.corners.get(&i);
            if inherited.is_none() && !physical_vertices[i] {
                continue;
            }
            let name = corners[i]
                .iter()
                .map(|o| o.role.as_str())
                .collect::<Vec<_>>()
                .join("|");
            let name = multicad_core::project_file::ContentHash::of(name.as_bytes());
            elements.push(TopologyElement {
                key: TopologyKey {
                    kind: TopologyKind::Vertex,
                    ordinal: vertex_index,
                },
                origins: BTreeSet::from([self.origin(
                    inherited
                        .cloned()
                        .unwrap_or_else(|| format!("entry.vertex/{name}")),
                )]),
                boundary: corners[i].clone(),
                frames: vec![],
                geometry: Some(AnalyticGeometry::Point {
                    position_mm: position(frame, point.nominal_mm),
                }),
            });
            vertex_index += 1;
        }
        let snapshot = TopologySnapshot {
            owner,
            revision,
            elements,
        };
        snapshot.validate().map_err(|e| {
            error(
                "Entry chamfer topology",
                GeometryErrorKind::InvalidResult,
                e,
            )
        })?;
        Ok(snapshot)
    }
}
