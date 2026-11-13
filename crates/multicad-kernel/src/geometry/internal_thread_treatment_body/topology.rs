//! Source-owned finite topology for the actual stock-minus-cavity material.
use super::*;
use multicad_core::references::{Origin, TopologyElement};
use std::collections::{BTreeMap, BTreeSet};
mod anchors;
mod graph;
mod identity;

fn unavailable(message: impl AsRef<str>) -> GeometryError {
    unsupported(message.as_ref())
}

fn face_role(key: source::FaceKey) -> String {
    use source::{FaceKey, ThreadBandRole};
    match key {
        FaceKey::StartCap => "cap.start".into(),
        FaceKey::EndCap => "cap.end".into(),
        FaceKey::EntryCone => "thread.entry-chamfer".into(),
        FaceKey::TerminalShoulder => "thread.end-relief/shoulder".into(),
        FaceKey::TerminalCylinder => "thread.end-relief/cylinder".into(),
        FaceKey::MaterialOuterWall => "stock.outer-wall".into(),
        FaceKey::Thread { start, band } => format!(
            "thread/start/{start}/{}",
            match band {
                ThreadBandRole::Lower => "lower-flank",
                ThreadBandRole::Crest => "crest",
                ThreadBandRole::Upper => "upper-flank",
                ThreadBandRole::Root => "root",
            }
        ),
    }
}

fn origin(feature: EntityId, role: String) -> Origin {
    Origin {
        source: feature,
        role,
    }
}

fn unique(identities: &BTreeMap<u32, Origin>) -> Result<(), GeometryError> {
    let mut used = BTreeSet::new();
    for value in identities.values() {
        if value.source.is_nil() || value.role.trim().is_empty() || value.role.len() > 256 {
            return Err(unavailable("Invalid source reference identity"));
        }
        if !used.insert(value) {
            return Err(unavailable(
                "Distinct source entities have ambiguous reference identities",
            ));
        }
    }
    Ok(())
}

impl InternalThreadTreatmentBody {
    pub(super) fn qualified_boundary(&self) -> Result<&source::Boundary, GeometryError> {
        if let Some(boundary) = self.boundary.get() {
            return Ok(boundary);
        }
        let boundary = self
            .material
            .boundary_with_interrupt(source::BoundaryOptions::default(), None)
            .map_err(failure)?;
        // Only a successful immutable source proof is retained. Placement and
        // per-call reference handoffs remain independently qualified.
        let _ = self.boundary.set(boundary);
        Ok(self
            .boundary
            .get()
            .expect("successful immutable boundary initialization"))
    }
    pub(crate) fn source_topology(
        &self,
        frame: PlaneFrame,
        owner: EntityId,
        revision: u64,
        policy: GeometryPolicy,
    ) -> Result<TopologySnapshot, GeometryError> {
        use multicad_core::references::{
            AnalyticGeometry, AnchoredFrame, TopologyKey, TopologyKind,
        };
        if owner.is_nil() || self.feature.is_nil() {
            return Err(unavailable(
                "Internal topology needs nonempty owner and source identities",
            ));
        }
        self.validate_frame(frame, policy)?;
        let boundary = self.qualified_boundary()?;
        let identities = identity::identities(boundary, self.feature)?;
        let transport = anchors::Transport::new(self.creation_frame, frame, policy)?;
        let mut faces = Vec::new();
        for index in 0..boundary.face_count() {
            let key = self
                .material
                .face_key(index)
                .expect("proved actual material face");
            let primary = origin(self.feature, face_role(key));
            let mut origins = BTreeSet::from([primary.clone()]);
            let mut frames = Vec::new();
            let station = match key {
                source::FaceKey::StartCap => Some(0.),
                source::FaceKey::EndCap => Some(self.material.cavity().full_length_mm()),
                source::FaceKey::TerminalShoulder => self.material.cavity().terminal_station_mm(),
                _ => None,
            };
            let plane = if let Some(z) = station {
                let zero = source::Enclosure {
                    lower: 0.,
                    upper: 0.,
                };
                let p = transport.local_point([
                    zero,
                    zero,
                    source::Enclosure { lower: z, upper: z },
                ])?;
                let support = multicad_core::references::ReferenceFrame {
                    origin_mm: p,
                    u: frame.u,
                    v: if matches!(key, source::FaceKey::StartCap) {
                        frame.v.map(|v| -v)
                    } else {
                        frame.v
                    },
                };
                support.validate().map_err(unavailable)?;
                frames.push(AnchoredFrame {
                    origin: primary.clone(),
                    frame: support,
                });
                Some(support)
            } else {
                None
            };
            let stock_role = match key {
                source::FaceKey::StartCap => Some("cap.start"),
                source::FaceKey::EndCap => Some("cap.end"),
                source::FaceKey::MaterialOuterWall => Some("side"),
                _ => None,
            };
            if let Some(role) = stock_role {
                let stock_id = EntityId::from_uuid_bits(self.material.stock().feature());
                let mut matches = self.stock_topology.elements.iter().filter(|element| {
                    element.key.kind == TopologyKind::Face
                        && element.origins.iter().any(|o| {
                            o.source == stock_id
                                && if role == "side" {
                                    o.role == "side.cylinder" || o.role.starts_with("side/")
                                } else {
                                    o.role == role
                                }
                        })
                });
                let original = matches
                    .next()
                    .ok_or_else(|| unavailable("Original stock face ancestry is unresolved"))?;
                if matches.next().is_some() {
                    return Err(unavailable("Original stock face ancestry is ambiguous"));
                }
                origins.extend(original.origins.iter().cloned());
                for anchor in &original.frames {
                    frames.push(AnchoredFrame {
                        origin: anchor.origin.clone(),
                        frame: transport.frame(anchor.frame)?,
                    });
                }
            }
            faces.push(TopologyElement {
                key: TopologyKey {
                    kind: TopologyKind::Face,
                    ordinal: index,
                },
                origins,
                boundary: BTreeSet::new(),
                frames,
                geometry: plane.map(|frame| AnalyticGeometry::Plane { frame }),
            });
        }
        let mut snapshot = graph::assemble(
            boundary,
            faces,
            &identities.edges,
            &identities.vertices,
            owner,
            revision,
        )?;
        let vertex_offset = boundary.face_count() as usize + boundary.topology_counts().1;
        for (ordinal, (_, point)) in boundary.vertices().enumerate() {
            let position_mm = transport.local_point(point.position_bounds())?;
            snapshot.elements[vertex_offset + ordinal].geometry =
                Some(AnalyticGeometry::Point { position_mm });
        }
        self.inherit_stock_rims(boundary, &mut snapshot)?;
        snapshot.validate().map_err(unavailable)?;
        Ok(snapshot)
    }
    fn inherit_stock_rims(
        &self,
        boundary: &source::Boundary,
        snapshot: &mut TopologySnapshot,
    ) -> Result<(), GeometryError> {
        use multicad_core::references::TopologyKind;
        let stock_id = EntityId::from_uuid_bits(self.material.stock().feature());
        let outer = self.material.face_count() - 1;
        for cap in 0..2 {
            let cap_origin = origin(
                stock_id,
                if cap == 0 { "cap.start" } else { "cap.end" }.into(),
            );
            let mut old =
                self.stock_topology.elements.iter().filter(|e| {
                    e.key.kind == TopologyKind::Edge && e.boundary.contains(&cap_origin)
                });
            let original = old
                .next()
                .ok_or_else(|| unavailable("Original full stock rim ancestry is unresolved"))?;
            if old.next().is_some() {
                return Err(unavailable("Original full stock rim ancestry is ambiguous"));
            }
            let old = original;
            // The actual primitive source has exactly one full outer rim per
            // cap. Internal treatments leave both circles unchanged.
            // The sealed source Graph::material constructor registers exactly
            // one additional full stock circle for each cap, attaches it to
            // the new outer face, then checks the closed material boundary.
            // No cavity curve can mention that new outer face. Reuse that
            // proved incidence rather than requalifying inaccessible states.
            let (ordinal, _) = boundary
                .curves()
                .enumerate()
                .find(|(_, (_, c))| c.faces() == [cap, outer] || c.faces() == [outer, cap])
                .expect("sealed material boundary retains each original full stock rim");
            let edge_index = boundary.face_count() as usize + ordinal;
            let primary = snapshot.elements[edge_index].origins.clone();
            snapshot.elements[edge_index]
                .origins
                .extend(old.origins.iter().cloned());
            for element in &mut snapshot.elements {
                if !element.boundary.is_disjoint(&primary) {
                    element.boundary.extend(old.origins.iter().cloned());
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use multicad_core::references::{TopologyKey, TopologyKind};
    use multicad_rust_kernel::Rational;
    #[test]
    fn actual_material_faces_curves_and_endpoint_incidence_survive_the_native_handoff() {
        for (left, starts) in [(false, 1), (true, 2)] {
            let frame = PlaneFrame::default();
            let stock = super::super::tests::fixture_stock(Rational::new(4, 1).unwrap(), frame);
            let solid = stock
                .rust_cut_authored_internal_thread_treatments_with_interrupt(
                    super::super::tests::cut(frame, left, starts),
                    EntityId::from_uuid_bits(8),
                    GeometryPolicy::default(),
                    None,
                )
                .unwrap();
            let (body, _) = solid.internal_thread_treatment_body().unwrap();
            let boundary = body
                .material
                .boundary_with_interrupt(source::BoundaryOptions::default(), None)
                .unwrap();
            let identities = identity::identities(&boundary, body.feature).unwrap();
            let faces = (0..boundary.face_count())
                .map(|index| TopologyElement {
                    key: TopologyKey {
                        kind: TopologyKind::Face,
                        ordinal: index,
                    },
                    origins: BTreeSet::from([origin(
                        body.feature,
                        face_role(body.material.face_key(index).unwrap()),
                    )]),
                    boundary: BTreeSet::new(),
                    frames: vec![],
                    geometry: None,
                })
                .collect();
            let snapshot = graph::assemble(
                &boundary,
                faces,
                &identities.edges,
                &identities.vertices,
                EntityId::from_uuid_bits(9),
                3,
            )
            .unwrap();
            snapshot.validate().unwrap();
            let (vertices, edges) = boundary.topology_counts();
            assert_eq!(
                snapshot.elements.len(),
                vertices + edges + boundary.face_count() as usize
            );
            for element in &snapshot.elements {
                assert!(!element.boundary.is_empty());
            }
            assert!(snapshot.elements.iter().all(|e| e.geometry.is_none()));
        }
    }
    #[test]
    fn stock_radius_edits_keep_source_roles_while_actual_outer_curves_change() {
        let frame = PlaneFrame::default();
        let stock = super::super::tests::fixture_stock(Rational::new(4, 1).unwrap(), frame);
        let solid = stock
            .rust_cut_authored_internal_thread_treatments_with_interrupt(
                super::super::tests::cut(frame, false, 1),
                EntityId::from_uuid_bits(8),
                GeometryPolicy::default(),
                None,
            )
            .unwrap();
        let (body, _) = solid.internal_thread_treatment_body().unwrap();
        let before = body
            .material
            .boundary_with_interrupt(source::BoundaryOptions::default(), None)
            .unwrap();
        let zero = Rational::new(0, 1).unwrap();
        let larger = multicad_rust_kernel::cylinder::AxialCylinder::new(
            7,
            multicad_rust_kernel::Circle::new([zero; 2], Rational::new(11, 4).unwrap()).unwrap(),
            zero,
            Rational::new(4, 1).unwrap(),
        )
        .unwrap();
        let changed = source::Material::new(larger, *body.material.cavity(), 0.1).unwrap();
        let after = changed
            .boundary_with_interrupt(source::BoundaryOptions::default(), None)
            .unwrap();
        let a = identity::identities(&before, body.feature).unwrap();
        let b = identity::identities(&after, body.feature).unwrap();
        assert_eq!(
            a.vertices.values().collect::<BTreeSet<_>>(),
            b.vertices.values().collect()
        );
        assert_eq!(
            a.edges.values().collect::<BTreeSet<_>>(),
            b.edges.values().collect()
        );
        let outer = changed.face_count() - 1;
        let source_curve = before
            .curves()
            .find(|(_, c)| c.faces().contains(&outer))
            .unwrap()
            .0;
        let edited_curve = after
            .curves()
            .find(|(_, c)| c.faces().contains(&outer))
            .unwrap()
            .0;
        let p = before.curve_point(source_curve, 0., None).unwrap();
        let q = after.curve_point(edited_curve, 0., None).unwrap();
        assert!((q.point_mm[0] - p.point_mm[0] - 0.25).abs() <= p.error_mm + q.error_mm);
    }
    #[test]
    fn cached_source_topology_preserves_stock_anchor_frames_after_actual_placement() {
        let frame = PlaneFrame {
            origin: [10., -20., 30.],
            u: [0.8, 0.6, 0.],
            v: [-0.36, 0.48, 0.8],
        };
        let stock = super::super::tests::fixture_stock(Rational::new(4, 1).unwrap(), frame);
        let solid = stock
            .rust_cut_authored_internal_thread_treatments_with_interrupt(
                super::super::tests::cut(frame, false, 1),
                EntityId::from_uuid_bits(8),
                GeometryPolicy::default(),
                None,
            )
            .unwrap();
        let (body, initial) = solid.internal_thread_treatment_body().unwrap();
        assert!(
            body.source_topology(
                initial,
                EntityId::from_uuid_bits(0),
                0,
                GeometryPolicy::default()
            )
            .is_err()
        );
        assert!(body.boundary.get().is_none());
        let first = body
            .source_topology(
                initial,
                EntityId::from_uuid_bits(9),
                3,
                GeometryPolicy::default(),
            )
            .unwrap();
        first.validate().unwrap();
        let cached = body.qualified_boundary().unwrap() as *const source::Boundary;
        let transform = RigidTransform {
            rotation: [[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]],
            translation_mm: [100., 200., 300.],
        };
        let moved = solid.transformed(transform).unwrap();
        let (same, placed) = moved.internal_thread_treatment_body().unwrap();
        assert_eq!(
            cached,
            same.qualified_boundary().unwrap() as *const source::Boundary
        );
        let snapshot = same
            .source_topology(
                placed,
                EntityId::from_uuid_bits(9),
                4,
                GeometryPolicy::default(),
            )
            .unwrap();
        snapshot.validate().unwrap();
        assert_eq!(snapshot.elements.len(), first.elements.len());
        for old in &body.stock_topology.elements {
            for old_anchor in &old.frames {
                let anchor = snapshot
                    .elements
                    .iter()
                    .flat_map(|e| &e.frames)
                    .find(|a| a.origin == old_anchor.origin)
                    .unwrap();
                for i in 0..3 {
                    let expected = transform.translation_mm[i]
                        + (0..3)
                            .map(|j| transform.rotation[i][j] * old_anchor.frame.origin_mm[j])
                            .sum::<f64>();
                    assert!((anchor.frame.origin_mm[i] - expected).abs() < 1e-9);
                    let u = (0..3)
                        .map(|j| transform.rotation[i][j] * old_anchor.frame.u[j])
                        .sum::<f64>();
                    let v = (0..3)
                        .map(|j| transform.rotation[i][j] * old_anchor.frame.v[j])
                        .sum::<f64>();
                    assert!((anchor.frame.u[i] - u).abs() < 1e-12);
                    assert!((anchor.frame.v[i] - v).abs() < 1e-12);
                }
            }
            assert!(
                snapshot
                    .elements
                    .iter()
                    .any(|e| e.origins.is_superset(&old.origins))
            );
        }
    }
    #[test]
    fn partial_or_stale_reference_transfers_do_not_publish_a_snapshot_or_change_the_source() {
        let frame = PlaneFrame::default();
        let policy = GeometryPolicy::default();
        let stock = super::super::tests::fixture_stock(Rational::new(4, 1).unwrap(), frame);
        let solid = stock
            .rust_cut_authored_internal_thread_treatments_with_interrupt(
                super::super::tests::cut(frame, false, 1),
                EntityId::from_uuid_bits(8),
                policy,
                None,
            )
            .unwrap();
        let (body, _) = solid.internal_thread_treatment_body().unwrap();
        let before = solid
            .topology_snapshot(EntityId::from_uuid_bits(9), 0)
            .unwrap();
        let boundary = body.qualified_boundary().unwrap();
        let names = identity::identities(boundary, body.feature).unwrap();
        let faces = before
            .elements
            .iter()
            .filter(|e| e.key.kind == TopologyKind::Face)
            .cloned()
            .map(|mut e| {
                e.boundary.clear();
                e
            })
            .collect::<Vec<_>>();
        let mut missing = faces.clone();
        missing.pop();
        assert!(
            graph::assemble(
                boundary,
                missing,
                &names.edges,
                &names.vertices,
                before.owner,
                0
            )
            .is_err()
        );
        let mut stale_edges = names.edges.clone();
        let id = *stale_edges.keys().next().unwrap();
        let name = stale_edges.remove(&id).unwrap();
        stale_edges.insert(u32::MAX, name);
        assert!(
            graph::assemble(
                boundary,
                faces.clone(),
                &stale_edges,
                &names.vertices,
                before.owner,
                0
            )
            .is_err()
        );
        let mut stale_vertices = names.vertices.clone();
        let id = *stale_vertices.keys().next().unwrap();
        let name = stale_vertices.remove(&id).unwrap();
        stale_vertices.insert(u32::MAX, name);
        assert!(
            graph::assemble(
                boundary,
                faces.clone(),
                &names.edges,
                &stale_vertices,
                before.owner,
                0
            )
            .is_err()
        );
        assert!(
            graph::assemble(
                boundary,
                faces,
                &names.edges,
                &names.vertices,
                EntityId::from_uuid_bits(0),
                0
            )
            .is_err()
        );
        let mut ambiguous = names.edges.clone();
        let ids = ambiguous.keys().copied().take(2).collect::<Vec<_>>();
        ambiguous.insert(ids[1], ambiguous[&ids[0]].clone());
        assert!(unique(&ambiguous).is_err());
        ambiguous.insert(ids[1], origin(body.feature, " ".into()));
        assert!(unique(&ambiguous).is_err());
        let cavity = body
            .material
            .cavity()
            .boundary_with_interrupt(source::BoundaryOptions::default(), None)
            .unwrap();
        assert!(identity::identities(&cavity, body.feature).is_err());
        assert_eq!(solid.topology_snapshot(before.owner, 0).unwrap(), before);
        solid
            .manufacturing_mesh_with_receipt(MeshOptions::default())
            .unwrap();
    }
    #[test]
    fn unresolved_and_ambiguous_original_stock_metadata_refuse_transfer() {
        let frame = PlaneFrame::default();
        let policy = GeometryPolicy::default();
        let stock = super::super::tests::fixture_stock(Rational::new(4, 1).unwrap(), frame);
        let solid = stock
            .rust_cut_authored_internal_thread_treatments_with_interrupt(
                super::super::tests::cut(frame, false, 1),
                EntityId::from_uuid_bits(8),
                policy,
                None,
            )
            .unwrap();
        let (original, _) = solid.internal_thread_treatment_body().unwrap();
        let stock_id = EntityId::from_uuid_bits(original.material.stock().feature());
        for kind in [TopologyKind::Face, TopologyKind::Edge] {
            let index = original
                .stock_topology
                .elements
                .iter()
                .position(|e| {
                    e.key.kind == kind
                        && if kind == TopologyKind::Face {
                            e.origins.contains(&origin(stock_id, "cap.start".into()))
                        } else {
                            e.boundary.contains(&origin(stock_id, "cap.start".into()))
                        }
                })
                .unwrap();
            for duplicate in [false, true] {
                let mut metadata = original.stock_topology.clone();
                if duplicate {
                    metadata.elements.push(metadata.elements[index].clone());
                } else {
                    metadata.elements.remove(index);
                }
                // Geometry, source boundary, placement and requested wall are
                // untouched. Only the incoming stock reference protocol is bad.
                let body = InternalThreadTreatmentBody {
                    material: original.material,
                    feature: original.feature,
                    requested_wall: original.requested_wall,
                    stock_topology: metadata,
                    creation_frame: original.creation_frame,
                    boundary: std::sync::OnceLock::new(),
                };
                let error = body
                    .source_topology(frame, EntityId::from_uuid_bits(9), 0, policy)
                    .unwrap_err();
                assert!(error.to_string().contains(if duplicate {
                    "ambiguous"
                } else {
                    "unresolved"
                }));
            }
        }
        let mut metadata = original.stock_topology.clone();
        metadata
            .elements
            .iter_mut()
            .find(|e| e.key.kind == TopologyKind::Face)
            .unwrap()
            .origins
            .insert(origin(
                EntityId::from_uuid_bits(0),
                "invalid-stock-source".into(),
            ));
        let invalid = InternalThreadTreatmentBody {
            material: original.material,
            feature: original.feature,
            requested_wall: original.requested_wall,
            stock_topology: metadata,
            creation_frame: original.creation_frame,
            boundary: std::sync::OnceLock::new(),
        };
        assert!(
            invalid
                .source_topology(frame, EntityId::from_uuid_bits(9), 0, policy)
                .is_err()
        );
        original
            .source_topology(frame, EntityId::from_uuid_bits(9), 0, policy)
            .unwrap();
        solid
            .manufacturing_mesh_with_receipt(MeshOptions::default())
            .unwrap();
    }
}
