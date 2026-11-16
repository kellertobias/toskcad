//! Actual thread-domain measures plus exact hex stock. Retained net measures
//! never acquire the separate original authored-law moment receipt.
use super::*;
fn planar_frame(origin: [f64; 3], u: [f64; 3], v: [f64; 3]) -> ReferenceFrame {
    ReferenceFrame {
        origin_mm: origin,
        u,
        v,
    }
}
pub(super) fn build(
    thread: Solid,
    intent: Material,
    minimum_wall: f64,
    feature: EntityId,
    policy: GeometryPolicy,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<Source, GeometryError> {
    checkpoint(interrupt)?;
    let spec = intent.thread().spec();
    let kind = intent.kind();
    let (af, head) = match kind {
        PipeHardwareKind::HexHead {
            across_flats_mm,
            height_mm,
        } => (across_flats_mm, Some(height_mm)),
        PipeHardwareKind::HexNut { across_flats_mm } => (across_flats_mm, None),
    };
    let mut local_topology = thread.topology_snapshot(feature, 0)?;
    let mut measured = BTreeMap::new();
    let (caps, thread_volume, radius) = if let Some(source) = thread.rust_clipped_npt_source() {
        let measures = thread
            .rust_clipped_npt_face_measures()
            .expect("actual clipped source measures");
        for (index, measure) in measures.iter().enumerate() {
            measured.insert(
                index as u32,
                (
                    interval(measure.area_mm2()),
                    measure.centroid_mm().map(interval),
                ),
            );
        }
        (
            [5u32, 6u32],
            interval(
                thread
                    .rust_clipped_npt_volume_bounds_mm3()
                    .expect("actual clipped volume"),
            ),
            I::point(source.stock_radius_limit_mm()),
        )
    } else {
        let (boundary, pose) = thread
            .rust_radial_boundary()
            .ok_or_else(|| failure(Cause::Unsupported))?;
        if pose != PlaneFrame::default() {
            return Err(failure(Cause::Unsupported));
        }
        let (provenance, construction) = thread
            .rust_authored_axial_pipe_thread_source()
            .ok_or_else(|| failure(Cause::Unsupported))?;
        if construction != PlaneFrame::default() || provenance.transport != *intent.thread() {
            return Err(failure(Cause::Unsupported));
        }
        let author = provenance
            .transport
            .radial_bounds([0., spec.pitch_mm], [0., spec.length_mm])
            .map_err(failure)?;
        let radius = I::point(author[1])
            + I::point(provenance.local_source_error_mm)
            + I::point(boundary.boundary_error_mm());
        if (I::point(af) / I::point(2.) - radius).lo <= minimum_wall {
            return Err(failure(Cause::Unsupported));
        }
        // Independent per-face relative area/centroid accuracy, not an
        // aggregate summed-area-width contract. Keep the actual intervals.
        // Full inertia is not part of this facade descriptor query.
        let properties = boundary
            .face_properties_with_interrupt(1e-4, interrupt)
            .map_err(radial_failure)?;
        for face in properties {
            measured.insert(
                face.face_index as u32,
                (interval(face.area_mm2), face.centroid_mm.map(interval)),
            );
        }
        (
            provenance.cap_faces.map(|i| i as u32),
            interval(
                boundary
                    .volume_with_interrupt(1e-4, interrupt)
                    .map_err(radial_failure)?
                    .volume_mm3,
            ),
            radius,
        )
    };
    let child_faces = local_topology
        .elements
        .iter()
        .filter(|e| e.key.kind == TopologyKind::Face)
        .map(|e| e.key.ordinal)
        .collect::<Vec<_>>();
    let clearance = I::point(af) / I::point(2.) - radius;
    if !clearance.valid() || clearance.lo <= minimum_wall {
        return Err(failure(Cause::Unsupported));
    }
    let removed = |index: u32| index == caps[1] || (head.is_none() && index == caps[0]);
    let mut faces = Vec::new();
    let mut face_map = BTreeMap::new();
    let mut loops = BTreeMap::new();
    let child_cap_loops = caps.map(|cap| crate::rust_display::planar_loops(&thread, cap));
    for &old_index in &child_faces {
        checkpoint(interrupt)?;
        if removed(old_index) {
            continue;
        }
        let index = faces.len() as u32;
        let (area, center) = *measured
            .get(&old_index)
            .expect("every actual finite source face measured");
        let frame = local_topology
            .elements
            .iter()
            .find(|e| e.key.kind == TopologyKind::Face && e.key.ordinal == old_index)
            .and_then(|e| match e.geometry {
                Some(AnalyticGeometry::Plane { frame }) => Some(frame),
                _ => None,
            });
        faces.push(LocalFace {
            area,
            center,
            frame,
            normal: frame.map(|f| f.normal().map(I::point)),
        });
        face_map.insert(old_index, index);
        if frame.is_some() {
            loops.insert(index, crate::rust_display::planar_loops(&thread, old_index));
        }
    }
    let root3 = I::point(3.).sqrt();
    let half = I::point(af) / I::point(2.);
    let r = I::point(af) / root3;
    let zero = I::ZERO;
    let xy = [
        [r, zero],
        [r / I::point(2.), half],
        [-r / I::point(2.), half],
        [-r, zero],
        [-r / I::point(2.), -half],
        [r / I::point(2.), -half],
    ];
    let lower = I::point(if head.is_some() { spec.length_mm } else { 0. });
    let upper = if let Some(h) = head {
        I::point(spec.length_mm) + I::point(h)
    } else {
        I::point(spec.length_mm)
    };
    if !upper.valid() || lower.hi >= upper.lo {
        return Err(failure(Cause::Unresolved));
    }
    let hex_vertices = std::array::from_fn(|i| {
        [
            xy[i % 6][0],
            xy[i % 6][1],
            if i < 6 { lower } else { upper },
        ]
    });
    for p in hex_vertices {
        let literal = midpoint(p);
        let error = (0..3).fold(I::ZERO, |s, k| {
            s + I::point((p[k] - I::point(literal[k])).magnitude())
        });
        if !error.valid() || error.hi > policy.linear_tolerance_mm / 8. {
            return Err(failure(Cause::Unresolved));
        }
        policy.point(literal)?;
    }
    let base = faces.len() as u32;
    let mesh_base = child_faces
        .iter()
        .copied()
        .max()
        .expect("actual finite thread faces")
        + 1;
    for i in 0..6 {
        let j = (i + 1) % 6;
        let dx = xy[j][0] - xy[i][0];
        let dy = xy[j][1] - xy[i][1];
        let length = (dx.square() + dy.square()).sqrt();
        let normal = [dy / length, -dx / length, zero];
        let center = [
            (xy[i][0] + xy[j][0]) / I::point(2.),
            (xy[i][1] + xy[j][1]) / I::point(2.),
            (lower + upper) / I::point(2.),
        ];
        let origin = midpoint(hex_vertices[i]);
        let u = midpoint([dx / length, dy / length, zero]);
        let frame = Some(planar_frame(origin, u, [0., 0., 1.]));
        faces.push(LocalFace {
            area: length * (upper - lower),
            center,
            frame,
            normal: Some(normal),
        });
        face_map.insert(mesh_base + i as u32, base + i as u32);
    }
    let hex_area = root3 * I::point(af).square() / I::point(2.);
    let cap_face = |cap: u32, z: I, positive: bool| -> Result<LocalFace, GeometryError> {
        let (cap_area, cap_center) = *measured.get(&cap).expect("actual cap source measured");
        let area = hex_area - cap_area;
        if !area.valid() || area.lo <= 0. {
            return Err(failure(Cause::Unsupported));
        }
        let center = [
            -cap_area * cap_center[0] / area,
            -cap_area * cap_center[1] / area,
            z,
        ];
        Ok(LocalFace {
            area,
            center,
            frame: Some(planar_frame(
                [0., 0., z.lo.midpoint(z.hi)],
                [1., 0., 0.],
                [0., if positive { 1. } else { -1. }, 0.],
            )),
            normal: Some([zero, zero, I::point(if positive { 1. } else { -1. })]),
        })
    };
    if head.is_some() {
        faces.push(cap_face(caps[1], lower, false)?);
        faces.push(LocalFace {
            area: hex_area,
            center: [zero, zero, upper],
            frame: Some(planar_frame(
                [0., 0., upper.lo.midpoint(upper.hi)],
                [1., 0., 0.],
                [0., 1., 0.],
            )),
            normal: Some([zero, zero, I::ONE]),
        });
    } else {
        faces.push(cap_face(caps[0], lower, false)?);
        faces.push(cap_face(caps[1], upper, true)?);
    }
    face_map.insert(mesh_base + 6, base + 6);
    face_map.insert(mesh_base + 7, base + 7);
    let hex_volume = hex_area * (upper - lower);
    let volume = if head.is_some() {
        thread_volume + hex_volume
    } else {
        hex_volume - thread_volume
    };
    if !volume.valid() || volume.lo <= 0. {
        return Err(failure(Cause::Unsupported));
    }
    let edge_offset = local_topology
        .elements
        .iter()
        .filter(|e| e.key.kind == TopologyKind::Edge)
        .map(|e| e.key.ordinal)
        .max()
        .map_or(0, |i| i + 1);
    let vertex_offset = local_topology
        .elements
        .iter()
        .filter(|e| e.key.kind == TopologyKind::Vertex)
        .map(|e| e.key.ordinal)
        .max()
        .map_or(0, |i| i + 1);
    local_topology.elements.retain_mut(|e| {
        if e.key.kind != TopologyKind::Face {
            return true;
        }
        if removed(e.key.ordinal) {
            return false;
        }
        e.key.ordinal = face_map[&e.key.ordinal];
        if head.is_none() {
            if let Some(AnalyticGeometry::Plane { frame }) = &mut e.geometry {
                frame.v = frame.v.map(|v| -v);
            }
            for f in &mut e.frames {
                f.frame.v = f.frame.v.map(|v| -v);
            }
        }
        true
    });
    if head.is_none() {
        for f in &mut faces[..base as usize] {
            if let Some(frame) = &mut f.frame {
                frame.v = frame.v.map(|v| -v);
            }
            if let Some(n) = &mut f.normal {
                *n = n.map(|v| -v);
            }
        }
    }
    let child_box = thread.bounds()?;
    let bounds = if head.is_some() {
        [
            [child_box.min[0].min(-r.hi), child_box.max[0].max(r.hi)],
            [
                child_box.min[1].min(-half.hi),
                child_box.max[1].max(half.hi),
            ],
            [
                child_box.min[2].min(lower.lo),
                child_box.max[2].max(upper.hi),
            ],
        ]
    } else {
        [[-r.hi, r.hi], [-half.hi, half.hi], [lower.lo, upper.hi]]
    };
    let mut result = Source {
        feature,
        intent,
        thread,
        faces,
        face_map,
        topology: local_topology,
        loops,
        hex_vertices,
        edge_offset,
        vertex_offset,
        bounds,
        volume,
        clearance,
        minimum_wall,
    };
    topology::extend(&mut result, base, child_cap_loops, head.is_some())?;
    checkpoint(interrupt)?;
    Ok(result)
}
