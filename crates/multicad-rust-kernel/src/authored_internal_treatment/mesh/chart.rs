//! Actual q/v facet candidates. Acceptance remains whole-source chord/normal.
use super::{owners::Owners, regular::Patch, *};
#[derive(Clone, Copy)]
pub(super) struct Edge {
    pub vertices: [usize; 2],
    pub v: [I; 2],
    pub slope: I,
    pub intercept: I,
    pub kind: Kind,
    pub errors: [f64; 2],
}
#[derive(Clone)]
pub(super) struct Shell {
    pub mesh: DisplayMesh,
    pub chord: f64,
    pub normal: f64,
    pub mouth: Vec<Edge>,
    pub top: Vec<Edge>,
}
#[cfg(test)]
pub(super) fn shell(
    body: &Cavity,
    options: Tessellation,
    subdivisions: usize,
    rows: u32,
    w: &mut Work<'_>,
) -> Result<Option<Shell>, Error> {
    Ok(match shell_attempt(body, options, subdivisions, rows, w)? {
        Attempt::Complete(shell) => Some(shell),
        Attempt::Refine | Attempt::Repartition => None,
    })
}
#[cfg(test)]
pub(super) fn shell_attempt(
    body: &Cavity,
    options: Tessellation,
    subdivisions: usize,
    rows: u32,
    w: &mut Work<'_>,
) -> Result<Attempt<Shell>, Error> {
    build(body, options, (subdivisions, rows, false), None, w)
}
pub(super) fn grid_attempt(
    body: &Cavity,
    options: Tessellation,
    grid: &Grid,
    w: &mut Work<'_>,
) -> Result<Attempt<Shell>, Error> {
    build(
        body,
        options,
        (grid.subdivisions, grid.rows(), grid.shifted),
        None,
        w,
    )
}
pub(super) fn traced_grid_attempt(
    body: &Cavity,
    options: Tessellation,
    grid: &Grid,
    w: &mut Work<'_>,
) -> Result<Attempt<(Shell, trace::Trace)>, Error> {
    traced(
        body,
        options,
        (grid.subdivisions, grid.rows(), grid.shifted),
        w,
    )
}
#[cfg(test)]
pub(super) fn traced_shell(
    body: &Cavity,
    options: Tessellation,
    subdivisions: usize,
    rows: u32,
    w: &mut Work<'_>,
) -> Result<Option<(Shell, trace::Trace)>, Error> {
    Ok(
        match traced_attempt(body, options, subdivisions, rows, w)? {
            Attempt::Complete(value) => Some(value),
            Attempt::Refine | Attempt::Repartition => None,
        },
    )
}
#[cfg(test)]
pub(super) fn traced_attempt(
    body: &Cavity,
    options: Tessellation,
    subdivisions: usize,
    rows: u32,
    w: &mut Work<'_>,
) -> Result<Attempt<(Shell, trace::Trace)>, Error> {
    traced(body, options, (subdivisions, rows, false), w)
}
fn traced(
    body: &Cavity,
    options: Tessellation,
    grid: (usize, u32, bool),
    w: &mut Work<'_>,
) -> Result<Attempt<(Shell, trace::Trace)>, Error> {
    let mut trace = trace::Trace::at_origin(
        options.triangle_limit.min(131_072),
        if grid.2 { sampling::TURN_OFFSET } else { 0. },
        w,
    )?;
    Ok(match build(body, options, grid, Some(&mut trace), w)? {
        Attempt::Complete(shell) => Attempt::Complete((shell, trace)),
        Attempt::Refine => Attempt::Refine,
        Attempt::Repartition => Attempt::Repartition,
    })
}
fn build(
    body: &Cavity,
    options: Tessellation,
    grid: (usize, u32, bool),
    mut trace: Option<&mut trace::Trace>,
    w: &mut Work<'_>,
) -> Result<Attempt<Shell>, Error> {
    w.step()?;
    let (subdivisions, rows, shifted) = grid;
    if subdivisions == 0
        || subdivisions > sampling::MAX_SUBDIVISIONS
        || rows == 0
        || rows > sampling::MAX_ROWS
    {
        return Err(Error::Limit);
    }
    let f = body.tail.map_or(body.full_length_mm(), |tail| tail.station);
    // Immutable Thread construction proves L/p + 2*starts <=128. This
    // bounds the physical phase span. The shifted sampling chart adds one
    // negative period below; all cells and allocations remain work charged.
    let periods = (f / body.thread.spec().pitch_mm).ceil();
    let maximum = options.triangle_limit.min(131_072);
    let mut owners = Owners::new(maximum, w)?;
    if let Some(trace) = &trace {
        trace.check_memory(owners.memory_bytes())?;
    }
    let mut mesh = DisplayMesh {
        vertices: vec![],
        indices: vec![],
        face_indices: vec![],
    };
    let mut mouths = vec![];
    let mut tops = vec![];
    // At most two rim receipts per retained point. Actual capacity is bounded
    // separately from compact exact-coordinate ownership (64MiB).
    mouths
        .try_reserve_exact(maximum.min(16_384))
        .map_err(|_| Error::Limit)?;
    tops.try_reserve_exact(maximum.min(16_384))
        .map_err(|_| Error::Limit)?;
    let allocations = [
        (mouths.capacity(), std::mem::size_of::<Edge>()),
        (tops.capacity(), std::mem::size_of::<Edge>()),
    ];
    storage::bytes(&allocations, 0, 16 * 1024 * 1024)?;
    let mut chord = 0_f64;
    let mut normal = 0_f64;
    // A positive angular seam shift needs one additional negative phase
    // period to retain the complete z=0 source clipping domain: construction
    // bounds starts by64, hence lead*TURN_OFFSET < pitch.
    let first = -i64::from(body.thread.spec().starts) - i64::from(shifted);
    for period in first..=periods as i64 {
        for band in 0..5 {
            let patch = Patch::with_offset(body, period, band, subdivisions, shifted, w)?;
            for index in 0..subdivisions {
                for row in 0..rows {
                    for kind in [Kind::Cone, Kind::Wall] {
                        if matches!(kind, Kind::Cone) && body.entry.is_none() {
                            continue;
                        }
                        let polygon = patch.cell(index, row, rows, kind, w)?;
                        if polygon.points.len() < 3 {
                            continue;
                        }
                        let field = patch.field(body, &polygon, w)?;
                        let mut ids = Vec::with_capacity(8);
                        let mut errors = Vec::with_capacity(8);
                        for p in &polygon.points {
                            let (key, encoded, rzv) =
                                patch.owned(p, &polygon, body.thread.spec().starts, w)?;
                            let proposed =
                                u32::try_from(mesh.vertices.len() / 3).map_err(|_| Error::Limit)?;
                            let reserved = trace.as_ref().map_or(0, |trace| trace.memory_bytes());
                            let vertex = if reserved == 0 {
                                owners.vertex(key, encoded, proposed, w)?
                            } else {
                                owners.vertex_reserved(key, encoded, proposed, reserved, w)?
                            } as usize;
                            let angle = rzv[2]
                                * I::new(
                                    std::f64::consts::TAU.next_down(),
                                    std::f64::consts::TAU.next_up(),
                                )
                                * I::point(if body.thread.spec().left_hand {
                                    -1.
                                } else {
                                    1.
                                });
                            let point = [rzv[0] * angle.cos(), rzv[0] * angle.sin(), rzv[1]];
                            if vertex == proposed as usize {
                                let literal = point.map(|p| p.lo.midpoint(p.hi));
                                // Checked exact r,z bounds are within1e7 and
                                // v within[0,1+TURN_OFFSET]. Angle intermediates are <7,
                                // trig within[-1,1], products/midpoints <2e7.
                                mesh.vertices.extend(literal);
                            }
                            let literal: [f64; 3] = mesh.vertices[vertex * 3..vertex * 3 + 3]
                                .try_into()
                                .unwrap();
                            let error = point
                                .into_iter()
                                .zip(literal)
                                .fold(I::point(0.), |sum, (p, l)| {
                                    sum + I::point((p - I::point(l)).magnitude())
                                })
                                .hi;
                            if !error.is_finite() || error > options.deflection / 4. {
                                return Err(Error::Unresolved);
                            }
                            ids.push(vertex);
                            errors.push(error);
                        }
                        for cut in 4..6 {
                            let positions: Vec<usize> = polygon
                                .points
                                .iter()
                                .enumerate()
                                .filter_map(|(i, p)| (p.active & (1 << cut) != 0).then_some(i))
                                .collect();
                            if positions.len() < 2 {
                                continue;
                            }
                            // Simple positive convex clipping has at most
                            // two distinct vertices on each physical support.
                            let [a, b] = [positions[0], positions[1]];
                            let edge = Edge {
                                vertices: [ids[a], ids[b]],
                                v: [
                                    polygon.points[a].v.enclosure(w)?.1,
                                    polygon.points[b].v.enclosure(w)?.1,
                                ],
                                slope: patch.slope.enclosure(w)?.1,
                                intercept: patch.intercept.enclosure(w)?.1,
                                kind,
                                errors: [errors[a], errors[b]],
                            };
                            let target = match polygon.stations[cut - 4] {
                                0 => &mut mouths,
                                2 => &mut tops,
                                _ => continue,
                            };
                            if let Some(trace) = trace.as_deref_mut() {
                                for position in [a, b] {
                                    let parameter = &polygon.points[position].v;
                                    let cache_bytes = owners.memory_bytes();
                                    trace.remember(ids[position], parameter, cache_bytes, w)?;
                                }
                            }
                            append_rim(target, edge, maximum.min(16_384))?;
                        }
                        let points: Vec<[f64; 3]> = ids
                            .iter()
                            .map(|&i| mesh.vertices[i * 3..i * 3 + 3].try_into().unwrap())
                            .collect();
                        // Convex exact-parameter polygon: fan triangles cover it
                        // without crossing clipping boundaries. Zero area exact
                        // clipping cells were removed by the polygon cardinality.
                        for i in 1..ids.len() - 1 {
                            let proof = match qualify_attempt(
                                &field,
                                &points,
                                &errors,
                                [0, i, i + 1],
                                options,
                                w,
                            )? {
                                Attempt::Complete(proof) => proof,
                                Attempt::Refine => return Ok(Attempt::Refine),
                                Attempt::Repartition => return Ok(Attempt::Repartition),
                            };
                            let node = Node {
                                band,
                                start: period.rem_euclid(i64::from(body.thread.spec().starts))
                                    as u32,
                            };
                            append_qualified(
                                &mut mesh,
                                &ids,
                                proof,
                                face(body, node, kind),
                                options,
                            )?;
                            chord = chord.max(proof.chord);
                            normal = normal.max(proof.angle);
                        }
                    }
                }
            }
        }
    }
    Ok(Attempt::Complete(Shell {
        mesh,
        chord,
        normal,
        mouth: mouths,
        top: tops,
    }))
}
fn append_rim(target: &mut Vec<Edge>, edge: Edge, maximum: usize) -> Result<(), Error> {
    if target.len() == maximum {
        return Err(Error::Limit);
    }
    target.push(edge);
    Ok(())
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn actual_boundary_receipts_obey_the_declared_rim_count_budget() {
        let body = super::super::tests::source(false, 1, Some(0.4), true);
        let source = shell(
            &body,
            super::super::tests::options(),
            1,
            64,
            &mut Work {
                remaining: 16_000_000,
                used: 0,
                interrupt: None,
            },
        )
        .unwrap()
        .unwrap();
        let mut target = vec![];
        append_rim(&mut target, source.top[0], 2).unwrap();
        append_rim(&mut target, source.top[1], 2).unwrap();
        assert_eq!(append_rim(&mut target, source.top[2], 2), Err(Error::Limit));
        assert_eq!(target.len(), 2);
        assert_eq!(target[0].vertices, source.top[0].vertices);
    }
    #[test]
    fn actual_regular_source_shell_retains_full_facet_certificates() {
        let body = super::super::tests::source(false, 1, Some(0.4), true);
        let options = super::super::tests::options();
        let w = &mut Work {
            remaining: 16_000_000,
            used: 0,
            interrupt: None,
        };
        let s = shell(&body, options, 1, 64, w)
            .unwrap()
            .expect("default actual source shell must qualify");
        assert!(s.chord <= options.deflection && s.normal <= options.angle);
        assert_eq!(s.mesh.face_indices.len(), 4160);
        assert_eq!(s.mouth.len(), 68);
        assert_eq!(s.top.len(), 68);
    }
}
#[cfg(test)]
mod grid_tests {
    use super::*;
    #[test]
    fn bounded_proposal_grids_reject_invalid_counts_before_allocating() {
        let body = super::super::tests::source(false, 1, Some(0.4), true);
        let w = &mut Work {
            remaining: 16_000_000,
            used: 0,
            interrupt: None,
        };
        for (sub, rows) in [(0, 64), (1, 0), (1025, 64), (1, 65537)] {
            assert!(matches!(
                shell(&body, super::super::tests::options(), sub, rows, w),
                Err(Error::Limit)
            ));
        }
    }
}
#[cfg(test)]
mod trace_tests {
    use super::*;
    #[test]
    fn optional_boundary_trace_keeps_default_cavity_coordinates_and_faces_identical() {
        let body = super::super::tests::source(false, 1, Some(0.4), true);
        let options = super::super::tests::options();
        let mut a = Work {
            remaining: 16_000_000,
            used: 0,
            interrupt: None,
        };
        let mut b = Work {
            remaining: 16_000_000,
            used: 0,
            interrupt: None,
        };
        let ordinary = shell(&body, options, 1, 64, &mut a).unwrap().unwrap();
        let (traced, trace) = traced_shell(&body, options, 1, 64, &mut b)
            .unwrap()
            .unwrap();
        assert_eq!(ordinary.mesh, traced.mesh);
        assert_eq!(ordinary.chord, traced.chord);
        assert_eq!(ordinary.normal, traced.normal);
        for edge in traced.mouth.iter().chain(&traced.top) {
            for &vertex in &edge.vertices {
                assert!(trace.bounds(vertex, &mut b).is_ok());
            }
        }
        let result = tessellate_candidate(&body, options, None).unwrap();
        assert_eq!(result.work_used, 3_484_305);
    }
}
