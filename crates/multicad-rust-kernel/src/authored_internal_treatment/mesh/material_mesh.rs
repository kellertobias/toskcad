//! Full actual stock boundary. Tool disks are omitted, rather than hidden
//! inside a disconnected stock cylinder. Every new facet has source receipts.
use super::{
    chart::Shell,
    close::{Disk, Prepared},
    jet::Jet,
    trace::Trace,
    *,
};
use std::cmp::Ordering;
mod radial;
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Link {
    from: usize,
    to: usize,
    angle: usize,
    next_angle: usize,
}
#[derive(Clone, Copy)]
struct Node {
    vertex: usize,
    angle: usize,
    error: f64,
}
fn tau(body: &Cavity) -> I {
    I::new(
        std::f64::consts::TAU.next_down(),
        std::f64::consts::TAU.next_up(),
    ) * I::point(if body.thread.spec().left_hand {
        -1.
    } else {
        1.
    })
}
fn curve(body: &Cavity, v: I, radius: I, z: I) -> source::Field {
    let v =
        Jet::constant(I::point(v.lo)).add(Jet::variable(0).scale(I::point(v.hi) - I::point(v.lo)));
    let theta = v.scale(tau(body));
    let r = Jet::constant(radius);
    source::Field {
        point: [r.mul(theta.cos()), r.mul(theta.sin()), Jet::constant(z)],
        normal: source::Normal::Radial {
            angle: theta.value,
            beta: I::point(0.),
            z: I::point(0.),
        },
    }
}
fn point(body: &Cavity, v: I, radius: I, z: I) -> [I; 3] {
    let theta = v * tau(body);
    [radius * theta.cos(), radius * theta.sin(), z]
}
fn span(edge: chart::Edge) -> I {
    I::new(
        edge.v[0].lo.min(edge.v[1].lo),
        edge.v[0].hi.max(edge.v[1].hi),
    )
}
fn annulus(body: &Cavity, disk: &Disk, radius: I) -> source::Field {
    let z = if disk.end == 0 {
        0.
    } else {
        body.full_length_mm()
    };
    let inner = close::field(body, disk.edge, z, disk.mode);
    let outer = curve(
        body,
        span(disk.edge),
        radius,
        if disk.end == 0 {
            I::point(0.)
        } else {
            body.end_interval()
        },
    );
    let t = Jet::variable(1);
    source::Field {
        point: std::array::from_fn(|i| {
            inner.point[i].add(outer.point[i].sub(inner.point[i]).mul(t))
        }),
        normal: source::Normal::Axial(I::point(if disk.end == 0 { -1. } else { 1. })),
    }
}
fn link(disk: &Disk, trace: &Trace, w: &mut Work<'_>) -> Result<Link, Error> {
    w.step()?;
    let ends = disk
        .edge
        .v
        .map(|v| v.lo == trace.end() && v.hi == trace.end());
    let order = match ends {
        [true, false] => Ordering::Greater,
        [false, true] => Ordering::Less,
        _ => trace.compare(disk.edge.vertices[0], disk.edge.vertices[1], w)?,
    };
    let order = match order {
        Ordering::Less => [0, 1],
        Ordering::Greater => [1, 0],
        Ordering::Equal => return Err(Error::Unresolved),
    };
    Ok(Link {
        from: disk.ids[order[0]],
        to: disk.ids[order[1]],
        angle: disk.edge.vertices[order[0]],
        next_angle: disk.edge.vertices[order[1]],
    })
}
fn ring(
    links: Vec<Link>,
    outer: &[usize],
    errors: &[f64],
    trace: &Trace,
    w: &mut Work<'_>,
) -> Result<Vec<Node>, Error> {
    w.charge(1 + outer.len() + links.len())?;
    if links.len() < 3 || links.len() > 16_384 {
        return Err(Error::Unresolved);
    }
    let mut next = Vec::new();
    next.try_reserve_exact(outer.len())
        .map_err(|_| Error::Limit)?;
    let allocation = [(next.capacity(), std::mem::size_of::<usize>())];
    storage::bytes(&allocation, 0, 4 * 1024 * 1024)?;
    next.resize(outer.len(), usize::MAX);
    let mut start = None;
    for (i, link) in links.iter().enumerate() {
        w.step()?;
        if link.from >= next.len() || link.to >= next.len() || next[link.from] != usize::MAX {
            return Err(Error::Unresolved);
        }
        next[link.from] = i;
        let v = trace.bounds(link.angle, w)?;
        if v.lo == trace.origin() && v.hi == trace.origin() && start.replace(link.from).is_some() {
            return Err(Error::Unresolved);
        }
    }
    let start = start.ok_or(Error::Unresolved)?;
    let mut at = start;
    let mut nodes = Vec::new();
    nodes
        .try_reserve_exact(links.len())
        .map_err(|_| Error::Limit)?;
    let allocations = [
        (next.capacity(), std::mem::size_of::<usize>()),
        (nodes.capacity(), std::mem::size_of::<Node>()),
    ];
    storage::bytes(&allocations, 0, 4 * 1024 * 1024)?;

    for _ in 0..links.len() {
        w.step()?;
        let i = *next.get(at).ok_or(Error::Unresolved)?;
        let edge = links.get(i).ok_or(Error::Unresolved)?;
        let vertex = *outer.get(at).ok_or(Error::Unresolved)?;
        if vertex == usize::MAX {
            return Err(Error::Unresolved);
        }
        nodes.push(Node {
            vertex,
            angle: edge.angle,
            error: errors[vertex],
        });
        at = edge.to;
        if at == start {
            break;
        }
        if trace.compare(edge.angle, edge.next_angle, w)? != Ordering::Less {
            return Err(Error::Unresolved);
        }
    }
    if at != start || nodes.len() != links.len() {
        return Err(Error::Unresolved);
    }
    Ok(nodes)
}
#[derive(Clone, Copy)]
struct StockSurface<'a> {
    body: &'a Cavity,
    radius: I,
    trace: &'a Trace,
}
fn boundary_point(
    source: StockSurface<'_>,
    disk: &Disk,
    i: usize,
    map: &mut [usize],
    errors: &mut Vec<f64>,
    mesh: &mut DisplayMesh,
    w: &mut Work<'_>,
) -> Result<(usize, f64), Error> {
    let StockSurface {
        body,
        radius,
        trace,
    } = source;
    w.step()?;
    let inner = disk.ids[i];
    let exact = point(
        body,
        trace.bounds(disk.edge.vertices[i], w)?,
        radius,
        if disk.end == 0 {
            I::point(0.)
        } else {
            body.end_interval()
        },
    );
    if map[inner] == usize::MAX {
        let id = mesh.vertices.len() / 3;
        mesh.vertices.extend(exact.map(|p| p.lo.midpoint(p.hi)));
        errors.push(0.);
        map[inner] = id;
    }
    let id = map[inner];
    let error = close::error(mesh, id, exact);
    errors[id] = errors[id].max(error);
    Ok((id, error))
}
fn outer_wall(
    source: StockSurface<'_>,
    shell: &mut Shell,
    rings: &[Vec<Node>; 2],
    face: u32,
    options: Tessellation,
    w: &mut Work<'_>,
) -> Result<bool, Error> {
    let StockSurface {
        body,
        radius,
        trace,
    } = source;
    let mut at = [0usize; 2];
    let counts = rings.each_ref().map(Vec::len);
    while at[0] < counts[0] || at[1] < counts[1] {
        w.step()?;
        let next = [at[0] + 1, at[1] + 1];
        let order = if next[0] >= counts[0] && next[1] >= counts[1] {
            Ordering::Equal
        } else if next[0] >= counts[0] {
            Ordering::Greater
        } else if next[1] >= counts[1] {
            Ordering::Less
        } else {
            trace.compare(rings[0][next[0]].angle, rings[1][next[1]].angle, w)?
        };
        let select = match order {
            Ordering::Less => vec![[[0, at[0]], [0, next[0]], [1, at[1]]]],
            Ordering::Greater => vec![[[0, at[0]], [1, next[1]], [1, at[1]]]],
            Ordering::Equal => vec![
                [[0, at[0]], [0, next[0]], [1, next[1]]],
                [[0, at[0]], [1, next[1]], [1, at[1]]],
            ],
        };
        for tri in select {
            let mut nodes = Vec::with_capacity(3);
            let mut ranges = Vec::with_capacity(3);
            for [end, index] in tri {
                let node = rings[end][index % counts[end]];
                let v = if index == counts[end] {
                    I::point(trace.end())
                } else {
                    trace.bounds(node.angle, w)?
                };
                nodes.push(node);
                ranges.push(v);
            }
            let domain = I::new(
                ranges.iter().map(|v| v.lo).fold(f64::INFINITY, f64::min),
                ranges
                    .iter()
                    .map(|v| v.hi)
                    .fold(f64::NEG_INFINITY, f64::max),
            );
            let mut field = curve(body, domain, radius, I::point(0.));
            field.point[2] = Jet::variable(1).scale(body.end_interval());
            let ids = nodes.iter().map(|n| n.vertex).collect::<Vec<_>>();
            let points = ids
                .iter()
                .map(|&i| shell.mesh.vertices[i * 3..i * 3 + 3].try_into().unwrap())
                .collect::<Vec<_>>();
            let errors = nodes.iter().map(|n| n.error).collect::<Vec<_>>();
            let Some(proof) = qualify(&field, &points, &errors, [0, 1, 2], options, w)? else {
                return Ok(false);
            };
            append_qualified(&mut shell.mesh, &ids, proof, face, options)?;
            shell.chord = shell.chord.max(proof.chord);
            shell.normal = shell.normal.max(proof.angle);
        }
        match order {
            Ordering::Less => at[0] += 1,
            Ordering::Greater => at[1] += 1,
            Ordering::Equal => {
                at[0] += 1;
                at[1] += 1;
            }
        }
    }
    Ok(true)
}
fn compact(mesh: &mut DisplayMesh, w: &mut Work<'_>) -> Result<(), Error> {
    let count = mesh.vertices.len() / 3;
    w.charge(count + mesh.indices.len())?;
    let mut used = Vec::new();
    used.try_reserve_exact(count).map_err(|_| Error::Limit)?;
    used.resize(count, false);
    for (i, &index) in mesh.indices.iter().enumerate() {
        if i % 1024 == 0 {
            w.step()?;
        }
        used[index as usize] = true;
    }
    let mut map = Vec::new();
    map.try_reserve_exact(count).map_err(|_| Error::Limit)?;
    map.resize(count, usize::MAX);
    let mut vertices = Vec::new();
    vertices
        .try_reserve_exact(3 * count)
        .map_err(|_| Error::Limit)?;
    let allocations = [
        (used.capacity(), std::mem::size_of::<bool>()),
        (map.capacity(), std::mem::size_of::<usize>()),
        (vertices.capacity(), std::mem::size_of::<f64>()),
    ];
    storage::bytes(&allocations, 0, 16 * 1024 * 1024)?;

    for (i, &active) in used.iter().enumerate() {
        if i % 256 == 0 {
            w.step()?;
        }
        if active {
            map[i] = vertices.len() / 3;
            vertices.extend_from_slice(&mesh.vertices[i * 3..i * 3 + 3]);
        }
    }
    for (i, index) in mesh.indices.iter_mut().enumerate() {
        if i % 1024 == 0 {
            w.step()?;
        }
        *index = map[*index as usize] as u32;
    }
    w.step()?;
    mesh.vertices = vertices;
    Ok(())
}
fn assemble(
    material: &Material,
    prepared: Prepared,
    trace: Trace,
    options: Tessellation,
    w: &mut Work<'_>,
) -> Result<Option<Candidate>, Error> {
    let Prepared { mut shell, caps } = prepared;
    w.charge(1 + shell.mesh.indices.len() / 3)?;
    for (i, triangle) in shell
        .mesh
        .indices
        .as_chunks_mut::<3>()
        .0
        .iter_mut()
        .enumerate()
    {
        if i % 1024 == 0 {
            w.step()?;
        }
        triangle.swap(1, 2);
    }
    let baseline = (
        shell.mesh.vertices.len(),
        shell.mesh.indices.len(),
        shell.mesh.face_indices.len(),
        shell.chord,
        shell.normal,
    );
    let mut strips = 1;
    loop {
        match stock_boundary(material, &mut shell, &caps, &trace, options, strips, w)? {
            radial::Verdict::Done => break,
            radial::Verdict::Radial if strips < radial::MAX_STRIPS => {
                w.step()?;
                shell.mesh.vertices.truncate(baseline.0);
                shell.mesh.indices.truncate(baseline.1);
                shell.mesh.face_indices.truncate(baseline.2);
                shell.chord = baseline.3;
                shell.normal = baseline.4;
                strips *= 2;
            }
            _ => return Ok(None),
        }
    }
    compact(&mut shell.mesh, w)?;
    Ok(Some(Candidate {
        mesh: shell.mesh,
        source_chord_bound_mm: shell.chord,
        source_normal_angle_bound_rad: shell.normal,
        work_used: w.used,
    }))
}
fn stock_boundary(
    material: &Material,
    shell: &mut Shell,
    caps: &[Disk],
    trace: &Trace,
    options: Tessellation,
    strips: usize,
    w: &mut Work<'_>,
) -> Result<radial::Verdict, Error> {
    let body = material.cavity();
    let radius = material.radius_interval();
    let count = shell.mesh.vertices.len() / 3;
    w.charge(count * 2 + caps.len())?;
    let mut outer = Vec::new();
    outer.try_reserve_exact(count).map_err(|_| Error::Limit)?;
    outer.resize(count, usize::MAX);
    let mut errors = Vec::new();
    errors
        .try_reserve_exact(count + 32_768)
        .map_err(|_| Error::Limit)?;
    errors.resize(count, 0.);
    let allocations = [
        (outer.capacity(), std::mem::size_of::<usize>()),
        (errors.capacity(), std::mem::size_of::<f64>()),
    ];
    storage::bytes(&allocations, 0, 8 * 1024 * 1024)?;
    let mut links: [Vec<Link>; 2] = std::array::from_fn(|_| Vec::new());
    for part in &mut links {
        part.try_reserve_exact(16_384).map_err(|_| Error::Limit)?;
    }
    let allocations = links
        .each_ref()
        .map(|part| (part.capacity(), std::mem::size_of::<Link>()));
    storage::bytes(&allocations, 0, 4 * 1024 * 1024)?;

    let mut spoke =
        radial::Spokes::new(count, caps.len(), strips, outer.capacity(), &mut errors, w)?;
    for &disk in caps {
        let [a, b] = [0, 1].map(|i| {
            boundary_point(
                StockSurface {
                    body,
                    radius,
                    trace,
                },
                &disk,
                i,
                &mut outer,
                &mut errors,
                &mut shell.mesh,
                w,
            )
        });
        let (a, ae) = a?;
        let (b, be) = b?;
        match radial::emit(
            &mut spoke,
            shell,
            StockSurface {
                body,
                radius,
                trace,
            },
            &disk,
            ([a, b], [ae, be]),
            options,
            w,
        )? {
            radial::Verdict::Done => {}
            other => return Ok(other),
        }
        w.charge(shell.mesh.vertices.len() / 3 - errors.len())?;
        errors.resize(shell.mesh.vertices.len() / 3, 0.);
        links[disk.end as usize].push(link(&disk, trace, w)?);
    }
    let rings = [
        ring(std::mem::take(&mut links[0]), &outer, &errors, trace, w)?,
        ring(std::mem::take(&mut links[1]), &outer, &errors, trace, w)?,
    ];
    let face = 2
        + 4 * body.thread.spec().starts
        + u32::from(body.entry.is_some())
        + 2 * u32::from(body.tail.is_some());
    if !outer_wall(
        StockSurface {
            body,
            radius,
            trace,
        },
        shell,
        &rings,
        face,
        options,
        w,
    )? {
        return Ok(radial::Verdict::Refine);
    }
    Ok(radial::Verdict::Done)
}
pub(super) fn tessellate(
    material: &Material,
    options: Tessellation,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<Candidate, Error> {
    options.check().map_err(|_| Error::Invalid)?;
    let body = material.cavity();
    let mut w = Work {
        remaining: 16_000_000,
        used: 0,
        interrupt,
    };
    w.step()?;
    if body
        .entry
        .is_some_and(|e| e.radius <= body.thread.maximum_radius_mm())
    {
        return Err(Error::Unresolved);
    }
    let mut grid = Grid::new();
    // A terminal station can lie arbitrarily close to a periodic knot at the
    // angular seam. Start with the existing exact shifted chart so that this
    // coincidence does not force a tiny literal seam edge. This changes only
    // the sampling chart: every cell still certifies the original source.
    grid.shifted = body.tail.is_some_and(|tail| {
        let phase = tail.station / body.thread.spec().pitch_mm;
        let knot = phase.round();
        // A conservative conditioning band only proposes another sampling
        // chart; it never proves geometric equality or separation. Keep the
        // existing exactly integral and ordinary fractional proposals.
        let distance = (phase - knot).abs();
        distance > 0. && distance <= 16. * f64::EPSILON * phase.abs().max(1.)
    });
    let mut angular_trial = body.tail.is_some();
    loop {
        let repartition = match chart::traced_grid_attempt(body, options, &grid, &mut w)? {
            Attempt::Complete((shell, trace)) => {
                if let Some(prepared) = close::prepare(shell, body, options, &mut w)? {
                    if let Some(candidate) = assemble(material, prepared, trace, options, &mut w)? {
                        return Ok(candidate);
                    }
                    if body.tail.is_some() && !grid.shifted {
                        // The two actual cap traces can produce distinct knots
                        // arbitrarily close to a dyadic angular grid vertex.
                        // Doubling that grid retains the ill-conditioned stock
                        // wall facet. Propose the existing exact shifted chart
                        // first; all failed work stays on this same ledger and
                        // every new source/closure/stock facet is re-certified.
                        grid.shifted = true;
                        continue;
                    }
                }
                false
            }
            Attempt::Refine => false,
            Attempt::Repartition => true,
        };
        if !repartition && angular_trial {
            // Angular curvature can require more rows while the same complete
            // axial phase band already suffices. Try that partition first;
            // failed work remains charged to this request's original ledger.
            angular_trial = false;
            grid.refine_angular()?;
        } else {
            grid.advance(repartition)?;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        Circle, Rational,
        authored_flat_thread::{AuthoredThread, SourceOptions},
        cylinder::AxialCylinder,
        helical_thread::Spec,
    };
    fn source(left_hand: bool, starts: u32, entry: bool, tail: bool) -> Material {
        let thread = AuthoredThread::new(
            Spec {
                major_radius_mm: 1.5,
                pitch_mm: 0.5,
                length_mm: 4.,
                depth_mm: 0.25,
                crest_width_mm: 0.125,
                lower_flank_deg: 30.,
                upper_flank_deg: 30.,
                starts,
                left_hand,
                root_radius_mm: None,
            },
            SourceOptions {
                position_limit_mm: 1e-7,
                normal_angle_limit_rad: 1e-7,
            },
        )
        .unwrap();
        let cavity = Cavity::new(
            thread,
            Intent {
                entry_depth_mm: entry.then_some(0.5),
                terminal: tail.then_some(Terminal {
                    width_mm: 0.5,
                    radial_depth_mm: 0.125,
                }),
            },
        )
        .unwrap();
        let zero = Rational::new(0, 1).unwrap();
        let stock = AxialCylinder::new(
            7,
            Circle::new([zero; 2], Rational::new(5, 2).unwrap()).unwrap(),
            zero,
            Rational::new(4, 1).unwrap(),
        )
        .unwrap();
        Material::new(stock, cavity, 0.5).unwrap()
    }
    #[test]
    fn actual_stock_minus_cavity_is_one_closed_boundary_with_all_source_facet_receipts() {
        let options = Tessellation {
            deflection: 0.1,
            angle: 0.3,
            triangle_limit: 1_000_000,
        };
        for (left, starts, entry, tail) in [
            (false, 1, true, true),
            (true, 2, true, true),
            (false, 1, false, true),
            (true, 1, true, false),
        ] {
            let material = source(left, starts, entry, tail);
            let timer = std::time::Instant::now();
            let result = material.tessellate_candidate(options, None).unwrap();
            let proof = crate::manufacturing::certify(&result.mesh).unwrap();
            let mass = material.volume_bounds(QueryOptions::default()).unwrap();
            eprintln!(
                "actualStockMFG {left}/{starts}/{entry}/{tail}: {}facets, {}work, {:?}, meshVolume{}, sourceMass[{:.9},{:.9}]",
                result.mesh.face_indices.len(),
                result.work_used,
                timer.elapsed(),
                proof.signed_volume_mm3,
                mass.lower,
                mass.upper
            );
            assert!(proof.edge_pairs > 0 && proof.vertex_links > 0 && proof.signed_volume_mm3 > 0.);
            assert!(
                result.source_chord_bound_mm <= options.deflection
                    && result.source_normal_angle_bound_rad <= options.angle
            );
            // Mesh signed volume is diagnostic; the independently certified
            // source mass interval remains separate, with no heuristic equality.
            assert!(mass.lower > 0. && mass.upper - mass.lower <= 1e-7);
            assert!(result.work_used <= 16_000_000);
        }
    }
    #[test]
    fn unequal_end_angular_grids_sew_by_exact_source_parameter_order() {
        let original = source(false, 1, true, true);
        let cavity = Cavity::new(
            *original.cavity().source_thread(),
            Intent {
                entry_depth_mm: Some(0.5),
                terminal: Some(Terminal {
                    width_mm: 0.3,
                    radial_depth_mm: 0.125,
                }),
            },
        )
        .unwrap();
        let material = Material::new(*original.stock(), cavity, 0.5).unwrap();
        let result = material
            .tessellate_candidate(
                Tessellation {
                    deflection: 0.1,
                    angle: 0.3,
                    triangle_limit: 1_000_000,
                },
                None,
            )
            .unwrap();
        assert!(
            crate::manufacturing::certify(&result.mesh)
                .unwrap()
                .signed_volume_mm3
                > 0.
        );
    }
    #[test]
    fn stock_mesh_cancellation_resource_and_quality_failures_preserve_actual_source() {
        let material = source(false, 1, true, true);
        let before = material;
        let options = Tessellation {
            deflection: 0.1,
            angle: 0.3,
            triangle_limit: 1_000_000,
        };
        assert_eq!(
            material
                .tessellate_candidate(
                    Tessellation {
                        deflection: 0.,
                        ..options
                    },
                    None
                )
                .unwrap_err(),
            Error::Invalid
        );
        assert_eq!(
            material
                .tessellate_candidate(
                    Tessellation {
                        triangle_limit: 1000,
                        ..options
                    },
                    None
                )
                .unwrap_err(),
            Error::Limit
        );
        assert_eq!(
            material
                .tessellate_candidate(
                    Tessellation {
                        deflection: 0.005,
                        ..options
                    },
                    None
                )
                .unwrap_err(),
            Error::Limit
        );
        let small = Cavity::new(
            *material.cavity().source_thread(),
            Intent {
                entry_depth_mm: Some(0.01),
                terminal: Some(Terminal {
                    width_mm: 0.5,
                    radial_depth_mm: 0.125,
                }),
            },
        )
        .unwrap();
        let small = Material::new(*material.stock(), small, 0.5).unwrap();
        assert_eq!(
            small.tessellate_candidate(options, None).unwrap_err(),
            Error::Unresolved
        );
        let calls = std::cell::Cell::new(0);
        let count = || {
            calls.set(calls.get() + 1);
            false
        };
        let expected = material
            .tessellate_candidate(options, Some(&count))
            .unwrap();
        let total = calls.get();
        for stop in [0, 100, 10_000, total - 1] {
            let calls = std::cell::Cell::new(0);
            let cancel = || {
                calls.set(calls.get() + 1);
                calls.get() > stop
            };
            assert_eq!(
                material
                    .tessellate_candidate(options, Some(&cancel))
                    .unwrap_err(),
                Error::Cancelled
            );
        }
        assert_eq!(material, before);
        assert_eq!(
            material.tessellate_candidate(options, None).unwrap().mesh,
            expected.mesh
        );
    }
    #[test]
    fn actual_closure_receipts_refuse_a_stricter_outer_wall_target() {
        let material = source(false, 1, true, true);
        let options = Tessellation {
            deflection: 0.1,
            angle: 0.3,
            triangle_limit: 1_000_000,
        };
        for target in [
            Tessellation {
                angle: 0.01,
                ..options
            },
            Tessellation {
                deflection: 1e-12,
                ..options
            },
        ] {
            let w = &mut Work {
                remaining: 16_000_000,
                used: 0,
                interrupt: None,
            };
            let (shell, trace) = chart::traced_shell(material.cavity(), options, 1, 64, w)
                .unwrap()
                .unwrap();
            let prepared = close::prepare(shell, material.cavity(), options, w)
                .unwrap()
                .unwrap();
            assert!(
                assemble(&material, prepared, trace, target, w)
                    .unwrap()
                    .is_none()
            );
        }
    }
    #[test]
    fn tiny_crest_unequal_grid_resource_refusal_does_not_drop_source_angles() {
        let original = source(false, 1, true, true);
        let mut spec = original.cavity().source_thread().spec();
        spec.crest_width_mm = 0.0001;
        let thread = AuthoredThread::new(
            spec,
            SourceOptions {
                position_limit_mm: 1e-7,
                normal_angle_limit_rad: 1e-7,
            },
        )
        .unwrap();
        let cavity = Cavity::new(
            thread,
            Intent {
                entry_depth_mm: Some(0.5),
                terminal: Some(Terminal {
                    width_mm: 0.3,
                    radial_depth_mm: 0.125,
                }),
            },
        )
        .unwrap();
        let material = Material::new(*original.stock(), cavity, 0.5).unwrap();
        assert_eq!(
            material
                .tessellate_candidate(
                    Tessellation {
                        deflection: 0.1,
                        angle: 0.3,
                        triangle_limit: 1_000_000
                    },
                    None
                )
                .unwrap_err(),
            Error::Limit
        );
        assert_eq!(
            material.cavity().source_thread().spec().crest_width_mm,
            0.0001
        );
    }
    #[test]
    fn actual_source_ring_receipts_reject_missing_mixed_and_misordered_boundary_incidence() {
        let material = source(false, 1, true, true);
        let options = Tessellation {
            deflection: 0.1,
            angle: 0.3,
            triangle_limit: 1_000_000,
        };
        let w = &mut Work {
            remaining: 16_000_000,
            used: 0,
            interrupt: None,
        };
        let (shell, trace) = chart::traced_shell(material.cavity(), options, 1, 64, w)
            .unwrap()
            .unwrap();
        let prepared = close::prepare(shell, material.cavity(), options, w)
            .unwrap()
            .unwrap();
        let count = prepared.shell.mesh.vertices.len() / 3;
        // This validator is topology-only: identity maps keep ACTUAL source
        // cap vertices. No stock-radius or surface-quality claim comes from it.
        let vertices: Vec<usize> = (0..count).collect();
        let errors = vec![0.; count];
        let links = prepared
            .caps
            .iter()
            .filter(|d| d.end == 0)
            .map(|d| link(d, &trace, w).unwrap())
            .collect::<Vec<_>>();
        let valid = ring(links.clone(), &vertices, &errors, &trace, w).unwrap();
        assert_eq!(valid.len(), links.len());
        assert!(matches!(
            ring(vec![], &vertices, &errors, &trace, w),
            Err(Error::Unresolved)
        ));
        let mut duplicate = links.clone();
        duplicate.push(links[0]);
        assert!(matches!(
            ring(duplicate, &vertices, &errors, &trace, w),
            Err(Error::Unresolved)
        ));
        let mut bad = links.clone();
        bad[0].from = count;
        assert!(matches!(
            ring(bad, &vertices, &errors, &trace, w),
            Err(Error::Unresolved)
        ));
        let top = prepared
            .caps
            .iter()
            .filter(|d| d.end == 1)
            .map(|d| link(d, &trace, w).unwrap())
            .collect::<Vec<_>>()
            .into_iter()
            .find(|l| {
                let b = trace.bounds(l.angle, w).unwrap();
                b.lo == 0. && b.hi == 0.
            })
            .unwrap();
        let mut mixed = links.clone();
        mixed.push(top);
        assert!(matches!(
            ring(mixed, &vertices, &errors, &trace, w),
            Err(Error::Unresolved)
        ));
        let start = valid[0].angle;
        let first = links.iter().position(|l| l.angle == start).unwrap();
        let mut missing = vertices.clone();
        missing[links[first].from] = usize::MAX;
        assert!(matches!(
            ring(links.clone(), &missing, &errors, &trace, w),
            Err(Error::Unresolved)
        ));
        let mut unordered = links.clone();
        unordered[first].next_angle = unordered[first].angle;
        assert!(matches!(
            ring(unordered, &vertices, &errors, &trace, w),
            Err(Error::Unresolved)
        ));
        let mut hole = links.clone();
        hole[first].to = hole[first].from;
        assert!(matches!(
            ring(hole, &vertices, &errors, &trace, w),
            Err(Error::Unresolved)
        ));
        let mut closing = *prepared
            .caps
            .iter()
            .filter(|d| d.end == 0)
            .find(|d| {
                d.edge
                    .v
                    .iter()
                    .any(|v| v.lo == trace.end() && v.hi == trace.end())
            })
            .unwrap();
        let expected = link(&closing, &trace, w).unwrap();
        closing.edge.v.swap(0, 1);
        closing.edge.vertices.swap(0, 1);
        closing.ids.swap(0, 1);
        closing.errors.swap(0, 1);
        assert_eq!(link(&closing, &trace, w).unwrap(), expected);
        closing.edge.v[1] = closing.edge.v[0];
        closing.edge.vertices[1] = closing.edge.vertices[0];
        closing.ids[1] = closing.ids[0];
        assert!(matches!(link(&closing, &trace, w), Err(Error::Unresolved)));
    }
    #[test]
    fn actual_stock_cylinder_accepts_either_end_ring_with_additional_angular_nodes() {
        let material = source(false, 1, true, true);
        let body = material.cavity();
        let radius = material.radius_interval();
        let options = Tessellation {
            deflection: 0.1,
            angle: 0.3,
            triangle_limit: 1_000_000,
        };
        for extra_end in [0, 1] {
            let w = &mut Work {
                remaining: 16_000_000,
                used: 0,
                interrupt: None,
            };
            let mut trace = Trace::new(130, w).unwrap();
            let mut shell = Shell {
                mesh: DisplayMesh {
                    vertices: vec![],
                    indices: vec![],
                    face_indices: vec![],
                },
                chord: 0.,
                normal: 0.,
                mouth: vec![],
                top: vec![],
            };
            let mut rings: [Vec<Node>; 2] = std::array::from_fn(|_| Vec::new());
            for (end, nodes) in rings.iter_mut().enumerate() {
                let mut angles = (0..64).map(|i| i as f64 / 64.).collect::<Vec<_>>();
                if end == extra_end {
                    angles.push(0.99);
                }
                for v in angles {
                    let id = shell.mesh.vertices.len() / 3;
                    let value = exact::Q::float(v, w).unwrap();
                    trace.remember(id, &value, 0, w).unwrap();
                    let source = point(
                        body,
                        I::point(v),
                        radius,
                        I::point(end as f64 * body.full_length_mm()),
                    );
                    shell
                        .mesh
                        .vertices
                        .extend(source.map(|a| a.lo.midpoint(a.hi)));
                    let error = close::error(&shell.mesh, id, source);
                    nodes.push(Node {
                        vertex: id,
                        angle: id,
                        error,
                    });
                }
            }
            assert!(
                outer_wall(
                    StockSurface {
                        body,
                        radius,
                        trace: &trace
                    },
                    &mut shell,
                    &rings,
                    10,
                    options,
                    w
                )
                .unwrap()
            );
            assert_eq!(shell.mesh.face_indices.len(), 129);
            assert!(shell.chord <= options.deflection && shell.normal <= options.angle);
        }
    }
    #[test]
    fn cancellation_inside_actual_annulus_proofs_never_publishes_partial_stock() {
        let material = source(false, 1, true, true);
        let original = material;
        let options = Tessellation {
            deflection: 0.1,
            angle: 0.3,
            triangle_limit: 1_000_000,
        };
        for stop in 10..15 {
            let mut setup = Work {
                remaining: 16_000_000,
                used: 0,
                interrupt: None,
            };
            let (shell, trace) = chart::traced_shell(material.cavity(), options, 1, 64, &mut setup)
                .unwrap()
                .unwrap();
            let prepared = close::prepare(shell, material.cavity(), options, &mut setup)
                .unwrap()
                .unwrap();
            let calls = std::cell::Cell::new(0);
            let cancel = || {
                calls.set(calls.get() + 1);
                calls.get() > stop
            };
            assert!(matches!(
                assemble(
                    &material,
                    prepared,
                    trace,
                    options,
                    &mut Work {
                        remaining: 16_000_000,
                        used: 0,
                        interrupt: Some(&cancel)
                    }
                ),
                Err(Error::Cancelled)
            ));
        }
        assert_eq!(material, original);
    }
    #[test]
    fn cancellation_in_outer_wall_and_final_compaction_never_publishes_stock() {
        let material = source(false, 1, true, true);
        let options = Tessellation {
            deflection: 0.1,
            angle: 0.3,
            triangle_limit: 1_000_000,
        };
        let setup = |subdivisions, rows| {
            let mut work = Work {
                remaining: 16_000_000,
                used: 0,
                interrupt: None,
            };
            let (shell, trace) =
                chart::traced_shell(material.cavity(), options, subdivisions, rows, &mut work)
                    .unwrap()?;
            let prepared = close::prepare(shell, material.cavity(), options, &mut work).unwrap()?;
            Some((prepared, trace))
        };
        let calls = std::cell::Cell::new(0usize);
        let count = || {
            calls.set(calls.get() + 1);
            false
        };
        let (subdivisions, rows) = (0..8)
            .find_map(|level| {
                let subdivisions = 1usize << level;
                let rows = 64u32 << level;
                let (prepared, trace) = setup(subdivisions, rows)?;
                calls.set(0);
                assemble(
                    &material,
                    prepared,
                    trace,
                    options,
                    &mut Work {
                        remaining: 16_000_000,
                        used: 0,
                        interrupt: Some(&count),
                    },
                )
                .unwrap()
                .map(|_| (subdivisions, rows))
            })
            .expect("actual admitted stock grid");
        let total = calls.get();
        assert!(total > 64);
        for stop in (total - 64)..total {
            let (prepared, trace) = setup(subdivisions, rows).unwrap();
            calls.set(0);
            let cancel = || {
                calls.set(calls.get() + 1);
                calls.get() > stop
            };
            assert!(matches!(
                assemble(
                    &material,
                    prepared,
                    trace,
                    options,
                    &mut Work {
                        remaining: 16_000_000,
                        used: 0,
                        interrupt: Some(&cancel)
                    }
                ),
                Err(Error::Cancelled)
            ));
        }
        assert_eq!(material, source(false, 1, true, true));
    }
    #[test]
    fn nonzero_exact_stock_origin_is_retained_while_candidates_are_stock_relative() {
        let original = source(false, 1, true, true);
        let origin = [
            Rational::new(33, 10).unwrap(),
            Rational::new(-7, 3).unwrap(),
            Rational::new(17, 10).unwrap(),
        ];
        let stock = AxialCylinder::new(
            9,
            Circle::new([origin[0], origin[1]], original.stock().circle().radius()).unwrap(),
            origin[2],
            original.stock().height(),
        )
        .unwrap();
        let moved = Material::new(stock, *original.cavity(), original.minimum_wall_mm()).unwrap();
        assert_eq!(moved.stock_chart_origin(), origin);
        assert_eq!(moved.stock(), &stock);
        assert_eq!(
            moved.volume_bounds(QueryOptions::default()).unwrap(),
            original.volume_bounds(QueryOptions::default()).unwrap()
        );
        let options = Tessellation {
            deflection: 0.1,
            angle: 0.3,
            triangle_limit: 1_000_000,
        };
        let local = moved.tessellate_candidate(options, None).unwrap();
        let baseline = original.tessellate_candidate(options, None).unwrap();
        assert_eq!(local.mesh, baseline.mesh);
        let points = local.mesh.vertices.as_chunks::<3>().0;
        assert_eq!(
            points.iter().map(|p| p[2]).fold(f64::INFINITY, f64::min),
            0.
        );
        assert_eq!(
            points
                .iter()
                .map(|p| p[2])
                .fold(f64::NEG_INFINITY, f64::max),
            4.
        );
        assert!(points.iter().any(|p| p[0] < 0.));
        assert!(origin[0].as_f64() - stock.circle().radius().as_f64() > 0.);
        assert_eq!(moved.stock_chart_origin(), origin);
    }
}

#[cfg(test)]
mod partition_closure_tests {
    use super::*;
    pub(super) fn material(radius: i64) -> Material {
        use crate::{
            Circle, Rational,
            authored_flat_thread::{AuthoredThread, SourceOptions},
            cylinder::AxialCylinder,
        };
        let thread = AuthoredThread::new(
            crate::helical_thread::Spec {
                major_radius_mm: 1.5,
                pitch_mm: 0.5,
                length_mm: 4.,
                depth_mm: 0.25,
                crest_width_mm: 0.125,
                lower_flank_deg: 30.,
                upper_flank_deg: 30.,
                starts: 1,
                left_hand: false,
                root_radius_mm: None,
            },
            SourceOptions {
                position_limit_mm: 1e-7,
                normal_angle_limit_rad: 1e-7,
            },
        )
        .unwrap();
        let cavity = Cavity::new(
            thread,
            Intent {
                entry_depth_mm: Some(0.5),
                terminal: Some(Terminal {
                    width_mm: 0.5,
                    radial_depth_mm: 0.125,
                }),
            },
        )
        .unwrap();
        let zero = Rational::new(0, 1).unwrap();
        let stock = AxialCylinder::new(
            7,
            Circle::new([zero; 2], Rational::new(radius, 1).unwrap()).unwrap(),
            zero,
            Rational::new(4, 1).unwrap(),
        )
        .unwrap();
        Material::new(stock, cavity, 0.5).unwrap()
    }
    #[test]
    fn actual_larger_stock_annulus_refines_after_qualified_inner_sources() {
        let material = material(8);
        let cavity = *material.cavity();
        let options = Tessellation {
            deflection: 0.1,
            angle: 0.3,
            triangle_limit: 1_000_000,
        };
        let mut w = Work {
            remaining: 16_000_000,
            used: 0,
            interrupt: None,
        };
        let (shell, trace) = chart::traced_shell(&cavity, options, 1, 64, &mut w)
            .unwrap()
            .unwrap();
        let prepared = close::prepare(shell, &cavity, options, &mut w)
            .unwrap()
            .unwrap();
        let direct = assemble(&material, prepared, trace, options, &mut w)
            .unwrap()
            .unwrap();
        let before = material;
        let candidate = material.tessellate_candidate(options, None).unwrap();
        assert_eq!(direct.mesh, candidate.mesh);
        assert!(candidate.source_chord_bound_mm <= options.deflection);
        assert!(candidate.source_normal_angle_bound_rad <= options.angle);
        assert!(candidate.work_used < 16_000_000);

        assert_eq!(material, before);
    }
    #[test]
    fn local_annulus_retries_cancel_atomically_and_keep_real_limits() {
        let actual = material(8);
        let body = actual.cavity();
        let options = Tessellation {
            deflection: 0.1,
            angle: 0.3,
            triangle_limit: 1_000_000,
        };
        let calls = std::cell::Cell::new(0usize);
        let count = || {
            calls.set(calls.get() + 1);
            false
        };
        let mut w = Work {
            remaining: 16_000_000,
            used: 0,
            interrupt: Some(&count),
        };
        let (shell, trace) = chart::traced_shell(body, options, 1, 64, &mut w)
            .unwrap()
            .unwrap();
        let mut prepared = close::prepare(shell, body, options, &mut w)
            .unwrap()
            .unwrap();
        let before_assembly = calls.get();
        assert!(matches!(
            stock_boundary(
                &actual,
                &mut prepared.shell,
                &prepared.caps,
                &trace,
                options,
                1,
                &mut w
            )
            .unwrap(),
            radial::Verdict::Radial
        ));
        let after_rejection = calls.get();
        calls.set(0);
        let complete = actual.tessellate_candidate(options, Some(&count)).unwrap();
        let total = calls.get();
        for stop in [
            0,
            before_assembly,
            after_rejection - 1,
            after_rejection,
            after_rejection + 1,
            total - 1,
        ] {
            calls.set(0);
            let cancel = || {
                calls.set(calls.get() + 1);
                calls.get() > stop
            };
            assert!(matches!(
                actual.tessellate_candidate(options, Some(&cancel)),
                Err(Error::Cancelled)
            ));
        }
        assert_eq!(
            actual.tessellate_candidate(options, None).unwrap().mesh,
            complete.mesh
        );
        let limited = Tessellation {
            triangle_limit: prepared.shell.mesh.face_indices.len() + 2,
            ..options
        };
        assert!(matches!(
            actual.tessellate_candidate(limited, None),
            Err(Error::Limit)
        ));
        let large = material(1_000_000);
        assert!(matches!(
            large.tessellate_candidate(options, None),
            Err(Error::Limit)
        ));
    }
    #[test]
    fn a_genuine_positive_single_ulp_wall_is_not_erased_to_publish_an_annulus() {
        use crate::{Circle, Rational, cylinder::AxialCylinder};
        let base = material(8);
        let zero = Rational::new(0, 1).unwrap();
        let den = 1_i64 << 49;
        let radius = Rational::new(8 * den + 1, den).unwrap();
        let stock = AxialCylinder::new(
            7,
            Circle::new([zero; 2], radius).unwrap(),
            zero,
            Rational::new(8, 1).unwrap(),
        )
        .unwrap();
        let thread = crate::authored_flat_thread::AuthoredThread::new(
            crate::helical_thread::Spec {
                length_mm: 8.,
                ..base.cavity().thread.spec()
            },
            crate::authored_flat_thread::SourceOptions {
                position_limit_mm: 1e-7,
                normal_angle_limit_rad: 1e-7,
            },
        )
        .unwrap();
        let cavity = Cavity::new(
            thread,
            Intent {
                entry_depth_mm: Some(6.75),
                ..base.cavity().intent()
            },
        )
        .unwrap();
        let actual = Material::new(stock, cavity, 1e-18).unwrap();
        let options = Tessellation {
            deflection: 0.1,
            angle: 0.3,
            triangle_limit: 1_000_000,
        };
        let mut w = Work {
            remaining: 16_000_000,
            used: 0,
            interrupt: None,
        };
        let (shell, trace) = chart::traced_shell(actual.cavity(), options, 1, 64, &mut w)
            .unwrap()
            .unwrap();
        let mut prepared = close::prepare(shell, actual.cavity(), options, &mut w)
            .unwrap()
            .unwrap();
        assert!(matches!(
            stock_boundary(
                &actual,
                &mut prepared.shell,
                &prepared.caps,
                &trace,
                options,
                4,
                &mut w
            )
            .unwrap(),
            radial::Verdict::Refine
        ));
        assert_eq!(actual.stock().circle().radius(), radius);
    }
}

#[cfg(test)]
mod conditioning_tests;
