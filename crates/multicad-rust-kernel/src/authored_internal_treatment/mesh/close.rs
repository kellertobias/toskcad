//! End disks and actual expanded-terminal annulus/cylinder, sharing source rim
//! indices. No floating nearest-point or tolerance-based sewing is used.
use super::{
    chart::{Edge, Shell},
    jet::Jet,
    *,
};
fn span(edge: Edge) -> I {
    I::new(
        edge.v[0].lo.min(edge.v[1].lo),
        edge.v[0].hi.max(edge.v[1].hi),
    )
}
fn angle(body: &Cavity, v: Jet) -> Jet {
    v.scale(
        I::new(
            std::f64::consts::TAU.next_down(),
            std::f64::consts::TAU.next_up(),
        ) * I::point(if body.thread.spec().left_hand {
            -1.
        } else {
            1.
        }),
    )
}
#[derive(Clone, Copy)]
enum Surface {
    InnerCap,
    Shoulder,
    Terminal,
    OuterCap,
}
fn field(body: &Cavity, edge: Edge, station: f64, mode: Surface) -> source::Field {
    let bound = span(edge);
    let v = Jet::constant(I::point(bound.lo))
        .add(Jet::variable(0).scale(I::point(bound.hi) - I::point(bound.lo)));
    let z = Jet::constant(I::point(station));
    let rho = if matches!(edge.kind, Kind::Cone) {
        Jet::constant(I::point(body.entry.unwrap().radius))
    } else {
        z.sub(v.scale(
            I::point(body.thread.spec().pitch_mm) * I::point(f64::from(body.thread.spec().starts)),
        ))
        .scale(edge.slope)
        .add(Jet::constant(edge.intercept))
    };
    let (r, z, normal) = match mode {
        Surface::InnerCap => (
            rho,
            z,
            source::Normal::Axial(I::point(if station == 0. { -1. } else { 1. })),
        ),
        Surface::Shoulder => {
            let t = Jet::variable(1);
            let outer = Jet::constant(I::point(body.tail.unwrap().radius));
            (
                rho.add(outer.sub(rho).mul(t)),
                z,
                source::Normal::Axial(I::point(-1.)),
            )
        }
        Surface::Terminal => {
            let r = Jet::constant(I::point(body.tail.unwrap().radius));
            let z =
                z.add(Jet::variable(1).scale(I::point(body.full_length_mm()) - I::point(station)));
            (
                r,
                z,
                source::Normal::Radial {
                    angle: angle(body, v).value,
                    beta: I::point(0.),
                    z: I::point(0.),
                },
            )
        }
        Surface::OuterCap => (
            Jet::constant(I::point(body.tail.unwrap().radius)),
            z,
            source::Normal::Axial(I::point(1.)),
        ),
    };
    let theta = angle(body, v);
    source::Field {
        point: [r.mul(theta.cos()), r.mul(theta.sin()), z],
        normal,
    }
}
fn source_point(body: &Cavity, v: I, r: f64, z: f64) -> [I; 3] {
    // A true complete turn is the SAME authored endpoint; evaluate its one
    // shared vertex at zero while retaining this exact periodic identity.
    let v = if v.lo == 1. && v.hi == 1. {
        I::point(0.)
    } else {
        v
    };
    let theta =
        v * I::new(
            std::f64::consts::TAU.next_down(),
            std::f64::consts::TAU.next_up(),
        ) * I::point(if body.thread.spec().left_hand {
            -1.
        } else {
            1.
        });
    [
        I::point(r) * theta.cos(),
        I::point(r) * theta.sin(),
        I::point(z),
    ]
}
fn error(mesh: &DisplayMesh, id: usize, p: [I; 3]) -> f64 {
    p.into_iter()
        .zip(&mesh.vertices[id * 3..id * 3 + 3])
        .fold(I::point(0.), |s, (p, l)| {
            s + I::point((p - I::point(*l)).magnitude())
        })
        .hi
}
fn outer(
    body: &Cavity,
    edge: Edge,
    z: f64,
    map: &mut [usize],
    mesh: &mut DisplayMesh,
    w: &mut Work<'_>,
) -> Result<([usize; 2], [f64; 2]), Error> {
    let mut ids = [0; 2];
    let mut errors = [0.; 2];
    for i in 0..2 {
        w.step()?;
        let original = edge.vertices[i];
        let source = source_point(body, edge.v[i], body.tail.unwrap().radius, z);
        if map[original] == usize::MAX {
            let id = mesh.vertices.len() / 3;
            let nominal = source.map(|v| v.lo.midpoint(v.hi));
            // Radius and z are bounded source dimensions <=1e7; v is
            // within[0,1], so angle<7 and every trig/product/midpoint is finite.
            mesh.vertices.extend(nominal);
            map[original] = id;
        }
        ids[i] = map[original];
        errors[i] = error(mesh, ids[i], source);
    }
    Ok((ids, errors))
}
fn quad(
    shell: &mut Shell,
    field: source::Field,
    ids: [usize; 4],
    errors: [f64; 4],
    face: u32,
    options: Tessellation,
    w: &mut Work<'_>,
) -> Result<bool, Error> {
    let points: Vec<[f64; 3]> = ids
        .iter()
        .map(|&i| shell.mesh.vertices[i * 3..i * 3 + 3].try_into().unwrap())
        .collect();
    for indices in [[0, 1, 2], [0, 2, 3]] {
        let Some(proof) = qualify(&field, &points, &errors, indices, options, w)? else {
            return Ok(false);
        };
        append_qualified(&mut shell.mesh, &ids, proof, face, options)?;
        shell.chord = shell.chord.max(proof.chord);
        shell.normal = shell.normal.max(proof.angle);
    }
    Ok(true)
}
struct Disk {
    edge: Edge,
    ids: [usize; 2],
    errors: [f64; 2],
    center: usize,
    end: u32,
    mode: Surface,
}
fn cap(
    shell: &mut Shell,
    body: &Cavity,
    disk: Disk,
    options: Tessellation,
    w: &mut Work<'_>,
) -> Result<bool, Error> {
    let Disk {
        edge,
        ids,
        errors,
        center,
        end,
        mode,
    } = disk;
    w.step()?;
    let station = if end == 0 { 0. } else { body.full_length_mm() };
    let field = field(body, edge, station, mode);
    // The exact disk is a radial fan of its actual boundary; at any radial
    // fraction the boundary chord discrepancy scales by <=1. Its plane normal
    // is exact throughout. No curved boundary is replaced by a source polygon.
    let chord = certificate::chord(&field, [errors[0], errors[1], 0.]);
    if !chord.is_finite() || chord > options.deflection {
        return Ok(false);
    }
    let ids = [center, ids[0], ids[1]];
    let points = ids.map(|i| shell.mesh.vertices[i * 3..i * 3 + 3].try_into().unwrap());
    let Some((reverse, normal)) =
        certificate::oriented_triangle(points, field.normal, options.angle)
    else {
        return Ok(false);
    };
    if shell.mesh.face_indices.len() == options.triangle_limit {
        return Err(Error::Limit);
    }
    shell.mesh.triangle(
        if reverse {
            [ids[0], ids[2], ids[1]]
        } else {
            ids
        },
        end,
    );
    shell.chord = shell.chord.max(chord);
    shell.normal = shell.normal.max(normal);
    Ok(true)
}
pub(super) fn close(
    mut shell: Shell,
    body: &Cavity,
    options: Tessellation,
    w: &mut Work<'_>,
) -> Result<Option<Candidate>, Error> {
    w.step()?;
    let count = shell.mesh.vertices.len() / 3;
    // Both maps are fixed by already bounded exact source vertices; no retained
    // rational heap and no geometric proximity search is introduced.
    w.charge(2 * count)?;
    let mut bottom = Vec::new();
    let mut top = Vec::new();
    bottom.try_reserve_exact(count).map_err(|_| Error::Limit)?;
    top.try_reserve_exact(count).map_err(|_| Error::Limit)?;
    let allocations = [
        (bottom.capacity(), std::mem::size_of::<usize>()),
        (top.capacity(), std::mem::size_of::<usize>()),
    ];
    storage::bytes(&allocations, 0, 4 * 1024 * 1024)?;
    bottom.resize(count, usize::MAX);
    top.resize(count, usize::MAX);
    let centers = [count, count + 1];
    shell
        .mesh
        .vertices
        .extend([0., 0., 0., 0., 0., body.full_length_mm()]);
    // One typed cap queue and one quality gate serve all physical disk kinds.
    // Rims are already bounded16k each; retain at most32k fixed-size receipts.
    w.charge(shell.mouth.len() + shell.top.len())?;
    let mut caps = Vec::new();
    caps.try_reserve_exact(shell.mouth.len() + shell.top.len())
        .map_err(|_| Error::Limit)?;
    let allocation = [(caps.capacity(), std::mem::size_of::<Disk>())];
    storage::bytes(&allocation, 0, 16 * 1024 * 1024)?;
    for edge in std::mem::take(&mut shell.mouth) {
        caps.push(Disk {
            edge,
            ids: edge.vertices,
            errors: edge.errors,
            center: centers[0],
            end: 0,
            mode: Surface::InnerCap,
        });
    }
    let base = 2 + 4 * body.thread.spec().starts + u32::from(body.entry.is_some());
    for edge in std::mem::take(&mut shell.top) {
        let (ids, errors, mode) = if let Some(tail) = body.tail {
            let (low, le) = outer(body, edge, tail.station, &mut bottom, &mut shell.mesh, w)?;
            let (high, he) = outer(
                body,
                edge,
                body.full_length_mm(),
                &mut top,
                &mut shell.mesh,
                w,
            )?;
            for (kind, ids, errors, face) in [
                (
                    Surface::Shoulder,
                    [edge.vertices[0], edge.vertices[1], low[1], low[0]],
                    [edge.errors[0], edge.errors[1], le[1], le[0]],
                    base,
                ),
                (
                    Surface::Terminal,
                    [low[0], low[1], high[1], high[0]],
                    [le[0], le[1], he[1], he[0]],
                    base + 1,
                ),
            ] {
                if !quad(
                    &mut shell,
                    field(body, edge, tail.station, kind),
                    ids,
                    errors,
                    face,
                    options,
                    w,
                )? {
                    return Ok(None);
                }
            }
            (high, he, Surface::OuterCap)
        } else {
            (edge.vertices, edge.errors, Surface::InnerCap)
        };
        caps.push(Disk {
            edge,
            ids,
            errors,
            center: centers[1],
            end: 1,
            mode,
        });
    }
    for disk in caps {
        if !cap(&mut shell, body, disk, options, w)? {
            return Ok(None);
        }
    }
    Ok(Some(Candidate {
        mesh: shell.mesh,
        source_chord_bound_mm: shell.chord,
        source_normal_angle_bound_rad: shell.normal,
        work_used: w.used,
    }))
}
#[cfg(test)]
mod tests {
    use super::*;
    fn work() -> Work<'static> {
        Work {
            remaining: 16_000_000,
            used: 0,
            interrupt: None,
        }
    }
    fn shell(body: &Cavity) -> Shell {
        chart::shell(body, super::super::tests::options(), 1, 64, &mut work())
            .unwrap()
            .unwrap()
    }
    #[test]
    fn actual_source_closure_receipts_independently_refuse_tighter_quality_and_collapsed_output() {
        let body = super::super::tests::source(false, 1, Some(0.4), true);
        let options = super::super::tests::options();
        let mut candidate = shell(&body);
        let edge = candidate.top[0];
        let mut map = vec![usize::MAX; candidate.mesh.vertices.len() / 3];
        let tail = body.tail.unwrap();
        let (low, le) = outer(
            &body,
            edge,
            tail.station,
            &mut map,
            &mut candidate.mesh,
            &mut work(),
        )
        .unwrap();
        let points = [edge.vertices[0], edge.vertices[1], low[1], low[0]];
        assert!(
            !quad(
                &mut candidate,
                field(&body, edge, tail.station, Surface::Shoulder),
                points,
                [edge.errors[0], edge.errors[1], le[1], le[0]],
                10,
                Tessellation {
                    deflection: 1e-12,
                    ..options
                },
                &mut work()
            )
            .unwrap()
        );
        // The independent disk gate receives an actual source boundary. Asking
        // for much tighter quality must not reuse its previous looser receipt.
        let edge = candidate.mouth[0];
        let center = candidate.mesh.vertices.len() / 3;
        candidate.mesh.vertices.extend([0., 0., 0.]);
        assert!(
            !cap(
                &mut candidate,
                &body,
                Disk {
                    edge,
                    ids: edge.vertices,
                    errors: edge.errors,
                    center,
                    end: 0,
                    mode: Surface::InnerCap
                },
                Tessellation {
                    deflection: 1e-12,
                    ..options
                },
                &mut work()
            )
            .unwrap()
        );
        // A degenerate emitted face is a bad candidate, never bad source data.
        assert!(
            !cap(
                &mut candidate,
                &body,
                Disk {
                    edge,
                    ids: [edge.vertices[0]; 2],
                    errors: edge.errors,
                    center,
                    end: 0,
                    mode: Surface::InnerCap
                },
                options,
                &mut work()
            )
            .unwrap()
        );
        assert!(
            close(
                shell(&body),
                &body,
                Tessellation {
                    deflection: 1e-12,
                    ..options
                },
                &mut work()
            )
            .unwrap()
            .is_none()
        );
        let entry_only = super::super::tests::source(false, 1, Some(0.4), false);
        assert!(
            close(
                shell(&entry_only),
                &entry_only,
                Tessellation {
                    deflection: 1e-12,
                    ..options
                },
                &mut work()
            )
            .unwrap()
            .is_none()
        );
        let large = Cavity::new(
            body.thread,
            Intent {
                entry_depth_mm: Some(0.4),
                terminal: Some(Terminal {
                    width_mm: 0.3,
                    radial_depth_mm: 100.,
                }),
            },
        )
        .unwrap();
        assert!(
            close(shell(&large), &large, options, &mut work())
                .unwrap()
                .is_none()
        );
    }
    #[test]
    fn every_closure_checkpoint_cancels_and_preserves_the_actual_source_receipt() {
        let body = super::super::tests::source(true, 2, Some(0.4), true);
        let unchanged = body;
        let original = shell(&body);
        let expected = original.mesh.clone();
        let options = super::super::tests::options();
        let calls = std::cell::Cell::new(0);
        let count = || {
            calls.set(calls.get() + 1);
            false
        };
        close(
            original.clone(),
            &body,
            options,
            &mut Work {
                remaining: 16_000_000,
                used: 0,
                interrupt: Some(&count),
            },
        )
        .unwrap()
        .unwrap();
        for stop in 0..calls.get() {
            let calls = std::cell::Cell::new(0);
            let cancel = || {
                calls.set(calls.get() + 1);
                calls.get() > stop
            };
            assert!(matches!(
                close(
                    original.clone(),
                    &body,
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
        assert_eq!(original.mesh, expected);
        assert_eq!(body, unchanged);
        assert!(
            close(original, &body, options, &mut work())
                .unwrap()
                .is_some()
        );
    }
}
