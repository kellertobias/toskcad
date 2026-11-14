//! Whole actual cap sets, using source-curve correspondence and independently
//! certified literal radial graphs. No literal angle is declared a source knot.
use super::*;
use ring::{Angle, Vertex};
pub(super) fn circle_count(radius: f64, options: Tessellation) -> R<usize> {
    let angle = (options.angle * 0.5).min((2. * options.deflection / radius).sqrt());
    let rows = (std::f64::consts::TAU / angle).ceil();
    if !rows.is_finite() || rows > 65536. {
        return Err(E::Limit);
    }
    Ok((rows as usize).next_power_of_two().max(8))
}
pub(super) fn plain(
    mesh: &mut DisplayMesh,
    source: &Hanger,
    count: usize,
    options: Tessellation,
    m: &mut Memory,
    w: &mut Work<'_>,
) -> R<(Ring, Ring, f64)> {
    let radius = source.spec().plain_radius_mm;
    let height = source.spec().plain_length_mm;
    let mut errors = m.allocate(2 * count, w)?;
    let mut rings = [
        Ring {
            vertices: m.allocate(count, w)?,
            curve_error: 0.,
            source_min: radius,
        },
        Ring {
            vertices: m.allocate(count, w)?,
            curve_error: 0.,
            source_min: radius,
        },
    ];
    let denominator = Q::integer(count as i64, w)?;
    for (end, ring) in rings.iter_mut().enumerate() {
        let z = if end == 0 { -height } else { 0. };
        for i in 0..count {
            w.charge(32)?;
            let q = Q::integer(i as i64, w)?.div(&denominator, w)?;
            let angle = q.enclosure(w)?.1 * tau();
            let actual = [
                I::point(radius) * angle.cos(),
                I::point(radius) * angle.sin(),
                I::point(z),
            ];
            let literal = actual.map(I::midpoint);
            let error = actual
                .into_iter()
                .zip(literal)
                .fold(I::point(0.), |s, (a, b)| {
                    s + I::point((a - I::point(b)).magnitude())
                })
                .hi;
            if !error.is_finite() {
                return Err(E::Unresolved);
            }
            let id = mesh.vertices.len() / 3;
            if mesh.vertices.len() + 3 > mesh.vertices.capacity() {
                return Err(E::Limit);
            }
            mesh.vertices.extend(literal);
            ring.vertices.push(Vertex {
                id,
                angle: Angle::new(&q, m, w)?,
            });
            errors.push(error);
        }
        let width = tau() / I::point(count as f64);
        ring.curve_error = (I::point(radius) * width.square() / I::point(8.)
            + I::point(
                errors[end * count..(end + 1) * count]
                    .iter()
                    .copied()
                    .fold(0., f64::max),
            ))
        .hi;
    }
    let mut chord = 0_f64;
    for i in 0..count {
        w.charge(32)?;
        let angle = jet::Jet::constant(tau() * I::point(i as f64 / count as f64))
            .add(jet::Jet::variable(0).scale(tau() / I::point(count as f64)));
        let z = jet::Jet::constant(I::point(-height))
            .add(jet::Jet::variable(1).scale(I::point(height)));
        let field = source::Field {
            point: [
                angle.cos().scale(I::point(radius)),
                angle.sin().scale(I::point(radius)),
                z,
            ],
            normal: source::Normal::Radial {
                angle: angle.value,
                beta: I::point(0.),
                z: I::point(0.),
            },
        };
        let ids = [
            rings[0].vertices[i].id,
            rings[0].vertices[(i + 1) % count].id,
            rings[1].vertices[(i + 1) % count].id,
            rings[1].vertices[i].id,
        ];
        let points = ids.map(|id| point(mesh, id));
        let err = [
            errors[i],
            errors[(i + 1) % count],
            errors[count + (i + 1) % count],
            errors[count + i],
        ];
        for tri in [[0, 1, 2], [0, 2, 3]] {
            let proof = match qualify_attempt(&field, &points, &err, tri, options, w)? {
                Attempt::Complete(p) => p,
                _ => return Err(E::Unresolved),
            };
            chord = chord.max(proof.chord);
            let ids = if proof.reverse {
                [ids[tri[0]], ids[tri[2]], ids[tri[1]]]
            } else {
                tri.map(|i| ids[i])
            };
            add_triangle(
                mesh,
                ids,
                source.face_index(FaceKey::PlainWall).unwrap(),
                options,
            )?;
        }
    }
    let [bottom, top] = rings;
    Ok((bottom, top, chord))
}
fn emit(
    mesh: &mut DisplayMesh,
    ids: [usize; 3],
    key: FaceKey,
    z: f64,
    source: &Hanger,
    options: Tessellation,
    w: &mut Work<'_>,
) -> R<()> {
    w.charge(12)?;
    let p = ids.map(|id| point(mesh, id));
    if p.iter().any(|v| v[2] != z) || ring::area(p, w)? != Ordering::Greater {
        return Err(E::Unresolved);
    }
    // Exact literal equality to the same z plane plus positive exact projected
    // area proves the actual normal, without normalizing an underflowed cross.
    let ids = if key == FaceKey::WoodTip {
        ids
    } else {
        [ids[0], ids[2], ids[1]]
    };
    add_triangle(
        mesh,
        ids,
        source.face_index(key).ok_or(E::Unresolved)?,
        options,
    )
}
fn next_angle(ring: &Ring, index: usize, w: &mut Work<'_>) -> R<Q> {
    let n = ring.vertices.len();
    let q = ring.vertices[index % n].angle.exact(w)?;
    if index >= n {
        Ok(q.add(&Q::integer(1, w)?, w)?)
    } else {
        Ok(q)
    }
}
pub(super) fn fill(
    mesh: &mut DisplayMesh,
    outer: &Ring,
    inner: Option<&Ring>,
    key: FaceKey,
    z: f64,
    plane_error: f64,
    source: &Hanger,
    options: Tessellation,
    w: &mut Work<'_>,
) -> R<f64> {
    let outside = ring::graph(outer, mesh, w)?;
    let error = if let Some(inner) = inner {
        let inside = ring::graph(inner, mesh, w)?;
        if inside.maximum >= outside.minimum {
            return Err(E::Unresolved);
        }
        outside.radial_error.max(inside.radial_error)
    } else {
        outside.radial_error
    };
    let error = (I::point(error) + I::point(plane_error)).hi;
    if !error.is_finite() || error > options.deflection {
        return Err(E::Unresolved);
    }
    if let Some(inner) = inner {
        let (mut i, mut j) = (0, 0);
        let (on, inn) = (outer.vertices.len(), inner.vertices.len());
        while i < on || j < inn {
            w.step()?;
            let advance_outer = if i == on {
                false
            } else if j == inn {
                true
            } else {
                next_angle(outer, i + 1, w)?.cmp(&next_angle(inner, j + 1, w)?, w)?
                    == Ordering::Less
            };
            let ids = if advance_outer {
                let ids = [
                    outer.vertices[i % on].id,
                    outer.vertices[(i + 1) % on].id,
                    inner.vertices[j % inn].id,
                ];
                i += 1;
                ids
            } else {
                let ids = [
                    outer.vertices[i % on].id,
                    inner.vertices[(j + 1) % inn].id,
                    inner.vertices[j % inn].id,
                ];
                j += 1;
                ids
            };
            emit(mesh, ids, key, z, source, options, w)?;
        }
    } else {
        if mesh.vertices.len() + 3 > mesh.vertices.capacity() {
            return Err(E::Limit);
        }
        let center = mesh.vertices.len() / 3;
        mesh.vertices.extend([0., 0., z]);
        for i in 0..outer.vertices.len() {
            emit(
                mesh,
                [
                    center,
                    outer.vertices[i].id,
                    outer.vertices[(i + 1) % outer.vertices.len()].id,
                ],
                key,
                z,
                source,
                options,
                w,
            )?;
        }
    }
    Ok(error)
}
