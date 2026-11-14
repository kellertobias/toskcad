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
                });
            let error = upper(error)?;
            let id = mesh.vertices.len() / 3;
            room(mesh.vertices.len(), 3, mesh.vertices.capacity())?;
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
            let ids = oriented(ids, tri, proof.reverse);
            let face = source.face_index(FaceKey::PlainWall).unwrap();
            add_triangle(mesh, ids, face, options)?;
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
    if !matches!(
        key,
        FaceKey::MetricTip | FaceKey::MetricShoulder | FaceKey::WoodShoulder | FaceKey::WoodTip
    ) {
        return Err(E::Unresolved);
    }
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
pub(super) struct Target {
    pub key: FaceKey,
    pub z: f64,
    pub plane_error: f64,
}
pub(super) fn fill(
    mesh: &mut DisplayMesh,
    outer: &Ring,
    inner: Option<&Ring>,
    target: Target,
    source: &Hanger,
    options: Tessellation,
    w: &mut Work<'_>,
) -> R<f64> {
    let Target {
        key,
        z,
        plane_error,
    } = target;
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
            let tri = [
                center,
                outer.vertices[i].id,
                outer.vertices[(i + 1) % outer.vertices.len()].id,
            ];
            emit(mesh, tri, key, z, source, options, w)?;
        }
    }
    Ok(error)
}
fn oriented(ids: [usize; 4], tri: [usize; 3], reverse: bool) -> [usize; 3] {
    if reverse {
        [ids[tri[0]], ids[tri[2]], ids[tri[1]]]
    } else {
        tri.map(|i| ids[i])
    }
}

#[cfg(test)]
pub(super) mod tests {
    use super::*;
    pub(in super::super) fn circles() -> (Hanger, DisplayMesh, Ring, Ring) {
        let source = Hanger::new(
            crate::authored_hanger_bolt::tests::fixture(false, false),
            None,
        )
        .unwrap();
        let mut w = Work {
            remaining: 16_000_000,
            used: 0,
            interrupt: None,
        };
        let mut m = Memory { used: SCRATCH };
        let mut mesh = DisplayMesh {
            vertices: Vec::with_capacity(3 * 34),
            indices: Vec::with_capacity(3 * 64),
            face_indices: Vec::with_capacity(64),
        };
        let options = Tessellation {
            deflection: 0.1,
            angle: 1.,
            triangle_limit: 100,
        };
        let (a, b, _) = plain(&mut mesh, &source, 16, options, &mut m, &mut w).unwrap();
        (source, mesh, a, b)
    }
    #[test]
    fn actual_plain_source_refuses_underresolved_facets_and_exhausted_output() {
        let (source, _, _, _) = circles();
        let options = Tessellation {
            deflection: 0.1,
            angle: 0.01,
            triangle_limit: 100,
        };
        let mut w = Work {
            remaining: 16_000_000,
            used: 0,
            interrupt: None,
        };
        let mut m = Memory { used: SCRATCH };
        let mut mesh = DisplayMesh {
            vertices: Vec::with_capacity(96),
            indices: Vec::with_capacity(96),
            face_indices: Vec::with_capacity(32),
        };
        assert!(matches!(
            plain(&mut mesh, &source, 16, options, &mut m, &mut w),
            Err(E::Unresolved)
        ));
        let mut mesh = DisplayMesh {
            vertices: Vec::with_capacity(96),
            indices: Vec::new(),
            face_indices: Vec::new(),
        };
        assert!(matches!(
            plain(
                &mut mesh,
                &source,
                16,
                Tessellation {
                    angle: 1.,
                    ..options
                },
                &mut m,
                &mut w
            ),
            Err(E::Limit)
        ));
    }
    #[test]
    fn exact_angular_zipper_handles_independently_shifted_circular_rings() {
        // Actual analytic circles exercise the primitive annulus fill without
        // claiming they are the Hanger's threaded shoulder boundary.
        let (source, _, _, _) = circles();
        let mut w = Work {
            remaining: 16_000_000,
            used: 0,
            interrupt: None,
        };
        let mut m = Memory { used: SCRATCH };
        let mut mesh = DisplayMesh {
            vertices: Vec::with_capacity(3 * 24),
            indices: Vec::with_capacity(3 * 24),
            face_indices: Vec::with_capacity(24),
        };
        let mut circle = |radius: f64, count: usize, shift: f64| {
            let mut ring = Ring {
                vertices: Vec::with_capacity(count),
                curve_error: radius * (std::f64::consts::TAU / count as f64).powi(2) / 8. + 1e-12,
                source_min: radius,
            };
            for i in 0..count {
                let q = Q::float(shift + i as f64 / count as f64, &mut w).unwrap();
                let theta = q.enclosure(&mut w).unwrap().1 * tau();
                let id = mesh.vertices.len() / 3;
                mesh.vertices.extend([
                    (I::point(radius) * theta.cos()).midpoint(),
                    (I::point(radius) * theta.sin()).midpoint(),
                    0.,
                ]);
                ring.vertices.push(Vertex {
                    id,
                    angle: Angle::new(&q, &mut m, &mut w).unwrap(),
                });
            }
            ring
        };
        let outer = circle(3., 8, 0.);
        let inner = circle(1., 16, 1. / 32.);
        let error = fill(
            &mut mesh,
            &outer,
            Some(&inner),
            Target {
                key: FaceKey::WoodShoulder,
                z: 0.,
                plane_error: 0.,
            },
            &source,
            Tessellation {
                deflection: 1.,
                angle: 1.,
                triangle_limit: 100,
            },
            &mut w,
        )
        .unwrap();
        assert!(error < 1.);
        assert_eq!(mesh.face_indices.len(), 24);
    }
    #[test]
    fn actual_plane_orientation_and_nested_ring_proof_refuse_wrong_geometry() {
        let (source, mesh, outer, inner) = circles();
        // Reverse an actual plain-wall facet proposal: the unchanged whole
        // source proof earns the sign which the emitted indices must retain.
        let angle = jet::Jet::variable(0).scale(tau() / I::point(16.));
        let z = jet::Jet::constant(I::point(-2.)).add(jet::Jet::variable(1).scale(I::point(2.)));
        let field = source::Field {
            point: [
                angle.cos().scale(I::point(1.5)),
                angle.sin().scale(I::point(1.5)),
                z,
            ],
            normal: source::Normal::Radial {
                angle: angle.value,
                beta: I::point(0.),
                z: I::point(0.),
            },
        };
        let ids = [0, 1, 17, 16];
        let points = ids.map(|id| point(&mesh, id));
        let mut proof_work = Work {
            remaining: 16_000_000,
            used: 0,
            interrupt: None,
        };
        let proof = qualify_attempt(
            &field,
            &points,
            &[1e-12; 4],
            [0, 2, 1],
            Tessellation {
                deflection: 0.1,
                angle: 1.,
                triangle_limit: 100,
            },
            &mut proof_work,
        )
        .unwrap();
        assert!(matches!(proof, Attempt::Complete(p) if p.reverse));
        let outward = oriented(ids, [0, 2, 1], true);
        assert_eq!(outward, [0, 1, 17]);
        let mut w = Work {
            remaining: 16_000_000,
            used: 0,
            interrupt: None,
        };
        let options = Tessellation {
            deflection: 0.1,
            angle: 1.,
            triangle_limit: 100,
        };
        assert!(
            circle_count(
                1.5,
                Tessellation {
                    deflection: f64::from_bits(1),
                    ..options
                }
            )
            .is_err()
        );
        let mut same = mesh.clone();
        assert!(
            fill(
                &mut same,
                &outer,
                Some(&inner),
                Target {
                    key: FaceKey::MetricShoulder,
                    z: -2.,
                    plane_error: 0.
                },
                &source,
                options,
                &mut w
            )
            .is_err()
        );
        let mut excess = mesh.clone();
        assert!(
            fill(
                &mut excess,
                &outer,
                None,
                Target {
                    key: FaceKey::MetricTip,
                    z: -2.,
                    plane_error: 1.
                },
                &source,
                options,
                &mut w
            )
            .is_err()
        );
        let mut tiny = mesh.clone();
        tiny.vertices.shrink_to_fit();
        assert!(
            fill(
                &mut tiny,
                &outer,
                None,
                Target {
                    key: FaceKey::MetricTip,
                    z: -2.,
                    plane_error: 0.
                },
                &source,
                options,
                &mut w
            )
            .is_err()
        );
        let mut triangle = DisplayMesh {
            vertices: vec![0., 0., 12., 0.1, 0., 12., 0., 0.1, 12.],
            indices: Vec::with_capacity(3),
            face_indices: Vec::with_capacity(1),
        };
        assert!(
            emit(
                &mut triangle,
                [0, 1, 2],
                FaceKey::PlainWall,
                12.,
                &source,
                options,
                &mut w
            )
            .is_err()
        );
        assert!(
            emit(
                &mut triangle,
                [0, 2, 1],
                FaceKey::WoodTip,
                12.,
                &source,
                options,
                &mut w
            )
            .is_err()
        );
        assert!(
            emit(
                &mut triangle,
                [0, 1, 2],
                FaceKey::WoodTip,
                11.,
                &source,
                options,
                &mut w
            )
            .is_err()
        );
        emit(
            &mut triangle,
            [0, 1, 2],
            FaceKey::WoodTip,
            12.,
            &source,
            options,
            &mut w,
        )
        .unwrap();
    }
}
