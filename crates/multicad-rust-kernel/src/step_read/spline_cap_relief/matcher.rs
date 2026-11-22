//! Bounded complete world-boundary equivalence of two retained descriptors.
use super::*;
use crate::spline_cap_relief::{Boundary, Curve as C, Surface as S};
fn curve_equal(
    a: &C,
    af: Frame,
    b: &C,
    bf: Frame,
    forward: bool,
    tol: f64,
    work: &mut usize,
) -> Result<bool> {
    charge(work, 1)?;
    let endpoints = |c: &C, f: Frame, reverse: bool| {
        [0., 1.].map(|t| f.world(c.point(if reverse { 1. - t } else { t })))
    };
    let ends_a = endpoints(a, af, false);
    let ends_b = endpoints(b, bf, !forward);
    if !ends_a.into_iter().zip(ends_b).all(|(a, b)| near(a, b, tol)) {
        return Ok(false);
    }
    match (a, b) {
        (
            C::Profile {
                edge: Edge::Cubic(a),
                z_mm: za,
            },
            C::Profile {
                edge: Edge::Cubic(b),
                z_mm: zb,
            },
        ) => {
            let a = a.bezier_spans().collect::<Vec<_>>();
            let mut b = b.bezier_spans().collect::<Vec<_>>();
            if a.len() != b.len() {
                return Ok(false);
            }
            charge(work, a.len() * 4)?;
            if !forward {
                b.reverse();
                for span in &mut b {
                    span.reverse();
                }
            }
            Ok(a.iter().flatten().zip(b.iter().flatten()).all(|(&a, &b)| {
                near(
                    af.world([a[0], a[1], *za]),
                    bf.world([b[0], b[1], *zb]),
                    tol,
                )
            }))
        }
        (
            C::Profile {
                edge: Edge::Analytic(a),
                z_mm: za,
            },
            C::Profile {
                edge: Edge::Analytic(b),
                z_mm: zb,
            },
        ) => {
            use crate::normal_sweep::Edge as A;
            match (a, b) {
                (A::Line { .. }, A::Line { .. }) => Ok(true),
                (
                    A::Circle {
                        center: ca,
                        clockwise: wa,
                        ..
                    },
                    A::Circle {
                        center: cb,
                        clockwise: wb,
                        ..
                    },
                )
                | (
                    A::Arc {
                        center: ca,
                        clockwise: wa,
                        ..
                    },
                    A::Arc {
                        center: cb,
                        clockwise: wb,
                        ..
                    },
                ) => {
                    let r = a.radius_mm().expect("circular profile");
                    Ok(near(
                        af.world([ca[0], ca[1], *za]),
                        bf.world([cb[0], cb[1], *zb]),
                        tol,
                    ) && (r - b.radius_mm().expect("circular profile")).abs() <= tol
                        && (a.signed_sweep_rad().expect("circle").abs()
                            - b.signed_sweep_rad().expect("circle").abs())
                        .abs()
                            * r
                            <= tol
                        && aligned(
                            af.axes[2].map(|x| if *wa { -x } else { x }),
                            bf.axes[2].map(|x| if *wb == forward { -x } else { x }),
                            r,
                            tol,
                        ))
                }
                _ => Ok(false),
            }
        }
        (C::Line { .. }, C::Line { .. }) => Ok(true),
        _ => Ok(false),
    }
}
fn support_equal(a: &S, af: Frame, b: &S, bf: Frame, extent: f64, tol: f64) -> bool {
    match (a, b) {
        (
            S::Plane {
                z_mm: a,
                upward: ua,
            },
            S::Plane {
                z_mm: b,
                upward: ub,
            },
        ) => {
            let na = af.axes[2].map(|x| if *ua { x } else { -x });
            let nb = bf.axes[2].map(|x| if *ub { x } else { -x });
            aligned(na, nb, extent, tol)
                && dot(sub(af.world([0., 0., *a]), bf.world([0., 0., *b])), na).abs() <= tol
        }
        (S::Extrusion { .. }, S::Extrusion { .. }) => aligned(af.axes[2], bf.axes[2], extent, tol),
        _ => false,
    }
}
fn verify(a: &Boundary, af: Frame, b: &Boundary, bf: Frame, tol: f64) -> Result<()> {
    require(
        a.vertices.len() == b.vertices.len()
            && a.edges.len() == b.edges.len()
            && a.faces.len() == b.faces.len(),
    )?;
    let mut work = MAX_WORK;
    let mut actual = b
        .vertices
        .iter()
        .enumerate()
        .map(|(i, &p)| (bf.world(p), i))
        .collect::<Vec<_>>();
    actual.sort_unstable_by(|a, b| a.0[0].total_cmp(&b.0[0]));
    let mut used = BTreeSet::new();
    let mut vertices = Vec::new();
    for &p in &a.vertices {
        let p = af.world(p);
        let start = actual.partition_point(|(q, _)| q[0] < p[0] - tol);
        let end = actual.partition_point(|(q, _)| q[0] <= p[0] + tol);
        charge(&mut work, 1 + end - start)?;
        let matches = actual[start..end]
            .iter()
            .filter(|(q, _)| near(p, *q, tol))
            .collect::<Vec<_>>();
        let &[&(_, id)] = matches.as_slice() else {
            return Err(ReadError::Invalid);
        };
        require(used.insert(id))?;
        vertices.push(id);
    }
    let mut by_ends: BTreeMap<[usize; 2], Vec<usize>> = BTreeMap::new();
    for (i, e) in b.edges.iter().enumerate() {
        charge(&mut work, 1)?;
        let mut k = e.vertices;
        k.sort_unstable();
        by_ends.entry(k).or_default().push(i);
    }
    used.clear();
    let mut edges = Vec::new();
    for e in &a.edges {
        let ends = e.vertices.map(|i| vertices[i]);
        let mut key = ends;
        key.sort_unstable();
        let candidates = by_ends.get(&key).ok_or(ReadError::Invalid)?;
        let mut matches = vec![];
        for &i in candidates {
            for forward in [true, false] {
                let actual = &b.edges[i];
                charge(&mut work, 1)?;
                if (if forward {
                    actual.vertices
                } else {
                    [actual.vertices[1], actual.vertices[0]]
                }) == ends
                    && curve_equal(&e.curve, af, &actual.curve, bf, forward, tol, &mut work)?
                {
                    matches.push((i, forward));
                }
            }
        }
        let [(i, forward)] = matches[..] else {
            return Err(ReadError::Invalid);
        };
        require(used.insert(i))?;
        edges.push((i, forward));
    }
    let cycle = |wire: Vec<(usize, bool)>| {
        let mut wire = wire;
        if let Some((i, _)) = wire.iter().enumerate().min_by_key(|(_, e)| **e) {
            wire.rotate_left(i);
        }
        wire
    };
    let mut faces = BTreeMap::new();
    for (i, f) in b.faces.iter().enumerate() {
        let outer = cycle(f.loops[0].iter().map(|e| (e.edge, e.forward)).collect());
        require(faces.insert(outer, i).is_none())?;
    }
    let extent = a
        .vertices
        .iter()
        .flatten()
        .map(|x| x.abs())
        .fold(1., f64::max);
    for f in &a.faces {
        let mut rings = f
            .loops
            .iter()
            .map(|r| {
                cycle(
                    r.iter()
                        .map(|e| {
                            let (i, forward) = edges[e.edge];
                            (i, e.forward == forward)
                        })
                        .collect(),
                )
            })
            .collect::<Vec<_>>();
        charge(&mut work, rings.iter().map(Vec::len).sum())?;
        let i = faces.remove(&rings[0]).ok_or(ReadError::Invalid)?;
        let actual = &b.faces[i];
        let mut holes = actual
            .loops
            .iter()
            .skip(1)
            .map(|r| cycle(r.iter().map(|e| (e.edge, e.forward)).collect()))
            .collect::<Vec<_>>();
        rings.remove(0);
        rings.sort();
        holes.sort();
        require(rings == holes && support_equal(&f.surface, af, &actual.surface, bf, extent, tol))?;
    }
    require(faces.is_empty())
}
pub fn spline_cap_relief_import_matches(
    imported: &ImportSolid,
    expected: &Relief,
    frame: [Point; 3],
    tolerance: f64,
) -> bool {
    if !tolerance.is_finite() || tolerance <= 0. {
        return false;
    }
    let ImportProfile::SplineCapRelief { relief } = &imported.profile else {
        return false;
    };
    let Ok(af) = Frame::new(frame[0], frame[1], frame[2]) else {
        return false;
    };
    let Ok(bf) = Frame::new(imported.origin_mm, imported.u, imported.v) else {
        return false;
    };
    if metric_frame(af, expected.bounds(), tolerance / 32.).is_err()
        || metric_frame(bf, relief.bounds(), tolerance / 32.).is_err()
    {
        return false;
    }
    verify(
        expected.boundary(),
        af,
        relief.boundary(),
        bf,
        tolerance / 16.,
    )
    .is_ok()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        normal_sweep::Edge as A,
        spline_cap_relief::{Kind, tests::fixture},
        spline_sweep::Cubic,
    };
    #[test]
    fn actual_curve_support_equivalence_distinguishes_chords_arcs_and_span_partitions() {
        let frame = Frame::new([0.; 3], [1., 0., 0.], [0., 1., 0.]).unwrap();
        let poles = vec![[0., 0.], [0.5, 1.], [1.5, 1.], [2., 0.]];
        let a = C::Profile {
            edge: Edge::Cubic(Cubic::new(poles.clone(), vec![0., 1.], false).unwrap()),
            z_mm: 0.,
        };
        let mut poles = poles;
        poles.reverse();
        let reversed = C::Profile {
            edge: Edge::Cubic(Cubic::new(poles, vec![0., 1.], false).unwrap()),
            z_mm: 0.,
        };
        let mut work = 100;
        assert!(curve_equal(&a, frame, &reversed, frame, false, 1e-7, &mut work).unwrap());
        let split = C::Profile {
            edge: Edge::Cubic(
                Cubic::new(
                    vec![
                        [0., 0.],
                        [0.2, 0.5],
                        [0.7, 0.5],
                        [1.3, 0.5],
                        [1.7, 0.5],
                        [2., 0.],
                    ],
                    vec![0., 0.5, 1.],
                    false,
                )
                .unwrap(),
            ),
            z_mm: 0.,
        };
        assert!(!curve_equal(&a, frame, &split, frame, true, 1e-7, &mut work).unwrap());
        let line = C::Profile {
            edge: Edge::Analytic(A::Line {
                a: [1., 0.],
                b: [0., 1.],
            }),
            z_mm: 0.,
        };
        let arc = C::Profile {
            edge: Edge::Analytic(A::Arc {
                a: [1., 0.],
                b: [0., 1.],
                center: [0.; 2],
                clockwise: false,
            }),
            z_mm: 0.,
        };
        assert!(!curve_equal(&line, frame, &arc, frame, true, 1e-7, &mut work).unwrap());
        let rail = C::Line {
            a_mm: [1., 0., 0.],
            b_mm: [0., 1., 0.],
        };
        assert!(!curve_equal(&rail, frame, &arc, frame, true, 1e-7, &mut work).unwrap());
        let other = C::Line {
            a_mm: [9., 0., 0.],
            b_mm: [0., 1., 0.],
        };
        assert!(!curve_equal(&rail, frame, &other, frame, true, 1e-7, &mut work).unwrap());
        let plane = S::Plane {
            z_mm: 0.,
            upward: true,
        };
        let wall = S::Extrusion {
            edge: Edge::Analytic(A::Line {
                a: [1., 0.],
                b: [0., 1.],
            }),
            span_mm: [0., 1.],
            inward: false,
        };
        assert!(!support_equal(&plane, frame, &wall, frame, 2., 1e-7));
        let one = fixture(Kind::Raised);
        let p = [[12., 2.], [14., 2.], [14., 4.], [12., 4.]];
        let extra = crate::spline_sweep::Sweep::new(
            (0..4)
                .map(|i| {
                    Edge::Analytic(A::Line {
                        a: p[i],
                        b: p[(i + 1) % 4],
                    })
                })
                .collect(),
            None,
            1.,
            1e-7,
        )
        .unwrap();
        let two = Relief::new(
            one.base().clone(),
            vec![one.motifs()[0].clone(), extra],
            Kind::Raised,
            1e-7,
            1e-4,
        )
        .unwrap();
        assert_eq!(
            verify(one.boundary(), frame, two.boundary(), frame, 1e-7),
            Err(ReadError::Invalid)
        );
    }
}
