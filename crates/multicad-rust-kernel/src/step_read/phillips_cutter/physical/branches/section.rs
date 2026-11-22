//! Uniform implicit-function contraction on a moving literal cubic section.
use super::*;
use functional::Functional;
use supports::chart::Chart;
#[derive(Clone, Copy, Debug)]
pub(in super::super) struct Section {
    pub parameter: [f64; 2],
    pub source_box: [[f64; 2]; 3],
    pub correction_mm: [f64; 3],
    /// Per unit parameter of this ORIGINAL Bezier span, not the subcell.
    pub derivative_correction_mm: [f64; 3],
    pub contraction: f64,
}
pub(super) fn split(p: [[I; 3]; 4]) -> [[[I; 3]; 4]; 2] {
    let mean = |a: [I; 3], b: [I; 3]| std::array::from_fn(|k| (a[k] + b[k]) / I::point(2.));
    let a = mean(p[0], p[1]);
    let b = mean(p[1], p[2]);
    let c = mean(p[2], p[3]);
    let d = mean(a, b);
    let e = mean(b, c);
    let f = mean(d, e);
    [[p[0], a, d, f], [f, e, c, p[3]]]
}
fn hull(p: [[I; 3]; 4]) -> [I; 3] {
    std::array::from_fn(|k| {
        I::new(
            p.iter().map(|p| p[k].lo).fold(f64::INFINITY, f64::min),
            p.iter().map(|p| p[k].hi).fold(f64::NEG_INFINITY, f64::max),
        )
    })
}
fn bound(p: &[I]) -> I {
    I::new(
        p.iter().map(|p| p.lo).fold(f64::INFINITY, f64::min),
        p.iter().map(|p| p.hi).fold(f64::NEG_INFINITY, f64::max),
    )
}
fn dot(a: [I; 3], b: [I; 3]) -> I {
    (0..3).fold(I::point(0.), |s, k| s + a[k] * b[k])
}
fn apply(m: [[f64; 2]; 2], v: [I; 2]) -> [I; 2] {
    m.map(|r| (0..2).fold(I::point(0.), |s, k| s + I::point(r[k]) * v[k]))
}
fn norm(v: [I; 2]) -> f64 {
    v.into_iter().map(I::magnitude).fold(0., f64::max)
}
fn physical(surface: &Surface, chart: &Chart, world: [I; 3]) -> bool {
    if let Surface::Cone {
        radius_mm,
        semi_angle_rad,
        ..
    } = surface
    {
        let a = support::adjacent(*semi_angle_rad);
        let c = a.cos();
        // The paired cached Chart::new already proved cos(a)>0 for this
        // SAME immutable parsed surface. Only nappe isolation remains here.
        let r = support::adjacent(*radius_mm) + (a.sin() / c) * chart.project(world)[2];
        r.valid() && r.lo > 0.
    } else {
        true
    }
}
pub(super) fn certify(
    surfaces: [&Surface; 2],
    charts: [&Chart; 2],
    f: Functional,
    p: [[I; 3]; 4],
    parameter: [f64; 2],
    position_limit: f64,
) -> Option<Section> {
    let base = hull(p);
    let basis = f.basis();
    let center = split(p)[0][3];
    let jacobian = |point: [I; 3]| {
        charts.map(|c| {
            let g = c.world_gradient(c.project(point));
            basis.map(|b| dot(g, b))
        })
    };
    let middle = jacobian(center);
    let Ok(m) = matrix::preconditioner(middle) else {
        return None;
    };
    // Preserve polynomial correlation F(C(t)), including all interval
    // intermediates. Interval operations cannot panic on a numerical refusal.
    let residual_raw = charts.map(|c| bound(&c.implicit_curve(p.map(|p| c.project(p)))));
    let residual = apply(m, residual_raw);
    let b = norm(residual);
    // Free unknowns are actual world-coordinate corrections. Any SAME-SECTION
    // root within Euclidean limit is in this search square; that limit is not
    // assigned as a measured source error.
    let trial = position_limit;
    let offset: [I; 3] = std::array::from_fn(|k| {
        (0..2).fold(I::point(0.), |s, j| s + basis[j][k] * I::new(-trial, trial))
    });
    let world = std::array::from_fn(|k| base[k] + offset[k]);
    let j = jacobian(world);
    let q = matrix::contraction(m, j);
    let row = q.map(|r| {
        r.into_iter()
            .fold(I::point(0.), |s, v| s + I::point(v.magnitude()))
    });
    let contraction = row.into_iter().map(|r| r.hi).fold(0., f64::max);
    let inclusion = I::point(b) + I::point(contraction) * I::point(trial);
    // Check EVERY intermediate, not only final coordinates. Grouping these
    // predicates changes neither numerical acceptance nor refusal semantics.
    let valid = base
        .into_iter()
        .chain(center)
        .chain(middle.into_iter().flatten())
        .chain(residual_raw)
        .chain(residual)
        .chain(offset)
        .chain(world)
        .chain(j.into_iter().flatten())
        .chain(q.into_iter().flatten())
        .chain(row)
        .chain([inclusion])
        .all(I::valid);
    if !valid
        || (0..2).any(|i| !physical(surfaces[i], charts[i], world))
        || !contraction.is_finite()
        || contraction >= 0.5
        || inclusion.hi >= trial
    {
        return None;
    }
    // G_t(y)=y-MF(C(t)+By) is a contraction mapping the CLOSED common
    // entire caller-admissible search square strictly into itself, uniformly
    // in t and literal input bounds.
    // Its unique fixed point is an actual root since finite M is invertible
    // (||I-MJ||<1). Adjacent source cells have the SAME C(t) and the SAME
    // fixed search square at their shared endpoint, so identify the same
    // root by uniqueness, never by welding independently evaluated points.
    let inverse_factor = I::point(1.) / (I::point(1.) - I::point(contraction));
    let radius_bound = I::point(b) * inverse_factor;
    let radius = radius_bound.hi;
    let correction_i: [I; 3] = std::array::from_fn(|k| {
        (I::point(basis[0][k].magnitude()) + I::point(basis[1][k].magnitude())) * I::point(radius)
    });
    let correction = correction_i.map(|x| x.hi);
    let distance = correction
        .into_iter()
        .fold(I::point(0.), |s, x| s + I::point(x).square())
        .sqrt();
    let root_box = std::array::from_fn(|k| {
        base[k]
            + (0..2).fold(I::point(0.), |s, j| {
                s + basis[j][k] * I::new(-radius, radius)
            })
    });
    let width = I::point(parameter[1]) - I::point(parameter[0]);
    let derivative: [I; 3] = std::array::from_fn(|k| {
        bound(&std::array::from_fn::<_, 3, _>(|i| {
            I::point(3.) * (p[i + 1][k] - p[i][k]) / width
        }))
    });
    let rhs_raw = charts.map(|c| dot(c.world_gradient(c.project(root_box)), derivative));
    let rhs = apply(m, rhs_raw);
    let derivative_radius = I::point(norm(rhs)) * inverse_factor;
    let derivative_i: [I; 3] = std::array::from_fn(|k| {
        (I::point(basis[0][k].magnitude()) + I::point(basis[1][k].magnitude()))
            * I::point(derivative_radius.hi)
    });
    let derivative_correction = derivative_i.map(|x| x.hi);
    let valid = [
        inverse_factor,
        radius_bound,
        distance,
        width,
        derivative_radius,
    ]
    .into_iter()
    .chain(correction_i)
    .chain(root_box)
    .chain(derivative)
    .chain(rhs_raw)
    .chain(rhs)
    .chain(derivative_i)
    .all(I::valid);
    if !valid || distance.hi > position_limit {
        return None;
    }
    Some(Section {
        parameter,
        source_box: root_box.map(|x| [x.lo, x.hi]),
        correction_mm: correction,
        derivative_correction_mm: derivative_correction,
        contraction,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    fn axis() -> Axis {
        Axis {
            origin: [0.; 3],
            normal: [0., 0., 1.],
            reference: [1., 0., 0.],
        }
    }
    fn point(p: [[I; 3]; 4], t: f64) -> [I; 3] {
        let mut q = p;
        for level in (1..4).rev() {
            for i in 0..level {
                q[i] =
                    std::array::from_fn(|k| q[i][k] * I::point(1. - t) + q[i + 1][k] * I::point(t));
            }
        }
        q[0]
    }
    fn quarter() -> [[f64; 3]; 4] {
        let k = 4. * (2_f64.sqrt() - 1.) / 3.;
        [[1., 0., 0.], [1., k, 0.], [k, 1., 0.], [0., 1., 0.]]
    }
    #[test]
    fn an_independent_circle_plane_oracle_is_enclosed_by_the_full_implicit_section_proof() {
        let plane = Surface::Plane { axis: axis() };
        let cylinder = Surface::Cylinder {
            axis: axis(),
            radius_mm: 1.,
        };
        let charts = [Chart::new(&plane).unwrap(), Chart::new(&cylinder).unwrap()];
        let raw = quarter();
        let span = CubicTrimSpan {
            parameters: [0., 1.],
            poles_mm: raw,
        };
        let f = Functional::new(&[span]).unwrap();
        let p = split(raw.map(|p| p.map(support::adjacent)))[0];
        let proof = certify(
            [&plane, &cylinder],
            [&charts[0], &charts[1]],
            f,
            p,
            [0., 0.5],
            0.01,
        )
        .unwrap();
        for t in [0., 0.25, 0.5, 0.75, 1.] {
            let c = point(p, t);
            let h = c[1] - c[0];
            let root = (I::point(2.) - h.square()).sqrt();
            let actual = [
                (-h + root) / I::point(2.),
                (h + root) / I::point(2.),
                I::point(0.),
            ];
            for k in 0..3 {
                assert!(
                    proof.source_box[k][0] <= actual[k].lo
                        && actual[k].hi <= proof.source_box[k][1]
                );
                assert!((actual[k] - c[k]).magnitude() <= proof.correction_mm[k]);
            }
        }
        let tightened = proof.correction_mm[0] * 1.1;
        assert!(
            certify(
                [&plane, &cylinder],
                [&charts[0], &charts[1]],
                f,
                p,
                [0., 0.5],
                tightened
            )
            .is_none()
        );
        assert!(proof.contraction < 0.5 && proof.correction_mm.into_iter().any(|v| v > 1e-5));
    }
    #[test]
    fn multiple_nearby_actual_intersection_branches_and_true_tangency_are_refused() {
        let plane = Surface::Plane { axis: axis() };
        let torus = Surface::Torus {
            axis: axis(),
            major_radius_mm: 1.,
            minor_radius_mm: 1e-6,
        };
        let charts = [Chart::new(&plane).unwrap(), Chart::new(&torus).unwrap()];
        let raw = quarter();
        let span = CubicTrimSpan {
            parameters: [0., 1.],
            poles_mm: raw,
        };
        let f = Functional::new(&[span]).unwrap();
        assert!(
            certify(
                [&plane, &torus],
                [&charts[0], &charts[1]],
                f,
                split(raw.map(|p| p.map(support::adjacent)))[0],
                [0., 0.5],
                0.01
            )
            .is_none()
        );
        let tangent = Surface::Plane {
            axis: Axis {
                origin: [1., 0., 0.],
                normal: [1., 0., 0.],
                reference: [0., 0., 1.],
            },
        };
        let cylinder = Surface::Cylinder {
            axis: axis(),
            radius_mm: 1.,
        };
        let charts = [
            Chart::new(&tangent).unwrap(),
            Chart::new(&cylinder).unwrap(),
        ];
        let raw = [
            [1., 0., 0.],
            [1., 0., 1. / 3.],
            [1., 0., 2. / 3.],
            [1., 0., 1.],
        ];
        let f = Functional::new(&[CubicTrimSpan {
            parameters: [0., 1.],
            poles_mm: raw,
        }])
        .unwrap();
        assert!(
            certify(
                [&tangent, &cylinder],
                [&charts[0], &charts[1]],
                f,
                raw.map(|p| p.map(support::adjacent)),
                [0., 1.],
                1e-7
            )
            .is_none()
        );
    }
}
