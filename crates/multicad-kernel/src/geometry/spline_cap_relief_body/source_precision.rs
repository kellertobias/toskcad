//! Outward source-map displacement and complete derivative-cone receipts.
use super::*;
#[derive(Clone, Copy, Default)]
pub(super) struct Receipt {
    pub distance_mm: f64,
    pub angle_rad: f64,
}
fn sum(a: f64, b: f64) -> (f64, f64) {
    let v = a + b;
    let p = v - a;
    (v, (a - (v - p)) + (b - p))
}
pub(super) fn add_up(a: f64, b: f64) -> f64 {
    let (v, e) = sum(a, b);
    if e > 0. { v.next_up() } else { v }
}
pub(super) fn remaining(requested: f64, spent: f64) -> f64 {
    let (value, error) = sum(requested, -spent);
    if error < 0. { value.next_down() } else { value }
}

fn mul_up(a: f64, b: f64) -> f64 {
    let v = a * b;
    if a != 0. && b != 0. && v.abs() < 1e-280 || a.mul_add(b, -v) > 0. {
        v.next_up()
    } else {
        v
    }
}
fn subtract(a: f64, b: f64) -> (f64, f64) {
    sum(a, -b)
}
fn mapped(a: [[f64; 2]; 2], p: [f64; 2]) -> [f64; 2] {
    a.map(|row| row[0] * p[0] + row[1] * p[1])
}
fn ratio_angle(r: f64) -> Result<f64, GeometryError> {
    if !r.is_finite() || r >= 0.25 {
        return Err(fail(
            GeometryErrorKind::Unsupported,
            "The source tangent map cannot be certified",
        ));
    }
    if r == 0. {
        return Ok(0.);
    }
    Ok((mul_up(2., r) / (1. - r).next_down()).next_up())
}
fn derivative_ratio(
    p: [[f64; 2]; 2],
    q: [[f64; 2]; 2],
    axis: [f64; 2],
) -> Result<f64, GeometryError> {
    let mut error = 0.;
    let mut projection = 0.;
    let mut projection_error = 0.;
    for i in 0..2 {
        let (source, se) = subtract(p[1][i], p[0][i]);
        let (target, te) = subtract(q[1][i], q[0][i]);
        let (difference, de) = subtract(target, source);
        error = add_up(
            error,
            add_up(
                difference.abs(),
                add_up(se.abs(), add_up(te.abs(), de.abs())),
            ),
        );
        let term = source * axis[i];
        let mut rounding = source.mul_add(axis[i], -term).abs();
        if source != 0. && axis[i] != 0. && term.abs() < 1e-280 {
            rounding = rounding.next_up();
        }
        let (next, loss) = sum(projection, term);
        projection = next;
        projection_error = add_up(
            projection_error,
            add_up(
                rounding,
                add_up(loss.abs(), mul_up(se.abs(), axis[i].abs())),
            ),
        );
    }
    if error == 0. {
        return Ok(0.);
    }
    let lower = (projection - projection_error).next_down();
    if lower <= 0. {
        return Err(fail(
            GeometryErrorKind::Unsupported,
            "The source derivative cone is unresolved",
        ));
    }
    let norm = axis[0].hypot(axis[1]).next_up();
    Ok((mul_up(error, norm) / lower).next_up())
}
pub(super) fn certify(
    source: &Edge,
    target: &Edge,
    a: [[f64; 2]; 2],
    point_error: f64,
) -> Result<Receipt, GeometryError> {
    let mut ratio = 0_f64;
    let angle = match (source, target) {
        (
            Edge::Analytic(AnalyticEdge::Line { a: p, b: r }),
            Edge::Analytic(AnalyticEdge::Line { a: q, b: s }),
        ) => {
            let p = [mapped(a, *p), mapped(a, *r)];
            let axis = std::array::from_fn(|i| p[1][i] - p[0][i]);
            ratio_angle(derivative_ratio(p, [*q, *s], axis)?)?
        }
        (Edge::Cubic(p), Edge::Cubic(q)) => {
            if p.knots() != q.knots()
                || p.poles().len() != q.poles().len()
                || p.reversed() != q.reversed()
            {
                return Err(fail(
                    GeometryErrorKind::InvalidResult,
                    "The source spline basis changed during mapping",
                ));
            }
            for (p, q) in p.bezier_spans().zip(q.bezier_spans()) {
                let p = p.map(|p| mapped(a, p));
                let axis = std::array::from_fn(|i| p[3][i] - p[0][i]);
                for i in 0..3 {
                    ratio = ratio.max(derivative_ratio([p[i], p[i + 1]], [q[i], q[i + 1]], axis)?);
                }
            }
            ratio_angle(ratio)?
        }
        (
            Edge::Analytic(AnalyticEdge::Circle { .. }),
            Edge::Analytic(AnalyticEdge::Circle { .. }),
        ) => 0.,
        (
            Edge::Analytic(AnalyticEdge::Arc {
                a: p,
                b: r,
                center: c,
                ..
            }),
            Edge::Analytic(AnalyticEdge::Arc {
                a: q,
                b: s,
                center: d,
                ..
            }),
        ) => {
            // Exact two-sum radial transport avoids inferring angular fidelity
            // from a small position error for a tiny, far-translated arc.
            for (p, q) in [(*p, *q), (*r, *s)] {
                let source = [subtract(p[0], c[0]), subtract(p[1], c[1])];
                let target = [subtract(q[0], d[0]), subtract(q[1], d[1])];
                for i in 0..2 {
                    let j = if a[i][0] != 0. { 0 } else { 1 };
                    let sign = a[i][j];
                    if target[i] != (sign * source[j].0, sign * source[j].1) {
                        return Err(fail(
                            GeometryErrorKind::Unsupported,
                            "The source arc radial map cannot be certified",
                        ));
                    }
                }
            }
            128. * f64::EPSILON
        }
        _ => {
            return Err(fail(
                GeometryErrorKind::InvalidResult,
                "A source curve changed class during cardinal mapping",
            ));
        }
    };
    let distance_mm = if let Edge::Analytic(edge @ AnalyticEdge::Arc { .. }) = source {
        add_up(
            point_error,
            mul_up(edge.radius_mm().expect("arc radius"), angle),
        )
    } else {
        point_error
    };
    Ok(Receipt {
        distance_mm,
        angle_rad: angle,
    })
}
#[cfg(test)]
mod tests {
    use super::*;
    use multicad_rust_kernel::spline_sweep::Cubic;
    const I: [[f64; 2]; 2] = [[1., 0.], [0., 1.]];
    fn tangent(p: [[f64; 2]; 4], t: f64) -> [f64; 2] {
        std::array::from_fn(|i| {
            3. * ((1. - t).powi(2) * (p[1][i] - p[0][i])
                + 2. * t * (1. - t) * (p[2][i] - p[1][i])
                + t * t * (p[3][i] - p[2][i]))
        })
    }
    #[test]
    fn full_derivative_receipt_bounds_dense_angles_and_refuses_far_uv_conditioning() {
        assert_eq!(remaining(1., 0.1), 0.9_f64.next_down());
        assert_eq!(remaining(1., 0.25), 0.75);
        let p = [[0., 0.], [0.1, 0.2], [0.4, 0.3], [0.5, 0.4]];
        let q = p.map(|p| [p[0] + 5., p[1] + 7.]);
        let edge =
            |p: [[f64; 2]; 4]| Edge::Cubic(Cubic::new(p.to_vec(), vec![0., 1.], false).unwrap());
        let receipt = certify(&edge(p), &edge(q), I, 2e-15).unwrap();
        assert!(receipt.angle_rad > 0.);
        for i in 0..=1000 {
            let t = i as f64 / 1000.;
            let a = tangent(p, t);
            let b = tangent(q, t);
            let angle = (a[0] * b[1] - a[1] * b[0])
                .abs()
                .atan2(a[0] * b[0] + a[1] * b[1]);
            assert!(angle <= receipt.angle_rad);
        }
        let small = p.map(|p| p.map(|x| x * 1e-5));
        let far = small.map(|p| [p[0] + 1e5, p[1] + 1e5]);
        let receipt = certify(&edge(small), &edge(far), I, 1e-9).unwrap();
        assert!(receipt.angle_rad > 1e-6);
        let flat = [[0., 0.], [0., 0.], [1., 1.], [2., 1.]];
        assert!(certify(&edge(flat), &edge(flat), I, 0.).is_ok());
        assert_eq!(
            derivative_ratio([[0., 0.], [0., 1.]], [[0., 0.], [0.1, 1.]], [1., 0.])
                .err()
                .unwrap()
                .kind,
            GeometryErrorKind::Unsupported
        );
        assert_eq!(
            ratio_angle(0.25).err().unwrap().kind,
            GeometryErrorKind::Unsupported
        );
        let arc = Edge::Analytic(AnalyticEdge::Arc {
            a: [1., 0.],
            b: [0., 1.],
            center: [0.; 2],
            clockwise: false,
        });
        assert!(certify(&arc, &arc, I, 0.).is_ok());
        let shifted = Edge::Analytic(AnalyticEdge::Arc {
            a: [1e7 + 1., 1e7],
            b: [1e7, 1e7 + 1.000001],
            center: [1e7, 1e7],
            clockwise: false,
        });
        assert!(certify(&arc, &shifted, I, 1e-9).is_err());
    }
}
