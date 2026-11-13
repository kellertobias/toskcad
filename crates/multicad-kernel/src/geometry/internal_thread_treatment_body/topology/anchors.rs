//! Transport retained literal anchors with the actual two placement matrices.
//! No fitted inverse, projection onto a nominal isometry or normalization.
use super::super::interval::Interval as B;
use super::*;
use multicad_core::references::ReferenceFrame;

type Matrix = [[B; 3]; 3];
fn dot(a: [B; 3], b: [B; 3]) -> B {
    (0..3).fold(B::ZERO, |s, i| s + a[i] * b[i])
}
fn cross(a: [B; 3], b: [B; 3]) -> [B; 3] {
    [(1, 2), (2, 0), (0, 1)].map(|(j, k)| a[j] * b[k] - a[k] * b[j])
}
fn axes(frame: PlaneFrame) -> Matrix {
    [frame.u, frame.v, frame.normal()].map(|a| a.map(B::point))
}
fn inverse(columns: Matrix) -> Result<Matrix, GeometryError> {
    let cofactors = [
        cross(columns[1], columns[2]),
        cross(columns[2], columns[0]),
        cross(columns[0], columns[1]),
    ];
    let determinant = dot(columns[0], cofactors[0]);
    if !determinant.valid() || determinant.lo <= 0. {
        return Err(unavailable(
            "Retained stock anchor frame has no proved proper inverse",
        ));
    }
    let result = cofactors.map(|row| row.map(|v| v / determinant));
    if result.into_iter().flatten().all(B::valid) {
        Ok(result)
    } else {
        Err(unavailable(
            "Retained stock anchor inverse exceeds numerical limits",
        ))
    }
}
fn midpoint(bounds: [B; 3], tolerance: f64) -> Result<[f64; 3], GeometryError> {
    if !tolerance.is_finite() || tolerance <= 0. || bounds.into_iter().any(|b| !b.valid()) {
        return Err(unavailable(
            "Retained anchor transport has no finite error bound",
        ));
    }
    let value = bounds.map(|b| b.lo.midpoint(b.hi));
    let error = bounds.into_iter().zip(value).fold(0_f64, |e, (b, v)| {
        (e + (v - b.lo).next_up().max((b.hi - v).next_up())).next_up()
    });
    if error > tolerance {
        Err(unavailable(
            "Retained anchor transport reaches the precision floor",
        ))
    } else {
        Ok(value)
    }
}
pub(super) struct Transport {
    original: PlaneFrame,
    current: PlaneFrame,
    inverse: Matrix,
    columns: Matrix,
    policy: GeometryPolicy,
}
impl Transport {
    pub fn new(
        original: PlaneFrame,
        current: PlaneFrame,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        policy.validate()?;
        original.checked(policy)?;
        current.checked(policy)?;
        Ok(Self {
            original,
            current,
            inverse: inverse(axes(original))?,
            columns: axes(current),
            policy,
        })
    }
    fn vector(&self, vector: [f64; 3]) -> [B; 3] {
        let local = self.inverse.map(|row| dot(row, vector.map(B::point)));
        std::array::from_fn(|i| dot(self.columns.map(|a| a[i]), local))
    }
    pub fn local_point(&self, point: [source::Enclosure; 3]) -> Result<[f64; 3], GeometryError> {
        let point = point.map(B::from);
        let bounds = std::array::from_fn(|i| {
            B::point(self.current.origin[i]) + dot(self.columns.map(|a| a[i]), point)
        });
        let result = midpoint(bounds, self.policy.linear_tolerance_mm / 8.)?;
        self.policy.point(result)?;
        Ok(result)
    }
    pub fn frame(&self, anchor: ReferenceFrame) -> Result<ReferenceFrame, GeometryError> {
        anchor.validate().map_err(unavailable)?;
        self.policy.point(anchor.origin_mm)?;
        if self.original == self.current {
            return Ok(anchor);
        }
        let delta = std::array::from_fn(|i| {
            B::point(anchor.origin_mm[i]) - B::point(self.original.origin[i])
        });
        let local = self.inverse.map(|row| dot(row, delta));
        let p = std::array::from_fn(|i| {
            B::point(self.current.origin[i]) + dot(self.columns.map(|a| a[i]), local)
        });
        let result = ReferenceFrame {
            origin_mm: midpoint(p, self.policy.linear_tolerance_mm / 8.)?,
            u: midpoint(self.vector(anchor.u), 1e-12)?,
            v: midpoint(self.vector(anchor.v), 1e-12)?,
        };
        self.policy.point(result.origin_mm)?;
        result.validate().map_err(unavailable)?;
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn at(frame: PlaneFrame, p: [f64; 3]) -> [f64; 3] {
        let n = frame.normal();
        std::array::from_fn(|i| {
            frame.origin[i] + p[0] * frame.u[i] + p[1] * frame.v[i] + p[2] * n[i]
        })
    }
    #[test]
    fn literal_stock_anchors_keep_offset_orientation_and_proper_axis_reversal() {
        let old = PlaneFrame {
            origin: [32., -16., 8.],
            u: [0.8, 0.6, 0.],
            v: [-0.36, 0.48, 0.8],
        };
        let local = [0.25, -0.5, 1.5];
        let anchor = ReferenceFrame {
            origin_mm: at(old, local),
            u: old.u,
            v: old.v,
        };
        let identity = Transport::new(old, old, GeometryPolicy::default()).unwrap();
        assert_eq!(identity.frame(anchor).unwrap(), anchor);
        for current in [
            PlaneFrame {
                origin: [-16., 32., 64.],
                u: [0., 1., 0.],
                v: [0., 0., 1.],
            },
            PlaneFrame {
                origin: [1., 2., 3.],
                u: [-1., 0., 0.],
                v: [0., 1., 0.],
            },
            PlaneFrame {
                u: [1. + 2_f64.powi(-44), 0., 0.],
                ..PlaneFrame::default()
            },
        ] {
            let transported = Transport::new(old, current, GeometryPolicy::default())
                .unwrap()
                .frame(anchor)
                .unwrap();
            let expected = at(current, local);
            for (i, coordinate) in expected.into_iter().enumerate() {
                assert!((transported.origin_mm[i] - coordinate).abs() < 1e-10);
                assert!((transported.u[i] - current.u[i]).abs() < 1e-12);
                assert!((transported.v[i] - current.v[i]).abs() < 1e-12);
            }
            transported.validate().unwrap();
        }
    }
    #[test]
    fn invalid_singular_and_unrepresentable_anchor_transports_refuse() {
        let p = GeometryPolicy::default();
        let good = PlaneFrame::default();
        assert!(Transport::new(PlaneFrame { u: [0.; 3], ..good }, good, p).is_err());
        assert!(
            Transport::new(
                good,
                PlaneFrame {
                    origin: [f64::INFINITY, 0., 0.],
                    ..good
                },
                p
            )
            .is_err()
        );
        assert!(inverse([[B::ZERO; 3]; 3]).is_err());
        assert!(
            inverse([
                [B::ONE, B::ZERO, B::ZERO],
                [B::ZERO, B::ONE, B::ZERO],
                [B::ZERO, B::ZERO, B::point(f64::from_bits(1))]
            ])
            .is_err()
        );
        assert!(inverse([[B::new(f64::NEG_INFINITY, f64::INFINITY); 3]; 3]).is_err());
        assert!(midpoint([B::new(0., 1.); 3], 1e-7).is_err());
        assert!(midpoint([B::ZERO; 3], 0.).is_err());
        let t = Transport::new(good, good, p).unwrap();
        assert!(
            t.frame(ReferenceFrame {
                origin_mm: [0.; 3],
                u: [0.; 3],
                v: [0., 1., 0.]
            })
            .is_err()
        );
    }
}
