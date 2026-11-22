//! Derivative transfer to the ACTUAL normalized literal analytic support.
use super::*;

impl Support {
    /// Bounds derivative differences in the caller's original source charts.
    /// `speeds` encloses the two unplaced source derivative norms, not a leaf
    /// displacement. This calculation does not establish finite trim filling.
    pub(crate) fn chart_derivative_bounds(
        &self,
        pose: &Pose,
        raw: &LiteralSupport,
        speeds: [f64; 2],
    ) -> Result<[f64; 2], Error> {
        if speeds.into_iter().any(|s| !s.is_finite() || s < 0.) {
            return Err(Error::InvalidGeometry);
        }
        if raw.forward != self.forward {
            return Err(Error::Unsupported);
        }
        let literal = Frame::literal(&raw.axis)?;
        let a = self.frame.axes.map(|v| map_vector(pose, v));
        let scaled = |v: BPoint, s: I| v.map(|x| x * s);
        let error = |a, b| norm(difference(a, b));
        let two = I::point(2.);
        let factor = match (self.kind, raw.kind) {
            (SourceKind::Plane, Kind::Plane) => {
                let dot = |v: BPoint| {
                    v.into_iter()
                        .zip(literal.axes[2])
                        .fold(I::point(0.), |s, (x, n)| s + x * n)
                };
                if dot(a[2]).lo <= 0. {
                    return Err(Error::Unsupported);
                }
                // Orthogonal projection onto the actual plane. Original source
                // tangent components each have magnitude <= source speed.
                absolute(dot(a[0])) + absolute(dot(a[1]))
            }
            (SourceKind::Cylinder(r), Kind::Cylinder(actual)) if r.lo > 0. && actual > 0. => {
                let rr = adjacent(actual);
                (error(scaled(a[0], r), scaled(literal.axes[0], rr))
                    + error(scaled(a[1], r), scaled(literal.axes[1], rr)))
                    / r
                    + error(a[2], literal.axes[2])
            }
            (SourceKind::Torus(major, minor), Kind::Torus(rmajor, rminor))
                if (major - minor).lo > 0. && minor.lo > 0. && rmajor > rminor && rminor > 0. =>
            {
                let em = error(
                    scaled(a[0], major),
                    scaled(literal.axes[0], adjacent(rmajor)),
                ) + error(
                    scaled(a[1], major),
                    scaled(literal.axes[1], adjacent(rmajor)),
                );
                let er = error(
                    scaled(a[0], minor),
                    scaled(literal.axes[0], adjacent(rminor)),
                ) + error(
                    scaled(a[1], minor),
                    scaled(literal.axes[1], adjacent(rminor)),
                ) + error(
                    scaled(a[2], minor),
                    scaled(literal.axes[2], adjacent(rminor)),
                );
                // theta speed <= S/(major-minor), minor-circle speed <= S/minor.
                (em + er) / (major - minor) + er / minor
            }
            (SourceKind::Cone(_, k), Kind::Cone(reference, angle))
                if k.lo > 0.
                    && reference >= 0.
                    && angle > 0.
                    && angle < std::f64::consts::FRAC_PI_2 =>
            {
                let angle = adjacent(angle);
                let cot = angle.cos() / angle.sin();
                // Use SAME radial rho from the actual apex. |rho'|<=S and
                // |rho theta'|<=S, including the limiting apex; no 1/rho.
                two * (error(a[0], literal.axes[0]) + error(a[1], literal.axes[1]))
                    + error(scaled(a[2], I::point(1.) / k), scaled(literal.axes[2], cot))
            }
            _ => return Err(Error::Unsupported),
        };
        if !factor.valid() {
            return Err(Error::Unsupported);
        }
        let bounds = speeds.map(|s| (I::point(s) * factor).hi);
        if bounds.into_iter().any(|s| !s.is_finite() || s < 0.) {
            return Err(Error::Unsupported);
        }
        Ok(bounds)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn analytic_coefficients_enclose_original_chart_derivative_transfer() {
        let pose = Pose::new([0.; 3], [1., 0., 0.], [0., 1., 0.]).unwrap();
        let sources = [
            SourceKind::Plane,
            SourceKind::Cylinder(I::point(2.)),
            SourceKind::Torus(I::point(4.), I::point(1.)),
            SourceKind::Cone(I::point(2.), I::point(1.)),
        ];
        for kind in sources {
            let source = Support {
                frame: Frame::new(
                    [I::point(0.); 3],
                    [I::point(0.), I::point(0.), I::point(1.)],
                    [I::point(1.), I::point(0.), I::point(0.)],
                )
                .unwrap(),
                kind,
                forward: true,
            };
            let raw = source.literal(&pose);
            let same = source
                .chart_derivative_bounds(&pose, &raw, [2., 3.])
                .unwrap();
            assert!(same.into_iter().all(|x| (0. ..1e-10).contains(&x)));
            let mut moved = pose;
            moved.axes = moved.axes.map(|v| v.map(|x| x * 2.));
            let doubled = source
                .chart_derivative_bounds(&moved, &raw, [2., 3.])
                .unwrap();
            // The source affine changes by identity. Each actual source tangent
            // therefore changes by exactly its norm (plane z-projection excepted).
            if matches!(kind, SourceKind::Plane) {
                assert!(doubled.into_iter().all(|x| x < 1e-10));
            } else {
                assert!(doubled[0] >= 2. && doubled[1] >= 3.);
            }
            assert_eq!(
                source.chart_derivative_bounds(&pose, &raw, [f64::NAN, 1.]),
                Err(Error::InvalidGeometry)
            );
            assert_eq!(
                source.chart_derivative_bounds(&pose, &raw, [-1., 1.]),
                Err(Error::InvalidGeometry)
            );
            let mut reversed = raw;
            reversed.forward = false;
            assert_eq!(
                source.chart_derivative_bounds(&pose, &reversed, [1.; 2]),
                Err(Error::Unsupported)
            );
            let incompatible = LiteralSupport {
                kind: Kind::Cylinder(-1.),
                ..raw
            };
            assert_eq!(
                source.chart_derivative_bounds(&pose, &incompatible, [1.; 2]),
                Err(Error::Unsupported)
            );
        }
    }
}
