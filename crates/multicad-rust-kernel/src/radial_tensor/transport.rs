//! Requalify actual rounded source values after an occurrence pose.
use super::*;
#[cfg(test)]
mod tests;

impl RadialBoundary {
    /// Apply a finite proper occurrence frame to every actual source value,
    /// then repeat cold ownership, trim and occupied-material admission.
    /// A successful result owns the rounded world values that later mesh and
    /// exchange proofs inspect; this is not an inferred construction recipe.
    pub fn transformed_with_interrupt(
        &self,
        frame: [Point; 3],
        tolerance: f64,
        minimum: f64,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, Error> {
        if !tolerance.is_finite() || tolerance <= 0. || !minimum.is_finite() || minimum <= 0. {
            return Err(Error::Invalid);
        }
        let [origin, u, v] = frame;
        checked_point(origin)?;
        let dot = |a: Point, b: Point| a.into_iter().zip(b).map(|(a, b)| a * b).sum::<f64>();
        if u.into_iter().chain(v).any(|x| !x.is_finite())
            || (dot(u, u) - 1.).abs() > 1e-12
            || (dot(v, v) - 1.).abs() > 1e-12
            || dot(u, v).abs() > 1e-12
        {
            return Err(Error::Invalid);
        }
        let n = [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ];
        let axes = [u, v, n];
        let vector =
            |p: Point| std::array::from_fn(|k| (0..3).map(|i| axes[i][k] * p[i]).sum::<f64>());
        let point = |p: Point| {
            let p = vector(p);
            let p = std::array::from_fn(|k| origin[k] + p[k]);
            checked_point(p)?;
            Ok(p)
        };
        let axis = |a: Axis| -> Result<Axis, Error> {
            Ok(Axis {
                origin: point(a.origin)?,
                normal: vector(a.normal),
                reference: vector(a.reference),
            })
        };
        let mut work = Work {
            remaining: 4_000_000,
            interrupt,
        };
        work.step()?;
        let mut input = self.input().clone();
        let mut transform_points = |points: &mut [Point]| -> Result<(), Error> {
            for chunk in points.chunks_mut(128) {
                work.charge(chunk.len())?;
                for p in chunk {
                    *p = point(*p)?;
                }
            }
            Ok(())
        };
        transform_points(&mut input.vertices)?;
        for curve in &mut input.curves {
            match curve {
                Curve::Spline { poles, .. } => transform_points(poles)?,
                Curve::Circle { axis: a, .. } => *a = axis(*a)?,
            }
        }
        for surface in &mut input.surfaces {
            match surface {
                Surface::Tensor(t) => {
                    let mut poles = t.poles().to_vec();
                    transform_points(&mut poles)?;
                    *t = RetainedTensor::new(t.spans(), poles)?;
                }
                Surface::Plane(a) | Surface::Cone { axis: a, .. } => *a = axis(*a)?,
            }
        }
        BoundaryOwnership::new(input, tolerance, interrupt)?
            .bind_surfaces(tolerance, interrupt)?
            .certify_axial(tolerance, minimum, interrupt)?
            .certify_caps(interrupt)?
            .certify_ordering(interrupt)?
            .certify_material(tolerance, interrupt)
    }
}
