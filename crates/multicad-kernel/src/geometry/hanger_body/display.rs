//! Whole-source C2 chord bounds with separately bounded emitted point error.
use super::*;
use crate::geometry::internal_thread_treatment_body::interval::Interval as I;
impl HangerBody {
    pub(in crate::geometry) fn display_edges(
        &self,
        frame: PlaneFrame,
        deflection: f64,
        policy: GeometryPolicy,
    ) -> Result<Vec<DisplayEdge>, GeometryError> {
        self.material(frame, policy)?;
        if !deflection.is_finite() || deflection <= 0. {
            return Err(input("Invalid hanger curve display quality"));
        }
        let mut remaining = 500_000usize;
        let mut output = Vec::with_capacity(self.boundary.curves().len());
        for (index, curve) in self.boundary.curves().iter().enumerate() {
            let full = self
                .boundary
                .curve_bounds(curve.key, [0., 1.], None)
                .map_err(super::super::hanger_material::failure)?;
            let m2 = full
                .second_derivative_mm
                .map(|[lo, hi]| I::new(lo, hi))
                .into_iter()
                .fold(I::ZERO, |s, x| s + I::point(x.magnitude()));
            let mut segments = 1usize;
            loop {
                // Admitted raw affine axes have Frobenius norm below3.
                // M2*h²/8 bounds a complete C2 source curve against its
                // endpoint chord. Half the allowance remains for literal
                // source evaluation and world-placement rounding.
                let n = I::point(segments as f64);
                let chord = I::point(3.) * m2 / (I::point(8.) * n * n);
                if chord.hi <= deflection / 2. {
                    break;
                }
                // Every previous iteration has segments < remaining <=500k;
                // doubling therefore fits the supported usize widths.
                segments *= 2;
                if segments >= remaining {
                    return Err(topology::refuse(
                        "Hanger display subdivision exceeds its point budget",
                    ));
                }
            }
            if segments >= remaining {
                return Err(topology::refuse("Hanger display exceeds its point budget"));
            }
            remaining -= segments + 1;
            let mut points = Vec::with_capacity(segments + 1);
            for i in 0..=segments {
                // A power-of-two denominator keeps every sampled parameter
                // exactly representable within this bounded grid.
                let t = i as f64 / segments as f64;
                let p = self
                    .boundary
                    .curve_bounds(curve.key, [t, t], None)
                    .map_err(super::super::hanger_material::failure)?
                    .point_mm;
                points.push(topology::world_point(frame, p, deflection / 2., policy)?);
            }
            output.push(DisplayEdge {
                index: index as u32,
                points_mm: points,
            });
        }
        Ok(output)
    }
}
