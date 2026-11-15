//! Whole-source comparison, including independent endpoint-ray corrections.
//! C0 and C1 use the same normalized parameter per paired piece, not an
//! endpoint fit or a sampled curve. Rotation about the common axis preserves
//! these Euclidean bounds; world placement requires a separate proof.
use super::*;
#[derive(Clone, Copy, Debug)]
pub struct RepairedHeadCorrespondence {
    pub position_mm: f64,
    /// Same normalized meridian parameter in each corresponding source piece.
    pub derivative_mm: f64,
    pub work_used: usize,
}
fn magnitude(v: [I; 2]) -> I {
    I::point(v[0].magnitude()) + I::point(v[1].magnitude())
}
impl RepairedHeadMeridian<'_> {
    pub fn translated_correspondence(
        &self,
        other: &RepairedHeadMeridian<'_>,
        source_z_shift_mm: f64,
        target_z_shift_mm: f64,
        maximum_work: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<RepairedHeadCorrespondence, RepairedHeadError> {
        poll(interrupt)?;
        if maximum_work > WORK
            || [source_z_shift_mm, target_z_shift_mm]
                .into_iter()
                .any(|v| !v.is_finite() || v.abs() > 1e7)
        {
            return Err(RepairedHeadError::Invalid);
        }
        // Fixed <=5-piece source, no allocation/subdivision. Per piece this
        // bounds <128 interval scalar operations (each <=16 scalar steps),
        // fixed-size array traversal and the role/kind/finiteness checks.
        // Charge a conservative 4096 primitive steps even for line pieces.
        let required = 1 + self.count * 4096;
        if maximum_work < required {
            return Err(RepairedHeadError::Limit);
        }
        ensure(self.count == other.count)?;
        let shift = I::point(source_z_shift_mm) - I::point(target_z_shift_mm);
        let mut result = RepairedHeadCorrespondence {
            position_mm: 0.,
            derivative_mm: 0.,
            work_used: required,
        };
        for index in 0..self.count {
            poll(interrupt)?;
            ensure(
                self.certificates[index].unwrap().role == other.certificates[index].unwrap().role,
            )?;
            let [position, derivative] =
                match (self.parts[index].unwrap(), other.parts[index].unwrap()) {
                    (Part::Line(a), Part::Line(b)) => {
                        let delta = |a: [f64; 2], b: [f64; 2]| {
                            [
                                I::point(a[0]) - I::point(b[0]),
                                I::point(a[1]) - I::point(b[1]) - shift,
                            ]
                        };
                        let first = delta(a.start, b.start);
                        let last = delta(a.end, b.end);
                        [
                            I::point(magnitude(first).hi.max(magnitude(last).hi)),
                            magnitude(std::array::from_fn(|i| last[i] - first[i])),
                        ]
                    }
                    (Part::Arc(a), Part::Arc(b)) => a.correspondence(b, shift)?,
                    _ => return Err(RepairedHeadError::Unresolved),
                };
            ensure([position, derivative].into_iter().all(finite))?;
            result.position_mm = result.position_mm.max(position.hi);
            result.derivative_mm = result.derivative_mm.max(derivative.hi);
        }
        poll(interrupt)?;
        Ok(result)
    }
}

#[cfg(test)]
mod tests;
