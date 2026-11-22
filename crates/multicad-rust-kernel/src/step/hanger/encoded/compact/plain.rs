//! Actual plain-wall tensor charts, using the SAME literal angular graph.
use super::*;
pub struct PlainWall<'b, 'c, 'd, 'f, 's> {
    source: &'b CompactBoundary<'c, 'd, 'f, 's>,
    face: usize,
    points: Vec<[EncodedScalar; 3]>,
    knots: Vec<f64>,
    error: [f64; 2],
    work: usize,
    bytes: usize,
}
impl PlainWall<'_, '_, '_, '_, '_> {
    pub fn source(&self) -> &CompactBoundary<'_, '_, '_, '_> {
        self.source
    }
    pub fn face(&self) -> usize {
        self.face
    }
    pub fn points(&self) -> &[[EncodedScalar; 3]] {
        &self.points
    }
    pub fn dimensions(&self) -> [usize; 2] {
        [2, self.points.len() / 2]
    }
    pub fn degrees(&self) -> [usize; 2] {
        [1, 10]
    }
    /// Global positive physical turns, agreeing with all cardinal graph owners.
    pub fn knots(&self) -> &[f64] {
        &self.knots
    }
    pub fn source_bounds(&self) -> [f64; 2] {
        self.error
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_bytes(&self) -> usize {
        self.bytes
    }
}
impl<'c, 'd, 'f, 's> CompactBoundary<'c, 'd, 'f, 's> {
    pub fn plain_wall<'b>(
        &'b self,
        face: usize,
        options: DomainOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Option<PlainWall<'b, 'c, 'd, 'f, 's>>, Error> {
        let owner = self.faces().get(face).ok_or(Error::Invalid)?;
        let ChartFaceKind::PlainWall { quarter } = owner.kind else {
            return Ok(None);
        };
        let field = self.surfaces().field();
        let inherited = self.retained_bytes()
            + self.surfaces().retained_bytes()
            + self.source().retained_bytes();
        let (mut w, mut m) = crate::step::hanger::domain::ledgers(options, inherited, interrupt)?;
        w.charge(self.work_used())?;
        m.add(SCRATCH)?;
        let base = m.used;
        m.add(size_of::<PlainWall<'_, '_, '_, '_, '_>>())?;
        let angular = &field.angular;
        let a = f64::from(quarter) / 4.;
        let b = a + 0.25;
        w.charge(angular.knots().len() * 2)?;
        let first = angular
            .knots()
            .iter()
            .position(|&k| k == a)
            .ok_or(Error::Unresolved)?;
        let last = angular
            .knots()
            .iter()
            .position(|&k| k == b)
            .ok_or(Error::Unresolved)?;
        if last <= first {
            return Err(Error::Unresolved);
        }
        let mut knots = m.allocate(last - first + 1, &mut w)?;
        knots.extend_from_slice(&angular.knots()[first..=last]);
        let width = (last - first) * 10 + 1;
        let mut points = m.allocate(2 * width, &mut w)?;
        let source = field.source.source_boundary().source();
        let radius = binary(source.spec().plain_radius_mm, &mut w)?;
        let bottom = N::integer(0).sub(&binary(source.spec().plain_length_mm, &mut w)?, &mut w)?;
        for end in 0..2 {
            let z = if end == 0 { &bottom } else { &N::integer(0) };
            for piece in first..last {
                let controls = &angular.coefficients()[piece];
                for q in controls
                    .iter()
                    .take(if piece + 1 == last { 11 } else { 10 })
                {
                    let x = radius.mul(&N(q[0].clone()), &mut w)?;
                    let y = radius.mul(&N(q[1].clone()), &mut w)?;
                    points.push([
                        EncodedScalar::retain(&x, &mut m, &mut w)?,
                        EncodedScalar::retain(&y, &mut m, &mut w)?,
                        EncodedScalar::retain(z, &mut m, &mut w)?,
                    ]);
                }
            }
        }
        let r = radius.enclosure(&mut w)?;
        let error = [
            (r * p(angular.position_error())).hi,
            (r * p(angular.derivative_error())).hi,
        ];
        w.charge(1)?;
        Ok(Some(PlainWall {
            source: self,
            face,
            points,
            knots,
            error,
            work: w.used,
            bytes: m.used - base,
        }))
    }
}
#[cfg(test)]
mod tests;
