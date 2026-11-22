//! Independently measured C0 coupling of retained polynomial rims to their
//! actual Plane/Circle supports, not a zero-gap sewing or tangent certificate.
use super::*;
use crate::step::helical::DEGREE;
use crate::step_read::phillips_cutter::physical::supports::chart::Chart;
pub(super) struct Caps {
    tip: Chart,
    crests: Vec<Circle>,
}
struct Circle {
    frame: support::Frame,
    radius: I,
}
impl Circle {
    fn local(&self, p: Point) -> [I; 3] {
        let d = std::array::from_fn::<_, 3, _>(|i| support::adjacent(p[i]) - self.frame.origin[i]);
        self.frame.axes.map(|a| {
            a.into_iter()
                .zip(d)
                .fold(I::point(0.), |sum, (a, d)| sum + a * d)
        })
    }
    fn distance(&self, points: &[Point], work: &mut Work<'_>) -> Result<f64> {
        require(points.len() > DEGREE && (points.len() - 1).is_multiple_of(DEGREE))?;
        let mut bound = 0_f64;
        for poles in points.windows(DEGREE + 1).step_by(DEGREE) {
            let p: Vec<_> = poles.iter().map(|&p| self.local(p)).collect();
            let mut radial = 0_f64;
            for k in 0..=2 * DEGREE {
                let mut squared = I::point(0.);
                for i in k.saturating_sub(DEGREE)..=DEGREE.min(k) {
                    work.step(1)?;
                    let weight = I::point(choose(DEGREE, i)) * I::point(choose(DEGREE, k - i))
                        / I::point(choose(2 * DEGREE, k));
                    squared = squared + (p[i][0] * p[k - i][0] + p[i][1] * p[k - i][1]) * weight;
                }
                radial = radial.max(((squared - self.radius.square()) / self.radius).magnitude());
            }
            let axial = p.iter().map(|p| p[2].magnitude()).fold(0_f64, f64::max);
            // |rho-r|=|rho²-r²|/(rho+r), and rho+r>=r>0.
            let distance = (I::point(radial).square() + I::point(axial).square()).sqrt();
            require(distance.valid() && distance.hi.is_finite())?;
            bound = bound.max(distance.hi);
        }
        Ok(bound)
    }
}
// All callers use n<=22, so these exact integer binomials and products remain
// below 2^53; interval division encloses the rational Bernstein weight.
fn choose(n: usize, k: usize) -> f64 {
    let mut result = 1_u64;
    for i in 0..k {
        result = result * (n - i) as u64 / (i + 1) as u64;
    }
    result as f64
}
impl Caps {
    pub(super) fn new(body: &Body<'_>, shell: &Shell<'_, '_>) -> Result<Self> {
        let (tip, forward, _) =
            helical::face(shell, body.material_faces[body.tensor_face_count + 1])?;
        let tip = Chart::new(&literal::surface(shell, tip, forward)?)?;
        let mut crests = Vec::with_capacity(body.tensor_face_count / 4);
        for &face in body.tensor_face_entities().iter().skip(3).step_by(4) {
            let (_, forward, bounds) = helical::face(shell, face)?;
            let wire = helical::orient(&helical::ring(shell, bounds[0], true)?, forward);
            let edge = wire.first().ok_or(ReadError::Invalid)?.0;
            let args = shell.args(edge, "EDGE_CURVE")?;
            let curve = reference(args.get(3).ok_or(ReadError::Invalid)?)?;
            let [Value::String(_), axis, Value::Number(radius)] = shell.args(curve, "CIRCLE")?
            else {
                return Err(ReadError::Invalid);
            };
            let frame = support::Frame::literal(&geometry::axis(shell, reference(axis)?)?)
                .map_err(source_error)?;
            let radius = support::adjacent(*radius);
            require(radius.valid() && radius.lo > 0.)?;
            crests.push(Circle { frame, radius });
        }
        Ok(Self { tip, crests })
    }
    pub(super) fn measure(
        &self,
        index: usize,
        controls: &[Vec<Point>],
        visible_start: usize,
        work: &mut Work<'_>,
    ) -> Result<(f64, Option<f64>)> {
        let mut tip = 0_f64;
        for row in controls {
            work.step(1)?;
            let p = *row.last().ok_or(ReadError::Invalid)?;
            let residual = self.tip.local(p)[2].magnitude();
            require(residual.is_finite())?;
            tip = tip.max(residual);
        }
        let circle = if index % 4 == 3 {
            let points: Vec<_> = controls.iter().map(|row| row[visible_start]).collect();
            Some(self.crests[index / 4].distance(&points, work)?)
        } else {
            None
        };
        Ok((tip, circle))
    }
}

#[cfg(test)]
mod tests;
