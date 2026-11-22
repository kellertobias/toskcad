//! Actual affine-radius/angle source curves, with a whole-curve C2 chord proof.
use super::*;
pub struct BoundaryPolyline {
    pub points_mm: Vec<[f64; 3]>,
    pub maximum_error_mm: f64,
}
impl<S> Boundary<S> {
    /// Stock-relative source curve. Every actual segment includes the numerical
    /// endpoint error; the source owns the exact periodic endpoint identity.
    pub fn curve_polyline(
        &self,
        id: u32,
        deflection_mm: f64,
        point_limit: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<BoundaryPolyline, Error> {
        self.curve_polyline_accounted(
            (id, deflection_mm, point_limit),
            self.graph.limits.work,
            0,
            interrupt,
        )
        .map(|(line, _)| line)
    }
    /// Compound callers charge all already-live source packets and output arenas.
    pub(crate) fn curve_polyline_accounted(
        &self,
        request: (u32, f64, usize),
        maximum_work: usize,
        reserved_bytes: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<(BoundaryPolyline, usize), Error> {
        let (id, deflection_mm, point_limit) = request;
        let mut work = Work {
            remaining: maximum_work.min(16_000_000),
            used: 0,
            interrupt,
        };
        work.step()?;
        if !deflection_mm.is_finite()
            || deflection_mm <= 0.
            || !(2..=500_000).contains(&point_limit)
        {
            return Err(Error::Invalid);
        }
        let curve = self
            .graph
            .curves
            .get(id as usize)
            .filter(|c| c.physical())
            .ok_or(Error::Invalid)?;
        let tau = I::new(
            std::f64::consts::TAU.next_down(),
            std::f64::consts::TAU.next_up(),
        );
        let radius = (curve.base[0] + curve.delta[0] * I::new(0., 1.)).magnitude();
        let omega = tau * I::point(curve.delta[2].magnitude());
        // P'' = (r'' - r theta'^2) er + (2 r' theta') et, z''=0.
        // Triangle inequality bounds Euclidean curvature without presuming the
        // nominal trig values or endpoint evaluations are exact.
        let second = (I::point(radius) * omega.square()
            + I::point(2. * curve.delta[0].magnitude()) * omega)
            .hi;
        let proposal = (I::point(second) / (I::point(4.) * I::point(deflection_mm)))
            .hi
            .sqrt()
            .ceil()
            .max(1.);
        if proposal > (point_limit - 1) as f64 {
            return Err(Error::Limit);
        }
        let segments = proposal as usize;
        let mut points = self.graph.buffer::<[f64; 3]>(
            segments + 1,
            self.query_headers()
                .checked_add(size_of::<BoundaryPolyline>())
                .and_then(|n| n.checked_add(reserved_bytes))
                .ok_or(Error::Limit)?,
            &mut work,
        )?;
        let mut previous = self.curve_point(id, 0., interrupt)?;
        points.push(previous.point_mm);
        let mut last_t = 0.;
        let mut maximum = 0_f64;
        for i in 1..=segments {
            work.step()?;
            let t = i as f64 / segments as f64;
            let point = self.curve_point(id, t, interrupt)?;
            let span = I::point(t) - I::point(last_t);
            let error = (I::point(second) * span.square() / I::point(8.)
                + I::point(previous.error_mm.max(point.error_mm)))
            .hi;
            if error > deflection_mm {
                return Err(Error::Unresolved);
            }
            maximum = maximum.max(error);
            points.push(point.point_mm);
            previous = point;
            last_t = t;
        }
        work.step()?;
        Ok((
            BoundaryPolyline {
                points_mm: points,
                maximum_error_mm: maximum,
            },
            work.used,
        ))
    }
}
