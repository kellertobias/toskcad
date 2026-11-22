//! Bounded quadrature of the nominal literal supports, separate from transport
//! and the authored tangent witness. These are not decoded-file mass receipts.
use super::*;
#[cfg(test)]
mod tests;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NominalSurfaceEvaluation {
    pub area_mm2: [f64; 2],
    /// XY in the contour chart; Z is positive depth below the open mouth.
    pub center_mm: [[f64; 2]; 3],
}
#[derive(Clone, Debug, PartialEq)]
pub struct NominalSurfaceEvaluations {
    pub mouth: NominalSurfaceEvaluation,
    pub walls: [NominalSurfaceEvaluation; 12],
    pub floor: NominalSurfaceEvaluation,
    /// Kept separate: quadrature bounds evaluate nominal literal supports.
    pub source_error_mm: f64,
}

fn bounds(value: I) -> [f64; 2] {
    [value.lo, value.hi]
}
fn evaluation(area: I, moments: [I; 3], scale: f64) -> Result<NominalSurfaceEvaluation, Error> {
    if !area.valid() || area.lo <= 0. {
        return Err(Error::Unresolved);
    }
    let center = moments.map(|v| v / area * I::point(scale));
    if !center.into_iter().all(I::valid) {
        return Err(Error::Unresolved);
    }
    Ok(NominalSurfaceEvaluation {
        area_mm2: bounds(area),
        center_mm: center.map(bounds),
    })
}

fn integrate(
    arc: Arc,
    width: f64,
    remaining: &mut usize,
    cancelled: &dyn Fn() -> bool,
    value: &dyn Fn(I) -> [Jet; 4],
) -> Result<[I; 4], Error> {
    let a = arc.start_angle_rad;
    let b = a + arc.sweep_rad;
    let end = I::point(a) + I::point(arc.sweep_rad);
    let tail_width = (end - I::point(b)).magnitude();
    let tail = value(I::new(end.lo.min(b), end.hi.max(b)));
    let mut result = tail.map(|v| {
        let error = (I::point(tail_width) * I::point(v.0[0].magnitude())).hi;
        I::new(-error, error)
    });
    let mut pending = vec![(a.min(b), a.max(b))];
    while let Some((lo, hi)) = pending.pop() {
        if cancelled() {
            return Err(Error::Cancelled);
        }
        if *remaining == 0 {
            return Err(Error::Limit);
        }
        *remaining -= 1;
        let span = I::point(hi) - I::point(lo);
        let (half, [first, last]) = gauss_nodes(lo, hi);
        let first = value(first);
        let last = value(last);
        let whole = value(I::new(lo, hi));
        let integral = std::array::from_fn::<_, 4, _>(|i| {
            let remainder = (I::point(whole[i].0[4].magnitude()) * span.square().square() * span
                / I::point(180.))
            .hi;
            half * (first[i].0[0] + last[i].0[0]) + I::new(-remainder, remainder)
        });
        let budget = width * ((hi - lo) / arc.sweep_rad.abs());
        if integral
            .into_iter()
            .all(|i| i.valid() && i.width() <= budget)
        {
            for i in 0..4 {
                result[i] = result[i] + integral[i];
            }
        } else {
            let middle = lo + (hi - lo) * 0.5;
            if middle <= lo || middle >= hi {
                return Err(Error::Unresolved);
            }
            pending.extend([(lo, middle), (middle, hi)]);
        }
    }
    Ok(result)
}

impl Drive {
    /// Evaluate actual nominal circular walls and the conical floor. Returned
    /// intervals bound numerical integration, not all geometry in a positional
    /// source-error tube. They do not certify exact support endpoint incidence.
    /// The source reserve/provenance remains separate.
    pub fn nominal_surface_evaluations(
        &self,
        integral_width_mm2: f64,
        work_limit: usize,
    ) -> Result<NominalSurfaceEvaluations, Error> {
        self.nominal_surface_evaluations_with_cancel(integral_width_mm2, work_limit, &|| false)
    }
    pub fn nominal_surface_evaluations_with_cancel(
        &self,
        width: f64,
        work_limit: usize,
        cancelled: &dyn Fn() -> bool,
    ) -> Result<NominalSurfaceEvaluations, Error> {
        if !width.is_finite() || width <= 0. {
            return Err(Error::Invalid);
        }
        if work_limit == 0 || work_limit > 1_000_000 {
            return Err(Error::Limit);
        }
        let scale = self.radius.max(self.cone.tip_depth_mm);
        let mut remaining = work_limit;
        let mut walls = [NominalSurfaceEvaluation {
            area_mm2: [0.; 2],
            center_mm: [[0.; 2]; 3],
        }; 12];
        let mut mouth = [I::point(0.); 4];
        for (index, &arc) in self.arcs.iter().enumerate() {
            let coordinate = |angle: I| {
                let cos = Jet::cos(angle);
                let sin = Jet::sin(angle);
                let x = Jet::constant(I::point(arc.center_mm[0])) + cos.scale(arc.radius_mm);
                let y = Jet::constant(I::point(arc.center_mm[1])) + sin.scale(arc.radius_mm);
                (x, y, cos, sin)
            };
            let planar = integrate(arc, width / 24., &mut remaining, cancelled, &|angle| {
                let (x, y, cos, sin) = coordinate(angle);
                [
                    (x * cos + y * sin).scale_interval(I::point(arc.radius_mm) * I::point(0.5)),
                    (x * x * cos)
                        .scale_interval(I::point(arc.radius_mm) * I::point(0.5) / I::point(scale)),
                    (y * y * sin)
                        .scale_interval(I::point(arc.radius_mm) * I::point(0.5) / I::point(scale)),
                    arc.moment(angle, self.radial_min_squared[index])
                        .scale_interval(I::point(1.) / I::point(scale)),
                ]
            })?;
            for i in 0..4 {
                mouth[i] = mouth[i] + planar[i] * I::point(arc.sweep_rad.signum());
            }
            let wall = integrate(arc, width / 2., &mut remaining, cancelled, &|angle| {
                let (x, y, _, _) = coordinate(angle);
                let floor = Jet::constant(I::point(self.cone.tip_depth_mm))
                    + arc
                        .radial(angle, self.radial_min_squared[index])
                        .scale(-self.cone.depth_per_radius);
                [
                    floor.scale(arc.radius_mm),
                    (x * floor).scale_interval(I::point(arc.radius_mm) / I::point(scale)),
                    (y * floor).scale_interval(I::point(arc.radius_mm) / I::point(scale)),
                    (floor * floor)
                        .scale_interval(I::point(arc.radius_mm) / (I::point(2.) * I::point(scale))),
                ]
            })?;
            walls[index] = evaluation(wall[0], [wall[1], wall[2], wall[3]], scale)?;
        }
        let area = mouth[0];
        let opening = evaluation(area, [mouth[1], mouth[2], I::point(0.)], scale)?;
        let slope = (I::point(1.) + I::point(self.cone.depth_per_radius).square()).sqrt();
        let depth_moment = I::point(self.cone.tip_depth_mm) / I::point(scale) * area
            - I::point(self.cone.depth_per_radius) * mouth[3];
        let floor_area = area * slope;
        let floor_moments = [mouth[1] * slope, mouth[2] * slope, depth_moment * slope];
        let mut floor = evaluation(floor_area, floor_moments, scale)?;
        // Constant cone area multiplier cancels from the XY centroid exactly.
        floor.center_mm[..2].copy_from_slice(&opening.center_mm[..2]);
        Ok(NominalSurfaceEvaluations {
            mouth: opening,
            walls,
            floor,
            source_error_mm: self.error,
        })
    }
}
