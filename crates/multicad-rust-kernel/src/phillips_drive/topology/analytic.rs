//! Actual finite analytic edge supports; intersection trims remain separate.
use super::*;
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum AnalyticCurve {
    Line(LineSupport),
    /// The circle is centered on the cutter axis in the horizontal plane.
    CircleArc {
        radius_mm_bounds: Bounds,
        depth_mm_bounds: Bounds,
        start_angle_rad_bounds: Bounds,
        sweep_rad_bounds: Bounds,
    },
    /// The ellipse lies in a horizontal plane, in its exact quarter-turn frame.
    EllipseArc {
        center_mm_bounds: [Bounds; 3],
        semiaxes_mm_bounds: [Bounds; 2],
        quarter_turn: u8,
        start_angle_rad_bounds: Bounds,
        sweep_rad_bounds: Bounds,
    },
}
pub(super) fn meridian_data(
    topology: &BoundaryTopology<'_>,
    station: MeridianStation,
    quarter: u8,
) -> [I; 4] {
    let p = topology
        .events
        .event(station.event())
        .point_mm_bounds()
        .map(|[a, b]| I::new(a, b));
    let phi = (p[1] / p[0]).atan();
    let pi = I::new(PI.next_down(), PI.next_up());
    let source = topology.events.source();
    let radius = match station {
        MeridianStation::Entry => {
            source.envelope.entry - I::point(source.spec.mouth_radius_mm) * source.envelope.cos_a
        }
        MeridianStation::Floor => I::point(source.spec.g_diameter_mm) / I::point(2.),
    };
    [
        radius,
        p[2],
        phi + I::point(quarter as f64) * pi / I::point(2.),
        pi / I::point(2.) - I::point(2.) * phi,
    ]
}
impl BoundaryTopology<'_> {
    /// Return the authored analytic support and finite angular domain. None
    /// denotes a genuine dual-support intersection, not a fitted analytic edge.
    /// These source bounds do not certify independently transmitted file data.
    pub fn curve_analytic_support(&self, index: usize) -> Result<Option<AnalyticCurve>, Error> {
        let key = self.edges.get(index).ok_or(Error::Invalid)?.key;
        let curve = match key {
            CurveKey::Trim { .. } => return Ok(None),
            CurveKey::Mouth { edge } => match self.events.mouth().edges()[edge].curve() {
                SectionCurve::Line { .. } => {
                    return self
                        .curve_line_support(index)
                        .map(|line| line.map(AnalyticCurve::Line));
                }
                SectionCurve::CircleArc {
                    radius_mm_bounds,
                    start_angle_rad_bounds,
                    sweep_rad_bounds,
                } => AnalyticCurve::CircleArc {
                    radius_mm_bounds,
                    depth_mm_bounds: bounds(I::point(0.)),
                    start_angle_rad_bounds,
                    sweep_rad_bounds,
                },
                SectionCurve::EllipseArc {
                    center_mm_bounds,
                    semiaxes_mm_bounds,
                    quarter_turn,
                    start_angle_rad_bounds,
                    sweep_rad_bounds,
                } => AnalyticCurve::EllipseArc {
                    center_mm_bounds: [
                        center_mm_bounds[0],
                        center_mm_bounds[1],
                        bounds(I::point(0.)),
                    ],
                    semiaxes_mm_bounds,
                    quarter_turn,
                    start_angle_rad_bounds,
                    sweep_rad_bounds,
                },
            },
            CurveKey::MeridianCircle { station, quarter } => {
                let [radius, depth, start, sweep] = meridian_data(self, station, quarter);
                AnalyticCurve::CircleArc {
                    radius_mm_bounds: bounds(radius),
                    depth_mm_bounds: bounds(depth),
                    start_angle_rad_bounds: bounds(start),
                    sweep_rad_bounds: bounds(sweep),
                }
            }
            _ => {
                return self
                    .curve_line_support(index)
                    .map(|line| line.map(AnalyticCurve::Line));
            }
        };
        Ok(Some(curve))
    }
}
#[cfg(test)]
mod tests;
