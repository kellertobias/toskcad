//! Actual finite cap curves and source-plane fan correspondence.
use super::*;
use properties::source::{Circle, Source};
use spline::Spline;
enum Segment {
    Spline(Spline),
    Circle(Circle),
}
pub(super) struct Cap {
    pub face: usize,
    pub normal: [I; 3],
    pub center: [I; 3],
    sides: Vec<Vec<Segment>>,
}
impl Cap {
    pub fn new(body: &RadialBoundary, end: usize, work: &mut Work<'_>) -> Result<Self, Error> {
        let caps = &body.ordered.caps;
        let chart = &caps.axial.chart;
        let face = caps.axial.caps[end];
        let axis = caps.axial.cap_axes[end];
        let normal = checked_axis(axis)?[2].map(|v| {
            if body.input().faces[face].forward {
                v
            } else {
                -v
            }
        });
        let delta =
            std::array::from_fn(|k| I::point(axis.origin[k]) - I::point(chart.axis_origin[k]));
        let station = dot(delta, chart.axes[2]);
        let center =
            std::array::from_fn(|k| I::point(chart.axis_origin[k]) + chart.axes[2][k] * station);
        let mut sides = Vec::with_capacity(caps.sides.len());
        for group in &caps.rims[end] {
            let mut segments = Vec::with_capacity(group.len());
            for &(index, direction) in group {
                work.step()?;
                let edge = body.input().edges[index];
                let segment = match &body.input().curves[edge.curve] {
                    Curve::Spline { poles, .. } => {
                        let mut controls = Vec::with_capacity(poles.len());
                        for chunk in poles.chunks(128) {
                            work.charge(chunk.len())?;
                            controls.extend(chunk.iter().map(|p| p.map(I::point)));
                        }
                        if edge.forward != direction {
                            controls.reverse();
                        }
                        Segment::Spline(Spline::new(controls))
                    }
                    Curve::Circle { axis, radius_mm } => {
                        work.charge(64)?;
                        let basis = checked_axis(*axis)?;
                        let (a, b) = caps::ends(body.input(), (index, direction));
                        let ray = |index: usize| {
                            let p = std::array::from_fn(|k| {
                                I::point(body.input().vertices[index][k]) - I::point(axis.origin[k])
                            });
                            [dot(p, basis[0]), dot(p, basis[1]), I::point(0.)]
                        };
                        let (a, b) = (ray(a), ray(b));
                        let sweep = if edge.forward == direction {
                            caps::phase(a, b)
                        } else {
                            caps::phase(b, a).map(|v| -v)
                        }
                        .ok_or(Error::Unresolved)?;
                        Segment::Circle(Circle {
                            center: axis.origin.map(I::point),
                            u: basis[0],
                            v: basis[1],
                            radius: I::point(*radius_mm),
                            start: properties::faces::argument([a[0], a[1]])?,
                            sweep,
                            range: [0., 1.],
                        })
                    }
                };
                segments.push(segment);
            }
            sides.push(segments);
        }
        Ok(Self {
            face,
            normal,
            center,
            sides,
        })
    }
    pub fn point(
        &self,
        body: &RadialBoundary,
        end: usize,
        side: usize,
        u: f64,
        work: &mut Work<'_>,
    ) -> Result<[I; 3], Error> {
        let group = &body.ordered.caps.rims[end][side];
        let vertex = if u == 0. {
            Some(caps::ends(body.input(), group[0]).0)
        } else if u == 1. {
            Some(caps::ends(body.input(), group[group.len() - 1]).1)
        } else if u == 0.5 && group.len() == 2 {
            Some(caps::ends(body.input(), group[0]).1)
        } else {
            None
        };
        if let Some(vertex) = vertex {
            work.step()?;
            return Ok(body.input().vertices[vertex].map(I::point));
        }
        let (segment, t) = self.segment(side, I::point(u));
        match segment {
            Segment::Spline(spline) => Ok(spline.bounds([t.lo, t.hi], work)?[0]),
            Segment::Circle(circle) => {
                work.charge(64)?;
                Ok(Source::Circle(circle.clone())
                    .at([t, I::point(0.)])
                    .p
                    .map(|v| v.c[0]))
            }
        }
    }
    fn segment(&self, side: usize, u: I) -> (&Segment, I) {
        let segments = &self.sides[side];
        if segments.len() == 1 {
            (&segments[0], u)
        } else if u.lo >= 0.5 {
            (&segments[1], I::point(2.) * u - I::point(1.))
        } else {
            (&segments[0], I::point(2.) * u)
        }
    }
    pub fn curve_error(
        &self,
        side: usize,
        u: [f64; 2],
        vertices: [Point; 2],
        work: &mut Work<'_>,
    ) -> Result<f64, Error> {
        if self.sides[side].len() == 2 && u[0] < 0.5 && u[1] > 0.5 {
            return Err(Error::Invalid);
        }
        let (segment, a) = self.segment(side, I::point(u[0]));
        let (_, b) = self.segment(side, I::point(u[1]));
        match segment {
            Segment::Spline(spline) => spline.chord_error(u, vertices, work),
            Segment::Circle(circle) => {
                // The caller's grid includes the actual two-arc seam, so
                // both endpoints use this one retained circle segment.
                let b = if u[1] == 0.5 && u[0] < 0.5 {
                    I::point(1.)
                } else {
                    b
                };
                work.charge(128)?;
                let error = |t: I, p: Point| {
                    Source::Circle(circle.clone())
                        .at([t, I::point(0.)])
                        .p
                        .into_iter()
                        .zip(p)
                        .fold(I::point(0.), |s, (v, p)| {
                            s + I::point((v.c[0] - I::point(p)).magnitude())
                        })
                };
                let endpoint = I::point(error(a, vertices[0]).hi.max(error(b, vertices[1]).hi));
                let curvature =
                    circle.radius * circle.sweep.square() * (b - a).square() / I::point(8.);
                Ok((endpoint + curvature).checked()?.hi)
            }
        }
    }
}
