//! Actual retained curve display, with no mesh-derived edge reconstruction.
use super::*;
use properties::source::{Circle, Source};
#[cfg(test)]
mod tests;

enum CurveSource {
    Spline(Spline),
    Circle(Circle),
}
impl CurveSource {
    fn at(&self, t: f64, work: &mut Work<'_>) -> Result<Point, Error> {
        work.charge(64)?;
        let point = match self {
            Self::Spline(spline) => spline.bounds([t, t], work)?[0],
            Self::Circle(circle) => Source::Circle(circle.clone())
                .at([I::point(t), I::point(0.)])
                .p
                .map(|p| p.c[0]),
        };
        for value in point {
            value.checked()?;
        }
        Ok(point.map(|p| p.midpoint()))
    }
    fn error(
        &self,
        range: [f64; 2],
        points: [Point; 2],
        work: &mut Work<'_>,
    ) -> Result<f64, Error> {
        match self {
            Self::Spline(spline) => spline.chord_error(range, points, work),
            Self::Circle(circle) => {
                work.charge(128)?;
                let mut endpoint: f64 = 0.;
                for (t, p) in range.into_iter().zip(points) {
                    let exact = Source::Circle(circle.clone())
                        .at([I::point(t), I::point(0.)])
                        .p;
                    let error = (0..3).fold(I::point(0.), |sum, k| {
                        sum + I::point((exact[k].c[0] - I::point(p[k])).magnitude())
                    });
                    endpoint = endpoint.max(error.checked()?.hi);
                }
                let width = I::point(range[1]) - I::point(range[0]);
                let curvature =
                    circle.radius * circle.sweep.square() * width.square() / I::point(8.);
                Ok((I::point(endpoint) + curvature).checked()?.hi)
            }
        }
    }
}

pub(super) fn sample(
    body: &RadialBoundary,
    edge: usize,
    deflection: f64,
    point_limit: usize,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<Vec<Point>, Error> {
    if !deflection.is_finite() || deflection <= 0. || point_limit < 2 {
        return Err(Error::Invalid);
    }
    if point_limit > 1_000_000 {
        return Err(Error::Limit);
    }
    let input = body.input();
    let edge = input.edges.get(edge).ok_or(Error::Invalid)?;
    let mut work = Work {
        remaining: 4_000_000,
        interrupt,
    };
    work.step()?;
    let points = edge.vertices.map(|v| input.vertices[v]);
    let curve = match &input.curves[edge.curve] {
        Curve::Spline { poles, .. } => {
            let mut controls = Vec::with_capacity(poles.len());
            for chunk in poles.chunks(128) {
                work.charge(chunk.len())?;
                controls.extend(chunk.iter().map(|p| p.map(I::point)));
            }
            if !edge.forward {
                controls.reverse();
            }
            CurveSource::Spline(Spline::new(controls))
        }
        Curve::Circle { axis, radius_mm } => {
            work.charge(64)?;
            let basis = checked_axis(*axis)?;
            let ray = |p: Point| {
                let p = std::array::from_fn(|k| I::point(p[k]) - I::point(axis.origin[k]));
                [dot(p, basis[0]), dot(p, basis[1]), I::point(0.)]
            };
            let [a, b] = points.map(ray);
            let sweep = if edge.forward {
                caps::phase(a, b)
            } else {
                caps::phase(b, a).map(|v| -v)
            }
            .ok_or(Error::Unresolved)?;
            CurveSource::Circle(Circle {
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
    let mut output = vec![points[0]];
    let mut pending = vec![([0., 1.], points)];
    while let Some((range, points)) = pending.pop() {
        work.step()?;
        if curve.error(range, points, &mut work)? <= deflection {
            output.push(points[1]);
        } else {
            if output.len() + pending.len() + 2 > point_limit {
                return Err(Error::Limit);
            }
            let middle = split_parameter(range)?;
            let p = curve.at(middle, &mut work)?;
            pending.push(([middle, range[1]], [p, points[1]]));
            pending.push(([range[0], middle], [points[0], p]));
        }
    }
    work.step()?;
    Ok(output)
}

fn split_parameter(range: [f64; 2]) -> Result<f64, Error> {
    let middle = range[0] + (range[1] - range[0]) * 0.5;
    if middle <= range[0] || middle >= range[1] {
        return Err(Error::Unresolved);
    }
    Ok(middle)
}
