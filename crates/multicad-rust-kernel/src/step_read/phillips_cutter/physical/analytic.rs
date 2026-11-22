//! Analytical carrier/numeric coupling: complete circle/ellipse supports and
//! the finite original-line parameter interval. Finite arc lifts, oriented face
//! filling, common intersection vertices, embedding and mass remain separate.
use super::*;
use supports::chart::Chart;
#[cfg(test)]
mod tests;

#[derive(Clone, Copy, Debug)]
pub struct AnalyticSupportOptions {
    pub maximum_distance_mm: f64,
    pub work_limit: usize,
}
#[derive(Debug)]
pub struct AnalyticSupportBounds<'a> {
    boundary: &'a Boundary,
    edges: Vec<EdgeSupportBound>,
    work_used: usize,
}
impl<'a> AnalyticSupportBounds<'a> {
    pub fn boundary(&self) -> &'a Boundary {
        self.boundary
    }
    pub fn edges(&self) -> &[EdgeSupportBound] {
        &self.edges
    }
    /// New work only; the prerequisite endpoint receipt is charged separately.
    pub fn work_used(&self) -> usize {
        self.work_used
    }
    pub fn retained_capacity_bytes(&self) -> usize {
        std::mem::size_of::<Self>().saturating_add(
            self.edges
                .capacity()
                .saturating_mul(std::mem::size_of::<EdgeSupportBound>()),
        )
    }
}
fn zero(coefficients: &[I]) -> Result<()> {
    require(
        coefficients
            .iter()
            .all(|i| i.valid() && i.lo <= 0. && i.hi >= 0.),
    )
}
fn hull(values: &[I]) -> I {
    if values.iter().any(|i| !i.valid()) {
        return I::new(f64::NAN, f64::NAN);
    }
    I::new(
        values.iter().map(|i| i.lo).fold(f64::INFINITY, f64::min),
        values
            .iter()
            .map(|i| i.hi)
            .fold(f64::NEG_INFINITY, f64::max),
    )
}
fn rotated([x, y]: [I; 2], quarter: u8) -> [I; 2] {
    match quarter {
        0 => [x, y],
        1 => [-y, x],
        2 => [-x, -y],
        _ => [y, -x],
    }
}
fn arc(axis: Axis, radii: [f64; 2], chart: &Chart, work: &mut Work<'_>) -> Result<f64> {
    let frame = support::Frame::literal(&LiteralAxis {
        origin: axis.origin,
        normal: axis.normal,
        reference: axis.reference,
    })
    .map_err(source_error)?;
    let radii = radii.map(support::adjacent);
    let one = I::point(1.);
    let zero_i = I::point(0.);
    let c = I::point(0.5).sqrt();
    let weights = [one, c, one];
    let w = hull(&weights);
    // Cancel translation before multiplying by the irrational weight. The
    // actual interval transport error remains; an origin*weight dependency
    // explosion is avoided at large world coordinates.
    let center = chart.project(frame.origin);
    let axes = [frame.axes[0], frame.axes[1]].map(|axis| chart.project_vector(axis));
    let mut upper = 0_f64;
    for quarter in 0..4 {
        work.step(512)?;
        let controls = [[one, zero_i], [c, c], [zero_i, one]].map(|p| rotated(p, quarter));
        let local = std::array::from_fn(|j| {
            std::array::from_fn(|i| {
                center[i] * weights[j]
                    + axes[0][i] * radii[0] * controls[j][0]
                    + axes[1][i] * radii[1] * controls[j][1]
            })
        });
        let coefficients = chart.implicit_homogeneous(local, weights);
        zero(&coefficients)?;
        let denominator = (0..chart.homogeneous_degree()).fold(one, |p, _| p * w);
        let residual = hull(&coefficients) / denominator;
        // Positive rational Bernstein weights put the whole curve inside the
        // convex hull of its Euclidean controls N_i/W_i. Dividing separate
        // numerator/weight hulls would destroy constant-plane correlation.
        let local_box = std::array::from_fn(|i| {
            hull(&std::array::from_fn::<_, 3, _>(|j| {
                local[j][i] / weights[j]
            }))
        });
        upper = upper.max(chart.uniform_distance(residual, local_box)?);
    }
    Ok(upper)
}
fn line(
    origin: Point,
    direction: Point,
    magnitude: f64,
    points: [Point; 2],
    chart: &Chart,
    work: &mut Work<'_>,
    cone: bool,
) -> Result<f64> {
    work.step(512)?;
    let seed = (0..3)
        .min_by(|&a, &b| direction[a].abs().total_cmp(&direction[b].abs()))
        .ok_or(ReadError::Invalid)?;
    let frame = support::Frame::literal(&LiteralAxis {
        origin,
        normal: direction,
        reference: std::array::from_fn(|i| if i == seed { 1. } else { 0. }),
    })
    .map_err(source_error)?;
    let magnitude = support::adjacent(magnitude);
    let vector = frame.axes[2].map(|v| v * magnitude);
    let length = vector.into_iter().fold(I::point(0.), |s, v| s + v.square());
    if !length.valid() || length.lo <= 0. {
        return Err(ReadError::Unsupported);
    }
    let local = points.map(|point| {
        let delta =
            std::array::from_fn::<_, 3, _>(|i| support::adjacent(point[i]) - frame.origin[i]);
        let parameter = delta
            .into_iter()
            .zip(vector)
            .fold(I::point(0.), |s, (a, b)| s + a * b)
            / length;
        chart.project(std::array::from_fn(|i| {
            frame.origin[i] + vector[i] * parameter
        }))
    });
    let [a, b] = local;
    let controls = [
        a,
        std::array::from_fn(|i| (I::point(2.) * a[i] + b[i]) / I::point(3.)),
        std::array::from_fn(|i| (a[i] + I::point(2.) * b[i]) / I::point(3.)),
        b,
    ];
    let coefficients = chart.implicit_curve(controls);
    zero(&coefficients)?;
    if cone {
        chart.cone_line_distance(local)
    } else {
        chart.uniform_distance(
            hull(&coefficients),
            std::array::from_fn(|i| hull(&[a[i], b[i]])),
        )
    }
}
#[derive(Clone, Copy)]
enum Primitive {
    Line(Point, Point, f64),
    Arc(Axis, [f64; 2]),
}
impl Boundary {
    pub fn analytic_support_bounds(
        &self,
        endpoints: &EndpointCouplingBounds<'_>,
        options: AnalyticSupportOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> std::result::Result<AnalyticSupportBounds<'_>, PhillipsReadError> {
        if !std::ptr::eq(self, endpoints.boundary())
            || !options.maximum_distance_mm.is_finite()
            || options.maximum_distance_mm <= 0.
        {
            return Err(PhillipsReadError::Invalid);
        }
        let mut work = Work {
            remaining: options.work_limit,
            interrupt,
            cancelled: Cell::new(false),
        };
        let result: Result<AnalyticSupportBounds<'_>> = (|| {
            let mut edges = Vec::new();
            for (edge_index, edge) in self.edges().iter().enumerate() {
                let primitive = match edge.curve {
                    Curve::CubicIntersection { .. } => continue,
                    Curve::Line {
                        origin,
                        direction,
                        magnitude,
                        ..
                    } => Primitive::Line(origin, direction, magnitude),
                    Curve::Circle { axis, radius_mm } => Primitive::Arc(axis, [radius_mm; 2]),
                    Curve::Ellipse { axis, semiaxes_mm } => Primitive::Arc(axis, semiaxes_mm),
                };
                let mut owner_distances_mm = [0.; 2];
                for (side, owner) in edge.owners.into_iter().enumerate() {
                    let surface = &self.faces()[owner].support;
                    let chart = Chart::new(surface)?;
                    owner_distances_mm[side] = match primitive {
                        Primitive::Line(origin, direction, magnitude) => line(
                            origin,
                            direction,
                            magnitude,
                            edge.vertices.map(|v| self.vertices()[v].point),
                            &chart,
                            &mut work,
                            matches!(surface, Surface::Cone { .. }),
                        )?,
                        Primitive::Arc(axis, radii) => arc(axis, radii, &chart, &mut work)?,
                    };
                    require(owner_distances_mm[side] <= options.maximum_distance_mm)?;
                }
                edges.push(EdgeSupportBound {
                    edge_index,
                    owner_distances_mm,
                });
            }
            Ok(AnalyticSupportBounds {
                boundary: self,
                edges,
                work_used: options.work_limit - work.remaining,
            })
        })();
        if work.cancelled.get() {
            Err(PhillipsReadError::Cancelled)
        } else {
            result.map_err(PhillipsReadError::from)
        }
    }
}
