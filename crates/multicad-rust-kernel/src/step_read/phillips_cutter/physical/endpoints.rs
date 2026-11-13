//! Literal vertex/curve endpoint coupling under independently derived numeric
//! enclosures. A caller's requested distance never expands a residual to zero.
use super::*;
#[cfg(test)]
mod tests;

#[derive(Clone, Copy, Debug)]
pub struct EndpointCouplingOptions {
    pub maximum_distance_mm: f64,
    pub work_limit: usize,
}
#[derive(Debug)]
pub struct EndpointCouplingBounds<'a> {
    boundary: &'a Boundary,
    edge_distances_mm: Vec<[f64; 2]>,
    work_used: usize,
}
impl<'a> EndpointCouplingBounds<'a> {
    pub fn boundary(&self) -> &'a Boundary {
        self.boundary
    }
    pub fn edge_distances_mm(&self) -> &[[f64; 2]] {
        &self.edge_distances_mm
    }
    pub fn work_used(&self) -> usize {
        self.work_used
    }
    pub fn retained_capacity_bytes(&self) -> usize {
        std::mem::size_of::<Self>().saturating_add(
            self.edge_distances_mm
                .capacity()
                .saturating_mul(std::mem::size_of::<[f64; 2]>()),
        )
    }
}
fn zero(i: I) -> Result<()> {
    require(i.valid() && i.lo <= 0. && i.hi >= 0.)
}
fn frame(axis: Axis) -> Result<support::Frame> {
    support::Frame::literal(&LiteralAxis {
        origin: axis.origin,
        normal: axis.normal,
        reference: axis.reference,
    })
    .map_err(source_error)
}
fn local(frame: support::Frame, point: Point) -> [I; 3] {
    let delta = std::array::from_fn::<_, 3, _>(|i| support::adjacent(point[i]) - frame.origin[i]);
    frame.axes.map(|axis| {
        axis.into_iter()
            .zip(delta)
            .fold(I::point(0.), |s, (a, b)| s + a * b)
    })
}
fn analytic_endpoint(curve: &Curve, point: Point) -> Result<f64> {
    let error = match curve {
        Curve::Line {
            origin, direction, ..
        } => {
            let delta = std::array::from_fn::<_, 3, _>(|i| {
                support::adjacent(point[i]) - support::adjacent(origin[i])
            });
            let direction = direction.map(support::adjacent);
            let cross = support::cross(delta, direction);
            for c in cross {
                zero(c)?;
            }
            let length = support::norm(direction);
            if length.lo <= 0. {
                return Err(ReadError::Unsupported);
            }
            support::norm(cross) / length
        }
        Curve::Circle { axis, radius_mm } => {
            let [x, y, z] = local(frame(*axis)?, point);
            let r = support::adjacent(*radius_mm);
            let square = x.square() + y.square();
            zero(z)?;
            zero(square - r.square())?;
            ((square.sqrt() - r).square() + z.square()).sqrt()
        }
        Curve::Ellipse {
            axis,
            semiaxes_mm: [a, b],
        } => {
            let [x, y, z] = local(frame(*axis)?, point);
            let (a, b) = (support::adjacent(*a), support::adjacent(*b));
            if a.lo <= 0. || b.lo <= 0. {
                return Err(ReadError::Unsupported);
            }
            let square = (x / a).square() + (y / b).square();
            zero(z)?;
            zero(square - I::point(1.))?;
            let radial = I::point(a.hi.max(b.hi)) * (square.sqrt() - I::point(1.));
            (radial.square() + z.square()).sqrt()
        }
        Curve::CubicIntersection { .. } => return Err(ReadError::Unsupported),
    };
    if !error.valid() {
        return Err(ReadError::Limit);
    }
    Ok(error.hi)
}
impl Boundary {
    pub fn endpoint_coupling_bounds(
        &self,
        options: EndpointCouplingOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> std::result::Result<EndpointCouplingBounds<'_>, PhillipsReadError> {
        if !options.maximum_distance_mm.is_finite() || options.maximum_distance_mm <= 0. {
            return Err(PhillipsReadError::Invalid);
        }
        let mut work = Work {
            remaining: options.work_limit,
            interrupt,
            cancelled: Cell::new(false),
        };
        let result: Result<EndpointCouplingBounds<'_>> = (|| {
            let mut edge_distances_mm = Vec::new();
            for edge in self.edges() {
                work.step(8)?;
                let mut distance = [0.; 2];
                // This zero reserve belongs only to the sealed arena's exact
                // shared point references, never to an arbitrary cubic query.
                // Parsing also proves the reversed EDGE_CURVE case.
                if !matches!(edge.curve, Curve::CubicIntersection { .. }) {
                    for (i, vertex) in edge.vertices.into_iter().enumerate() {
                        distance[i] =
                            analytic_endpoint(&edge.curve, self.vertices()[vertex].point)?;
                        require(distance[i] <= options.maximum_distance_mm)?;
                    }
                }
                edge_distances_mm.push(distance);
            }
            Ok(EndpointCouplingBounds {
                boundary: self,
                edge_distances_mm,
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
