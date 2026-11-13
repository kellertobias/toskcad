//! Finite declared-curve parameter bounds, not common true vertices or face filling.
use super::*;
#[derive(Clone, Copy)]
pub struct LiftOptions {
    pub maximum_work: usize,
    pub maximum_bytes: usize,
}
impl Default for LiftOptions {
    fn default() -> Self {
        Self {
            maximum_work: WORK,
            maximum_bytes: BYTES,
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum LiftKind {
    /// Parameters at EDGE_CURVE start/end on its original VECTOR, not length.
    Line { parameters: [[f64; 2]; 2] },
    /// A connected start-angle chart and one directed sweep on the actual
    /// circle/ellipse. Both are independent numeric enclosures, never snaps.
    Arc {
        start_angle_rad: [f64; 2],
        signed_sweep_rad: [f64; 2],
    },
}
#[derive(Clone, Copy, Debug)]
pub struct FiniteLift {
    edge: usize,
    kind: LiftKind,
}
impl FiniteLift {
    pub fn edge_index(&self) -> usize {
        self.edge
    }
    pub fn kind(&self) -> LiftKind {
        self.kind
    }
}
/// Retains original numeric endpoint/support receipts. This parameter receipt
/// does not identify independently rounded intersection endpoints as equal.
pub struct FiniteLifts<'a> {
    analytic: AnalyticSupportBounds<'a>,
    endpoints: EndpointCouplingBounds<'a>,
    lifts: Vec<FiniteLift>,
    work_used: usize,
}
fn pi() -> I {
    I::new(
        std::f64::consts::PI.next_down(),
        std::f64::consts::PI.next_up(),
    )
}
fn tau() -> I {
    I::point(2.) * pi()
}
fn bounds(x: I) -> [f64; 2] {
    [x.lo, x.hi]
}
fn local(frame: support::Frame, p: Point) -> [I; 3] {
    let delta = std::array::from_fn::<_, 3, _>(|k| support::adjacent(p[k]) - frame.origin[k]);
    frame
        .axes
        .map(|axis| (0..3).fold(I::point(0.), |s, k| s + axis[k] * delta[k]))
}
fn frame(axis: Axis, budget: &mut Budget<'_>) -> BranchResult<support::Frame> {
    budget.charge(1)?;
    support::Frame::literal(&LiteralAxis {
        origin: axis.origin,
        normal: axis.normal,
        reference: axis.reference,
    })
    .map_err(source_error)
    .map_err(PhillipsReadError::from)
}
fn scaled(p: [I; 2]) -> BranchResult<[I; 2]> {
    if p.iter().any(|x| !x.valid()) {
        return Err(PhillipsReadError::Unsupported);
    }
    let scale = p.into_iter().map(I::magnitude).fold(0., f64::max);
    if scale == 0. {
        return Err(PhillipsReadError::Unsupported);
    }
    // Positive homogeneous rescaling changes no angle or curve radius; direct
    // endpoint quotients avoid overflowing the reciprocal of a tiny scale.
    Ok(p.map(|x| I::new((x.lo / scale).next_down(), (x.hi / scale).next_up())))
}
fn angle(y: I, x: I) -> BranchResult<I> {
    if !x.valid() || !y.valid() || !(x.lo > 0. || x.hi < 0. || y.lo > 0. || y.hi < 0.) {
        return Err(PhillipsReadError::Unsupported);
    }
    // Each chosen rectangle chart is continuous and avoids the origin. Angle
    // extrema lie at corners; this also avoids any unbounded ratio/0*infinity.
    let negative_x = x.hi < 0.;
    let mut out = I::new(f64::INFINITY, f64::NEG_INFINITY);
    for a in [y.lo, y.hi] {
        for b in [x.lo, x.hi] {
            let raw = a.atan2(b);
            let mut v = I::point(raw) + I::new(-16. * f64::EPSILON, 16. * f64::EPSILON);
            if negative_x && raw < 0. {
                v = v + tau();
            }
            out.lo = out.lo.min(v.lo);
            out.hi = out.hi.max(v.hi);
        }
    }
    Ok(out)
}
fn directed_sweep(first: [I; 2], last: [I; 2], forward: bool) -> BranchResult<I> {
    let cross = first[0] * last[1] - first[1] * last[0];
    let dot = first[0] * last[0] + first[1] * last[1];
    let relative = angle(cross, dot)?;
    let sweep = if forward {
        if relative.lo > 0. {
            relative
        } else if relative.hi < 0. {
            relative + tau()
        } else {
            return Err(PhillipsReadError::Unsupported);
        }
    } else if relative.hi < 0. {
        relative
    } else if relative.lo > 0. {
        relative - tau()
    } else {
        return Err(PhillipsReadError::Unsupported);
    };
    // Distinct uncertain near-coincident endpoints cannot select between a
    // positive tiny arc and a nearly complete turn by numeric proximity.
    if (forward && !(sweep.lo > 0. && sweep.hi < tau().lo))
        || (!forward && !(sweep.hi < 0. && sweep.lo > -tau().lo))
    {
        return Err(PhillipsReadError::Unsupported);
    }
    Ok(sweep)
}
fn arc(
    axis: Axis,
    radii: Option<[f64; 2]>,
    points: [Point; 2],
    forward: bool,
    budget: &mut Budget<'_>,
) -> BranchResult<LiftKind> {
    let frame = frame(axis, budget)?;
    let directions = points.map(|p| {
        let p = local(frame, p);
        if let Some(r) = radii {
            [
                p[0] / support::adjacent(r[0]),
                p[1] / support::adjacent(r[1]),
            ]
        } else {
            [p[0], p[1]]
        }
    });
    let first = scaled(directions[0])?;
    let last = scaled(directions[1])?;
    let start = angle(first[1], first[0])?;
    let sweep = directed_sweep(first, last, forward)?;
    Ok(LiftKind::Arc {
        start_angle_rad: bounds(start),
        signed_sweep_rad: bounds(sweep),
    })
}
fn line(
    origin: Point,
    direction: Point,
    magnitude: f64,
    points: [Point; 2],
    forward: bool,
    budget: &mut Budget<'_>,
) -> BranchResult<LiftKind> {
    let seed = (0..3)
        .min_by(|&a, &b| direction[a].abs().total_cmp(&direction[b].abs()))
        .unwrap();
    let frame = frame(
        Axis {
            origin,
            normal: direction,
            reference: std::array::from_fn(|k| if k == seed { 1. } else { 0. }),
        },
        budget,
    )?;
    let vector = frame.axes[2].map(|x| x * support::adjacent(magnitude));
    let length = vector.into_iter().fold(I::point(0.), |s, x| s + x.square());
    // The SAME immutable earned analytic receipt already proves this original
    // VECTOR's positive finite squared length and numeric endpoint coupling.
    let parameters = points.map(|p| {
        let delta = std::array::from_fn::<_, 3, _>(|k| support::adjacent(p[k]) - frame.origin[k]);
        (0..3).fold(I::point(0.), |s, k| s + delta[k] * vector[k]) / length
    });
    let d = parameters[1] - parameters[0];
    if !d.valid()
        || parameters.iter().any(|p| !p.valid())
        || (forward && d.lo <= 0.)
        || (!forward && d.hi >= 0.)
    {
        return Err(PhillipsReadError::Unsupported);
    }
    Ok(LiftKind::Line {
        parameters: parameters.map(bounds),
    })
}
fn footprint(
    base: usize,
    capacity: usize,
    limit: usize,
    budget: &mut Budget<'_>,
) -> BranchResult<()> {
    budget.charge(1)?;
    if add(
        base,
        capacity.saturating_mul(std::mem::size_of::<FiniteLift>()),
    ) > limit
    {
        Err(PhillipsReadError::Limit)
    } else {
        Ok(())
    }
}
impl<'a> AnalyticSupportBounds<'a> {
    pub fn finite_edge_lifts(
        self,
        endpoints: EndpointCouplingBounds<'a>,
        options: LiftOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> BranchResult<FiniteLifts<'a>> {
        if options.maximum_work == 0
            || options.maximum_work > WORK
            || options.maximum_bytes == 0
            || options.maximum_bytes > BYTES
            || !std::ptr::eq(self.boundary(), endpoints.boundary())
        {
            return Err(PhillipsReadError::Invalid);
        }
        if interrupt.is_some_and(|f| f()) {
            return Err(PhillipsReadError::Cancelled);
        }
        let inherited = add(self.work_used(), endpoints.work_used());
        let mut budget = Budget {
            remaining: options
                .maximum_work
                .checked_sub(inherited)
                .ok_or(PhillipsReadError::Limit)?,
            interrupt,
        };
        let boundary = self.boundary();
        let base = add(
            add(
                add(arena(boundary), self.retained_capacity_bytes()),
                endpoints.retained_capacity_bytes(),
            ),
            add(std::mem::size_of::<FiniteLifts<'_>>(), 8192),
        );
        footprint(base, self.edges().len(), options.maximum_bytes, &mut budget)?;
        let mut lifts = Vec::new();
        lifts
            .try_reserve_exact(self.edges().len())
            .map_err(|_| PhillipsReadError::Limit)?;
        footprint(base, lifts.capacity(), options.maximum_bytes, &mut budget)?;
        for (edge_index, edge) in boundary.edges().iter().enumerate() {
            let points = edge.vertices.map(|i| boundary.vertices()[i].point);
            let kind = match edge.curve {
                Curve::Line {
                    origin,
                    direction,
                    magnitude,
                    ..
                } => {
                    budget.charge(1024)?;
                    line(
                        origin,
                        direction,
                        magnitude,
                        points,
                        edge.same_sense,
                        &mut budget,
                    )?
                }
                Curve::Circle { axis, .. } => {
                    budget.charge(1024)?;
                    arc(axis, None, points, edge.same_sense, &mut budget)?
                }
                Curve::Ellipse { axis, semiaxes_mm } => {
                    budget.charge(1024)?;
                    arc(
                        axis,
                        Some(semiaxes_mm),
                        points,
                        edge.same_sense,
                        &mut budget,
                    )?
                }
                Curve::CubicIntersection { .. } => continue,
            };
            lifts.push(FiniteLift {
                edge: edge_index,
                kind,
            });
        }
        budget.charge(1)?;
        Ok(FiniteLifts {
            analytic: self,
            endpoints,
            lifts,
            work_used: options.maximum_work - budget.remaining,
        })
    }
}
impl<'a> FiniteLifts<'a> {
    pub fn boundary(&self) -> &'a Boundary {
        self.analytic.boundary()
    }
    pub fn edges(&self) -> &[FiniteLift] {
        &self.lifts
    }
    pub fn endpoint_bounds(&self) -> &EndpointCouplingBounds<'a> {
        &self.endpoints
    }
    pub fn work_used(&self) -> usize {
        self.work_used
    }
    pub fn retained_capacity_bytes(&self) -> usize {
        add(
            std::mem::size_of::<Self>(),
            add(
                allocation(&self.lifts),
                add(
                    self.analytic
                        .retained_capacity_bytes()
                        .saturating_sub(std::mem::size_of::<AnalyticSupportBounds<'_>>()),
                    self.endpoints
                        .retained_capacity_bytes()
                        .saturating_sub(std::mem::size_of::<EndpointCouplingBounds<'_>>()),
                ),
            ),
        )
    }
}
#[cfg(test)]
mod tests;
