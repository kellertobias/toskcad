//! Shared-variable interval geometry. This is a family of possible planar
//! loops, not a midpoint curve and not an independently parsed file.
use super::*;
#[derive(Clone, Copy)]
pub struct IntervalCurve<'a> {
    pub curve_owner: usize,
    pub spans: usize,
    /// Every coefficient names ONE variable in IntervalInput::points.
    pub point_indices: &'a [usize],
}
#[derive(Clone, Copy)]
pub struct IntervalInput<'a> {
    /// [x interval,y interval] for each distinct defining point variable.
    pub points: &'a [[[f64; 2]; 2]],
    pub curves: &'a [IntervalCurve<'a>],
    pub counter_clockwise: bool,
}
pub struct IntervalCertificate<'a> {
    input: IntervalInput<'a>,
    geometry: Geometry,
    work: usize,
    input_bytes: usize,
}
impl<'a> IntervalCertificate<'a> {
    pub fn input(&self) -> IntervalInput<'a> {
        self.input
    }
    pub fn angular_advance(&self) -> [f64; 2] {
        self.geometry.angle
    }
    pub fn winding(&self) -> i8 {
        if self.input.counter_clockwise { 1 } else { -1 }
    }
    pub fn nodes(&self) -> usize {
        self.geometry.nodes
    }
    pub fn leaves(&self) -> usize {
        self.geometry.leaves
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn input_budget_bytes(&self) -> usize {
        self.input_bytes
    }
    pub fn scratch_budget_bytes(&self) -> usize {
        SCRATCH
    }
}
/// Certify EVERY assignment inside the supplied boxes with the declared shared
/// point variables. Joins require the same registry index, not overlapping
/// intervals, midpoint equality or a caller's unsupported endpoint assertion.
/// A file consumer must separately bind these indices to its actual numeric
/// leaves, knot structure, oriented owning graph and physical plane.
pub fn certify_intervals<'a>(
    input: IntervalInput<'a>,
    options: Options,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<IntervalCertificate<'a>, Error> {
    let mut budget = Budget {
        remaining: options.maximum_work,
        nodes: 0,
        interrupt,
    };
    budget.charge(1)?;
    if options.maximum_work > WORK
        || options.maximum_nodes > NODES
        || options.maximum_bytes > BYTES
        || options.maximum_depth > 32
        || input.curves.is_empty()
        || input.curves.len() > CURVES
        || input.points.is_empty()
    {
        return Err(Error::Invalid);
    }
    if input.points.len() > POINTS {
        return Err(Error::Limit);
    }
    let mut references = 0usize;
    for curve in input.curves {
        if curve.point_indices.len() > POINTS {
            return Err(Error::Limit);
        }
        references += curve.point_indices.len();
    }
    if references > POINTS {
        return Err(Error::Limit);
    }
    let input_bytes = std::mem::size_of_val(input.points)
        + references * std::mem::size_of::<usize>()
        + std::mem::size_of_val(input.curves);
    if input_bytes + SCRATCH > options.maximum_bytes {
        return Err(Error::Limit);
    }
    for point in input.points {
        budget.charge(8)?;
        if point.iter().any(|b| {
            !b[0].is_finite()
                || !b[1].is_finite()
                || b[0] > b[1]
                || b[0].abs() > 1e7
                || b[1].abs() > 1e7
        }) {
            return Err(Error::Invalid);
        }
    }
    budget.charge(input.curves.len())?;
    for (i, c) in input.curves.iter().enumerate() {
        budget.charge(1 + input.curves.len() + c.point_indices.len())?;
        if c.spans == 0
            || c.spans.checked_mul(11).and_then(|n| n.checked_add(1)) != Some(c.point_indices.len())
            || c.point_indices
                .iter()
                .any(|index| *index >= input.points.len())
        {
            return Err(Error::Invalid);
        }
        if input.curves[..i].iter().any(|previous| {
            previous.curve_owner == c.curve_owner
                || previous.point_indices.first() == c.point_indices.first()
        }) {
            return Err(Error::Invalid);
        }
        let next = &input.curves[(i + 1) % input.curves.len()];
        if c.point_indices.last() != next.point_indices.first() {
            return Err(Error::Invalid);
        }
    }
    let geometry = geometry(
        input.curves.len(),
        |i| input.curves[i].spans,
        |curve, span| {
            std::array::from_fn(|i| {
                let bounds = input.points[input.curves[curve].point_indices[11 * span + i]];
                [
                    I::new(bounds[0][0], bounds[0][1]),
                    I::new(bounds[1][0], bounds[1][1]),
                    I::point(0.),
                ]
            })
        },
        input.counter_clockwise,
        options,
        &mut budget,
    )?;
    Ok(IntervalCertificate {
        input,
        geometry,
        work: options.maximum_work - budget.remaining,
        input_bytes,
    })
}
#[cfg(test)]
mod tests;
