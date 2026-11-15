//! Standalone literal degree11 planar loops. Strict angular monotonicity,
//! exact numeric closure and one integer winding certify a radial Jordan loop.
//! No authoring recipe, file ownership, disk face or solid admission is inferred;
//! source decoders must bind this SAME input to their actual owning graph.
use crate::tapered_thread::interval::I;
mod interval;
pub use interval::{IntervalCertificate, IntervalCurve, IntervalInput, certify_intervals};
const POINTS: usize = 524_288;
const CURVES: usize = 64;
const NODES: usize = 65_536;
const WORK: usize = 64_000_000;
const BYTES: usize = 16 * 1024 * 1024;
const SCRATCH: usize = 128 * 1024;
type P = [f64; 3];
type B = [I; 3];
#[derive(Clone, Copy)]
pub struct Curve<'a> {
    pub curve_owner: usize,
    pub start_owner: usize,
    pub end_owner: usize,
    pub spans: usize,
    pub poles: &'a [P],
}
#[derive(Clone, Copy)]
pub struct Input<'a> {
    pub curves: &'a [Curve<'a>],
    pub counter_clockwise: bool,
}
#[derive(Clone, Copy)]
pub struct Options {
    pub maximum_work: usize,
    pub maximum_nodes: usize,
    pub maximum_bytes: usize,
    pub maximum_depth: u8,
}
impl Default for Options {
    fn default() -> Self {
        Self {
            maximum_work: WORK,
            maximum_nodes: NODES,
            maximum_bytes: BYTES,
            maximum_depth: 32,
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Error {
    Invalid,
    Limit,
    Unresolved,
    Cancelled,
}
pub struct Certificate<'a> {
    input: Input<'a>,
    plane_z: f64,
    angle: [f64; 2],
    work: usize,
    nodes: usize,
    leaves: usize,
    input_bytes: usize,
}
impl<'a> Certificate<'a> {
    pub fn input(&self) -> Input<'a> {
        self.input
    }
    pub fn plane_z(&self) -> f64 {
        self.plane_z
    }
    /// Positive angular advance in the requested orientation, in radians.
    pub fn angular_advance(&self) -> [f64; 2] {
        self.angle
    }
    pub fn winding(&self) -> i8 {
        if self.input.counter_clockwise { 1 } else { -1 }
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn nodes(&self) -> usize {
        self.nodes
    }
    pub fn leaves(&self) -> usize {
        self.leaves
    }
    pub fn input_budget_bytes(&self) -> usize {
        self.input_bytes
    }
    pub fn scratch_budget_bytes(&self) -> usize {
        SCRATCH
    }
}
struct Budget<'a> {
    remaining: usize,
    nodes: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl Budget<'_> {
    fn charge(&mut self, n: usize) -> Result<(), Error> {
        if self.interrupt.is_some_and(|f| f()) {
            return Err(Error::Cancelled);
        }
        self.remaining = self.remaining.checked_sub(n).ok_or(Error::Limit)?;
        Ok(())
    }
}
fn hull(points: impl Iterator<Item = B>) -> B {
    let mut h = [I::new(f64::INFINITY, f64::NEG_INFINITY); 3];
    for p in points {
        for a in 0..3 {
            h[a].lo = h[a].lo.min(p[a].lo);
            h[a].hi = h[a].hi.max(p[a].hi);
        }
    }
    h
}
pub(crate) fn split_half_controls(mut p: [B; 12]) -> [[B; 12]; 2] {
    let mut left = [[I::point(0.); 3]; 12];
    let mut right = left;
    for level in 0..12 {
        left[level] = p[0];
        right[11 - level] = p[11 - level];
        for i in 0..11 - level {
            p[i] = std::array::from_fn(|a| (p[i][a] + p[i + 1][a]) * I::point(0.5));
        }
    }
    [left, right]
}
// For0<r<=1/2, successive alternating terms decrease. The11-term partial
// sum through+r^21/21 is an upper bound; subtracting r^23/23 is a lower bound.
// Every arithmetic operation is directed. No libm/error-ULP assumption enters.
fn atan_small(r: I) -> I {
    let square = r * r;
    let mut power = r;
    let mut sum = r;
    for k in 1..=10 {
        power = power * square;
        let term = power / I::point((2 * k + 1) as f64);
        sum = if k % 2 == 0 { sum + term } else { sum - term };
    }
    let remainder = (power * square / I::point(23.)).hi;
    I::new((sum - I::point(remainder)).lo, sum.hi)
}
fn positive_ratio(numerator: I, denominator: I) -> I {
    I::new(
        (numerator.lo / denominator.hi).next_down(),
        (numerator.hi / denominator.lo).next_up(),
    )
}
fn increment(p: &[B; 12], counter_clockwise: bool) -> Option<I> {
    let position = hull(p.iter().copied());
    let derivative = hull(
        p.windows(2)
            .map(|q| std::array::from_fn(|a| I::point(11.) * (q[1][a] - q[0][a]))),
    );
    let angular = position[0] * derivative[1] - position[1] * derivative[0];
    let angular = if counter_clockwise { angular } else { -angular };
    if angular.lo <= 0. {
        return None;
    }
    // A chosen literal reference is only an untrusted half-plane proposal.
    // Whole coefficient convex hull positivity earns the actual containment,
    // excludes the origin and prevents a hidden extra turn in this leaf.
    let reference = p[0].map(|x| x.lo / 2. + x.hi / 2.);
    if p.iter()
        .any(|q| (I::point(reference[0]) * q[0] + I::point(reference[1]) * q[1]).lo <= 0.)
    {
        return None;
    }
    positive_chord_angle([p[0][0], p[0][1]], [p[11][0], p[11][1]], counter_clockwise)
}
/// Rigorous positive angle of one literal straight chord. Positive cross/dot
/// and ratio<=1/2 earn an origin-free sector narrower thanpi; callers charge
/// their own bounded endpoint arithmetic on their continuing ledger.
pub(crate) fn positive_chord_angle(a: [I; 2], b: [I; 2], counter_clockwise: bool) -> Option<I> {
    let dot = a[0] * b[0] + a[1] * b[1];
    let cross = a[0] * b[1] - a[1] * b[0];
    let cross = if counter_clockwise { cross } else { -cross };
    if dot.lo <= 0. || cross.lo <= 0. {
        return None;
    }
    // Both intervals are strictly positive here. Divide their endpoints
    // directly: forming reciprocal(dot) first can overflow for a valid tiny
    // literal loop even though its angle ratio is finite and small.
    let ratio = positive_ratio(cross, dot);
    if !ratio.valid() || ratio.hi > 0.5 {
        return None;
    }
    Some(atan_small(ratio))
}
/// Source-coordinate work counts bounded interval/scalar primitives, not CPU
/// instructions. Each leaf precharges512 for hull/derivative/half-plane/angle
/// arithmetic. A split additionally precharges1024 for its66*3 midpoint
/// operations and both arrays. The <=33-leaf DFS frontier is <24KiB;128KiB
/// covers its vector, current/split/derivative arrays and fixed metadata.
pub fn certify<'a>(
    input: Input<'a>,
    options: Options,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<Certificate<'a>, Error> {
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
    {
        return Err(Error::Invalid);
    }
    budget.charge(input.curves.len())?;
    let mut count = 0usize;
    for curve in input.curves {
        if curve.poles.len() > POINTS {
            return Err(Error::Limit);
        }
        // At most64*524288 controls, which fits even usize32.
        count += curve.poles.len();
    }
    if count > POINTS {
        return Err(Error::Limit);
    }
    let input_bytes = count * std::mem::size_of::<P>() + std::mem::size_of_val(input.curves);
    if input_bytes + SCRATCH > options.maximum_bytes {
        return Err(Error::Limit);
    }
    let plane = input.curves[0].poles.first().ok_or(Error::Invalid)?[2];
    for (i, c) in input.curves.iter().enumerate() {
        budget.charge(1 + input.curves.len())?;
        if c.spans == 0
            || c.spans.checked_mul(11).and_then(|n| n.checked_add(1)) != Some(c.poles.len())
            || input.curves[..i]
                .iter()
                .any(|p| p.curve_owner == c.curve_owner || p.start_owner == c.start_owner)
        {
            return Err(Error::Invalid);
        }
        for p in c.poles {
            budget.charge(8)?;
            if p.iter().any(|x| !x.is_finite() || x.abs() > 1e7) || p[2] != plane {
                return Err(Error::Invalid);
            }
        }
        let next = &input.curves[(i + 1) % input.curves.len()];
        if c.end_owner != next.start_owner || c.poles.last() != next.poles.first() {
            return Err(Error::Invalid);
        }
    }
    let geometry = geometry(
        input.curves.len(),
        |i| input.curves[i].spans,
        |curve, span| {
            std::array::from_fn(|i| input.curves[curve].poles[11 * span + i].map(I::point))
        },
        input.counter_clockwise,
        options,
        &mut budget,
    )?;
    Ok(Certificate {
        input,
        plane_z: plane,
        angle: geometry.angle,
        work: options.maximum_work - budget.remaining,
        nodes: geometry.nodes,
        leaves: geometry.leaves,
        input_bytes,
    })
}
struct Geometry {
    angle: [f64; 2],
    nodes: usize,
    leaves: usize,
}
fn geometry(
    curve_count: usize,
    spans: impl Fn(usize) -> usize,
    controls: impl Fn(usize, usize) -> [B; 12],
    counter_clockwise: bool,
    options: Options,
    budget: &mut Budget<'_>,
) -> Result<Geometry, Error> {
    let mut total = I::point(0.);
    let mut leaves = 0usize;
    for curve in 0..curve_count {
        for span in 0..spans(curve) {
            let p = controls(curve, span);
            let mut pending = vec![(p, 0u8)];
            while let Some((p, depth)) = pending.pop() {
                budget.charge(512)?;
                budget.nodes += 1;
                if budget.nodes > options.maximum_nodes {
                    return Err(Error::Limit);
                }
                if let Some(angle) = increment(&p, counter_clockwise) {
                    total = total + angle;
                    leaves += 1;
                    continue;
                }
                if depth == options.maximum_depth {
                    return Err(Error::Unresolved);
                }
                budget.charge(1024)?;
                let [a, b] = split_half_controls(p);
                pending.push((b, depth + 1));
                pending.push((a, depth + 1));
            }
        }
    }
    // Exact numeric endpoint closure makes the total argument advance an
    // integer multiple of2pi. Strict angular progress and the half-plane
    // enclosures earned each local lift.4<total<8 admits only winding1:
    //3<pi<4 implies2pi in(6,8), whereas4pi>12. No inferred rounded turn count.
    if !total.valid() || total.lo <= 4. || total.hi >= 8. {
        return Err(Error::Unresolved);
    }
    Ok(Geometry {
        angle: [total.lo, total.hi],
        nodes: budget.nodes,
        leaves,
    })
}
#[cfg(test)]
mod tests;
