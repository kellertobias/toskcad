use super::*;
fn split(span: BezierSpan, work: &mut Work<'_>) -> Result<[BezierSpan; 2], Error> {
    work.charge(2048)?;
    let mut net = span.poles.map(|p| p.map(I::point));
    let mut left = [[I::point(0.); 3]; 12];
    let mut right = left;
    left[0] = net[0];
    right[11] = net[11];
    for level in (1..12).rev() {
        for i in 0..level {
            let next = net[i + 1];
            for (axis, value) in net[i].iter_mut().enumerate() {
                *value = (*value + next[axis]) * I::point(0.5);
            }
        }
        left[12 - level] = net[0];
        right[level - 1] = net[level - 1];
    }
    let center = span.domain[0].midpoint(span.domain[1]);
    let mut children = [span; 2];
    for (index, ideal) in [left, right].into_iter().enumerate() {
        let child = &mut children[index];
        child.domain = if index == 0 {
            [span.domain[0], center]
        } else {
            [center, span.domain[1]]
        };
        child.poles = ideal.map(|p| p.map(I::midpoint));
        child.poles[if index == 0 { 0 } else { 11 }] = span.poles[if index == 0 { 0 } else { 11 }];
        let errors: [[I; 3]; 12] = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                I::new(child.poles[i][j].next_down(), child.poles[i][j].next_up()) - ideal[i][j]
            })
        });
        let norm = |p: [I; 3]| {
            p.into_iter()
                .fold(I::point(0.), |s, x| s + I::point(x.magnitude()))
        };
        let position = errors.into_iter().map(norm).fold(0f64, |a, b| a.max(b.hi));
        let derivative = (0..11)
            .map(|i| norm(std::array::from_fn(|j| errors[i + 1][j] - errors[i][j])))
            .fold(0f64, |a, b| a.max(b.hi));
        child.position_error = (I::point(span.position_error) + I::point(position)).hi;
        child.derivative_error = (I::point(span.derivative_error)
            + I::point(11.) * I::point(derivative) / I::point(child.domain[1] - child.domain[0]))
        .hi;
    }
    Ok(children)
}
pub(super) fn curve(
    curve: &mut CurveCarrier,
    count: usize,
    memory: &mut Memory,
    work: &mut Work<'_>,
) -> Result<(), Error> {
    while curve.spans.len() < count {
        let mut next = memory.allocate(curve.spans.len() * 2, work)?;
        for &span in &curve.spans {
            next.extend(split(span, work)?);
        }
        // Shared station IDs, not rounded-coordinate equality, own each join.
        let count = next.len();
        for (i, span) in next.iter_mut().enumerate() {
            span.owners = [
                curves::owner(curve.source, i, count),
                curves::owner(curve.source, i + 1, count),
            ];
        }
        memory.release(&curve.spans);
        curve.spans = next;
    }
    Ok(())
}
#[cfg(test)]
mod tests;
