//! Taylor's integral remainder uses the order12 derivative enclosure over the
//! entire cell. Closed face domains are partitioned, never sampled or truncated.
use super::{series::ORDER, source::Part, *};
pub(super) struct Work<'a> {
    remaining: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl<'a> Work<'a> {
    pub fn new(
        options: FacePropertyOptions,
        interrupt: Option<&'a dyn Fn() -> bool>,
    ) -> Result<Self, Error> {
        let w = Self {
            remaining: options.maximum_cells,
            interrupt,
        };
        w.checkpoint()?;
        if [
            options.area_tolerance_mm2,
            options.first_moment_tolerance_mm3,
            options.centroid_tolerance_mm,
        ]
        .into_iter()
        .any(|v| !v.is_finite() || v <= 0.)
            || options.maximum_cells == 0
            || options.maximum_cells > 65536
        {
            return Err(Error::Invalid);
        }
        Ok(w)
    }
    pub fn checkpoint(&self) -> Result<(), Error> {
        if self.interrupt.is_some_and(|f| f()) {
            Err(Error::Cancelled)
        } else {
            Ok(())
        }
    }
    fn cell(&mut self) -> Result<(), Error> {
        self.checkpoint()?;
        self.remaining = self.remaining.checked_sub(1).ok_or(Error::Limit)?;
        Ok(())
    }
}
fn power(v: I, n: usize) -> I {
    (0..n).fold(I::point(1.), |r, _| r * v)
}
fn cell(part: &Part, a: f64, b: f64) -> Result<[I; 4], Error> {
    let c = a + (b - a) * 0.5;
    let center = part.field(I::point(c))?;
    let whole = part.field(I::new(a, b))?;
    let lo = I::point(a) - I::point(c);
    let hi = I::point(b) - I::point(c);
    let moments: [I; ORDER + 1] =
        std::array::from_fn(|n| (power(hi, n + 1) - power(lo, n + 1)) / I::point((n + 1) as f64));
    let result = std::array::from_fn(|axis| {
        let polynomial =
            (0..ORDER).fold(I::point(0.), |sum, n| sum + center[axis].0[n] * moments[n]);
        polynomial + whole[axis].0[ORDER] * moments[ORDER]
    });
    checked_values(result)
}
pub(super) fn parts(
    parts: &[Part],
    extent: f64,
    options: FacePropertyOptions,
    work: &mut Work<'_>,
) -> Result<[I; 4], Error> {
    let mut minimum = I::point(0.);
    for part in parts {
        work.cell()?;
        minimum = minimum + part.field(I::new(0., 1.))?[0].0[0];
    }
    if minimum.lo <= 0. || !minimum.valid() {
        return Err(Error::Unresolved);
    }
    let centroid_budget =
        (I::point(options.centroid_tolerance_mm) * I::point(minimum.lo) / I::point(32.)).lo;
    let area_budget = options
        .area_tolerance_mm2
        .min(minimum.lo / 32.)
        .min(centroid_budget / (extent + 1.));
    let moment_budget = options.first_moment_tolerance_mm3.min(centroid_budget);
    if area_budget <= 0. || moment_budget <= 0. {
        return Err(Error::Unresolved);
    }
    let denominator = I::point(8. * parts.len() as f64);
    let per_part = [
        (I::point(area_budget) / denominator).lo,
        (I::point(moment_budget) / denominator).lo,
    ];
    if per_part.into_iter().any(|v| v <= 0. || !v.is_finite()) {
        return Err(Error::Unresolved);
    }
    let mut sum = [I::point(0.); 4];
    let mut pending = Vec::with_capacity(128);
    for (i, _) in parts.iter().enumerate() {
        pending.push((i, 0., 1., per_part));
    }
    while let Some((index, a, b, tolerance)) = pending.pop() {
        work.cell()?;
        let value = cell(&parts[index], a, b)?;
        if value[0].width() <= tolerance[0] && value[1..].iter().all(|v| v.width() <= tolerance[1])
        {
            for axis in 0..4 {
                sum[axis] = sum[axis] + value[axis]
            }
            continue;
        }
        let mid = a + (b - a) * 0.5;
        if mid <= a || mid >= b || tolerance.into_iter().any(|v| v * 0.5 == 0.) {
            return Err(Error::Unresolved);
        }
        push_children(&mut pending, index, a, mid, b, tolerance)?;
    }
    Ok(sum)
}

type Pending = (usize, f64, f64, [f64; 2]);
fn push_children(
    pending: &mut Vec<Pending>,
    index: usize,
    a: f64,
    mid: f64,
    b: f64,
    tolerance: [f64; 2],
) -> Result<(), Error> {
    if pending.len() + 2 > 2048 {
        return Err(Error::Limit);
    }
    let tolerance = tolerance.map(|v| v * 0.5);
    pending.push((index, a, mid, tolerance));
    pending.push((index, mid, b, tolerance));
    Ok(())
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn empty_domains_and_full_refinement_stack_never_publish_partial_integrals() {
        let opts = FacePropertyOptions::default();
        let mut work = Work::new(opts, None).unwrap();
        assert!(matches!(
            parts(&[], 1., opts, &mut work),
            Err(Error::Unresolved)
        ));
        let mut pending = vec![(0, 0., 1., [1.; 2]); 2048];
        assert_eq!(
            push_children(&mut pending, 0, 0., 0.5, 1., [1.; 2]),
            Err(Error::Limit)
        );
        assert_eq!(pending.len(), 2048);
        pending.clear();
        push_children(&mut pending, 0, 0., 0.5, 1., [1.; 2]).unwrap();
        assert_eq!(
            pending,
            vec![(0, 0., 0.5, [0.5; 2]), (0, 0.5, 1., [0.5; 2])]
        );
    }
}
