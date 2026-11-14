//! Actual band coefficients with whole-cell Gauss2/C4 remainder.
use super::*;
use crate::tapered_thread::{Scalar, Series};
type S = Series<5>;
struct Work<'a> {
    left: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl Work<'_> {
    fn cell(&mut self) -> Result<(), Error> {
        checkpoint(self.interrupt)?;
        self.left = self.left.checked_sub(1).ok_or(Error::Limit)?;
        Ok(())
    }
}
impl Hanger {
    pub fn volume(
        &self,
        options: QueryOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Volume, Error> {
        checkpoint(interrupt)?;
        if !options.volume_error_mm3.is_finite() || options.volume_error_mm3 <= 0. {
            return Err(Error::Invalid);
        }
        if options.maximum_cells == 0 || options.maximum_cells > 65536 {
            return Err(Error::Limit);
        }
        let mut work = Work {
            left: options.maximum_cells,
            interrupt,
        };
        let mut wood = I::point(0.);
        for b in self.wood().bands() {
            work.cell()?;
            let [a, c] = b.radius_mm.map(I::point);
            wood = wood
                + (I::point(b.phase_mm[1]) - I::point(b.phase_mm[0]))
                    * (a.square() + a * c + c.square())
                    / I::point(3.);
        }
        wood = wood * pi() * I::point(self.wood().full_length_mm())
            / I::point(self.wood().spec().pitch_mm);
        work.cell()?;
        let plain = pi()
            * I::point(self.spec.plain_radius_mm).square()
            * I::point(self.spec.plain_length_mm);
        let thread = *self.metric();
        let spec = thread.spec();
        let reserve = thread
            .phase_domain_volume_error_mm3([0., spec.length_mm])
            .map_err(|_| Error::Unresolved)?;
        let available = (I::point(options.volume_error_mm3) - I::point(2.) * I::point(reserve)).lo;
        if available <= 0. {
            return Err(Error::Unresolved);
        }
        let factor = pi() * I::point(spec.length_mm) / I::point(spec.pitch_mm);
        let mut metric = I::point(0.);
        for band in 0..5 {
            let mut pending = vec![(0., 1., available / 32., 0u8)];
            while let Some((a, b, budget, depth)) = pending.pop() {
                work.cell()?;
                if let Some(value) = cell(thread, band, a, b, factor)
                    && value.width() <= budget
                {
                    metric = metric + value;
                    continue;
                }
                pending.extend(children(a, b, budget, depth)?);
            }
        }
        metric = metric + I::new(-reserve, reserve);
        let total = wood + plain + metric;
        if !total.valid() || total.lo <= 0. || total.width() > options.volume_error_mm3 {
            return Err(Error::Unresolved);
        }
        let result = Volume {
            material_mm3: pair(total),
            wood_mm3: pair(wood),
            plain_mm3: pair(plain),
            metric_mm3: pair(metric),
            metric_phase_error_mm3: reserve,
        };
        checkpoint(interrupt)?;
        Ok(result)
    }
}
fn field(thread: Thread, band: usize, t: I) -> S {
    let knots = thread.profile_breaks();
    let span = I::point(knots[band + 1]) - I::point(knots[band]);
    let phase = S::constant(knots[band]) + S::variable(t) * S::interval(span);
    let major = S::constant(thread.spec().major_radius_mm);
    let r = match band {
        1 | 3 => {
            let base = if band == 1 { knots[1] } else { knots[4] };
            major + (phase - S::constant(base)) * S::constant(thread.slope_at_phase(base, band))
        }
        2 => {
            let root = thread.root_circle().expect("sealed rounded source");
            let x = phase - S::constant(root.phase);
            S::constant(root.center_radius)
                - (S::constant(root.radius) * S::constant(root.radius) - x * x).sqrt()
        }
        _ => major,
    };
    r * r * S::interval(span)
}
fn cell(thread: Thread, band: usize, a: f64, b: f64, factor: I) -> Option<I> {
    let width = I::point(b) - I::point(a);
    let mid = (I::point(a) + I::point(b)) / I::point(2.);
    let offset = width / I::point(12.).sqrt();
    let left = field(thread, band, mid - offset);
    let right = field(thread, band, mid + offset);
    let whole = field(thread, band, I::new(a, b));
    if left
        .c
        .into_iter()
        .chain(right.c)
        .chain(whole.c)
        .any(|v| !v.valid())
    {
        return None;
    }
    let e =
        (I::point(whole.c[4].magnitude()) * width.square().square() * width / I::point(180.)).hi;
    let value = ((left.c[0] + right.c[0]) * width / I::point(2.) + I::new(-e, e)) * factor;
    value.valid().then_some(value)
}

// One depth-first binary path, at most49 pending siblings. No source band or
// numeric cell can bypass the finite refinement/parameter-resolution guard.
type Pending = (f64, f64, f64, u8);
fn children(a: f64, b: f64, budget: f64, depth: u8) -> Result<[Pending; 2], Error> {
    let mid = a.midpoint(b);
    if depth >= 48 || mid <= a || mid >= b || budget * 0.5 == 0. {
        return Err(Error::Unresolved);
    }
    Ok([
        (mid, b, budget * 0.5, depth + 1),
        (a, mid, budget * 0.5, depth + 1),
    ])
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn actual_source_numeric_cells_and_finite_refinement_limits_never_admit_invalid_proofs() {
        let t = super::super::tests::fixture(false, false).metric;
        assert!(cell(t, 2, f64::MAX, f64::MAX, I::point(1.)).is_none());
        assert_eq!(children(0., 1., 1., 48), Err(Error::Unresolved));
        assert_eq!(children(1., 1_f64.next_up(), 1., 0), Err(Error::Unresolved));
        assert_eq!(
            children(0., 1., f64::from_bits(1), 0),
            Err(Error::Unresolved)
        );
        assert_eq!(
            children(0., 1., 1., 0).unwrap(),
            [(0.5, 1., 0.5, 1), (0., 0.5, 0.5, 1)]
        );
    }
}
