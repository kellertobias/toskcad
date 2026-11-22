//! Exact-source affine and circular section integrals, outward evaluated.
use super::*;
use crate::authored_metric_thread::Knot;
struct Work<'a> {
    remaining: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl Work<'_> {
    fn cell(&mut self) -> Result<(), Error> {
        checkpoint(self.interrupt)?;
        self.remaining = self.remaining.checked_sub(1).ok_or(Error::Limit)?;
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
        let mut w = Work {
            remaining: options.maximum_cells,
            interrupt,
        };
        let mut wood = I::point(0.);
        for b in self.wood().bands() {
            w.cell()?;
            let [a, b_radius] = b.radius_mm.map(I::point);
            let width = I::point(b.phase_mm[1]) - I::point(b.phase_mm[0]);
            wood = wood + width * (a.square() + a * b_radius + b_radius.square()) / I::point(3.);
        }
        wood = wood * pi() * I::point(self.wood().full_length_mm())
            / I::point(self.wood().spec().pitch_mm);
        w.cell()?;
        let plain = pi()
            * I::point(self.spec.plain_radius_mm).square()
            * I::point(self.spec.plain_length_mm);
        let metric = metric_integral(self.metric(), &mut w)?;
        let total = wood + plain + metric;
        if !total.valid() || total.lo <= 0. || total.width() > options.volume_error_mm3 {
            return Err(Error::Unresolved);
        }
        let result = Volume {
            material_mm3: pair(total),
            wood_mm3: pair(wood),
            plain_mm3: pair(plain),
            metric_mm3: pair(metric),
        };
        checkpoint(interrupt)?;
        Ok(result)
    }
}
fn metric_integral(source: &AuthoredMetricThread, w: &mut Work<'_>) -> Result<I, Error> {
    let knots = [
        Knot::Begin,
        Knot::CrestEnd,
        Knot::UpperRoot,
        Knot::LowerRoot,
        Knot::CrestStart,
        Knot::End,
    ]
    .map(|k| source.knot(k));
    let mut integral = I::point(0.);
    for band in 0..5 {
        w.cell()?;
        let value = if band == 2 {
            // The OWNED canonical circle has half trim width a=p/8, tangent
            // height R/2 and asin(a/R)=π/3 by its sharedQ(sqrt3) construction.
            // ∫(C−sqrt(R²−x²))²dx across[-a,a] is evaluated without fitting.
            let a = I::point(source.spec().pitch_mm) / I::point(8.);
            let radius = I::new(source.root_radius_mm()[0], source.root_radius_mm()[1]);
            let tangent = I::new(knots[2].radius_mm[0], knots[2].radius_mm[1]);
            let center = tangent + radius / I::point(2.);
            I::point(2.) * a * (center.square() + radius.square())
                - I::point(2.) * a * a * a / I::point(3.)
                - I::point(2.)
                    * center
                    * (a * radius / I::point(2.) + pi() * radius.square() / I::point(3.))
        } else {
            let width = I::new(knots[band + 1].phase_mm[0], knots[band + 1].phase_mm[1])
                - I::new(knots[band].phase_mm[0], knots[band].phase_mm[1]);
            let a = I::new(knots[band].radius_mm[0], knots[band].radius_mm[1]);
            let b = I::new(knots[band + 1].radius_mm[0], knots[band + 1].radius_mm[1]);
            width * (a.square() + a * b + b.square()) / I::point(3.)
        };
        integral = integral + value;
    }
    // At every axial station the actual one-lead chart covers N full periods;
    // its angular Jacobian true2π/(N*p) cancels that SAME source's start count.
    Ok(pi() * I::point(source.full_length_mm()) * integral / I::point(source.spec().pitch_mm))
}
