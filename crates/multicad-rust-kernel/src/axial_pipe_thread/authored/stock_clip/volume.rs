use super::*;
/// Bounds on the original authored clipped law, not decoded or mesh mass.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StockClipVolume {
    volume: [f64; 2],
    work: usize,
    bytes: usize,
}
impl StockClipVolume {
    pub fn source_volume_mm3(&self) -> [f64; 2] {
        self.volume
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.bytes
    }
}
fn mean(a: I, b: I) -> I {
    (a + b) / I::point(2.)
}
fn square_mean(a: I, b: I) -> I {
    (a.square() + a * b + b.square()) / I::point(3.)
}
fn product_mean(a: [I; 2], b: [I; 2]) -> I {
    (a[0] * b[0] + (a[0] * b[1] + a[1] * b[0]) / I::point(2.) + a[1] * b[1]) / I::point(3.)
}
fn triple_mean(c: [I; 2], b: [I; 2]) -> I {
    // Integral of the cubic Bernstein product c(t)*b(t)^2. Every
    // coefficient is nonnegative in this earned source-domain partition.
    (c[0] * b[0].square()
        + (I::point(2.) * c[0] * b[0] * b[1] + c[1] * b[0].square()) / I::point(3.)
        + (c[0] * b[1].square() + I::point(2.) * c[1] * b[0] * b[1]) / I::point(3.)
        + c[1] * b[1].square())
        / I::point(4.)
}
impl StockClipPlan {
    pub fn volume_bounds(
        &self,
        maximum_width_mm3: f64,
        options: StockClipOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<StockClipVolume, Error> {
        if !maximum_width_mm3.is_finite() || maximum_width_mm3 <= 0. {
            return Err(Error::Invalid);
        }
        let mut work = Meter::new(
            options,
            size_of::<Self>() + size_of::<StockClipVolume>(),
            interrupt,
        )?;
        work.charge(self.work)?;
        work.charge(512)?;
        let p = I::point(self.source.spec.pitch_mm);
        let k = I::point(self.source.spec.radial_taper);
        let d = I::point(self.source.depth());
        let limit = I::point(self.stock_radius_limit_mm());
        let slopes = self.slopes[0] + self.slopes[1];
        let b0 = self.root_radius;
        // The SAME symbolic contact identity fixes B(zC)=stock−depth.
        let bc = limit - d;
        let be = self.root_radius + k * I::point(self.source.spec.length_mm);
        let b2 = square_mean(b0, bc);
        let natural = self.root_width * b2
            + self.width * square_mean(b0 + d, limit)
            + d * slopes * (b2 + d * mean(b0, bc) + d.square() / I::point(3.));
        let c = [d, limit - be];
        let b = [bc, be];
        let plateau = mean(self.width, self.width + (d - c[1]) * slopes);
        let clipped = self.root_width * square_mean(b[0], b[1])
            + slopes
                * (triple_mean(c, b)
                    + limit * product_mean(c, b)
                    + limit.square() * mean(c[0], c[1]))
                / I::point(3.)
            + plateau * limit.square();
        // dtheta = -(2pi/p)dq. Half the squared-radius angular integral
        // gives pi/p times these phase integrals, independent of roll/phase.
        let raw = (natural * self.event
            + clipped * (I::point(self.source.spec.length_mm) - self.event))
            * pi()
            / p;
        // Each sealed dimension is <=1e7, p>=.001 and all fixed-degree
        // positive products here have magnitude <1e35. The enclosure is
        // finite. Intersect the separately earned positive material volume
        // with [0,+inf); a zero numeric lower bound is not empty material.
        let volume = I::new(raw.lo.max(0.), raw.hi);
        if volume.width() > maximum_width_mm3 {
            return Err(Error::Unresolved);
        }
        work.charge(1)?;
        Ok(StockClipVolume {
            volume: bounds(volume),
            work: work.used,
            bytes: size_of::<Self>() + size_of::<StockClipVolume>() + SCRATCH,
        })
    }
}
