//! Bounded source-curve evaluator retaining one exact rotation packet.
//! No display/manufacturing geometry is used to define these curves.
use super::*;
#[cfg(test)]
mod tests;
pub struct PeriodicCurveQuery<'s, 'i> {
    source: &'s OwnedPeriodicMaterial,
    rotation: rotation::Rotation,
    rate: I,
    caps: [I; 2],
    work: Work<'i>,
    initial: usize,
}
impl OwnedPeriodicMaterial {
    pub fn source_curve_query<'s, 'i>(
        &'s self,
        maximum_work: usize,
        interrupt: Option<&'i dyn Fn() -> bool>,
    ) -> FResult<PeriodicCurveQuery<'s, 'i>> {
        let mut work = Work {
            remaining: maximum_work,
            interrupt,
            cancelled: Cell::new(false),
        };
        let prepared = (|| {
            work.step(1)?;
            if maximum_work > WORK {
                return Err(PhillipsReadError::Invalid);
            }
            let rotation = rotation::Rotation::new(&mut work)?;
            let rate = I::point(2.) * rotation.pi(&mut work)? / I::point(self.thread.pitch_mm());
            Ok((rotation, rate))
        })();
        if work.cancelled.get() {
            return Err(PhillipsReadError::Cancelled);
        }
        let (rotation, rate) = prepared?;
        Ok(PeriodicCurveQuery {
            source: self,
            rotation,
            rate,
            caps: self.caps.map(expression),
            work,
            initial: maximum_work,
        })
    }
}
impl PeriodicCurveQuery<'_, '_> {
    fn finish<T>(&self, result: Result<T>) -> FResult<T> {
        if self.work.cancelled.get() {
            Err(PhillipsReadError::Cancelled)
        } else {
            result.map_err(Into::into)
        }
    }
    fn step(&mut self, n: usize) -> FResult<()> {
        let result = self.work.step(n);
        self.finish(result)
    }
    pub fn source(&self) -> &OwnedPeriodicMaterial {
        self.source
    }
    pub fn work_used(&self) -> usize {
        self.initial - self.work.remaining
    }
    pub fn remaining_work(&self) -> usize {
        self.work.remaining
    }
    /// Fixed2048-bit reduction and degree31 Taylor scratch, including the
    /// retained four735-bit integers. This is an operation-storage bound.
    pub fn scratch_upper_bound_bytes(&self) -> usize {
        64 * 1024
    }
    fn parameter(&mut self, t: f64) -> FResult<I> {
        self.step(256)?;
        if !t.is_finite() || !(0. ..=1.).contains(&t) {
            return Err(PhillipsReadError::Invalid);
        }
        Ok(I::point(t))
    }
    fn point(
        &mut self,
        band: usize,
        u: I,
        z: I,
        radius: Option<I>,
        physical_z: Option<I>,
    ) -> FResult<[[f64; 2]; 3]> {
        if band >= 4 {
            return Err(PhillipsReadError::Invalid);
        }
        let p = self.source.thread.phases();
        let phase = p[band] + u * (p[band + 1] - p[band]);
        let hand = I::point(if self.source.thread.left_hand() {
            -1.
        } else {
            1.
        });
        let result = self
            .rotation
            .sin_cos(hand * self.rate * (z - phase), &mut self.work);
        let [c, s] = self.finish(result)?;
        let r = radius.unwrap_or_else(|| self.source.thread.radii(band, u));
        let p = [r * c, r * s, physical_z.unwrap_or(z)];
        require(p.into_iter().all(I::valid))?;
        Ok(p.map(|v| [v.lo, v.hi]))
    }
    pub fn thread_rail_point_mm(&mut self, band: usize, t: f64) -> FResult<[[f64; 2]; 3]> {
        let t = self.parameter(t)?;
        let z = self.caps[0] + t * (self.caps[1] - self.caps[0]);
        self.point(band, I::point(0.), z, None, None)
    }
    pub fn thread_rim_point_mm(
        &mut self,
        cap: usize,
        band: usize,
        t: f64,
    ) -> FResult<[[f64; 2]; 3]> {
        let t = self.parameter(t)?;
        let z = *self.caps.get(cap).ok_or(PhillipsReadError::Invalid)?;
        self.point(band, t, z, None, None)
    }
    pub fn join_rim_point_mm(&mut self, band: usize, t: f64) -> FResult<[[f64; 2]; 3]> {
        let t = self.parameter(t)?;
        let radius = I::point(self.source.bottom[0]);
        self.point(band, t, self.caps[0], Some(radius), None)
    }
    /// Four original phase bands cover a head station circle exactly once.
    /// Each phase band is sampled independently so no chord spans a change
    /// in this convenient (nonuniform angular) parameter speed.
    pub fn head_rim_point_mm(
        &mut self,
        station: usize,
        phase_band: usize,
        t: f64,
    ) -> FResult<[[f64; 2]; 3]> {
        let t = self.parameter(t)?;
        let [r, z] = self
            .source
            .head
            .view()
            .station(station)
            .ok_or(PhillipsReadError::Invalid)?;
        self.point(
            phase_band,
            t,
            self.caps[0],
            Some(I::point(r)),
            Some(I::point(z) - I::point(self.source.head_shift)),
        )
    }
    /// Whole-parameter Euclidean second-derivative bounds for the actual
    /// finite curves. The values are with respect to each method's t in[0,1].
    pub fn thread_rail_second_bound_mm(&mut self, band: usize) -> FResult<f64> {
        self.step(256)?;
        if band >= 4 {
            return Err(PhillipsReadError::Invalid);
        }
        let speed = self.rate * (self.caps[1] - self.caps[0]);
        let bound = I::point(self.source.thread.radius_nodes_mm()[band]) * speed.square();
        require(bound.valid() && bound.lo >= 0.)?;
        Ok(bound.hi)
    }
    fn angular_second(&mut self, band: usize, r: I, dr: I) -> FResult<f64> {
        self.step(256)?;
        if band >= 4 {
            return Err(PhillipsReadError::Invalid);
        }
        let p = self.source.thread.phases();
        let speed = self.rate * (p[band + 1] - p[band]);
        let bound = r * speed.square() + I::point(2.) * I::point(dr.magnitude()) * speed;
        require(bound.valid() && bound.lo >= 0.)?;
        Ok(bound.hi)
    }
    pub fn thread_rim_second_bound_mm(&mut self, band: usize) -> FResult<f64> {
        if band >= 4 {
            return Err(PhillipsReadError::Invalid);
        }
        let r = self.source.thread.radius_nodes_mm();
        self.angular_second(
            band,
            I::point(r[band].max(r[band + 1])),
            I::point(r[band + 1]) - I::point(r[band]),
        )
    }
    pub fn join_rim_second_bound_mm(&mut self, band: usize) -> FResult<f64> {
        self.angular_second(band, I::point(self.source.bottom[0]), I::point(0.))
    }
    pub fn head_rim_second_bound_mm(&mut self, station: usize, band: usize) -> FResult<f64> {
        let p = self
            .source
            .head
            .view()
            .station(station)
            .ok_or(PhillipsReadError::Invalid)?;
        self.angular_second(band, I::point(p[0]), I::point(0.))
    }
}
