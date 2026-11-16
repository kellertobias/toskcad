//! Original retained source integrals; no polyhedral/file mass is inferred.
use super::super::authored::StockClipMeasureOptions;
use super::*;
fn p(v: f64) -> I {
    I::point(v)
}
fn interval(v: [f64; 2]) -> I {
    I::new(v[0], v[1])
}
fn pair(v: I) -> [f64; 2] {
    [v.lo, v.hi]
}
fn pi() -> I {
    I::new(
        std::f64::consts::PI.next_down(),
        std::f64::consts::PI.next_up(),
    )
}
#[derive(Clone, Copy, Debug)]
pub struct Volume {
    value: I,
    work: usize,
    bytes: usize,
}
impl Volume {
    pub fn source_volume_mm3(&self) -> [f64; 2] {
        pair(self.value)
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.bytes
    }
}
#[derive(Clone, Copy, Debug)]
pub struct FaceMeasure {
    face: Face,
    values: [I; 4],
    centroid: [I; 3],
    work: usize,
    bytes: usize,
}
impl FaceMeasure {
    pub fn face(&self) -> Face {
        self.face
    }
    pub fn area_mm2(&self) -> [f64; 2] {
        pair(self.values[0])
    }
    pub fn first_moments_mm3(&self) -> [[f64; 2]; 3] {
        [
            pair(self.values[1]),
            pair(self.values[2]),
            pair(self.values[3]),
        ]
    }
    pub fn centroid_mm(&self) -> [[f64; 2]; 3] {
        self.centroid.map(pair)
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.bytes
    }
}
impl Work<'_> {
    pub(super) fn child(&self, prefix: usize) -> Result<StockClipOptions, Error> {
        if self.left == 0 {
            return Err(Error::Limit);
        }
        Ok(StockClipOptions {
            maximum_work: 16_000_000.min(prefix.saturating_add(self.left)),
            maximum_bytes: self
                .limit
                .checked_sub(size_of::<ExternalStockCut>() + size_of::<FaceMeasure>())
                .filter(|n| *n > 0)
                .ok_or(Error::Limit)?,
        })
    }
}
impl ExternalStockCut {
    pub fn volume_bounds(
        &self,
        maximum_width_mm3: f64,
        options: StockClipOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Volume, Error> {
        if !maximum_width_mm3.is_finite() || maximum_width_mm3 <= 0. {
            return Err(Error::Invalid);
        }
        let bytes = size_of::<Self>() + size_of::<Volume>() + 16 * 1024;
        let mut work = Work::new(options, bytes, interrupt)?;
        work.charge(self.work)?;
        work.charge(128)?;
        let tail = pi() * self.radius.square() * self.tail;
        let radial_area = self.radial.square();
        if !radial_area.valid() || radial_area.lo <= 0. {
            return Err(Error::Unresolved);
        }
        let available = (p(maximum_width_mm3) - p(tail.width())).lo;
        if available <= 0. {
            return Err(Error::Unresolved);
        }
        let threaded = self.thread.volume_bounds(
            (p(available) / radial_area / p(4.)).lo,
            work.child(self.thread.work_used())?,
            interrupt,
        )?;
        work.charge(threaded.work_used() - self.thread.work_used())?;
        let total = interval(threaded.source_volume_mm3()) * radial_area + tail;
        if !total.valid() || total.lo <= 0. || total.width() > maximum_width_mm3 {
            return Err(Error::Unresolved);
        }
        work.charge(1)?;
        Ok(Volume {
            value: total,
            work: work.used,
            bytes,
        })
    }
    pub fn face_measure(
        &self,
        face: Face,
        options: StockClipMeasureOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<FaceMeasure, Error> {
        if [
            options.area_width_mm2,
            options.moment_width_mm3,
            options.centroid_width_mm,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v <= 0.)
        {
            return Err(Error::Invalid);
        }
        if face == Face::Thread(StockClipFace::EndCap) {
            return Err(Error::Invalid);
        }
        let bytes = size_of::<Self>() + size_of::<FaceMeasure>() + 16 * 1024;
        let mut work = Work::new(options.limits, bytes, interrupt)?;
        work.charge(self.work)?;
        work.charge(256)?;
        let r = self.radius;
        let lo = self.join;
        let hi = self.top;
        let lambda = self.radial;
        let area_scale = lambda.square();
        let moment_scale = area_scale * lambda;
        if !area_scale.valid() || !moment_scale.valid() {
            return Err(Error::Unresolved);
        }
        let disk = pi() * r.square();
        let mut fixed_z = None;
        let mut values = match face {
            Face::Thread(key) => {
                let measure = self.thread.face_measure(
                    key,
                    StockClipMeasureOptions {
                        area_width_mm2: options.area_width_mm2 / 4.,
                        moment_width_mm3: options.moment_width_mm3 / 4.,
                        centroid_width_mm: options.centroid_width_mm / 4.,
                        limits: work.child(self.thread.work_used())?,
                    },
                    interrupt,
                )?;
                work.charge(measure.work_used() - self.thread.work_used())?;
                let area = interval(measure.area_mm2());
                let moments = measure.first_moments_mm3().map(interval);
                if key == StockClipFace::StartCap {
                    fixed_z = Some(self.bottom);
                    [
                        area * area_scale,
                        moments[0] * moment_scale,
                        moments[1] * moment_scale,
                        p(0.),
                    ]
                } else {
                    // Singular values of diag(lambda,lambda,1) bound each
                    // local area factor between min(lambda,lambda²) and max.
                    // Signed moment cancellations require absolute perturbation
                    // bounds, not multiplication by an area-factor interval.
                    let factor = I::new(lambda.lo.min(area_scale.lo), lambda.hi.max(area_scale.hi));
                    let change = (factor - p(1.)).magnitude();
                    let xy = p(change) * area * p(self.thread.stock_radius_limit_mm()) * lambda;
                    let zz = p(change) * area * p(self.thread.source().spec().length_mm);
                    [
                        area * factor,
                        moments[0] * lambda + I::new(-xy.hi, xy.hi),
                        moments[1] * lambda + I::new(-xy.hi, xy.hi),
                        moments[2] + I::new(-zz.hi, zz.hi),
                    ]
                }
            }
            Face::Shoulder => {
                let measure = self.thread.face_measure(
                    StockClipFace::EndCap,
                    StockClipMeasureOptions {
                        area_width_mm2: options.area_width_mm2 / 4.,
                        moment_width_mm3: options.moment_width_mm3 / 4.,
                        centroid_width_mm: options.centroid_width_mm / 4.,
                        limits: work.child(self.thread.work_used())?,
                    },
                    interrupt,
                )?;
                work.charge(measure.work_used() - self.thread.work_used())?;
                let area = disk - interval(measure.area_mm2()) * area_scale;
                let m = measure.first_moments_mm3().map(interval);
                fixed_z = Some(lo);
                [
                    area,
                    -m[0] * moment_scale,
                    -m[1] * moment_scale,
                    area * p(self.thread.source().spec().length_mm),
                ]
            }
            Face::TailCylinder => {
                let area = p(2.) * pi() * r * self.tail;
                [
                    area,
                    p(0.),
                    p(0.),
                    area * ((self.join - self.bottom) + (self.top - self.bottom)) / p(2.),
                ]
            }
            Face::TailCap => {
                fixed_z = Some(hi);
                [disk, p(0.), p(0.), disk * (hi - self.bottom)]
            }
        };
        values[1] = values[1] + self.center[0] * values[0];
        values[2] = values[2] + self.center[1] * values[0];
        values[3] = values[3] + self.bottom * values[0];
        if values.iter().any(|v| !v.valid())
            || values[0].lo <= 0.
            || values[0].width() > options.area_width_mm2
            || values[1..]
                .iter()
                .any(|v| v.width() > options.moment_width_mm3)
        {
            return Err(Error::Unresolved);
        }
        let mut centroid = std::array::from_fn(|i| values[i + 1] / values[0]);
        if let Some(z) = fixed_z {
            centroid[2] = z;
        }
        if centroid
            .iter()
            .any(|v| !v.valid() || v.width() > options.centroid_width_mm)
        {
            return Err(Error::Unresolved);
        }
        work.charge(1)?;
        Ok(FaceMeasure {
            face,
            values,
            centroid,
            work: work.used,
            bytes,
        })
    }
}
