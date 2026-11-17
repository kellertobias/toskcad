//! Explicit immutable preparation for requested curved properties. No body,
//! decoder, query or selected-planar helper forces this plan's construction.
use super::*;
use std::sync::Arc;
#[derive(Clone, Copy, Default)]
pub struct HeadGeometryOptions {
    pub plan: LiteralGeometryOptions,
    pub additional_resident_bytes: usize,
}
enum Owner {
    Local(Arc<OwnedActualHead>),
    Mapped(Arc<transport::MappedActualHead>),
}
impl Owner {
    fn local(&self) -> &OwnedActualHead {
        match self {
            Self::Local(v) => v,
            Self::Mapped(v) => v.local_geometry(),
        }
    }
    fn retained(&self) -> usize {
        match self {
            Self::Local(v) => v.receipt().retained_bytes,
            Self::Mapped(v) => v.receipt().retained_bytes,
        }
    }
}
pub struct HeadGeometryPlan {
    owner: Owner,
    auxiliary: LiteralFaceGeometryPlan,
    ordinal: usize,
    preparation: usize,
    retained: usize,
    peak: usize,
    prior_external: usize,
}
impl HeadGeometryPlan {
    pub fn new(
        source: Arc<OwnedActualHead>,
        ordinal: usize,
        options: HeadGeometryOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self> {
        Self::prepare(Owner::Local(source), ordinal, options, interrupt)
    }
    /// Retain the original mapped file as well as its shared local owner;
    /// mapped-only bytes cannot disappear from subsequent local properties.
    pub fn new_mapped(
        source: Arc<transport::MappedActualHead>,
        ordinal: usize,
        options: HeadGeometryOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self> {
        Self::prepare(Owner::Mapped(source), ordinal, options, interrupt)
    }
    fn prepare(
        owner: Owner,
        ordinal: usize,
        options: HeadGeometryOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self> {
        if options.plan.maximum_work == 0
            || options.plan.maximum_work > MEASURE_WORK
            || options.plan.maximum_bytes == 0
            || options.plan.maximum_bytes > MAX_BYTES
        {
            return Err(Error::Invalid);
        }
        let source = owner.local();
        if ordinal >= 9 {
            return Err(if ordinal >= 17 {
                Error::Invalid
            } else {
                Error::Unsupported
            });
        }
        let entity = source.physical_roles().thread_face_entities[ordinal];
        let face = source
            .auxiliary_thread()
            .faces()
            .iter()
            .position(|f| f.entity() == entity)
            .ok_or(Error::Invalid)?;
        if !matches!(
            source.auxiliary_thread().faces()[face].support(),
            Support::Tensor(_)
        ) {
            return Err(Error::Unsupported);
        }

        let outside = owner
            .retained()
            .checked_sub(source.auxiliary_thread().retained_budget_bytes())
            .and_then(|n| n.checked_add(size_of::<Self>()))
            .and_then(|n| n.checked_add(MEASURE_SCRATCH))
            .and_then(|n| n.checked_add(options.additional_resident_bytes))
            .ok_or(Error::Limit)?;
        let bytes = owner
            .retained()
            .checked_add(size_of::<Self>())
            .and_then(|n| n.checked_add(MEASURE_SCRATCH))
            .and_then(|n| n.checked_add(options.additional_resident_bytes))
            .ok_or(Error::Limit)?;
        if bytes > options.plan.maximum_bytes {
            return Err(Error::Limit);
        }
        let mut b = Budget {
            used: 0,
            remaining: options.plan.maximum_work,
            bytes,
            peak: bytes,
            limit: options.plan.maximum_bytes,
            interrupt,
        };
        b.step(256)?;
        if b.remaining == 0 {
            return Err(Error::Limit);
        }
        let auxiliary = source.auxiliary_thread().face_geometry_plan(
            face,
            LiteralGeometryOptions {
                maximum_work: b.remaining,
                maximum_bytes: options
                    .plan
                    .maximum_bytes
                    .checked_sub(outside)
                    .ok_or(Error::Limit)?,
            },
            interrupt,
        )?;
        b.step(auxiliary.preparation_work_used())?;
        b.peak = b.peak.max(
            auxiliary
                .construction_peak_bytes()
                .checked_add(outside)
                .ok_or(Error::Limit)?,
        );
        b.step(0)?;
        // The child plan header moves into this inline field; only its actual
        // additional normal-net buffers remain in addition to our own header.
        let retained = owner
            .retained()
            .checked_add(auxiliary.retained_budget_bytes())
            .and_then(|n| n.checked_sub(source.auxiliary_thread().retained_budget_bytes()))
            .and_then(|n| n.checked_sub(size_of::<LiteralFaceGeometryPlan>()))
            .and_then(|n| n.checked_add(size_of::<Self>()))
            .ok_or(Error::Limit)?;
        Ok(Self {
            owner,
            auxiliary,
            ordinal,
            preparation: b.used,
            retained,
            peak: b.peak,
            prior_external: options.additional_resident_bytes,
        })
    }
    pub fn source(&self) -> &OwnedActualHead {
        self.owner.local()
    }
    pub fn mapped_source(&self) -> Option<&transport::MappedActualHead> {
        match &self.owner {
            Owner::Local(_) => None,
            Owner::Mapped(v) => Some(v),
        }
    }
    pub fn selected_face_ordinal(&self) -> usize {
        self.ordinal
    }
    pub fn preparation_work_used(&self) -> usize {
        self.preparation
    }
    /// Complete actual owner and normal nets, sharing local source bytes once.
    pub fn retained_budget_bytes(&self) -> usize {
        self.retained
    }
    pub fn construction_peak_bytes(&self) -> usize {
        self.peak
    }
    /// Measures the LOCAL sewn reference, even when the original mapped
    /// owner is retained. Actual file B and extra Native A are not applied here.
    /// Fresh integration and inherited preparation share the SAME 512M caller
    /// allowance. Already-earned immutable material construction is separate.
    pub fn face_measure(
        &self,
        ordinal: usize,
        options: HeadMeasureOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<HeadFaceMeasure<'_>> {
        checked(options)?;
        if ordinal >= 17 {
            return Err(Error::Invalid);
        }
        if ordinal != self.ordinal {
            return Err(Error::Unsupported);
        }
        let source = self.source();
        let bytes = self
            .retained
            .checked_add(options.additional_resident_bytes)
            .and_then(|n| n.checked_add(size_of::<HeadFaceMeasure<'_>>()))
            .and_then(|n| n.checked_add(MEASURE_SCRATCH))
            .ok_or(Error::Limit)?;
        let prior_peak = self
            .peak
            .checked_sub(self.prior_external)
            .and_then(|n| n.checked_add(options.additional_resident_bytes))
            .and_then(|n| n.checked_add(size_of::<HeadFaceMeasure<'_>>()))
            .ok_or(Error::Limit)?;
        if bytes.max(prior_peak) > options.measure.maximum_bytes {
            return Err(Error::Limit);
        }
        let mut b = Budget {
            used: self.preparation,
            remaining: options
                .measure
                .maximum_work
                .checked_sub(self.preparation)
                .ok_or(Error::Limit)?,
            bytes,
            peak: bytes.max(prior_peak),
            limit: options.measure.maximum_bytes,
            interrupt,
        };
        b.step(256)?;
        let entity = source.physical_roles().thread_face_entities[ordinal];
        b.step(512)?;
        let scale = source.axial_sewing_transfer().positive_axial_scale_bounds();
        let minimum = scale[0].min(1.);
        let maximum = scale[1].max(1.);
        let span = source.auxiliary_thread().local_enclosure_mm();
        let shift: [I; 3] = std::array::from_fn(|k| {
            (I::point(span[k][1]) - I::point(span[k][0])) * (I::point(maximum) - I::point(minimum))
                / (I::point(2.) * I::point(minimum))
        });
        let stretch = scale[1].max(1.);
        if shift.iter().any(|s| {
            !s.valid()
                || (I::point(2.) * I::point(stretch) * *s).hi >= options.measure.centroid_width_mm
        }) {
            return Err(Error::Unsupported);
        }
        if b.remaining == 0 {
            return Err(Error::Limit);
        }
        let outside = bytes
            .checked_sub(self.auxiliary.retained_budget_bytes())
            .ok_or(Error::Limit)?;
        let area_width = (I::point(options.measure.area_width_mm2) / I::point(maximum)).lo;
        let shift_width = I::point(2.)
            * I::point(stretch)
            * I::point(shift.iter().map(|s| s.hi).fold(0., f64::max));
        let centroid_width =
            ((I::point(options.measure.centroid_width_mm) - shift_width) / I::point(stretch)).lo;
        if !area_width.is_finite()
            || area_width <= 0.
            || !centroid_width.is_finite()
            || centroid_width <= 0.
        {
            return Err(Error::Unsupported);
        }
        let raw = self.auxiliary.face_measure(
            FaceMeasureOptions {
                area_width_mm2: area_width,
                centroid_width_mm: centroid_width,
                maximum_work: b.remaining,
                maximum_bytes: options
                    .measure
                    .maximum_bytes
                    .checked_sub(outside)
                    .ok_or(Error::Limit)?,
            },
            interrupt,
        )?;
        b.step(raw.work_used())?;
        b.peak = b.peak.max(
            raw.accounted_peak_bytes()
                .checked_add(outside)
                .ok_or(Error::Limit)?,
        );
        b.step(256)?;
        let area = raw.area_mm2();
        let area = I::new(area[0], area[1]) * I::new(minimum, maximum);
        let center = raw.centroid_mm();
        // For q=|cofactor(T)n|/|n| in[m,M], centered covariance gives
        // |mean_q(x)-mean(x)|<=span(x)*(M-m)/(2m). This weights the SAME
        // actual original face; it is not merely mapping its old centroid.
        let center = std::array::from_fn(|k| {
            let c = I::new(center[k][0], center[k][1]) + I::new(-shift[k].hi, shift[k].hi);
            [c.lo, c.hi]
        });
        let centroid = source
            .axial_sewing_transfer()
            .map_raw_box_mm(center)?
            .map(|v| I::new(v[0], v[1]));
        if !area.valid()
            || area.lo <= 0.
            || area.width() > options.measure.area_width_mm2
            || centroid
                .iter()
                .any(|v| !v.valid() || v.width() > options.measure.centroid_width_mm)
        {
            return Err(Error::Unsupported);
        }
        b.step(0)?;
        Ok(HeadFaceMeasure {
            source,
            ordinal,
            entity,
            area: [area.lo, area.hi],
            centroid: centroid.map(|v| [v.lo, v.hi]),
            work: b.used - self.preparation,
            inherited_plan_work: self.preparation,
            peak: b.peak,
        })
    }
}
#[cfg(test)]
mod tests;
