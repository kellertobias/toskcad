//! Strict whole-box queries of the named sewn reference: global T(actual
//! auxiliary thread) union the earned actual polygonal head extrusion. This
//! does not classify the unsewn raw physical shell or infer author intent.
use super::*;
use crate::tapered_thread::interval::I;
use std::mem::size_of;
const QUERY_WORK: usize = 16_000_000;
const QUERY_SCRATCH: usize = 16 * 1024;
#[derive(Clone, Copy, Debug, Default)]
pub struct HeadQueryOptions {
    pub query: QueryOptions,
    pub additional_resident_bytes: usize,
}

pub struct HeadMaterialQuery<'source> {
    source: &'source OwnedActualHead,
    input: [[f64; 2]; 3],
    location: Location,
    work: usize,
    peak: usize,
}
impl<'source> HeadMaterialQuery<'source> {
    pub fn source(&self) -> &'source OwnedActualHead {
        self.source
    }
    pub fn local_bounds_mm(&self) -> [[f64; 2]; 3] {
        self.input
    }
    pub fn location(&self) -> Location {
        self.location
    }
    /// Fresh query work; immutable construction is not repeated.
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn query_budget_bytes(&self) -> usize {
        self.peak
    }
}
fn checked(input: [[f64; 2]; 3], options: HeadQueryOptions) -> Result<[I; 3]> {
    let p = input.map(|v| I::new(v[0], v[1]));
    if p.iter().any(|v| !v.valid())
        || options.query.maximum_work == 0
        || options.query.maximum_work > QUERY_WORK
        || options.query.maximum_bytes == 0
        || options.query.maximum_bytes > MAX_BYTES
        || options.query.maximum_nodes > 65_536
        || options.query.maximum_depth > 64
    {
        return Err(Error::Invalid);
    }
    Ok(p)
}
impl AxialSewingTransfer {
    /// Enclose the inverse of the SAME positive global axial sewing map.
    pub fn inverse_sewn_box_mm(self, input: [[f64; 2]; 3]) -> Result<[[f64; 2]; 3]> {
        let mut p = input.map(|b| I::new(b[0], b[1]));
        if p.iter().any(|v| !v.valid()) {
            return Err(Error::Invalid);
        }
        let lower = self.lower_raw_height_mm();
        let lower = I::new(lower[0], lower[1]);
        let scale = self.positive_axial_scale_bounds();
        p[2] = lower + (p[2] - lower) / I::new(scale[0], scale[1]);
        if p.iter().any(|v| !v.valid()) {
            return Err(Error::Unsupported);
        }
        Ok(p.map(|v| [v.lo, v.hi]))
    }
}
impl OwnedActualHead {
    pub fn classify_local_box(
        &self,
        input: [[f64; 2]; 3],
        options: HeadQueryOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<HeadMaterialQuery<'_>> {
        let p = checked(input, options)?;
        let overhead = options
            .additional_resident_bytes
            .checked_add(size_of::<HeadMaterialQuery<'_>>())
            .and_then(|n| n.checked_add(QUERY_SCRATCH))
            .ok_or(Error::Limit)?;
        let bytes = self
            .receipt()
            .retained_bytes
            .checked_add(overhead)
            .ok_or(Error::Limit)?;
        if bytes > options.query.maximum_bytes {
            return Err(Error::Limit);
        }
        let mut b = Budget {
            used: 0,
            remaining: options.query.maximum_work,
            bytes,
            peak: bytes,
            limit: options.query.maximum_bytes,
            interrupt,
        };
        b.step(256)?;
        let stations = self.head_stations_mm();
        let collar = I::point(self.planar_sewing_mm());
        let mut inside = p[2].lo > (I::point(stations[0]) + collar).hi
            && p[2].hi < (I::point(stations[1]) - collar).lo;
        let mut outside = p[2].hi < (I::point(stations[0]) - collar).lo
            || p[2].lo > (I::point(stations[1]) + collar).hi;
        let polygon = self.head_polygon_mm();
        for i in 0..6 {
            b.step(128)?;
            let a = polygon[i];
            let z = polygon[(i + 1) % 6];
            let edge = [
                I::point(z[0]) - I::point(a[0]),
                I::point(z[1]) - I::point(a[1]),
            ];
            let cross = edge[0] * (p[1] - I::point(a[1])) - edge[1] * (p[0] - I::point(a[0]));
            let length = (edge[0].square() + edge[1].square()).sqrt();
            let margin = collar * length;
            if !cross.valid() || !margin.valid() {
                return Err(Error::Unsupported);
            }
            inside &= cross.lo > margin.hi;
            outside |= cross.hi < -margin.hi;
        }
        let location = if inside {
            Location::Inside
        } else {
            b.step(128)?;
            let raw = self.axial_sewing_transfer().inverse_sewn_box_mm(input)?;
            let extra = bytes
                .checked_sub(self.auxiliary_thread().retained_budget_bytes())
                .ok_or(Error::Limit)?;
            if b.remaining == 0 {
                return Err(Error::Limit);
            }
            let result = self.auxiliary_thread().classify_local_box(
                raw,
                QueryOptions {
                    maximum_work: b.remaining,
                    maximum_bytes: options
                        .query
                        .maximum_bytes
                        .checked_sub(extra)
                        .ok_or(Error::Limit)?,
                    ..options.query
                },
                interrupt,
            )?;
            b.step(result.work_used())?;
            b.peak = b.peak.max(
                result
                    .query_budget_bytes()
                    .checked_add(extra)
                    .ok_or(Error::Limit)?,
            );
            match result.location() {
                Location::Inside => Location::Inside,
                Location::Outside if outside => Location::Outside,
                _ => Location::Unresolved,
            }
        };
        b.step(0)?;
        Ok(HeadMaterialQuery {
            source: self,
            input,
            location,
            work: b.used,
            peak: b.peak,
        })
    }
}
impl transport::MappedHeadFrame {
    /// Actual EXPRESS positive uniform frame inverse, applied to a whole box.
    pub fn inverse_box_mm(self, input: [[f64; 2]; 3]) -> Result<[[f64; 2]; 3]> {
        let input = input.map(|v| I::new(v[0], v[1]));
        if input.iter().any(|v| !v.valid()) {
            return Err(Error::Invalid);
        }
        let origin = self.origin_mm().map(|v| I::new(v[0], v[1]));
        let axes = self.axes();
        let scale = self.scale_bounds();
        let result = std::array::from_fn::<_, 3, _>(|j| {
            (0..3).fold(I::point(0.), |sum, i| {
                sum + I::new(axes[j][i][0], axes[j][i][1]) * (input[i] - origin[i])
            }) / I::new(scale[0], scale[1])
        });
        if result.iter().any(|v| !v.valid()) {
            return Err(Error::Unsupported);
        }
        Ok(result.map(|v| [v.lo, v.hi]))
    }
}
pub struct MappedHeadMaterialQuery<'source> {
    source: &'source transport::MappedActualHead,
    occurrence: Option<&'source LiteralOccurrencePlacement>,
    local: HeadMaterialQuery<'source>,
    work: usize,
    peak: usize,
}
impl<'source> MappedHeadMaterialQuery<'source> {
    pub fn source(&self) -> &'source transport::MappedActualHead {
        self.source
    }
    pub fn occurrence(&self) -> Option<&'source LiteralOccurrencePlacement> {
        self.occurrence
    }
    pub fn local_bounds_mm(&self) -> [[f64; 2]; 3] {
        self.local.local_bounds_mm()
    }
    pub fn location(&self) -> Location {
        self.local.location()
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn query_budget_bytes(&self) -> usize {
        self.peak
    }
}
impl transport::MappedActualHead {
    pub fn classify_box(
        &self,
        input: [[f64; 2]; 3],
        options: HeadQueryOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<MappedHeadMaterialQuery<'_>> {
        self.query_placed(input, None, options, interrupt)
    }
    pub fn classify_placed_box<'source>(
        &'source self,
        input: [[f64; 2]; 3],
        occurrence: &'source LiteralOccurrencePlacement,
        options: HeadQueryOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<MappedHeadMaterialQuery<'source>> {
        self.query_placed(input, Some(occurrence), options, interrupt)
    }
    fn query_placed<'source>(
        &'source self,
        input: [[f64; 2]; 3],
        occurrence: Option<&'source LiteralOccurrencePlacement>,
        options: HeadQueryOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<MappedHeadMaterialQuery<'source>> {
        checked(input, options)?;
        let extra = self
            .receipt()
            .retained_bytes
            .checked_sub(self.local_geometry().receipt().retained_bytes)
            .and_then(|n| n.checked_add(occurrence.map_or(0, |a| a.retained_budget_bytes())))
            .and_then(|n| n.checked_add(size_of::<MappedHeadMaterialQuery<'_>>()))
            .and_then(|n| n.checked_add(QUERY_SCRATCH))
            .and_then(|n| n.checked_add(options.additional_resident_bytes))
            .ok_or(Error::Limit)?;
        let bytes = self
            .local_geometry()
            .receipt()
            .retained_bytes
            .checked_add(extra)
            .ok_or(Error::Limit)?;
        if bytes > options.query.maximum_bytes {
            return Err(Error::Limit);
        }
        let mut b = Budget {
            used: 0,
            remaining: options.query.maximum_work,
            bytes,
            peak: bytes,
            limit: options.query.maximum_bytes,
            interrupt,
        };
        b.step(256)?;
        let input = if let Some(a) = occurrence {
            a.inverse_box_mm(input)?
        } else {
            input
        };
        let input = self.actual_frame().inverse_box_mm(input)?;
        if b.remaining == 0 {
            return Err(Error::Limit);
        }
        let local = self.local_geometry().classify_local_box(
            input,
            HeadQueryOptions {
                query: QueryOptions {
                    maximum_work: b.remaining,
                    ..options.query
                },
                additional_resident_bytes: extra,
            },
            interrupt,
        )?;
        b.step(local.work_used())?;
        b.peak = b.peak.max(local.query_budget_bytes());
        b.step(0)?;
        Ok(MappedHeadMaterialQuery {
            source: self,
            occurrence,
            local,
            work: b.used,
            peak: b.peak,
        })
    }
}
#[cfg(test)]
mod tests;
