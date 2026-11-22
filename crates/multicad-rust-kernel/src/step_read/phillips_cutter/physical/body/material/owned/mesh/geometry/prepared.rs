//! Invocation-local exact angular query memoization, never a new geometry.
//! Each private record retains the ACTUAL curve view used to earn its values.
use super::*;
pub(super) struct PreparedAngular<'a, const N: usize> {
    _source: AngularCurve<'a, N>,
    start: B,
    end: B,
    middle: B,
    whole: B,
    sag: I,
}
pub(super) enum AngularInput<'s, 'a, const N: usize> {
    Original(&'s [AngularCurve<'a, N>]),
    Prepared(&'s [PreparedAngular<'a, N>]),
}
impl<const N: usize> AngularInput<'_, '_, N> {
    pub(super) fn len(&self) -> usize {
        match self {
            Self::Original(a) => a.len(),
            Self::Prepared(a) => a.len(),
        }
    }
    pub(super) fn position(
        &self,
        h: &crate::head_profile::RepairedHeadMeridian<'_>,
        part: usize,
        t: f64,
        index: usize,
        end: bool,
        work: &mut Work<'_>,
    ) -> R<B> {
        match self {
            Self::Original(curves) => {
                let mut curve = curves[index];
                if end {
                    curve.range = [curve.range[1]; 2];
                }
                position(h, part, t, curve, work)
            }
            Self::Prepared(curves) => {
                work.charge(32)?;
                let rz = meridian(h, part, [t; 2])[0];
                let d = if end {
                    curves[index].end
                } else {
                    curves[index].start
                };
                Ok([rz[0] * d[0], rz[0] * d[1], rz[1]])
            }
        }
    }
    pub(super) fn sag(&self, index: usize, work: &mut Work<'_>) -> R<(I, B)> {
        match self {
            Self::Original(a) => a[index].sag(work),
            Self::Prepared(a) => {
                work.charge(16)?;
                Ok((a[index].sag, a[index].whole))
            }
        }
    }
    pub(super) fn middle(&self, index: usize, work: &mut Work<'_>) -> R<B> {
        match self {
            Self::Original(a) => {
                let a = a[index];
                a.direction([a.range[0] + (a.range[1] - a.range[0]) / 2.; 2], work)
            }
            Self::Prepared(a) => {
                work.charge(16)?;
                Ok(a[index].middle)
            }
        }
    }
}
fn prepare<'a>(
    angular: &[AngularCurve<'a, 12>],
    work: &mut Work<'_>,
) -> R<Vec<PreparedAngular<'a, 12>>> {
    let bytes = angular
        .len()
        .checked_mul(std::mem::size_of::<PreparedAngular<'_, 12>>())
        .ok_or(PhillipsReadError::Limit)?;
    work.storage(bytes)?;
    work.charge(angular.len())?;
    let mut prepared = Vec::with_capacity(angular.len());
    work.storage(prepared.capacity() * std::mem::size_of::<PreparedAngular<'_, 12>>())?;
    for &source in angular {
        // Four actual gamma queries debit themselves. Normalization/sag,
        // literal metadata/copies/validation/control fit this extra1024.
        work.charge(1024)?;
        let start = source.direction([source.range[0]; 2], work)?;
        let end = source.direction([source.range[1]; 2], work)?;
        let (sag, whole) = minor_arc_sag(start, end, source.gamma(source.range, work)?)?;
        let middle = source.direction(
            [source.range[0] + (source.range[1] - source.range[0]) / 2.; 2],
            work,
        )?;
        prepared.push(PreparedAngular {
            _source: source,
            start,
            end,
            middle,
            whole,
            sag,
        });
    }
    Ok(prepared)
}
pub(in super::super) struct PreparedHeadAngular<'a> {
    values: Vec<PreparedAngular<'a, 12>>,
}
impl<'a> PreparedHeadAngular<'a> {
    pub(in super::super) fn new(angular: &[AngularCurve<'a, 12>], work: &mut Work<'_>) -> R<Self> {
        Ok(Self {
            values: prepare(angular, work)?,
        })
    }
    pub(in super::super) fn bytes(&self) -> usize {
        self.values.capacity() * std::mem::size_of::<PreparedAngular<'_, 12>>()
    }
    pub(in super::super) fn attempt(
        &self,
        head: &crate::head_profile::OwnedRepairedHead,
        counts: &[usize],
        options: Tessellation,
        work: &mut Work<'_>,
    ) -> R<std::result::Result<Head, Split>> {
        head_attempt_curve_impl(
            head,
            AngularInput::Prepared(&self.values),
            counts,
            options,
            work,
        )
    }
}
#[cfg(test)]
pub(in super::super) fn head_attempt_curve_cached(
    head: &crate::head_profile::OwnedRepairedHead,
    angular: &[AngularCurve<'_, 12>],
    counts: &[usize],
    options: Tessellation,
    work: &mut Work<'_>,
) -> R<std::result::Result<Head, Split>> {
    let prepared = match PreparedHeadAngular::new(angular, work) {
        Ok(prepared) => prepared,
        Err(error) => return angular_rejection(error),
    };
    let bytes = prepared.bytes();
    work.held += bytes;
    // The immutable actual curve views are retained inside the sealed records;
    // no separately supplied source list can mismatch cached interval values.
    let result = prepared.attempt(head, counts, options, work);
    drop(prepared);
    work.held -= bytes;
    result
}
