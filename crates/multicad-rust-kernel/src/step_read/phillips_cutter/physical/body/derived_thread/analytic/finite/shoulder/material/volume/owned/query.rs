//! Whole boxes against the actual owned periodic head/thread minus Cutter.
//! Analytic angular extrema, directed arithmetic and all intersecting real
//! phase bands enclose the whole query. No corner material vote or mesh proxy.
use super::*;
use crate::step_read::phillips_cutter::physical::body::{
    OwnedMaterialQueryError as E, OwnedMaterialQueryOptions as Options,
};
use crate::{
    head_profile::{RepairedHeadLocation, RepairedHeadQueryOptions},
    phillips_drive::Classification as C,
};
mod angle;
mod profile;
#[cfg(test)]
mod tests;
const MAX: usize = 18_008_192;
const THREAD_WORK: usize = 2_000_000;
type R<T> = std::result::Result<T, E>;
pub struct PeriodicOwnedMaterialQuery<'s> {
    source: &'s OwnedPeriodicMaterial,
    bounds: [[f64; 2]; 3],
    location: C,
    used: usize,
    remaining: usize,
}
impl PeriodicOwnedMaterialQuery<'_> {
    pub fn source(&self) -> &OwnedPeriodicMaterial {
        self.source
    }
    pub fn bounds_mm(&self) -> [[f64; 2]; 3] {
        self.bounds
    }
    pub fn location(&self) -> C {
        self.location
    }
    pub fn work_used(&self) -> usize {
        self.used
    }
    pub fn remaining_work(&self) -> usize {
        self.remaining
    }
}
fn charge(w: &mut Work<'_>, n: usize) -> R<()> {
    w.step(n).map_err(|_| {
        if w.cancelled.get() {
            E::Cancelled
        } else {
            E::Limit
        }
    })
}
fn thread_error(e: PhillipsReadError) -> E {
    use crate::tapered_thread::Error as T;
    match e {
        PhillipsReadError::Cancelled => E::Cancelled,
        PhillipsReadError::Limit => E::Thread(T::Limit),
        PhillipsReadError::Invalid => E::Thread(T::Invalid),
        _ => E::Thread(T::Unresolved),
    }
}
fn thread_location(
    source: &PeriodicThreadSource,
    caps: [[f64; 2]; 2],
    p: [I; 3],
    maximum_work: usize,
    interrupt: Option<&dyn Fn() -> bool>,
) -> R<(C, usize)> {
    let mut work = Work {
        remaining: maximum_work,
        interrupt,
        cancelled: Cell::new(false),
    };
    let result = (|| -> FResult<C> {
        work.step(512)?;
        let low = expression(caps[0]);
        let high = expression(caps[1]);
        if maximum_work > THREAD_WORK
            || p.into_iter().any(|v| !v.valid())
            || !low.valid()
            || !high.valid()
            || low.lo < 0.
            || (high - low).lo <= 0.
            || high.hi > source.length_mm()
        {
            return Err(PhillipsReadError::Invalid);
        }
        let major = source.radius_nodes_mm()[0];
        let root = source.radius_nodes_mm()[1];
        if p[2].hi < low.lo
            || p[2].lo > high.hi
            || p[..2].iter().any(|v| v.lo > major || v.hi < -major)
        {
            return Ok(C::Outside);
        }
        let radius = (p[0].square() + p[1].square()).sqrt();
        if radius.lo > major {
            return Ok(C::Outside);
        }
        // Any portion beyond this containing cylinder prevents whole-box Inside.
        // Refuse an enclosure rather than overflowing its squared-coordinate upper.
        if p[..2].iter().any(|v| v.lo < -major || v.hi > major) || !radius.valid() {
            return Ok(C::UnresolvedBoundary);
        }
        let interior = p[2].lo > low.hi && p[2].hi < high.lo;
        if interior && radius.hi < root {
            return Ok(C::Inside);
        }
        if radius.lo <= 0. {
            return Ok(C::UnresolvedBoundary);
        }
        // radius.lo>0 proves a coordinate interval is separated from zero.
        // All coordinates are within finite major<=1e7. Every corner's dominant
        // coordinate is >=sqrt(min-subnormal), so its reciprocal is finite;
        // ratio reduction and the fixed atan sum therefore cannot refuse.
        let (theta, pi) = angle::enclose([p[0], p[1]], &mut work)?
            .expect("positive radial lower bound gives a finite angle halfplane");
        let z = I::new(p[2].lo.max(low.lo), p[2].hi.min(high.hi));
        let hand = I::point(if source.left_hand() { -1. } else { 1. });
        let pitch = I::point(source.pitch_mm());
        let phase = z - hand * pitch * theta / (I::point(2.) * pi);
        let radial = profile::bounds(source, phase, &mut work)?;
        let (lower, upper) = (radial.lo, radial.hi);
        Ok(if radius.lo > upper {
            C::Outside
        } else if interior && radius.hi < lower {
            C::Inside
        } else {
            C::UnresolvedBoundary
        })
    })();
    if work.cancelled.get() {
        Err(E::Cancelled)
    } else {
        result
            .map(|c| (c, maximum_work - work.remaining))
            .map_err(thread_error)
    }
}
impl OwnedPeriodicMaterial {
    pub fn classify_local_box(
        &self,
        bounds: [[f64; 2]; 3],
        o: Options,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> R<PeriodicOwnedMaterialQuery<'_>> {
        let mut w = Work {
            remaining: o.maximum_work,
            interrupt,
            cancelled: Cell::new(false),
        };
        let result = (|| {
            charge(&mut w, 64)?;
            if o.maximum_work > MAX
                || o.maximum_head_cells == 0
                || o.maximum_head_cells > 1_000_000
                || bounds.into_iter().map(interval).any(|v| !v.valid())
            {
                return Err(E::Invalid);
            }
            let outside = bounds
                .iter()
                .zip(self.enclosure)
                .any(|(b, e)| b[1] < e[0] || b[0] > e[1]);
            let location = if outside {
                C::Outside
            } else {
                // Portions outside the earned source envelope are already Outside. Clip
                // only to keep numeric queries finite, and never upgrade a partial box to
                // Inside when that discarded portion exists.
                let extends = bounds
                    .iter()
                    .zip(self.enclosure)
                    .any(|(b, e)| b[0] < e[0] || b[1] > e[1]);
                let p = std::array::from_fn(|k| {
                    I::new(
                        bounds[k][0].max(self.enclosure[k][0]),
                        bounds[k][1].min(self.enclosure[k][1]),
                    )
                });
                let loc = self.locate(p, o.maximum_head_cells, &mut w)?;
                // The earned envelope contains both complete components, and
                // each provider's Inside is strict interior. A clipped box on
                // an envelope boundary cannot lie in either strict interior.
                assert!(
                    !extends || loc != C::Inside,
                    "source envelope boundary is not strict material interior"
                );
                loc
            };
            charge(&mut w, 1)?;
            Ok(PeriodicOwnedMaterialQuery {
                source: self,
                bounds,
                location,
                used: o.maximum_work - w.remaining,
                remaining: w.remaining,
            })
        })();
        if w.cancelled.get() {
            Err(E::Cancelled)
        } else {
            result
        }
    }
    fn locate(&self, p: [I; 3], cells: usize, w: &mut Work<'_>) -> R<C> {
        charge(w, 4096)?;
        let c = [p[0], p[1], p[2] - expression(self.mouth)];
        let cavity = self
            .cavity
            .classify_box_mm(c.map(|v| [v.lo, v.hi]))
            .map_err(E::Cavity)?;
        if cavity == C::Inside {
            return Ok(C::Outside);
        }
        let head = self.head.view();
        let hp = [p[0], p[1], p[2] + I::point(self.head_shift)];
        let h = head
            .classify_box(
                hp.map(|v| [v.lo, v.hi]),
                RepairedHeadQueryOptions {
                    maximum_work: w.remaining.min(16_000_000),
                    maximum_cells: cells,
                },
                w.interrupt,
            )
            .map_err(E::Head)?;
        charge(w, h.work_used())?;
        let union = if h.location() == RepairedHeadLocation::Inside {
            C::Inside
        } else {
            let (thread, used) = thread_location(
                &self.thread,
                self.caps,
                p,
                w.remaining.min(THREAD_WORK),
                w.interrupt,
            )?;
            charge(w, used)?;
            match (h.location(), thread) {
                (_, C::Inside) => C::Inside,
                (RepairedHeadLocation::Outside, C::Outside) => C::Outside,
                _ => C::UnresolvedBoundary,
            }
        };
        Ok(match (union, cavity) {
            (C::Outside, _) => C::Outside,
            (C::Inside, C::Outside) => C::Inside,
            _ => C::UnresolvedBoundary,
        })
    }
}
