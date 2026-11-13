//! One actual trim law, evaluated as intervals or normalized Taylor jets.
use super::*;
use crate::tapered_thread::interval::J;
use crate::tapered_thread::{Scalar, Series};
pub(super) trait Value: Scalar {
    fn bound(v: I) -> Self;
    fn parameter(v: I) -> Self;
    fn value(self) -> I;
    fn checked(self) -> Result<Self, Error>;
    fn sine(self) -> Self;
}
impl Value for I {
    fn bound(v: I) -> Self {
        v
    }
    fn parameter(v: I) -> Self {
        v
    }
    fn value(self) -> I {
        self
    }
    fn checked(self) -> Result<Self, Error> {
        I::checked(self).map_err(|_| Error::Unresolved)
    }
    fn sine(self) -> Self {
        self.sin()
    }
}
impl Value for Series<5> {
    fn bound(v: I) -> Self {
        Self::interval(v)
    }
    fn parameter(v: I) -> Self {
        Self::variable(v)
    }
    fn value(self) -> I {
        self.c[0]
    }
    fn checked(self) -> Result<Self, Error> {
        Series::checked(self).map_err(|_| Error::Unresolved)
    }
    fn sine(self) -> Self {
        self.sin_cos()[0]
    }
}
impl Value for Series<6> {
    fn bound(v: I) -> Self {
        Self::interval(v)
    }
    fn parameter(v: I) -> Self {
        Self::variable(v)
    }
    fn value(self) -> I {
        self.c[0]
    }
    fn checked(self) -> Result<Self, Error> {
        Series::checked(self).map_err(|_| Error::Unresolved)
    }
    fn sine(self) -> Self {
        self.sin_cos()[0]
    }
}
impl Value for J {
    fn bound(v: I) -> Self {
        Self {
            v,
            ..Self::constant(0.)
        }
    }
    fn parameter(v: I) -> Self {
        Self::d(v)
    }
    fn value(self) -> I {
        self.v
    }
    fn checked(self) -> Result<Self, Error> {
        J::checked(self).map_err(|_| Error::Unresolved)
    }
    fn sine(self) -> Self {
        self.sin()
    }
}
fn sq<T: Value>(x: T) -> T {
    x * x
}
pub(super) fn root<T: Value>(x: T) -> Result<T, Error> {
    if x.value().lo <= 0. {
        return Err(Error::Unresolved);
    }
    x.sqrt().checked()
}
pub(super) fn wall_offset<T: Value>(k: &Corner, wall: WallKind, y: T) -> Result<T, Error> {
    let c = T::bound;
    match wall {
        WallKind::Outer => (c(k.join_s) + (y - c(k.join_y)) * c(k.cos_theta) / c(k.sin_theta))
            .checked()
            .map(|v| v / c(k.cos_beta)),
        WallKind::Inner => match k.inner {
            SourceInner::ShallowV { sine, cosine } => {
                (y * c(cosine) / c(sine) / c(k.cos_beta)).checked()
            }
            SourceInner::Rounded { radius_mm } => {
                let r = c(I::point(radius_mm));
                (sq(y) / (r + root(sq(r) - sq(y))?) / c(k.cos_beta)).checked()
            }
        },
    }
}
pub(super) fn cone_point<T: Value>(
    base: T,
    beta: T,
    slope: T,
    intercept: T,
    y: T,
    offset: T,
) -> Result<[T; 3], Error> {
    let q = beta / slope;
    let a = base + offset - q * intercept;
    let den = T::constant(1.) - sq(q);
    if a.value().lo <= 0. || den.value().lo <= 0. {
        return Err(Error::Unresolved);
    }
    let r = (a * q + root(sq(a) + den * sq(y))?) / den;
    let d = (intercept - r) / slope;
    Ok([(a + q * r).checked()?, y.checked()?, d.checked()?])
}
fn line_xy<T: Value>(r: T, a: T, b: T) -> Result<[T; 2], Error> {
    let difference = (r - a) * (r + a);
    let y = difference / (root(sq(r * b) + difference)? + a * b);
    Ok([(a + b * y).checked()?, y.checked()?])
}
pub(super) fn trim_state<T: Value>(
    atlas: &BoundaryEvents<'_>,
    index: usize,
    t: I,
) -> Result<[T; 3], Error> {
    let trim = atlas.trims[index];
    let c = T::bound;
    let k = &atlas.source.corners[0];
    let parameter =
        c(trim.parameters[0]) + T::parameter(t) * c(trim.parameters[1] - trim.parameters[0]);
    if trim.envelope == EnvelopeRole::EntryRoundover {
        let r = c(I::point(atlas.source.spec.mouth_radius_mm));
        let radius = c(atlas.source.envelope.entry) - r * parameter.sine();
        let depth = T::constant(2.) * r * sq((parameter / T::constant(2.)).sine());
        let rb = c(atlas.base) - c(atlas.beta_tangent) * depth;
        let slope = c(k.cos_theta / k.sin_theta / k.cos_beta);
        let offset = rb + c((k.join_s - k.join_y * k.cos_theta / k.sin_theta) / k.cos_beta);
        let [x, y] = line_xy(radius, offset, slope)?;
        Ok([x, y, depth.checked()?])
    } else {
        let i = if trim.envelope == EnvelopeRole::OuterFlank {
            0
        } else {
            1
        };
        cone_point(
            c(atlas.base),
            c(atlas.beta_tangent),
            c(atlas.cone_slopes[i]),
            c(atlas.cone_intercepts[i]),
            parameter,
            wall_offset(k, trim.wall, parameter)?,
        )
    }
}
#[cfg(test)]
mod tests;
