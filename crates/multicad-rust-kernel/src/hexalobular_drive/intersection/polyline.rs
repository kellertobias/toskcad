//! Whole-curve display bounds through directed Bernstein subdivision.
use super::*;
#[cfg(test)]
mod tests;
type Poles = [[I; 3]; 4];
fn norm(p: [I; 3]) -> f64 {
    p.into_iter()
        .fold(I::point(0.), |a, b| a + b.square())
        .sqrt()
        .hi
}
fn point(p: [I; 3]) -> [f64; 3] {
    p.map(I::midpoint)
}
fn split(p: Poles) -> [Poles; 2] {
    let mid = |a: [I; 3], b: [I; 3]| std::array::from_fn(|i| (a[i] + b[i]) / I::point(2.));
    let a = mid(p[0], p[1]);
    let b = mid(p[1], p[2]);
    let c = mid(p[2], p[3]);
    let d = mid(a, b);
    let e = mid(b, c);
    let m = mid(d, e);
    [[p[0], a, d, m], [m, e, c, p[3]]]
}
fn chord_error(p: Poles) -> f64 {
    let a = point(p[0]);
    let b = point(p[3]);
    p.into_iter()
        .enumerate()
        .map(|(j, control)| {
            let t = I::point(j as f64) / I::point(3.);
            norm(std::array::from_fn(|i| {
                control[i] - (I::point(a[i]) + (I::point(b[i]) - I::point(a[i])) * t)
            }))
        })
        .fold(0., f64::max)
}
fn append(
    p: Poles,
    target: f64,
    points: &mut Vec<[f64; 3]>,
    remaining: &mut usize,
    cancelled: &dyn Fn() -> bool,
) -> Result<(), Error> {
    let mut pending = vec![p];
    while let Some(p) = pending.pop() {
        if cancelled() {
            return Err(Error::Cancelled);
        }
        if pending.len() + 1 > *remaining {
            return Err(Error::Limit);
        }
        if chord_error(p) <= target {
            points.push(point(p[3]));
            *remaining -= 1;
        } else {
            let [a, b] = split(p);
            // With s splits and a accepted leaves, pending=s-a+1 and
            // remaining=limit-a. The check above bounds s by limit-1 even
            // when finite precision prevents further geometric contraction.
            pending.extend([b, a]);
        }
    }
    Ok(())
}

impl Drive {
    /// Every point of the retained finite circle/cone curve lies within
    /// error_mm of this ordered polyline. Positional source error is included.
    pub fn intersection_polyline(
        &self,
        index: usize,
        error_mm: f64,
        point_limit: usize,
    ) -> Result<Vec<[f64; 3]>, Error> {
        self.intersection_polyline_with_cancel(index, error_mm, point_limit, &|| false)
    }
    pub fn intersection_polyline_with_cancel(
        &self,
        index: usize,
        error_mm: f64,
        point_limit: usize,
        cancelled: &dyn Fn() -> bool,
    ) -> Result<Vec<[f64; 3]>, Error> {
        if !error_mm.is_finite() || error_mm <= 0. {
            return Err(Error::Invalid);
        }
        if !(2..=1_000_000).contains(&point_limit) {
            return Err(Error::Limit);
        }
        let request = (I::point(error_mm) / I::point(4.)).lo;
        let spans =
            self.intersection_carrier_with_cancel(index, request, point_limit, cancelled)?;
        carriers_polyline(&spans, error_mm, point_limit, cancelled)
    }
}

fn carriers_polyline(
    spans: &[CubicIntersectionSpan],
    error_mm: f64,
    point_limit: usize,
    cancelled: &dyn Fn() -> bool,
) -> Result<Vec<[f64; 3]>, Error> {
    let mut points = Vec::new();
    let mut remaining = point_limit;
    for span in spans {
        let poles = span.poles_mm.map(|p| p.map(I::point));
        if points.is_empty() {
            points.push(point(poles[0]));
            remaining -= 1;
        }
        append(
            poles,
            (I::point(error_mm) / I::point(2.)).lo,
            &mut points,
            &mut remaining,
            cancelled,
        )?;
    }
    Ok(points)
}

impl FiniteBinding {
    /// Approximate the retained literal master curves, preserving their shared
    /// vertices. This does not regenerate nominal circle/cone carriers.
    pub fn intersection_polyline(
        &self,
        index: usize,
        error_mm: f64,
        point_limit: usize,
    ) -> Result<Vec<[f64; 3]>, Error> {
        self.intersection_polyline_with_cancel(index, error_mm, point_limit, &|| false)
    }
    pub fn intersection_polyline_with_cancel(
        &self,
        index: usize,
        error_mm: f64,
        point_limit: usize,
        cancelled: &dyn Fn() -> bool,
    ) -> Result<Vec<[f64; 3]>, Error> {
        if index >= 12 || !error_mm.is_finite() || error_mm <= 0. {
            return Err(Error::Invalid);
        }
        if !(2..=1_000_000).contains(&point_limit) {
            return Err(Error::Limit);
        }
        carriers_polyline(&self.carriers()[index], error_mm, point_limit, cancelled)
    }
}
