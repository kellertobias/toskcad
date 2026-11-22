//! Cached Taylor enclosures of an immutable polynomial normal field.
use super::*;
type Grid = Vec<Vec<[I; 3]>>;
#[cfg(test)]
mod tests;
pub(super) struct NormalJet {
    domain: [[f64; 2]; 2],
    center: [I; 3],
    gradient: [[I; 3]; 2],
    hessian: [[I; 3]; 3],
}
impl NormalJet {
    pub fn new(normal: &Grid, domain: [[f64; 2]; 2], work: &mut Work<'_>) -> Result<Self, Error> {
        work.charge(32768)?;
        let width = domain.map(|d| (I::point(d[1]) - I::point(d[0])).positive());
        let width = [width[0]?, width[1]?];
        let first = [
            derivative(normal, 0, width[0])?,
            derivative(normal, 1, width[1])?,
        ];
        let value = center(normal)?;
        let gradient = [center(&first[0])?, center(&first[1])?];
        let hessian = [
            hull(&derivative(&first[0], 0, width[0])?)?,
            hull(&derivative(&first[0], 1, width[1])?)?,
            hull(&derivative(&first[1], 1, width[1])?)?,
        ];
        Ok(Self {
            domain,
            center: value,
            gradient,
            hessian,
        })
    }
    pub fn enclose(&self, range: [[f64; 2]; 2], work: &mut Work<'_>) -> Result<[I; 3], Error> {
        work.charge(64)?;
        if range.into_iter().flatten().any(|v| !v.is_finite()) {
            return Err(Error::Invalid);
        }
        let mut delta = [I::point(0.); 2];
        for axis in 0..2 {
            let d = self.domain[axis];
            let center = I::point(d[0]) + (I::point(d[1]) - I::point(d[0])) * I::point(0.5);
            delta[axis] =
                I::new(range[axis][0].max(d[0]), range[axis][1].min(d[1])).checked()? - center;
        }
        // Taylor's integral remainder along the segment from the exact source
        // midpoint to x stays inside this convex source rectangle. Whole-source
        // Hessian hulls enclose every point of that segment, including both
        // mixed derivative contributions. No sampled normal is authoritative.
        let out = std::array::from_fn(|k| {
            self.center[k]
                + self.gradient[0][k] * delta[0]
                + self.gradient[1][k] * delta[1]
                + I::point(0.5) * self.hessian[0][k] * delta[0].square()
                + self.hessian[1][k] * delta[0] * delta[1]
                + I::point(0.5) * self.hessian[2][k] * delta[1].square()
        });
        for value in out {
            value.checked()?;
        }
        Ok(out)
    }
}
fn derivative(poles: &Grid, axis: usize, width: I) -> Result<Grid, Error> {
    let n = [poles.len() - 1, poles[0].len() - 1];
    (0..=n[0] - usize::from(axis == 0))
        .map(|u| {
            (0..=n[1] - usize::from(axis == 1))
                .map(|v| {
                    let next = poles[u + usize::from(axis == 0)][v + usize::from(axis == 1)];
                    let out = std::array::from_fn(|k| {
                        (next[k] - poles[u][v][k]) * I::point(n[axis] as f64) / width
                    });
                    for value in out {
                        value.checked()?;
                    }
                    Ok(out)
                })
                .collect()
        })
        .collect()
}
fn center(poles: &Grid) -> Result<[I; 3], Error> {
    let weight = |degree: usize, i: usize| {
        I::point(trims::binomial(degree, i)) * I::point(2_f64.powi(-(degree as i32)))
    };
    let anchor = poles[0][0];
    let mut out = anchor;
    for (u, row) in poles.iter().enumerate() {
        for (v, &point) in row.iter().enumerate() {
            let weight = weight(poles.len() - 1, u) * weight(row.len() - 1, v);
            for k in 0..3 {
                out[k] = out[k] + weight * (point[k] - anchor[k]);
            }
        }
    }
    for value in out {
        value.checked()?;
    }
    Ok(out)
}
fn hull(poles: &Grid) -> Result<[I; 3], Error> {
    let mut out = poles[0][0];
    for &point in poles.iter().flatten() {
        for k in 0..3 {
            point[k].checked()?;
            out[k].lo = out[k].lo.min(point[k].lo);
            out[k].hi = out[k].hi.max(point[k].hi);
        }
    }
    Ok(out)
}
