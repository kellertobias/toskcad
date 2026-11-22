//! Auxiliary normal-field proposal with an EXACT degree-elevation residual.
//! The original degree-eleven geometry and its exact N21 owner stay intact.
//! Sample/interpolation values never authorize a normal bound by themselves.
use super::*;
use std::sync::Arc;
mod mvt;
const SMALL: usize = 8;
const LARGE: usize = 22;
pub struct ReducedNormal11 {
    source: Arc<NormalPolynomial11>,
    controls: [[Point; SMALL]; SMALL],
    residual: [I; 3],
}
fn basis(t: f64) -> [f64; LARGE] {
    let mut b = [0.; LARGE];
    b[0] = 1.;
    let complement = 1. - t;
    for degree in 1..LARGE {
        for i in (1..=degree).rev() {
            b[i] = t * b[i - 1] + complement * b[i];
        }
        b[0] *= complement;
    }
    b
}
fn interpolation() -> [[f64; SMALL]; SMALL] {
    // L_k(u)=prod_{m!=k}(7u-m)/prod_{m!=k}(k-m). Convert this
    // integer power polynomial to Bernstein7. Every exact intermediate
    // numerator has magnitude<2^53 (7 fixed factors, no external input).
    let mut matrix = [[0.; SMALL]; SMALL];
    for node in 0..SMALL {
        let mut powers = [0i64; SMALL];
        powers[0] = 1;
        let mut degree = 0;
        let mut denominator = 1i64;
        for other in 0..SMALL {
            if other == node {
                continue;
            }
            denominator *= node as i64 - other as i64;
            let mut next = [0i64; SMALL];
            for p in 0..=degree {
                next[p] -= other as i64 * powers[p];
                next[p + 1] += 7 * powers[p];
            }
            powers = next;
            degree += 1;
        }
        for (j, row) in matrix.iter_mut().enumerate() {
            let numerator = (0..=j)
                .map(|p| powers[p] * choose(j, p) as i64 * (105 / choose(7, p) as i64))
                .sum::<i64>();
            row[node] = numerator as f64 / (105 * denominator) as f64;
        }
    }
    matrix
}
fn elevation() -> [[I; SMALL]; LARGE] {
    std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            if i < j || i - j > 14 {
                I::point(0.)
            } else {
                I::point((choose(7, j) * choose(14, i - j)) as f64) / I::point(choose(21, i) as f64)
            }
        })
    })
}
// Keep this cache in its own bounded frame. The parent constructor's full
// array frame remains conservatively reserved even before its later stages.
// Cache11616B +local/returned rows8448B +fixed metadata fit24KiB; the
// constructor's named arrays and bounded copies fit40KiB: total<=64KiB.
#[inline(never)]
fn proposal_rows(
    source: &NormalPolynomial11,
    weights: &[[f64; LARGE]; SMALL],
) -> [[Point; LARGE]; SMALL] {
    let mut midpoint = [[[0.; 3]; LARGE]; LARGE];
    for (i, row) in midpoint.iter_mut().enumerate() {
        for (j, value) in row.iter_mut().enumerate() {
            *value = source.polynomial.controls[i][j].map(I::midpoint);
        }
    }
    let mut rows = [[[0.; 3]; LARGE]; SMALL];
    for (node, values) in rows.iter_mut().enumerate() {
        for (j, value) in values.iter_mut().enumerate() {
            for (i, row) in midpoint.iter().enumerate() {
                for k in 0..3 {
                    value[k] += weights[node][i] * row[j][k];
                }
            }
        }
    }
    rows
}
impl ReducedNormal11 {
    /// Proposal contractions, exact elevation, complete coefficient residual,
    /// integer matrix setup, all validation/copies and loop bookkeeping.
    // One actual N21 midpoint per component (1452), not once for EACH of8
    // proposal nodes (11616):10164calls/30492 scalar arithmetic operations
    // removed. New fixed cache/copy/control overhead <=8192; starting with
    // the prior131072 bound gives108772, rounded conservatively to118784.
    pub const WORK_BOUND: usize = 118_784;
    /// Fresh constructor scratch, excluding the immutable source owner.
    pub const WORKSPACE_BYTES: usize = 64 * 1024;
    /// Two bounded degree7 tensor restrictions, hull/error arithmetic and
    /// callback/validation bookkeeping. Source owners are caller-resident.
    pub const QUERY_WORK_BOUND: usize = 16384;
    pub const QUERY_WORKSPACE_BYTES: usize = 32 * 1024;
    pub fn new(source: Arc<NormalPolynomial11>) -> Result<Self, Error> {
        // The proposal is ordinary finite arithmetic only. Its accuracy is
        // NOT trusted: actual N21 minus exact elevation below is authority.
        let weights = std::array::from_fn::<_, SMALL, _>(|i| basis(i as f64 / 7.));
        let rows = proposal_rows(&source, &weights);
        let samples: [[Point; SMALL]; SMALL] = std::array::from_fn(|u| {
            std::array::from_fn(|v| {
                std::array::from_fn(|k| (0..LARGE).map(|j| weights[v][j] * rows[u][j][k]).sum())
            })
        });
        let matrix = interpolation();
        let along_u: [[Point; SMALL]; SMALL] = std::array::from_fn(|i| {
            std::array::from_fn(|v| {
                std::array::from_fn(|k| (0..SMALL).map(|u| matrix[i][u] * samples[u][v][k]).sum())
            })
        });
        let controls: [[Point; SMALL]; SMALL] = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                std::array::from_fn(|k| (0..SMALL).map(|v| matrix[j][v] * along_u[i][v][k]).sum())
            })
        });
        if controls.iter().flatten().flatten().any(|x| !x.is_finite()) {
            return Err(Error::Unresolved);
        }
        let weights = elevation();
        let elevated_u: [[[I; 3]; SMALL]; LARGE] = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                std::array::from_fn(|k| {
                    (0..SMALL).fold(I::point(0.), |s, u| {
                        s + weights[i][u] * I::point(controls[u][j][k])
                    })
                })
            })
        });
        let mut residual = [I::new(f64::INFINITY, f64::NEG_INFINITY); 3];
        for (i, row) in source.polynomial.controls.iter().enumerate() {
            for (j, p) in row.iter().enumerate() {
                for k in 0..3 {
                    let elevated = (0..SMALL)
                        .fold(I::point(0.), |s, v| s + weights[j][v] * elevated_u[i][v][k]);
                    let d = (p[k] - elevated).checked()?;
                    residual[k].lo = residual[k].lo.min(d.lo);
                    residual[k].hi = residual[k].hi.max(d.hi);
                }
            }
        }
        if residual.iter().any(|x| !x.valid()) {
            return Err(Error::Unresolved);
        }
        Ok(Self {
            source,
            controls,
            residual,
        })
    }
    pub fn source(&self) -> &NormalPolynomial11 {
        &self.source
    }
    pub fn retained_budget_bytes(&self) -> usize {
        std::mem::size_of::<Self>() + self.source.retained_budget_bytes() + 64
    }
    /// Complete coefficient residual in ROOT-unit normal coordinates.
    /// This is not a source-geometry transfer or a measured angular error.
    pub fn coefficient_residual_bounds(&self) -> [[f64; 2]; 3] {
        self.residual.map(|x| [x.lo, x.hi])
    }
    /// Guaranteed original N21 field on the entire requested ROOT domain.
    pub fn bounds(&self, domain: [[f64; 2]; 2]) -> Result<[[f64; 2]; 3], Error> {
        if !NormalPolynomial11::valid(domain) {
            return Err(Error::Invalid);
        }
        let rows = self
            .controls
            .iter()
            .map(|row| row.iter().map(|p| p.map(I::point)).collect())
            .collect();
        let restricted = restricted_grid(rows, domain)?;
        let b = coefficient_hull(&restricted);
        let result = std::array::from_fn::<_, 3, _>(|k| b[k] + self.residual[k]);
        if result.iter().any(|x| !x.valid()) {
            return Err(Error::Unresolved);
        }
        Ok(result.map(|x| [x.lo, x.hi]))
    }
}
/// A restricted auxiliary field with the SAME actual N21 owner and unchanged
/// root residual. Subdivision never multiplies by a source-area Jacobian:
/// these are root-unit normals, used only for direction enclosures.
#[derive(Clone)]
pub struct ReducedNormalCell11 {
    source: Arc<ReducedNormal11>,
    controls: [[[I; 3]; 8]; 8],
    domain: [[f64; 2]; 2],
    depth: usize,
}
impl ReducedNormal11 {
    pub const CELL_ROOT_WORK: usize = 512;
    pub const CELL_SPLIT_WORK: usize = 3072;
    pub const CELL_BOUND_WORK: usize = 1024;
    pub const CELL_CONE_WORK: usize = 8192;
    pub const CELL_PROJECTED_WORK: usize = 3072;
    pub const CELL_PREPARE_WORK: usize = PreparedNormalCones11::PREPARE_WORK;
    pub const CELL_HULL_WORK: usize = PreparedNormalCones11::HULL_WORK;
    pub const CELL_GROUP_WORK: usize = PreparedNormalCones11::GROUP_WORK;
    pub const CELL_REFINE_WORK: usize = PreparedNormalCones11::COEFFICIENT_WORK;
    pub const CELL_BYTES: usize = 4096;
    pub fn root_cell(self: &Arc<Self>) -> ReducedNormalCell11 {
        ReducedNormalCell11 {
            source: Arc::clone(self),
            controls: self.controls.map(|r| r.map(|p| p.map(I::point))),
            domain: [[0., 1.]; 2],
            depth: 0,
        }
    }
}
fn split_curve(mut p: [[I; 3]; 8]) -> ([[I; 3]; 8], [[I; 3]; 8]) {
    let mut left = [[I::point(0.); 3]; 8];
    let mut right = left;
    left[0] = p[0];
    right[7] = p[7];
    for level in 1..8 {
        for i in 0..8 - level {
            p[i] = std::array::from_fn(|k| (p[i][k] + p[i + 1][k]) * I::point(0.5));
        }
        left[level] = p[0];
        right[7 - level] = p[7 - level];
    }
    (left, right)
}
/// Fixed4KiB numerical buffer; construction is private to the SAME retained
/// cell. Every value projects coefficient+UNCHANGED root residual exactly as
/// the original per-coefficient cone predicate. No geometry admission follows.
pub(crate) struct PreparedNormalCones11 {
    values: [[I; 4]; 64],
}
fn cone_tangent(angle: f64) -> Result<f64, Error> {
    if !(angle > 0. && angle < std::f64::consts::FRAC_PI_2) {
        return Err(Error::Invalid);
    }
    let tangent = (I::point(angle).sin() / I::point(angle).cos()).checked()?;
    if tangent.lo < 0. {
        return Err(Error::Unresolved);
    }
    Ok(tangent.lo)
}
fn projected_cone(value: [I; 4], tangent: f64) -> Result<bool, Error> {
    if value[0].lo <= 0. {
        return Ok(false);
    }
    let perpendicular = value[1..]
        .iter()
        .fold(I::point(0.), |sum, v| sum + v.square())
        .checked()?;
    let allowed = (I::point(value[0].lo) * I::point(tangent))
        .square()
        .checked()?;
    Ok(perpendicular.hi <= allowed.lo)
}
fn projected_hull(values: impl Iterator<Item = [I; 4]>) -> [I; 4] {
    values.fold(
        [I::new(f64::INFINITY, f64::NEG_INFINITY); 4],
        |mut hull, value| {
            for k in 0..4 {
                hull[k].lo = hull[k].lo.min(value[k].lo);
                hull[k].hi = hull[k].hi.max(value[k].hi);
            }
            hull
        },
    )
}
impl PreparedNormalCones11 {
    pub(crate) const PREPARE_WORK: usize = 3072;
    pub(crate) const HULL_WORK: usize = 1024;
    pub(crate) const GROUP_WORK: usize = 1536;
    pub(crate) const COEFFICIENT_WORK: usize = 3072;
    pub(crate) fn hull_cone(&self, angle: f64) -> Result<bool, Error> {
        projected_cone(
            projected_hull(self.values.iter().copied()),
            cone_tangent(angle)?,
        )
    }
    /// All16 disjoint2×2 coefficient groups cover ALL64 original vectors.
    /// A group hull in a convex circular cone proves every member is in it.
    pub(crate) fn grouped_cones(&self, angle: f64) -> Result<bool, Error> {
        let tangent = cone_tangent(angle)?;
        for u in 0..4 {
            for v in 0..4 {
                let values = [
                    [2 * u, 2 * v],
                    [2 * u + 1, 2 * v],
                    [2 * u, 2 * v + 1],
                    [2 * u + 1, 2 * v + 1],
                ]
                .into_iter()
                .map(|[u, v]| self.values[u * 8 + v]);
                if !projected_cone(projected_hull(values), tangent)? {
                    return Ok(false);
                }
            }
        }
        Ok(true)
    }
    pub(crate) fn coefficient_cones(&self, angle: f64) -> Result<bool, Error> {
        let tangent = cone_tangent(angle)?;
        for value in self.values {
            if !projected_cone(value, tangent)? {
                return Ok(false);
            }
        }
        Ok(true)
    }
}
impl ReducedNormalCell11 {
    pub fn source(&self) -> &ReducedNormal11 {
        &self.source
    }
    pub fn domain(&self) -> [[f64; 2]; 2] {
        self.domain
    }
    pub fn bounds(&self) -> Result<[[f64; 2]; 3], Error> {
        let mut b = [I::new(f64::INFINITY, f64::NEG_INFINITY); 3];
        for p in self.controls.iter().flatten() {
            for k in 0..3 {
                b[k].lo = b[k].lo.min(p[k].lo);
                b[k].hi = b[k].hi.max(p[k].hi);
            }
        }
        let result = std::array::from_fn::<_, 3, _>(|k| b[k] + self.source.residual[k]);
        if result.iter().any(|x| !x.valid()) {
            return Err(Error::Unresolved);
        }
        Ok(result.map(|x| [x.lo, x.hi]))
    }
    pub(crate) fn prepare_normal_cones(
        &self,
        direction: [[f64; 2]; 3],
    ) -> Result<PreparedNormalCones11, Error> {
        if direction
            .iter()
            .any(|d| !d[0].is_finite() || !d[1].is_finite() || d[0] > d[1])
        {
            return Err(Error::Invalid);
        }
        let n = direction.map(|[a, b]| I::new(a, b));
        let mut values = [[I::point(0.); 4]; 64];
        for (target, p) in values.iter_mut().zip(self.controls.iter().flatten()) {
            let p: [I; 3] = std::array::from_fn(|k| p[k] + self.source.residual[k]);
            let projection = (0..3).fold(I::point(0.), |sum, k| sum + n[k] * p[k]);
            let cross: [I; 3] = std::array::from_fn(|k| {
                n[(k + 1) % 3] * p[(k + 2) % 3] - n[(k + 2) % 3] * p[(k + 1) % 3]
            });
            *target = [projection, cross[0], cross[1], cross[2]];
            for v in *target {
                v.checked()?;
            }
        }
        Ok(PreparedNormalCones11 { values })
    }
    /// Project the complete restricted proposal BEFORE four scalar hulls.
    /// Actual N21=Q7+R uses this SAME unchanged global residual; projecting
    /// R separately preserves the full actual field, not sampled directions.
    pub fn certifies_projected_normal_cone(
        &self,
        direction: [[f64; 2]; 3],
        angle: f64,
    ) -> Result<bool, Error> {
        if !(angle > 0. && angle < std::f64::consts::FRAC_PI_2)
            || direction
                .iter()
                .any(|d| !d[0].is_finite() || !d[1].is_finite() || d[0] > d[1])
        {
            return Err(Error::Invalid);
        }
        let n = direction.map(|[a, b]| I::new(a, b));
        let project = |p: [I; 3]| -> [I; 4] {
            let projection = (0..3).fold(I::point(0.), |sum, k| sum + n[k] * p[k]);
            let cross: [I; 3] = std::array::from_fn(|k| {
                n[(k + 1) % 3] * p[(k + 2) % 3] - n[(k + 2) % 3] * p[(k + 1) % 3]
            });
            [projection, cross[0], cross[1], cross[2]]
        };
        let residual = project(self.source.residual);
        if residual.iter().any(|v| !v.valid()) {
            return Err(Error::Unresolved);
        }
        let mut hull = [I::new(f64::INFINITY, f64::NEG_INFINITY); 4];
        for p in self.controls.iter().flatten() {
            let value = project(*p);
            for k in 0..4 {
                value[k].checked()?;
                hull[k].lo = hull[k].lo.min(value[k].lo);
                hull[k].hi = hull[k].hi.max(value[k].hi);
            }
        }
        let value: [I; 4] = std::array::from_fn(|k| hull[k] + residual[k]);
        for v in value {
            v.checked()?;
        }
        if value[0].lo <= 0. {
            return Ok(false);
        }
        let tangent = (I::point(angle).sin() / I::point(angle).cos()).checked()?;
        if tangent.lo < 0. {
            return Err(Error::Unresolved);
        }
        let perpendicular = value[1..]
            .iter()
            .fold(I::point(0.), |sum, v| sum + v.square())
            .checked()?;
        let allowed = (I::point(value[0].lo) * I::point(tangent.lo))
            .square()
            .checked()?;
        Ok(perpendicular.hi <= allowed.lo)
    }
    /// Sufficient whole-field normal cone proof. Each restricted B7 vector
    /// coefficient plus the UNCHANGED root residual box must lie in the same
    /// positive circular cone. This cone is convex; positive Bernstein weights
    /// sum to one, so the complete actual N21 field lies in it. This preserves
    /// vector correlation lost by boxing three components independently.
    pub fn certifies_normal_cone(
        &self,
        direction: [[f64; 2]; 3],
        angle: f64,
    ) -> Result<bool, Error> {
        if !(angle > 0. && angle < std::f64::consts::FRAC_PI_2)
            || direction
                .iter()
                .any(|d| !d[0].is_finite() || !d[1].is_finite() || d[0] > d[1])
        {
            return Err(Error::Invalid);
        }
        let n = direction.map(|[a, b]| I::new(a, b));
        let tangent = (I::point(angle).sin() / I::point(angle).cos()).checked()?;
        // Squaring a negative directed lower bound would turn an uncertain
        // tiny angle into a spuriously positive acceptance threshold.
        if tangent.lo < 0. {
            return Err(Error::Unresolved);
        }
        for p in self.controls.iter().flatten() {
            let p: [I; 3] = std::array::from_fn(|k| p[k] + self.source.residual[k]);
            let projection = (0..3)
                .fold(I::point(0.), |s, k| s + n[k] * p[k])
                .checked()?;
            if projection.lo <= 0. {
                return Ok(false);
            }
            let c: [I; 3] = std::array::from_fn(|k| {
                n[(k + 1) % 3] * p[(k + 2) % 3] - n[(k + 2) % 3] * p[(k + 1) % 3]
            });
            let perpendicular = c
                .into_iter()
                .fold(I::point(0.), |s, x| s + x.square())
                .checked()?;
            let allowed = (I::point(projection.lo) * I::point(tangent.lo))
                .square()
                .checked()?;
            if perpendicular.hi > allowed.lo {
                return Ok(false);
            }
        }
        Ok(true)
    }
    pub fn subdivide(&self, axis: usize) -> Result<[Self; 2], Error> {
        if axis > 1 {
            return Err(Error::Invalid);
        }
        if self.depth >= 64 {
            return Err(Error::Limit);
        }
        let [lo, hi] = self.domain[axis];
        let mid = lo + (hi - lo) * 0.5;
        if mid <= lo || mid >= hi {
            return Err(Error::Unresolved);
        }
        let mut left = [[[I::point(0.); 3]; 8]; 8];
        let mut right = left;
        for other in 0..8 {
            let p = std::array::from_fn(|i| {
                if axis == 0 {
                    self.controls[i][other]
                } else {
                    self.controls[other][i]
                }
            });
            let (a, b) = split_curve(p);
            for i in 0..8 {
                if axis == 0 {
                    left[i][other] = a[i];
                    right[i][other] = b[i];
                } else {
                    left[other][i] = a[i];
                    right[other][i] = b[i];
                }
            }
        }
        let mut a = self.domain;
        let mut b = a;
        a[axis][1] = mid;
        b[axis][0] = mid;
        Ok([
            Self {
                source: Arc::clone(&self.source),
                controls: left,
                domain: a,
                depth: self.depth + 1,
            },
            Self {
                source: Arc::clone(&self.source),
                controls: right,
                domain: b,
                depth: self.depth + 1,
            },
        ])
    }
}
#[cfg(test)]
mod tests;
