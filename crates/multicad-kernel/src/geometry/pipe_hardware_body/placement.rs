//! Actual Native construction and subsequent occurrence factors. Original
//! coefficients are retained; no orthogonalization or rounded composition.
use super::*;
type P = [f64; 3];
type B = [I; 3];
fn dot(a: B, b: B) -> I {
    (0..3).fold(I::ZERO, |s, k| s + a[k] * b[k])
}
fn cross(a: B, b: B) -> B {
    std::array::from_fn(|k| a[(k + 1) % 3] * b[(k + 2) % 3] - a[(k + 2) % 3] * b[(k + 1) % 3])
}
fn difference(a: P, b: P) -> B {
    std::array::from_fn(|k| I::point(a[k]) - I::point(b[k]))
}
#[derive(Clone, Copy)]
struct Permutation {
    index: [usize; 3],
    negative: [bool; 3],
}
impl Permutation {
    fn new(columns: [P; 4]) -> Option<Self> {
        let mut index = [0; 3];
        let mut negative = [false; 3];
        for row in 0..3 {
            let mut found = None;
            for col in 0..3 {
                match columns[col + 1][row] {
                    1. | -1. => {
                        if found.is_some() {
                            return None;
                        }
                        found = Some(col);
                        negative[row] = columns[col + 1][row] < 0.;
                    }
                    0. => {}
                    _ => return None,
                }
            }
            index[row] = found?;
        }
        if index[0] == index[1] || index[0] == index[2] || index[1] == index[2] {
            return None;
        }
        let inversions = (0..3)
            .flat_map(|i| (i + 1..3).map(move |j| (i, j)))
            .filter(|&(i, j)| index[i] > index[j])
            .count();
        if (inversions + negative.iter().filter(|v| **v).count()) % 2 != 0 {
            return None;
        }
        Some(Self { index, negative })
    }
    fn interval(self, p: B) -> B {
        std::array::from_fn(|k| {
            if self.negative[k] {
                -p[self.index[k]]
            } else {
                p[self.index[k]]
            }
        })
    }
    fn point(self, p: P) -> P {
        std::array::from_fn(|k| {
            if self.negative[k] {
                -p[self.index[k]]
            } else {
                p[self.index[k]]
            }
        })
    }
}
// Knuth TwoSum: for finite inputs/result under IEEE nearest rounding and
// gradual underflow, sum+residual is the exact real sum. No matrix coefficient
// is replaced; exact signed permutations need only these three translations.
fn two_sum(a: f64, b: f64) -> Option<(f64, f64)> {
    let sum = a + b;
    if !sum.is_finite() {
        return None;
    }
    let virtual_b = sum - a;
    let virtual_a = sum - virtual_b;
    let residual = (a - virtual_a) + (b - virtual_b);
    residual.is_finite().then_some((sum, residual))
}
#[derive(Clone)]
struct Affine {
    columns: [P; 4],
    cofactor: [B; 3],
    determinant: I,
    sigma: [f64; 2],
    permutation: Option<Permutation>,
}
impl Affine {
    fn new(columns: [P; 4]) -> Result<Self, GeometryError> {
        if columns.iter().flatten().any(|x| !x.is_finite()) {
            return Err(failure(Cause::Invalid));
        }
        let a = [columns[1], columns[2], columns[3]].map(|p| p.map(I::point));
        if let Some(permutation) = Permutation::new(columns) {
            // Integer coefficient binding proves an EXACT proper isometry.
            // This branch is not an approximate rigid-matrix normalization.
            return Ok(Self {
                columns,
                cofactor: a,
                determinant: I::ONE,
                sigma: [1., 1.],
                permutation: Some(permutation),
            });
        }
        let mut gram = I::ZERO;
        for i in 0..3 {
            for j in 0..3 {
                gram = gram
                    + I::point(
                        (dot(a[i], a[j]) - I::point(if i == j { 1. } else { 0. })).magnitude(),
                    );
            }
        }
        let cofactor = [cross(a[1], a[2]), cross(a[2], a[0]), cross(a[0], a[1])];
        let determinant = dot(a[0], cofactor[0]);
        if !gram.valid() || gram.hi > 1e-6 || !determinant.valid() || determinant.lo <= 0. {
            return Err(failure(Cause::Unresolved));
        }
        let sigma = [
            (I::ONE - I::point(gram.hi)).sqrt().lo,
            (I::ONE + I::point(gram.hi)).sqrt().hi,
        ];
        Ok(Self {
            columns,
            cofactor,
            determinant,
            sigma,
            permutation: None,
        })
    }
    fn linear(&self, p: B) -> B {
        if let Some(permutation) = self.permutation {
            return permutation.interval(p);
        }
        std::array::from_fn(|k| {
            (0..3).fold(I::ZERO, |s, j| s + I::point(self.columns[j + 1][k]) * p[j])
        })
    }
    fn point(&self, p: B) -> B {
        let q = self.linear(p);
        std::array::from_fn(|k| q[k] + I::point(self.columns[0][k]))
    }
    fn rounded_point(&self, p: P) -> Result<(P, f64), GeometryError> {
        let mut out = [0.; 3];
        let mut error = 0f64;
        for (k, target) in out.iter_mut().enumerate() {
            let mut sum = self.columns[0][k];
            for (j, &coordinate) in p.iter().enumerate() {
                let coefficient = self.columns[j + 1][k];
                if coefficient == 0. || coordinate == 0. {
                    continue;
                }
                let product = coefficient * coordinate;
                if !product.is_finite() {
                    return Err(failure(Cause::Unresolved));
                }
                let residue = coefficient.mul_add(coordinate, -product);
                if !residue.is_finite() {
                    return Err(failure(Cause::Unresolved));
                }
                let (next, addition) =
                    two_sum(sum, product).ok_or_else(|| failure(Cause::Unresolved))?;
                // FMA is correctly rounded, so the next larger magnitude
                // encloses the exact product residual, including underflow.
                // TwoSum encloses no guessed ulp: its addition residue is exact.
                let contribution = (residue.abs().next_up() + addition.abs()).next_up();
                error = (error + contribution).next_up();
                sum = next;
            }
            *target = sum;
        }
        if !error.is_finite() {
            return Err(failure(Cause::Unresolved));
        }
        Ok((out, error))
    }
    fn inverse(&self, p: B) -> B {
        let q = std::array::from_fn(|k| p[k] - I::point(self.columns[0][k]));
        self.cofactor.map(|row| dot(row, q) / self.determinant)
    }
    fn normal(&self, n: B) -> B {
        std::array::from_fn(|k| (0..3).fold(I::ZERO, |s, j| s + self.cofactor[j][k] * n[j]))
    }
}
#[derive(Clone)]
pub(super) struct Placement {
    factors: Vec<Affine>,
    linear_columns: [B; 3],
    composed_origin: B,
    linear_permutation: Option<Permutation>,
    cofactor_columns: [B; 3],
    pub(super) determinant: I,
    pub(super) sigma: [f64; 2],
}
impl Placement {
    pub(super) fn new(
        frame: PlaneFrame,
        kind: PipeHardwareKind,
        length: f64,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        frame.checked(policy)?;
        let mut p = Self {
            factors: Vec::new(),
            composed_origin: [I::ZERO; 3],
            linear_columns: [
                [I::ONE, I::ZERO, I::ZERO],
                [I::ZERO, I::ONE, I::ZERO],
                [I::ZERO, I::ZERO, I::ONE],
            ],
            linear_permutation: Some(Permutation {
                index: [0, 1, 2],
                negative: [false; 3],
            }),
            cofactor_columns: [
                [I::ONE, I::ZERO, I::ZERO],
                [I::ZERO, I::ONE, I::ZERO],
                [I::ZERO, I::ZERO, I::ONE],
            ],
            determinant: I::ONE,
            sigma: [1., 1.],
        };
        if matches!(kind, PipeHardwareKind::HexHead { .. }) {
            p.append([[0., 0., length], [1., 0., 0.], [0., -1., 0.], [0., 0., -1.]])?;
        }
        if frame != PlaneFrame::default() {
            p.append([frame.origin, frame.u, frame.v, frame.normal()])?;
        }
        Ok(p)
    }
    fn append(&mut self, columns: [P; 4]) -> Result<(), GeometryError> {
        if self.factors.len() >= 256 {
            return Err(failure(Cause::Limit));
        }
        let a = Affine::new(columns)?;
        self.composed_origin = a.point(self.composed_origin);
        self.linear_columns = self.linear_columns.map(|p| a.linear(p));
        self.linear_permutation = if self.linear_columns.iter().flatten().all(|v| v.lo == v.hi) {
            Permutation::new([
                [0.; 3],
                self.linear_columns[0].map(|v| v.lo),
                self.linear_columns[1].map(|v| v.lo),
                self.linear_columns[2].map(|v| v.lo),
            ])
        } else {
            None
        };
        self.cofactor_columns = if self.linear_permutation.is_some() {
            self.linear_columns
        } else {
            [
                cross(self.linear_columns[1], self.linear_columns[2]),
                cross(self.linear_columns[2], self.linear_columns[0]),
                cross(self.linear_columns[0], self.linear_columns[1]),
            ]
        };
        self.determinant = self.determinant * a.determinant;
        self.sigma = [
            (I::point(self.sigma[0]) * I::point(a.sigma[0])).lo,
            (I::point(self.sigma[1]) * I::point(a.sigma[1])).hi,
        ];
        if !self.determinant.valid()
            || self.determinant.lo <= 0.
            || !self.sigma[1].is_finite()
            || self.sigma[0] <= 0.
        {
            return Err(failure(Cause::Unresolved));
        }
        self.factors.push(a);
        Ok(())
    }
    pub(super) fn transformed(
        &self,
        t: RigidTransform,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        t.validate(policy)?;
        let mut result = self.clone();
        result.append([
            t.translation_mm,
            std::array::from_fn(|k| t.rotation[k][0]),
            std::array::from_fn(|k| t.rotation[k][1]),
            std::array::from_fn(|k| t.rotation[k][2]),
        ])?;
        Ok(result)
    }
    /// Each factor as (row-major matrix, translation), in application order.
    pub(super) fn factor_rows(&self) -> Vec<([[f64; 3]; 3], [f64; 3])> {
        self.factors
            .iter()
            .map(|a| {
                let c = a.columns;
                (
                    std::array::from_fn(|i| std::array::from_fn(|j| c[j + 1][i])),
                    c[0],
                )
            })
            .collect()
    }
    /// Enclosures of the exact composed linear columns and origin image.
    pub(super) fn composed(&self) -> ([B; 3], B) {
        (self.linear_columns, self.composed_origin)
    }
    pub(super) fn is_identity(&self) -> bool {
        self.factors.is_empty()
    }
    pub(super) fn retained_bytes(&self) -> usize {
        std::mem::size_of::<Self>() + self.factors.capacity() * std::mem::size_of::<Affine>()
    }
    pub(super) fn condition(&self) -> f64 {
        (I::point(self.sigma[1]) / I::point(self.sigma[0])).hi
    }
    pub(super) fn point_bounds(&self, mut p: B) -> Result<B, GeometryError> {
        for a in &self.factors {
            p = a.point(p);
        }
        if p.iter().any(|x| !x.valid()) {
            return Err(failure(Cause::Unresolved));
        }
        Ok(p)
    }
    pub(super) fn inverse_point(&self, mut p: B) -> Result<B, GeometryError> {
        for a in self.factors.iter().rev() {
            p = a.inverse(p);
        }
        if p.iter().any(|x| !x.valid()) {
            return Err(failure(Cause::Unresolved));
        }
        Ok(p)
    }
    pub(super) fn linear_bounds(&self, p: B) -> B {
        if let Some(permutation) = self.linear_permutation {
            return permutation.interval(p);
        }
        // Constructed once from all ORIGINAL factors using directed interval
        // products. This encloses their exact linear composition; it is never
        // a rounded replacement matrix or a normalized pose.
        std::array::from_fn(|k| (0..3).fold(I::ZERO, |s, j| s + self.linear_columns[j][k] * p[j]))
    }
    pub(super) fn normal(&self, mut p: B) -> Result<B, GeometryError> {
        for a in &self.factors {
            p = a.normal(p);
        }
        let length = dot(p, p).sqrt();
        if !length.valid() || length.lo <= 0. {
            return Err(failure(Cause::Unresolved));
        }
        Ok(p.map(|x| x / length))
    }
    pub(super) fn normal_factor(&self, mut p: B) -> Result<I, GeometryError> {
        let original = dot(p, p).sqrt();
        for a in &self.factors {
            p = a.normal(p);
        }
        let mapped = dot(p, p).sqrt();
        if !original.valid() || original.lo <= 0. || !mapped.valid() || mapped.lo <= 0. {
            return Err(failure(Cause::Unresolved));
        }
        Ok(mapped / original)
    }
    pub(super) fn point(&self, p: P) -> Result<(P, f64), GeometryError> {
        if p.iter().any(|v| !v.is_finite()) {
            return Err(failure(Cause::Invalid));
        }
        if self.linear_permutation.is_none() {
            // ORIGINAL factors still define the exact composition. These
            // cached directed columns/origin enclose that product; independent
            // interval evaluation may widen it but cannot substitute a rounded
            // matrix. Publish a literal midpoint only with its complete L1
            // distance to the original composed point enclosure.
            let linear = self.linear_bounds(p.map(I::point));
            let exact: B = std::array::from_fn(|k| linear[k] + self.composed_origin[k]);
            if exact.iter().any(|v| !v.valid()) {
                return Err(failure(Cause::Unresolved));
            }
            let point = exact.map(|v| v.lo.midpoint(v.hi));
            let error = (0..3).fold(I::ZERO, |v, k| {
                v + I::point((I::point(point[k]) - exact[k]).magnitude())
            });
            if !error.valid() {
                return Err(failure(Cause::Unresolved));
            }
            return Ok((point, error.hi));
        }
        let mut point = p;
        let mut error = 0f64;
        for a in &self.factors {
            if let Some(permutation) = a.permutation {
                let linear = permutation.point(point);
                for k in 0..3 {
                    let (value, residual) = two_sum(linear[k], a.columns[0][k])
                        .ok_or_else(|| failure(Cause::Unresolved))?;
                    point[k] = value;
                    error = (I::point(error) + I::point(residual.abs())).hi;
                }
            } else {
                let (literal, rounding) = a.rounded_point(point)?;
                let total = I::point(a.sigma[1]) * I::point(error) + I::point(rounding);
                if !total.valid() {
                    return Err(failure(Cause::Unresolved));
                }
                error = total.hi;
                point = literal;
            }
        }
        // Exact permutations preserve prior Euclidean error; their L1
        // translation residual bounds Euclidean rounding too. Generic factors
        // multiply prior error by the actual norm bound before adding their
        // directed rounding. Induction preserves the original whole chain.
        if !error.is_finite() {
            return Err(failure(Cause::Unresolved));
        }
        Ok((point, error))
    }
    pub(super) fn point_work(&self) -> usize {
        if self.linear_permutation.is_none() {
            // Nine interval products, twelve affine additions, three midpoint
            // residuals plus L1 accumulation/finite checks: <40 interval/scalar
            // primitives plus bounded loads/branches, conservatively128 units.
            return 128;
        }
        // Three TwoSum translations plus L1 directed accumulation need fewer
        // than64 scalar operations/factor. Generic interval mapping retains
        // nine product/FMA/TwoSum terms and directed error accumulation need
        // <=192 logical operations PER generic factor. Original factors stay
        // separate, including long pose chains; no quota is reset or enlarged.
        self.factors
            .iter()
            .map(|a| if a.permutation.is_some() { 64 } else { 192 })
            .sum::<usize>()
            .max(64)
    }
    pub(super) fn normal_work(&self) -> usize {
        if self.linear_permutation.is_some() {
            128
        } else {
            160
        }
    }
    /// The sphere map n->C*n/|C*n| has differential norm <=cond(C).
    /// Integrating along the original short normal great-circle gives angular
    /// error <=condition*local_angle; this preserves the angular contract.
    /// For a cofactor C, cond(C)=cond(A), including every actual factor.
    pub(super) fn local_quality(&self, q: Tessellation) -> Result<Tessellation, GeometryError> {
        quality(q)?;
        if self.is_identity() {
            return Ok(q);
        }
        let local = Tessellation {
            deflection: (I::point(q.deflection) * I::point(0.99) / I::point(self.sigma[1])).lo,
            // Reserve only the independently checked world-rounding angle;
            // a broad1% allowance can unnecessarily double a coupled source
            // grid. The actual transport check consumes q/4096 below, hence
            // condition*local+reserve <=q*(1-3/4096), never a relaxed target.
            angle: (I::point(q.angle) * (I::ONE - I::ONE / I::point(1024.))
                / I::point(self.condition()))
            .lo,
            triangle_limit: q.triangle_limit,
        };
        quality(local)?;
        Ok(local)
    }
    pub(super) fn require_normal(
        &self,
        local: [P; 3],
        world: [P; 3],
        reserve: f64,
        mut charge_extra: impl FnMut(usize) -> Result<(), GeometryError>,
    ) -> Result<B, GeometryError> {
        let local_normal = cross(
            difference(local[1], local[0]),
            difference(local[2], local[0]),
        );
        // cross(Au,Av)=cofactor(A)*cross(u,v). The original factor product's
        // directed cofactor was earned once, avoiding two complete affine
        // edge maps at EVERY facet without replacing any source coefficient.
        let exact = if let Some(p) = self.linear_permutation {
            p.interval(local_normal)
        } else {
            std::array::from_fn(|k| {
                (0..3).fold(I::ZERO, |s, j| {
                    s + self.cofactor_columns[j][k] * local_normal[j]
                })
            })
        };
        let actual = cross(
            difference(world[1], world[0]),
            difference(world[2], world[0]),
        );
        let lower = exact
            .into_iter()
            .map(|v| {
                if v.lo > 0. {
                    v.lo
                } else if v.hi < 0. {
                    -v.hi
                } else {
                    0.
                }
            })
            .fold(0f64, f64::max);
        let error = (0..3).fold(I::ZERO, |s, k| {
            s + I::point((actual[k] - exact[k]).magnitude())
        });
        let threshold = I::point(lower) * I::point(reserve.min(1.) / 2.);
        // A ball of relative radius e<1 around a nonzero vector subtends at
        // most asin(e). For e<=1/2, asin(e)<=2e. This whole-interval bound
        // proves nonzero orientation and the unchanged angular reserve.
        // Cached-cofactor work uses <=56 interval operations (<=38 for an
        // exact permutation), plus bounded loads/branches:160/128 units.
        if lower > 0. && error.valid() && threshold.valid() && error.hi <= threshold.lo {
            return Ok(actual);
        }
        // Parallel magnitude error can make this bound conservative. Retain
        // the original cross/dot proof, charging BEFORE its extra arithmetic.
        charge_extra(64)?;
        let projection = dot(exact, actual);
        let perpendicular = cross(exact, actual)
            .into_iter()
            .fold(I::ZERO, |s, v| s + I::point(v.magnitude()));
        let ratio = perpendicular / projection;
        // atan(ratio)<=ratio; no libm angle value is a proof input.
        if !projection.valid() || projection.lo <= 0. || !ratio.valid() || ratio.hi > reserve {
            return Err(failure(Cause::Unresolved));
        }
        Ok(actual)
    }
    pub(super) fn area_scale(&self) -> I {
        I::new(
            (I::point(self.sigma[0]).square()).lo,
            (I::point(self.sigma[1]).square()).hi,
        )
    }
}
#[cfg(test)]
mod tests;
