//! Actual repaired meridian/source-ring proofs, independent of display normals.
use super::*;
mod prepared;
pub(super) use prepared::PreparedHeadAngular;
#[cfg(test)]
pub(super) use prepared::head_attempt_curve_cached;
pub(super) type P = [f64; 3];
pub(super) type B = [I; 3];
pub(super) type R<T> = std::result::Result<T, PhillipsReadError>;
pub(super) type Angular<'a> = AngularCurve<'a, 10>;
#[derive(Clone, Copy)]
pub(super) struct AngularCurve<'a, const N: usize> {
    pub(super) poles: &'a [P],
    pub(super) spans: usize,
    pub(super) reverse: bool,
    pub(super) range: [f64; 2],
}
pub(super) fn sum_norm(p: B) -> I {
    p.into_iter()
        .fold(I::point(0.), |s, v| s + I::point(v.magnitude()))
}
pub(super) fn dot(a: B, b: B) -> I {
    (0..3).fold(I::point(0.), |s, k| s + a[k] * b[k])
}
pub(super) fn cross(a: B, b: B) -> B {
    std::array::from_fn(|k| a[(k + 1) % 3] * b[(k + 2) % 3] - a[(k + 2) % 3] * b[(k + 1) % 3])
}
pub(super) fn hull(a: I, b: I) -> I {
    I::new(a.lo.min(b.lo), a.hi.max(b.hi))
}
fn curve_split<const N: usize>(mut p: [[I; 3]; N], t: I) -> ([[I; 3]; N], [[I; 3]; N]) {
    let mut left = [[I::point(0.); 3]; N];
    let mut right = left;
    left[0] = p[0];
    right[N - 1] = p[N - 1];
    let complement = I::point(1.) - t;
    for level in 1..N {
        for i in 0..N - level {
            p[i] = std::array::from_fn(|k| complement * p[i][k] + t * p[i + 1][k]);
        }
        left[level] = p[0];
        right[N - 1 - level] = p[N - 1 - level];
    }
    (left, right)
}
impl<const N: usize> AngularCurve<'_, N> {
    pub(super) fn gamma(&self, r: [f64; 2], work: &mut Work<'_>) -> R<B> {
        let mut out = [I::new(f64::INFINITY, f64::NEG_INFINITY); 3];
        for span in 0..self.spans {
            work.charge(1)?;
            // Cap::point -> mesh::Spline::new defines this trace's global
            // parameter with literal rounded binary64 knots, not exact i/n.
            // Match that same immutable source parameterization here.
            let a = span as f64 / self.spans as f64;
            let b = (span + 1) as f64 / self.spans as f64;
            let lo = r[0].max(a);
            let hi = r[1].min(b);
            if lo > hi {
                continue;
            }
            // Two de Casteljau restrictions perform at most 3*N*(N-1)
            // directed coordinate blends. Preserve the established N=10
            // charge; N=12 has a separate conservative full-operation bound.
            work.charge(if N == 10 { 1024 } else { 16 * N * N })?;
            let controls = std::array::from_fn::<_, N, _>(|j| {
                let i = (N - 1) * span + j;
                self.poles[if self.reverse {
                    self.poles.len() - 1 - i
                } else {
                    i
                }]
                .map(I::point)
            });
            let t = (I::new(lo, hi) - I::point(a)) / (I::point(b) - I::point(a));
            // Restrict the ACTUAL polynomial first; a coefficient hull
            // preserves dependency lost by natural interval Bernstein weights.
            let lower = t.lo.max(0.);
            let upper = t.hi.min(1.);
            let (left, _) = curve_split(controls, I::point(upper));
            let restricted = if lower > 0. {
                curve_split(left, I::point(lower) / I::point(upper)).1
            } else {
                left
            };
            let p: B =
                std::array::from_fn(|k| restricted.iter().map(|p| p[k]).reduce(hull).unwrap());

            for k in 0..3 {
                out[k] = hull(out[k], p[k]);
            }
        }
        if out.iter().any(|p| !p.valid()) {
            return Err(PhillipsReadError::Unsupported);
        }
        Ok(out)
    }
    pub(super) fn direction(&self, r: [f64; 2], work: &mut Work<'_>) -> R<B> {
        let g = self.gamma(r, work)?;
        unit_direction(g)
    }

    fn sag(&self, work: &mut Work<'_>) -> R<(I, B)> {
        let a = self.direction([self.range[0]; 2], work)?;
        let b = self.direction([self.range[1]; 2], work)?;
        minor_arc_sag(a, b, self.gamma(self.range, work)?)
    }
}
// Normalize a whole actual position enclosure. A broad box can reach the
// axis even when the source does not; that is a conservative refusal.
pub(super) fn unit_direction(g: B) -> R<B> {
    let norm = (g[0].square() + g[1].square()).sqrt();
    if norm.lo <= 0. || !norm.valid() {
        return Err(PhillipsReadError::Unsupported);
    }
    Ok([g[0] / norm, g[1] / norm, I::point(0.)])
}
// Caller earns the angular ordering from the immutable retained source (or
// the explicit positive-angle circle chart in the numerical oracle). The
// halfplane therefore chooses the minor branch; it is not a monotonicity test.
pub(super) fn minor_arc_sag(a: B, b: B, gamma: B) -> R<(I, B)> {
    let mid = std::array::from_fn(|k| a[k] + b[k]);
    let halfplane = dot(gamma, mid);
    if !halfplane.valid() || halfplane.lo <= 0. {
        return Err(PhillipsReadError::Unsupported);
    }
    let cosine = ((I::point(1.) + dot(a, b)) / I::point(2.)).sqrt();
    if !cosine.valid() {
        return Err(PhillipsReadError::Unsupported);
    }
    // No inverse of cosine is taken. A zero lower bound simply gives the
    // conservative sag<=1. Finite sqrt bounds make1-cosine finite as well.
    let sag = I::point(1.) - cosine;
    let direction = std::array::from_fn(|k| {
        if k == 2 {
            I::point(0.)
        } else {
            hull(a[k], b[k]) + I::new(-sag.hi, sag.hi)
        }
    });
    Ok((sag, direction))
}
fn meridian(
    h: &crate::head_profile::RepairedHeadMeridian<'_>,
    part: usize,
    r: [f64; 2],
) -> [[I; 2]; 2] {
    let b = h.certified_interval_bounds(part, r);
    [
        b.position_mm.map(|p| I::new(p[0], p[1])),
        b.derivative_mm.map(|p| I::new(p[0], p[1])),
    ]
}
fn position<const N: usize>(
    h: &crate::head_profile::RepairedHeadMeridian<'_>,
    part: usize,
    t: f64,
    a: AngularCurve<'_, N>,
    work: &mut Work<'_>,
) -> R<B> {
    let rz = meridian(h, part, [t; 2])[0];
    let d = a.direction([a.range[0]; 2], work)?;
    Ok([rz[0] * d[0], rz[0] * d[1], rz[1]])
}
fn normal(d: B, derivative: [I; 2]) -> B {
    [derivative[1] * d[0], derivative[1] * d[1], -derivative[0]]
}
pub(super) fn normal_match(points: [P; 3], source: B, angle: f64) -> bool {
    let p = points.map(|p| p.map(I::point));
    let n = cross(
        std::array::from_fn(|k| p[1][k] - p[0][k]),
        std::array::from_fn(|k| p[2][k] - p[0][k]),
    );
    let projection = dot(n, source);
    let c = cross(n, source);
    let tan = I::point(angle).sin() / I::point(angle).cos();
    let perpendicular = dot(c, c);
    let allowed = (I::point(projection.lo) * I::point(tan.lo)).square();
    n.iter()
        .chain(source.iter())
        .chain(c.iter())
        .all(|v| v.valid())
        && projection.valid()
        && perpendicular.valid()
        && allowed.valid()
        && projection.lo > 0.
        && perpendicular.hi <= allowed.lo
}
#[derive(Debug)]
pub(super) enum Split {
    U,
    T(usize),
}
// An inconclusive angular enclosure proposes a finer sampling partition.
// Resource and cancellation failures are never converted into geometric retry.
pub(super) fn angular_rejection(error: PhillipsReadError) -> R<std::result::Result<Head, Split>> {
    match error {
        PhillipsReadError::Unsupported => Ok(Err(Split::U)),
        error => Err(error),
    }
}
pub(super) struct Head {
    pub(super) mesh: crate::display::DisplayMesh,
    pub(super) top: Vec<u32>,
    pub(super) bottom: Vec<u32>,
    pub(super) distance: f64,
}
pub(super) fn head_attempt(
    source: &OwnedDerivedMaterial,
    angular: &[Angular<'_>],
    counts: &[usize],
    options: Tessellation,
    work: &mut Work<'_>,
) -> R<std::result::Result<Head, Split>> {
    head_attempt_curve(source.head(), angular, counts, options, work)
}
// The actual head is independent of the degree of the retained entry curve.
// Both callers retain their distinct sealed material/source authority.
pub(super) fn head_attempt_curve<const N: usize>(
    head: &crate::head_profile::OwnedRepairedHead,
    angular: &[AngularCurve<'_, N>],
    counts: &[usize],
    options: Tessellation,
    work: &mut Work<'_>,
) -> R<std::result::Result<Head, Split>> {
    head_attempt_curve_impl(
        head,
        prepared::AngularInput::Original(angular),
        counts,
        options,
        work,
    )
}
fn head_attempt_curve_impl<const N: usize>(
    head: &crate::head_profile::OwnedRepairedHead,
    angular: prepared::AngularInput<'_, '_, N>,
    counts: &[usize],
    options: Tessellation,
    work: &mut Work<'_>,
) -> R<std::result::Result<Head, Split>> {
    let h = &head.view();
    let n = angular.len();
    let vertices = (1 + counts.iter().sum::<usize>()) * n;
    if vertices * 2 > options.triangle_limit {
        return Err(PhillipsReadError::Limit);
    }
    work.storage(vertices * 3 * 8 + (vertices - n) * 2 * 16 + 4 * n * 4)?;
    let mut mesh = crate::display::DisplayMesh {
        vertices: Vec::with_capacity(vertices * 3),
        indices: Vec::with_capacity((vertices - n) * 6),
        face_indices: Vec::with_capacity((vertices - n) * 2),
    };
    let mut row = Vec::with_capacity(n);
    for i in 0..n {
        let p = angular.position(h, 0, 0., i, false, work)?;
        row.push((mesh.vertices.len() / 3) as u32);
        mesh.vertices.extend(p.map(I::midpoint));
    }
    let top = row.clone();
    work.storage(mesh_bytes(&mesh) + 4 * n * 4)?;
    let mut distance = 0f64;
    for (part, &count) in counts.iter().enumerate() {
        for k in 0..count {
            let ts = [k as f64 / count as f64, (k + 1) as f64 / count as f64];
            work.storage(mesh_bytes(&mesh) + (top.capacity() + row.capacity() + n) * 4)?;
            let mut next = Vec::with_capacity(n);
            work.storage(
                mesh_bytes(&mesh) + (top.capacity() + row.capacity() + next.capacity()) * 4,
            )?;
            for i in 0..n {
                let p = angular.position(h, part, ts[1], i, false, work)?;
                next.push((mesh.vertices.len() / 3) as u32);
                mesh.vertices.extend(p.map(I::midpoint));
            }
            let ends = ts.map(|t| meridian(h, part, [t; 2])[0]);
            let tmid = (ts[0] + ts[1]) / 2.;
            let center = meridian(h, part, [tmid; 2])[0];
            let field = meridian(h, part, ts);
            let width = I::point(ts[1]) - I::point(ts[0]);
            let em = (0..2).fold(I::point(0.), |s, j| {
                s + I::point((center[j] - (ends[0][j] + ends[1][j]) / I::point(2.)).magnitude())
                    + width / I::point(2.)
                        * I::point((field[1][j] - (ends[1][j] - ends[0][j]) / width).magnitude())
            });
            for i in 0..n {
                work.charge(4096)?;
                let (sag, d) = match angular.sag(i, work) {
                    Ok(v) => v,
                    Err(error) => return angular_rejection(error),
                };
                let ideal = [
                    angular.position(h, part, ts[0], i, false, work)?,
                    angular.position(h, part, ts[0], i, true, work)?,
                    angular.position(h, part, ts[1], i, true, work)?,
                    angular.position(h, part, ts[1], i, false, work)?,
                ];
                let ids = [row[i], row[(i + 1) % n], next[(i + 1) % n], next[i]];
                let points = ids.map(|id| {
                    std::array::from_fn::<_, 3, _>(|j| mesh.vertices[3 * id as usize + j])
                });
                let ec = (0..4)
                    .map(|j| {
                        sum_norm(std::array::from_fn(|k| {
                            ideal[j][k] - I::point(points[j][k])
                        }))
                        .hi
                    })
                    .fold(0., f64::max);
                let ea = I::point(ends[0][0].hi.max(ends[1][0].hi)) * sag;
                let error = (em + ea + I::point(ec)).hi;
                if error > options.deflection {
                    return Ok(Err(if em.hi > ea.hi {
                        Split::T(part)
                    } else {
                        Split::U
                    }));
                }
                let nf = normal(d, field[1]);
                for tri in [[0, 1, 2], [0, 2, 3]] {
                    if !normal_match(tri.map(|j| points[j]), nf, options.angle) {
                        let du = angular.middle(i, work)?;
                        let nu = normal(d, meridian(h, part, [tmid; 2])[1]);
                        let nt = normal(du, field[1]);
                        let wu = nu.iter().map(|i| i.hi - i.lo).sum::<f64>();
                        let wt = nt.iter().map(|i| i.hi - i.lo).sum::<f64>();
                        return Ok(Err(if wt > wu { Split::T(part) } else { Split::U }));
                    }
                    mesh.indices.extend(tri.map(|j| ids[j]));
                    mesh.face_indices.push(part as u32);
                }
                distance = distance.max(error);
            }
            row = next;
            work.storage(mesh_bytes(&mesh) + (top.capacity() + row.capacity()) * 4)?;
        }
    }
    Ok(Ok(Head {
        mesh,
        top,
        bottom: row,
        distance,
    }))
}
// The source translation is exact real addition. Its literal binary64 image
// is a new facet candidate, so certify positional and normal transfer rather
// than treating translation as a bitwise rigid operation.
pub(super) fn translation_transfer(
    original: &crate::display::DisplayMesh,
    target: &crate::display::DisplayMesh,
    offset: u32,
    z: f64,
    positional_budget: f64,
    angular_budget: f64,
    work: &mut Work<'_>,
) -> R<f64> {
    let mut displacement = 0f64;
    for (i, p) in original.vertices.as_chunks::<3>().0.iter().enumerate() {
        work.charge(1)?;
        let actual: P = std::array::from_fn(|k| target.vertices[3 * (offset as usize + i) + k]);
        let ideal = [I::point(p[0]), I::point(p[1]), I::point(p[2]) + I::point(z)];
        displacement =
            displacement.max(sum_norm(std::array::from_fn(|k| ideal[k] - I::point(actual[k]))).hi);
    }
    if displacement > positional_budget {
        return Err(PhillipsReadError::Unsupported);
    }
    for ids in original.indices.as_chunks::<3>().0 {
        work.charge(1)?;
        let before = ids.map(|id| {
            std::array::from_fn::<_, 3, _>(|k| I::point(original.vertices[3 * id as usize + k]))
        });
        let n = cross(
            std::array::from_fn(|k| before[1][k] - before[0][k]),
            std::array::from_fn(|k| before[2][k] - before[0][k]),
        );
        let after =
            ids.map(|id| std::array::from_fn(|k| target.vertices[3 * (offset + id) as usize + k]));
        if !normal_match(after, n, angular_budget) {
            return Err(PhillipsReadError::Unsupported);
        }
    }
    Ok(displacement)
}

#[cfg(test)]
mod curve_tests;
