//! Occupancy of the actual measured D11 carrier, without degree reduction.
//!
//! Exact rails sew one connected parameter annulus. Its distinct constant-z
//! caps, positive radius and consistently signed full (theta,z) Jacobian make
//! the cylindrical-coordinate map piecewise regular and a proper local
//! homeomorphism. At exact shared rails the two consistent coorientations
//! occupy opposite sides; differentiability across knots is not assumed.
//! Complete cap curves have degree one. The covering therefore has
//! one sheet: every physical height has a single ordered radial graph. This
//! proves embedding and cyclic ordering, not merely equal-v rail separation.
//! Adding the two actual planar star disks closes that graph around the axis.
//! The guaranteed inner cylinder is occupied. No head/cavity join, raw-exact
//! closure, mass, manufacturing mesh or STEP admission follows from this type.
use super::*;
use crate::tapered_thread::TensorCell11;

const MAX_WORK: usize = 64_000_000;
const MAX_TOTAL: usize = 64_000_000;
const MAX_CELLS: usize = 65_536;
const MAX_DEPTH: usize = 24;
// D11 split rows grow from one to twelve poles and can retain capacity16.
// 16KiB covers16×16 interval triples,16 row Vec headers and a cell header.
// At most depth+1
// DFS siblings live. Eight extra cells cover current input, both split outputs,
// de Casteljau transposes/temporaries and either first-derivative grid. No
// retained per-leaf vectors exist. The degree20×22 product's483 intervals,
// two12×12 interval weight matrices, and144 literal triples fit within two
// of the eight extra cells (7728+4608+3456 bytes). All coexistence is reserved
// before construction. Bounds are allocation accounting, not allocator RSS.
const CELL_BYTES: usize = 16384;
type Pending = [Option<(TensorCell11, usize)>; MAX_DEPTH + 2];

#[derive(Clone, Copy)]
pub struct Degree11OccupiedOptions {
    pub minimum_radius_mm: f64,
    pub minimum_height_mm: f64,
    /// New occupied-stage work only. The aggregate cap also includes both
    /// measured carrier construction and its inherited closed-head work.
    pub maximum_work: usize,
    pub maximum_total_work: usize,
    pub maximum_bytes: usize,
    /// Complete visited leaf/subdivision nodes, across every original cell.
    pub maximum_cells: usize,
    pub maximum_depth: usize,
}
impl Default for Degree11OccupiedOptions {
    fn default() -> Self {
        Self {
            minimum_radius_mm: 1e-5,
            minimum_height_mm: 1e-5,
            maximum_work: MAX_WORK,
            maximum_total_work: MAX_TOTAL,
            maximum_bytes: BYTES,
            maximum_cells: MAX_CELLS,
            maximum_depth: MAX_DEPTH,
        }
    }
}
#[derive(Clone, Copy)]
struct Proof {
    enclosure: [[f64; 2]; 3],
    radius: f64,
    jacobian: f64,
    axial: f64,
    winding: [[f64; 2]; 2],
    cells: usize,
}
/// Sealed occupied radial graph retaining the SAME defining D11 poles, face
/// ownership, original-control provenance and measured transfer records.
pub struct OccupiedDegree11Thread<'head, 'a, 'input> {
    carrier: CoherentDegree11ThreadCarrier<'head, 'a, 'input>,
    proof: Proof,
    new_work: usize,
    peak: usize,
}
impl<'head, 'a, 'input> OccupiedDegree11Thread<'head, 'a, 'input> {
    pub fn carrier(&self) -> &CoherentDegree11ThreadCarrier<'head, 'a, 'input> {
        &self.carrier
    }
    pub fn body(&self) -> &'a Body<'input> {
        self.carrier.body()
    }
    pub fn chart(&self) -> CommonMouthChart {
        self.carrier.chart()
    }
    pub fn cap_stations_mm(&self) -> [f64; 2] {
        self.carrier.cap_stations_mm()
    }
    /// A full occupied cylinder around the local axis between the cap stations.
    pub fn occupied_core_radius_mm(&self) -> f64 {
        self.proof.radius
    }
    pub fn local_enclosure_bounds_mm(&self) -> [[f64; 2]; 3] {
        self.proof.enclosure
    }
    /// Minimum certified signed cylindrical Jacobian over final leaf charts.
    /// Derivatives use each leaf's normalized square, not global STEP units.
    pub fn minimum_leaf_jacobian_mm3(&self) -> f64 {
        self.proof.jacobian
    }
    pub fn minimum_leaf_axial_derivative_mm(&self) -> f64 {
        self.proof.axial
    }
    pub fn cap_winding_angle_bounds_rad(&self) -> [[f64; 2]; 2] {
        self.proof.winding
    }
    pub fn visited_cells(&self) -> usize {
        self.proof.cells
    }
    pub fn new_work_used(&self) -> usize {
        self.new_work
    }
    pub fn total_work_used(&self) -> usize {
        self.carrier.total_work_used() + self.new_work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.carrier.retained_budget_bytes() + additional_bytes()
    }
    pub fn accounted_peak_bytes(&self) -> usize {
        self.peak
    }
}
fn additional_bytes() -> usize {
    std::mem::size_of::<OccupiedDegree11Thread<'_, '_, '_>>()
        - std::mem::size_of::<CoherentDegree11ThreadCarrier<'_, '_, '_>>()
}
fn valid(o: Degree11OccupiedOptions) -> bool {
    o.minimum_radius_mm.is_finite()
        && o.minimum_radius_mm > 0.
        && o.minimum_height_mm.is_finite()
        && o.minimum_height_mm > 0.
        && o.maximum_work > 0
        && o.maximum_work <= MAX_WORK
        && o.maximum_total_work > 0
        && o.maximum_total_work <= MAX_TOTAL
        && o.maximum_bytes > 0
        && o.maximum_bytes <= BYTES
        && o.maximum_cells > 0
        && o.maximum_cells <= MAX_CELLS
        && o.maximum_depth <= MAX_DEPTH
}
impl<'head, 'a, 'input> CoherentDegree11ThreadCarrier<'head, 'a, 'input> {
    pub fn certify_occupied(
        self,
        options: Degree11OccupiedOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> std::result::Result<OccupiedDegree11Thread<'head, 'a, 'input>, PhillipsReadError> {
        if !valid(options) {
            return Err(PhillipsReadError::Invalid);
        }
        let available = options
            .maximum_total_work
            .checked_sub(self.total_work_used())
            .ok_or(PhillipsReadError::Limit)?
            .min(options.maximum_work);
        let mut work = Work {
            remaining: available,
            interrupt,
            cancelled: Cell::new(false),
        };
        let result = (|| {
            work.step(1)?;
            let scratch = (options.maximum_depth + 8)
                .checked_mul(CELL_BYTES)
                .and_then(|n| {
                    // Fixed metadata is reserved even for a lower depth.
                    n.checked_add(std::mem::size_of::<Pending>())
                })
                .ok_or(ReadError::Limit)?;
            let peak = self
                .retained_budget_bytes()
                .checked_add(additional_bytes())
                .and_then(|n| n.checked_add(scratch))
                .ok_or(ReadError::Limit)?;
            if peak > options.maximum_bytes {
                return Err(ReadError::Limit);
            }
            let proof = certify(self.faces(), self.cap_stations_mm(), options, &mut work)?;
            Ok((proof, peak))
        })();
        let (proof, peak) = result.map_err(|e| bound_error(e, &work))?;
        Ok(OccupiedDegree11Thread {
            carrier: self,
            proof,
            new_work: available - work.remaining,
            peak,
        })
    }
}
fn positive(i: I) -> bool {
    i.valid() && i.lo.is_finite() && i.hi.is_finite() && i.lo > 0.
}
fn refusal() -> ReadError {
    ReadError::Unsupported
}
fn cross(a: [I; 2], b: [I; 2]) -> I {
    a[0] * b[1] - a[1] * b[0]
}
fn dot(a: [I; 2], b: [I; 2]) -> I {
    a[0] * b[0] + a[1] * b[1]
}
fn pi() -> I {
    I::new(
        std::f64::consts::PI.next_down(),
        std::f64::consts::PI.next_up(),
    )
}
fn choose(n: usize, k: usize) -> usize {
    (0..k.min(n - k)).fold(1, |value, i| value * (n - i) / (i + 1))
}
struct Weights {
    u: [[I; 12]; 12],
    v: [[I; 12]; 12],
}
impl Weights {
    fn new(work: &mut Work<'_>) -> Result<Self> {
        // Bounded integer binomials (<2^53), directed ratios, and arrays.
        work.step(32768)?;
        Ok(Self {
            u: std::array::from_fn(|i| {
                std::array::from_fn(|k| {
                    if k <= i {
                        return I::point(0.);
                    }
                    I::point(((k - i) * choose(11, i) * choose(11, k)) as f64)
                        / I::point(choose(20, i + k - 1) as f64)
                })
            }),
            v: std::array::from_fn(|j| {
                std::array::from_fn(|b| {
                    I::point((choose(11, j) * choose(11, b)) as f64)
                        / I::point(choose(22, j + b) as f64)
                })
            }),
        })
    }
}
fn angular_coefficients(
    points: &[[Point; 12]; 12],
    weights: &Weights,
    work: &mut Work<'_>,
) -> Result<Vec<I>> {
    work.step(21 * 23)?;
    let mut coefficients = vec![I::point(0.); 21 * 23];
    // Antisymmetry cancels the leading U terms exactly:
    // B_i B_k' - B_k B_i' = (k-i) B_i B_k/[u(1-u)].
    // Thus A_u=x*y_u-y*x_u has degree20,22. Every ordered V pair
    // remains present; no separability, analytic law or fitted recipe is used.
    for i in 0..12 {
        for k in i + 1..12 {
            for j in 0..12 {
                for b in 0..12 {
                    // <=8 directed interval scalar operations per term,
                    // including exact-pole transport, product and accumulation.
                    work.step(8)?;
                    let a = points[i][j];
                    let p = points[k][b];
                    let c = cross(
                        [I::point(a[0]), I::point(a[1])],
                        [I::point(p[0]), I::point(p[1])],
                    );
                    let index = (i + k - 1) * 23 + j + b;
                    coefficients[index] =
                        coefficients[index] + c * weights.u[i][k] * weights.v[j][b];
                }
            }
        }
    }
    Ok(coefficients)
}
fn certify(
    faces: &[Degree11ThreadFace],
    stations: [f64; 2],
    options: Degree11OccupiedOptions,
    work: &mut Work<'_>,
) -> Result<Proof> {
    if !(3..=32).contains(&faces.len())
        || !stations.iter().all(|z| z.is_finite())
        || (I::point(stations[1]) - I::point(stations[0])).lo <= options.minimum_height_mm
    {
        return Err(refusal());
    }
    let mut proof = Proof {
        enclosure: [[f64::INFINITY, f64::NEG_INFINITY]; 3],
        radius: f64::INFINITY,
        jacobian: f64::INFINITY,
        axial: f64::INFINITY,
        winding: [[0.; 2]; 2],
        cells: 0,
    };
    let forward = faces[0].forward();
    let sign = I::point(if forward { 1. } else { -1. });
    let mut angles = [I::point(0.); 2];
    let weights = Weights::new(work)?;
    work.step(MAX_DEPTH + 2)?;
    let mut pending: Pending = std::array::from_fn(|_| None);
    for (index, face) in faces.iter().enumerate() {
        let source = face.source();
        let next = &faces[(index + 1) % faces.len()];
        let [width, height] = source.dimensions();
        if face.forward() != forward
            || height != next.source().dimensions()[1]
            || face.transfer().original_rail_entities[1]
                != next.transfer().original_rail_entities[0]
        {
            return Err(refusal());
        }
        work.step(height * 3)?;
        if source.poles()[source.poles().len() - height..] != next.source().poles()[..height] {
            return Err(refusal());
        }
        for row in source.poles().chunks_exact(height) {
            work.step(4)?;
            if row[0][2] != stations[0] || row[height - 1][2] != stations[1] {
                return Err(refusal());
            }
        }
        for p in source.poles() {
            work.step(6)?;
            for (bound, &value) in proof.enclosure.iter_mut().zip(p) {
                bound[0] = bound[0].min(value);
                bound[1] = bound[1].max(value);
            }
        }
        for (end, angle) in angles.iter_mut().enumerate() {
            for u in 0..source.spans()[0] {
                // Fixed degree11 whole-curve hull, not endpoint sampling:
                // every pole must lie in one strict open half-plane. Together
                // with the independently proved signed angular Jacobian this
                // bounds the complete curve's angle increment to (0,pi).
                work.step(256)?;
                let points: [[I; 2]; 12] = std::array::from_fn(|i| {
                    let p = source.poles()[(11 * u + i) * height + end * (height - 1)];
                    [I::point(p[0]), I::point(p[1])]
                });
                let a = points[0];
                let b = points[11];
                let center = [a[0] + b[0], a[1] + b[1]];
                if points.iter().any(|&p| !positive(dot(center, p))) {
                    return Err(refusal());
                }
                let c = sign * cross(a, b);
                if !positive(c) {
                    return Err(refusal());
                }
                let increment = pi() / I::point(2.) - (dot(a, b) / c).atan();
                if !positive(increment) || increment.hi >= pi().lo {
                    return Err(refusal());
                }
                *angle = *angle + increment;
            }
        }
        for u in 0..source.spans()[0] {
            for v in 0..source.spans()[1] {
                let poles = std::array::from_fn(|i| {
                    std::array::from_fn(|j| source.poles()[(11 * u + i) * height + 11 * v + j])
                });
                let angular = angular_coefficients(&poles, &weights, work)?;
                work.step(angular.len() * 2)?;
                let minimum = angular
                    .into_iter()
                    .map(|a| (sign * a).lo)
                    .fold(f64::INFINITY, f64::min);
                work.step(1024)?;
                // RetainedTensor11 validated nonzero spans, exact144-cell
                // grids and finite |pole|<=1e7. Both indices are within those
                // immutable spans: cell construction has no fallible input.
                pending[0] = Some((source.cell([u, v]).expect("validated original D11 cell"), 0));
                let mut pending_len = 1;
                while pending_len > 0 {
                    pending_len -= 1;
                    let (cell, depth) = pending[pending_len]
                        .take()
                        .expect("nonempty bounded D11 DFS");
                    work.step(8192)?;
                    proof.cells += 1;
                    if proof.cells > options.maximum_cells {
                        return Err(ReadError::Limit);
                    }
                    // Every cell stays12×12. <=24 midpoint restrictions keep
                    // |control|<2e7;11× differences stay<5e8. Position/first
                    // derivative hulls are therefore nonempty and finite.
                    let d = cell
                        .first_derivative_bounds()
                        .expect("finite bounded D11 derivative hulls");
                    let p = d.position.map(|x| I::new(x[0], x[1]));
                    let du = d.du.map(|x| I::new(x[0], x[1]));
                    let dv = d.dv.map(|x| I::new(x[0], x[1]));
                    let xy = [p[0], p[1]];
                    let angular_u = sign * cross(xy, [du[0], du[1]]);
                    let angular_v = sign * cross(xy, [dv[0], dv[1]]);
                    let mut jacobian = angular_u * dv[2] - angular_v * du[2];
                    if minimum > 0. && positive(dv[2]) {
                        // Original-cell angular derivative becomes a leaf
                        // derivative with its actual U restriction factor.
                        // All V/axial factors below already use that leaf.
                        let domain = cell.domain();
                        let factor = I::point(domain[0][1]) - I::point(domain[0][0]);
                        let lower = I::point(minimum) * factor * I::point(dv[2].lo)
                            - I::point(angular_v.magnitude()) * I::point(du[2].magnitude());
                        if positive(lower) {
                            jacobian = lower;
                        }
                    }
                    let radius2 = p[0].square() + p[1].square();
                    let mut radius = radius2.sqrt();
                    work.step(32)?;
                    if positive(jacobian) {
                        // J=P_xy dot(N_y,-N_x), where N=P_u×P_v.
                        // Cauchy gives |P_xy| >= J_lower/|N_xy|_upper.
                        // This earns a positive core even when the Cartesian
                        // hull of a long rotating cell contains the axis.
                        let normal = [du[1] * dv[2] - du[2] * dv[1], du[2] * dv[0] - du[0] * dv[2]];
                        let length = (normal[0].square() + normal[1].square()).sqrt();
                        // Derivatives<5e8 bound Nxy<5e17 and length<1e18,
                        // far below overflow. I::square rounds its upper end
                        // upward even when powi underflows to0; sqrt likewise
                        // rounds upward. Thus the enclosing hi remains finite
                        // and strictly positive even for subnormal normals.
                        assert!(
                            length.valid() && length.hi > 0. && length.hi.is_finite(),
                            "bounded directed normal norm"
                        );
                        let lower = I::point(jacobian.lo) / I::point(length.hi);
                        if positive(lower) {
                            radius.lo = radius.lo.max(lower.lo);
                        }
                    }
                    if positive(dv[2])
                        && positive(jacobian)
                        && radius.valid()
                        && radius.lo > options.minimum_radius_mm
                    {
                        proof.radius = proof.radius.min(radius.lo);
                        proof.jacobian = proof.jacobian.min(jacobian.lo);
                        proof.axial = proof.axial.min(dv[2].lo);
                        continue;
                    }
                    if dv[2].hi <= 0.
                        || jacobian.hi <= 0.
                        || (radius.valid() && radius.hi <= options.minimum_radius_mm)
                    {
                        return Err(refusal());
                    }
                    if depth >= options.maximum_depth {
                        return Err(refusal());
                    }
                    // Heuristic only; the SAME complete D11 polynomial is
                    // restricted by directed de Casteljau in both children.
                    let extent = |d: [I; 3]| d.into_iter().map(I::magnitude).sum::<f64>();
                    let axis = usize::from(extent(dv) > extent(du));
                    work.step(16384)?;
                    // Axis is exactly0/1; depth<24<helper64; every domain is
                    // dyadic in[0,1] with width>=2^-24, so its midpoint cannot
                    // collapse. All coefficients remain finite as above.
                    let [a, b] = cell
                        .subdivide(axis)
                        .expect("bounded dyadic D11 restriction");
                    // DFS has <=depth+1 siblings, so26 fixed slots suffice.
                    pending[pending_len] = Some((b, depth + 1));
                    pending[pending_len + 1] = Some((a, depth + 1));
                    pending_len += 2;
                }
            }
        }
        debug_assert_eq!(width, 11 * source.spans()[0] + 1);
    }
    for (end, angle) in angles.into_iter().enumerate() {
        // Exact closed seams force an integer winding. Every full cap curve
        // is regular with consistent orientation, and this strict interval
        // contains only the integer multiple2*pi. No approximate2*pi equality.
        if !angle.valid() || angle.lo <= pi().hi || angle.hi >= (I::point(3.) * pi()).lo {
            return Err(refusal());
        }
        proof.winding[end] = [angle.lo, angle.hi];
    }
    Ok(proof)
}

#[cfg(test)]
mod tests;
