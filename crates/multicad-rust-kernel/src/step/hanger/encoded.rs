//! Shared scalar owners for a FACTORIZED encoded Hanger support field.
//! This is distinct from BodyField's independently rounded f64 tensor poles.
//! Finite trim/closed-shell and decoded material certificates remain separate.
use super::*;
use crate::{
    authored_hanger_bolt::{Boundary, Component, CurveKey, End, FaceKey},
    authored_internal_treatment::ThreadBandRole,
    step::insert::{ReducedUnitAngular, domains::N, reduced_unit::fit_reduced_unit},
};
use num_bigint::Sign;
use num_rational::BigRational;
mod cap;
pub(super) mod decimal;
pub use decimal::FiniteDecimal;
mod curve;
pub use curve::{
    CapCurvePiece, LiteralCapCurve, LiteralVertex, LiteralVertices, RestrictedCapCurve,
    RestrictionBounds, SidePole, SideSupport, SideSupportBatch, SideSupportKind,
};
mod closed;
mod finite;
pub use cap::{CapEndpoint, CapTrace};
pub use closed::{
    ClosedCoedge, ClosedCurve, ClosedEdge, ClosedFace, ClosedFaceKind, ClosedField, ClosedVertex,
    FieldVolume,
};
pub use finite::{FiniteCell, FiniteCorner, FiniteDomains, FiniteEdge, FiniteVertex};
const SCRATCH: usize = 16 * 1024 * 1024;
/// Canonical exact coefficient, owned compactly; no retained bigint arena.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct EncodedScalar {
    sign: i8,
    n: Box<[u8]>,
    d: Box<[u8]>,
}
impl EncodedScalar {
    pub fn signum(&self) -> i8 {
        self.sign
    }
    pub fn numerator_magnitude_le(&self) -> &[u8] {
        &self.n
    }
    pub fn denominator_magnitude_le(&self) -> &[u8] {
        &self.d
    }
    fn retain(q: &N, m: &mut Memory, w: &mut Work<'_>) -> Result<Self, Error> {
        w.charge(
            3 + (q.0.numer().bits() as usize).div_ceil(8)
                + (q.0.denom().bits() as usize).div_ceil(8),
        )?;
        let (sign, n) = q.0.numer().to_bytes_le();
        let (_, d) = q.0.denom().to_bytes_le();
        // Bigint->byte Vec and boxing shrink transients lie in SCRATCH;
        // compact retained payloads + their descriptor Vecs count exactly.
        m.add(n.len() + d.len())?;
        Ok(Self {
            sign: match sign {
                Sign::Minus => -1,
                Sign::NoSign => 0,
                Sign::Plus => 1,
            },
            n: n.into_boxed_slice(),
            d: d.into_boxed_slice(),
        })
    }
    fn number(&self, w: &mut Work<'_>) -> Result<N, Error> {
        N::parts(self.sign, &self.n, &self.d, w)
    }
}
pub struct EncodedMeridian {
    component: Component,
    band: ThreadBandRole,
    piece: BandPiece,
    phase: Vec<EncodedScalar>,
    radius: Vec<EncodedScalar>,
    weights: Vec<EncodedScalar>,
    bounds: SupportBounds,
}
impl EncodedMeridian {
    pub fn component(&self) -> Component {
        self.component
    }
    pub fn band(&self) -> ThreadBandRole {
        self.band
    }
    pub fn piece(&self) -> BandPiece {
        self.piece
    }
    pub fn degree(&self) -> usize {
        self.weights.len() - 1
    }
    pub fn phase_coefficients(&self) -> &[EncodedScalar] {
        &self.phase
    }
    pub fn radial_coefficients(&self) -> &[EncodedScalar] {
        &self.radius
    }
    pub fn weights(&self) -> &[EncodedScalar] {
        &self.weights
    }
    pub fn bounds(&self) -> SupportBounds {
        self.bounds
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EncodedCell {
    pub source_key: MeridianRowKey,
    pub meridian: usize,
}
pub struct EncodedField<'a> {
    source: SourceField<'a>,
    angular: ReducedUnitAngular,
    angular_jacobian: f64,
    meridians: Vec<EncodedMeridian>,
    cells: Vec<EncodedCell>,
    work: usize,
    bytes: usize,
}
pub struct EncodedControl {
    pub source_key: MeridianRowKey,
    /// Exact homogeneous [x*w,y*w,z*w,w], NOT Cartesian control points.
    pub homogeneous: [EncodedScalar; 4],
    pub work_used: usize,
}
fn binary(f: f64, w: &mut Work<'_>) -> Result<N, Error> {
    // Every finite binary64 has<=1075-bit exact numerator/denominator.
    w.charge(20 * 20)?;
    Ok(N(BigRational::from_float(f).ok_or(Error::Invalid)?))
}
fn query_options(options: CurveOptions, w: &Work<'_>) -> Result<DomainOptions, Error> {
    if w.remaining() == 0 {
        return Err(Error::Limit);
    }
    Ok(DomainOptions {
        maximum_work: w.remaining(),
        maximum_bytes: options.bytes,
    })
}
impl<'a> EncodedField<'a> {
    pub fn prepare(
        boundary: &'a Boundary,
        options: CurveOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, Error> {
        options.validate()?;
        let mut w = Work::new(options.work, interrupt)?;
        w.charge(1)?;
        let mut memory = Memory::new(boundary.retained_bytes() + size_of::<Self>(), options.bytes)?;
        // Same audited exact arithmetic contract as insert domains:<=32 live
        // <=16384-bit values + bigint normalization/byte-boxing temporaries.
        memory.add(SCRATCH)?;
        let source = SourceField::prepare(boundary, query_options(options, &w)?, interrupt)?;
        w.charge(source.work_used())?;
        let s = boundary.source();
        let metric = s.metric().spec();
        let wood = s.wood().spec();
        let scale = (p(s
            .wood()
            .maximum_radius_mm()
            .max(s.metric().maximum_radius_mm()))
            + p(wood.pitch_mm.max(metric.pitch_mm)) * p(f64::from(wood.starts.max(metric.starts)))
            + p(1.))
        .hi;
        let ao = CurveOptions {
            position_error_mm: options.position_error_mm / (4. * scale),
            derivative_error_mm: options.derivative_error_mm / (4. * scale),
            ..options
        };
        if ao.position_error_mm == 0. || ao.derivative_error_mm == 0. {
            return Err(Error::Unresolved);
        }
        // ONE positive unit angular graph. Left-hand supports reverse this
        // exact coefficient sequence, rather than fitting another graph.
        let angular = fit_reduced_unit(false, ao, &mut memory, &mut w)?;
        w.charge(64)?;
        let angular_jacobian =
            super::periodic::radial_proof(angular.position_error(), angular.derivative_error())?;
        let mut cells = memory.allocate(
            [Component::Metric, Component::Wood]
                .into_iter()
                .map(|c| {
                    let r = source.cell_periods(c);
                    (r[1] - r[0]) as usize * 5
                })
                .sum(),
            &mut w,
        )?;
        let mut meridians = memory.allocate(10, &mut w)?;
        for c in [Component::Metric, Component::Wood] {
            let starts = if c == Component::Metric {
                metric.starts
            } else {
                wood.starts
            };
            let range = source.cell_periods(c);
            for (role, piece) in [
                (ThreadBandRole::Crest, BandPiece::CrestLowPhase),
                (ThreadBandRole::Upper, BandPiece::Whole),
                (ThreadBandRole::Root, BandPiece::Whole),
                (ThreadBandRole::Lower, BandPiece::Whole),
                (ThreadBandRole::Crest, BandPiece::CrestHighPhase),
            ] {
                let index = meridians.len();
                let mut first = None;
                for curve in boundary.curves() {
                    w.charge(1)?;
                    let CurveKey::Cap {
                        component: owner,
                        end: End::Start,
                        period,
                        band: actual,
                    } = curve.key
                    else {
                        continue;
                    };
                    if owner != c || actual != role {
                        continue;
                    }
                    let face = *curve
                        .faces
                        .iter()
                        .find(|f| {
                            matches!(f, FaceKey::MetricThread { .. } | FaceKey::WoodThread { .. })
                        })
                        .ok_or(Error::Invalid)?;
                    for j in range[0]..range[1] {
                        w.charge(1)?;
                        let reference = j - i32::from(piece == BandPiece::CrestLowPhase);
                        if reference.rem_euclid(starts as i32) != period as i32 {
                            continue;
                        }
                        let key = MeridianRowKey {
                            face,
                            total_period: j,
                            piece,
                        };
                        first.get_or_insert(key);
                        cells.push(EncodedCell {
                            source_key: key,
                            meridian: index,
                        });
                    }
                }
                let key = first.ok_or(Error::Invalid)?;
                let row = source.meridian_row(
                    key.face,
                    key.piece,
                    CellSeam {
                        period: key.total_period,
                        upper_turn: false,
                    },
                    query_options(options, &w)?,
                    interrupt,
                )?;
                w.charge(row.work_used())?;
                meridians.push(profile(
                    &source,
                    (c, role, piece),
                    &row,
                    (&angular, options),
                    &mut memory,
                    &mut w,
                )?);
            }
        }
        w.charge(1)?;
        memory.used -= SCRATCH;
        Ok(Self {
            source,
            angular,
            angular_jacobian,
            meridians,
            cells,
            work: w.used,
            bytes: memory.used,
        })
    }
    pub fn source_field(&self) -> &SourceField<'a> {
        &self.source
    }
    /// Always positive RH physical-v; LH evaluation uses exact reversal.
    pub fn angular(&self) -> &ReducedUnitAngular {
        &self.angular
    }
    pub fn angular_jacobian_lower(&self) -> f64 {
        self.angular_jacobian
    }
    pub fn meridians(&self) -> &[EncodedMeridian] {
        &self.meridians
    }
    pub fn cells(&self) -> &[EncodedCell] {
        &self.cells
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_bytes(&self) -> usize {
        self.bytes
    }
    // Caller owns the shared exact-arithmetic scratch/ledger. Numeric controls
    // are constructed once from the SAME retained scalar owners; the public
    // detached query below separately accounts its compact result storage.
    fn homogeneous_control(
        &self,
        cell: usize,
        row: usize,
        piece: usize,
        column: usize,
        w: &mut Work<'_>,
    ) -> Result<(MeridianRowKey, [N; 4]), Error> {
        let key = self.cells.get(cell).ok_or(Error::Invalid)?;
        let meridian = &self.meridians[key.meridian];
        if row > meridian.degree() || piece >= self.angular.coefficients().len() || column > 10 {
            return Err(Error::Invalid);
        }
        let source = self.source.source_boundary().source();
        let (pitch, starts, left) = if meridian.component == Component::Metric {
            let s = source.metric().spec();
            (s.pitch_mm, s.starts, s.left_hand)
        } else {
            let s = source.wood().spec();
            (s.pitch_mm, s.starts, s.left_hand)
        };
        let r = meridian.radius[row].number(w)?;
        let q = meridian.phase[row].number(w)?;
        let weight = meridian.weights[row].number(w)?;
        let pitch = binary(pitch, w)?;
        let mut offset = pitch.mul(&N::integer(i64::from(key.source_key.total_period)), w)?;
        if meridian.component == Component::Metric {
            offset = offset
                .sub(&binary(source.metric().full_length_mm(), w)?, w)?
                .sub(&binary(source.spec().plain_length_mm, w)?, w)?;
        }
        let lead = pitch.mul(&N::integer(i64::from(starts)), w)?;
        let a = binary(self.angular.knots()[piece], w)?;
        let b = binary(self.angular.knots()[piece + 1], w)?;
        let v = a.add(
            &b.sub(&a, w)?
                .mul(&N::integer(column as i64), w)?
                .div(&N::integer(10), w)?,
            w,
        )?;
        let (ap, ac) = if left {
            (self.angular.coefficients().len() - 1 - piece, 10 - column)
        } else {
            (piece, column)
        };
        let xy = &self.angular.coefficients()[ap][ac];
        let x = r.mul(&N(xy[0].clone()), w)?;
        let y = r.mul(&N(xy[1].clone()), w)?;
        let z = q.add(&weight.mul(&offset.add(&lead.mul(&v, w)?, w)?, w)?, w)?;
        Ok((key.source_key, [x, y, z, weight]))
    }
    pub fn control(
        &self,
        cell: usize,
        row: usize,
        piece: usize,
        column: usize,
        options: DomainOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<EncodedControl, Error> {
        let mut w = super::domain::work(options, self.bytes, interrupt)?;
        let mut memory = Memory::new(self.bytes, options.maximum_bytes)?;
        memory.add(SCRATCH)?;
        let (source_key, [x, y, z, weight]) =
            self.homogeneous_control(cell, row, piece, column, &mut w)?;
        let output = [
            EncodedScalar::retain(&x, &mut memory, &mut w)?,
            EncodedScalar::retain(&y, &mut memory, &mut w)?,
            EncodedScalar::retain(&z, &mut memory, &mut w)?,
            EncodedScalar::retain(&weight, &mut memory, &mut w)?,
        ];
        w.charge(1)?;
        Ok(EncodedControl {
            source_key,
            homogeneous: output,
            work_used: w.used,
        })
    }
}
fn phase_row(
    source: &SourceField<'_>,
    c: Component,
    index: usize,
    w: &mut Work<'_>,
) -> Result<Vec<N>, Error> {
    let s = source.source_boundary().source();
    if c == Component::Wood {
        let band = s.wood().bands()[index];
        return Ok(vec![
            binary(band.phase_mm[0], w)?,
            binary(band.phase_mm[1], w)?,
        ]);
    }
    use crate::authored_metric_thread::{Knot, phase_coefficients};
    let p = binary(s.metric().spec().pitch_mm, w)?;
    let crest = binary(s.metric().spec().crest_width_mm, w)?;
    let knots = [
        Knot::Begin,
        Knot::CrestEnd,
        Knot::UpperRoot,
        Knot::LowerRoot,
        Knot::CrestStart,
        Knot::End,
    ];
    let mut at = |i: usize| -> Result<N, Error> {
        let [(pn, pd), (cn, cd)] = phase_coefficients(knots[i]);
        p.mul(&N::integer(pn), w)?
            .div(&N::integer(pd), w)?
            .add(&crest.mul(&N::integer(cn), w)?.div(&N::integer(cd), w)?, w)
    };
    let first = at(index)?;
    let last = at(index + 1)?;
    if index == 2 {
        Ok(vec![first, p.div(&N::integer(4), w)?, last])
    } else {
        Ok(vec![first, last])
    }
}
fn profile(
    source: &SourceField<'_>,
    key: (Component, ThreadBandRole, BandPiece),
    row: &MeridianRow,
    target: (&ReducedUnitAngular, CurveOptions),
    m: &mut Memory,
    w: &mut Work<'_>,
) -> Result<EncodedMeridian, Error> {
    let (component, band, piece) = key;
    let (angular, options) = target;
    let index = match band {
        ThreadBandRole::Upper => 1,
        ThreadBandRole::Root => 2,
        ThreadBandRole::Lower => 3,
        ThreadBandRole::Crest => {
            if piece == BandPiece::CrestLowPhase {
                0
            } else {
                4
            }
        }
    };
    let q = phase_row(source, component, index, w)?;
    let degree = row.degree();
    let root = component == Component::Metric && band == ThreadBandRole::Root;
    // SAME exact [1,1] or [1,1/2,1] constants as the sealed source's
    // RationalMeridian. q uses its actual phase_coefficients and p/4 middle;
    // weight/phase are NOT independently fitted approximation inputs.
    let minimum_weight = p(if root { 0.5 } else { 1. });
    let dw = p(if root { 1. } else { 0. });
    let delta = row
        .controls()
        .iter()
        .map(|c| hull(c.source_bounds[0]) - p(c.literal[0]))
        .collect::<Vec<_>>();
    let er = p(delta.iter().map(|x| x.magnitude()).fold(0_f64, f64::max));
    let der = p(delta
        .windows(2)
        .map(|x| (p(degree as f64) * (x[1] - x[0])).magnitude())
        .fold(0_f64, f64::max));
    let radial_error = er / minimum_weight;
    let radial_du_error = der / minimum_weight + er * dw / minimum_weight.square();
    let s = source.source_boundary().source();
    let (radius, rmin, pitch, starts) = if component == Component::Metric {
        let x = s.metric().spec();
        (
            s.metric().maximum_radius_mm(),
            s.metric().minimum_radius_mm()[0],
            x.pitch_mm,
            x.starts,
        )
    } else {
        let x = s.wood().spec();
        (
            s.wood().maximum_radius_mm(),
            s.wood().minimum_radius_mm(),
            x.pitch_mm,
            x.starts,
        )
    };
    w.charge(512)?;
    let (dq, radial_du) = if root {
        let b = s.metric().root_meridian().whole_bounds();
        (
            hull(b.first_derivative_mm[0]),
            p(hull(b.first_derivative_mm[1]).magnitude()),
        )
    } else {
        let mut dr = 0_f64;
        for pair in row.controls().windows(2) {
            dr = dr
                .max((hull(pair[1].source_bounds[0]) - hull(pair[0].source_bounds[0])).magnitude());
        }
        (q[degree].sub(&q[0], w)?.enclosure(w)?, p(dr))
    };
    let tau = I::new(
        std::f64::consts::TAU.next_down(),
        std::f64::consts::TAU.next_up(),
    );
    let proxy_radius = p(radius) + radial_error;
    let ep = p(angular.position_error());
    let ed = p(angular.derivative_error());
    let position = radial_error + proxy_radius * ep;
    let du = radial_du_error + (radial_du + radial_du_error) * ep;
    let dv = tau * radial_error + proxy_radius * ed;
    let normal_lower = tau * p(rmin) * dq;
    let pu = radial_du + p(dq.magnitude());
    let pv = tau * p(radius) + p(pitch) * p(f64::from(starts));
    let cross_error = du * pv + dv * pu + du * dv;
    let normal = (p(2.) * cross_error / normal_lower).hi;
    let bounds = SupportBounds {
        position_error_mm: position.hi,
        derivative_error_mm: [du.hi, dv.hi],
        normal_angle_rad: normal,
    };
    if (p(rmin) - radial_error).lo <= 0.
        || dq.lo <= 0.
        || (normal_lower - p(2.) * cross_error).lo <= 0.
        || !normal.is_finite()
        || !position.valid()
        || position.hi > options.position_error_mm
        || !du.valid()
        || du.hi > options.derivative_error_mm
        || !dv.valid()
        || dv.hi > options.derivative_error_mm
    {
        return Err(Error::Unresolved);
    }
    let mut phase = m.allocate(q.len(), w)?;
    let mut radius = m.allocate(q.len(), w)?;
    let mut weights = m.allocate(q.len(), w)?;
    for (i, c) in row.controls().iter().enumerate() {
        phase.push(EncodedScalar::retain(&q[i], m, w)?);
        radius.push(EncodedScalar::retain(&binary(c.literal[0], w)?, m, w)?);
        weights.push(EncodedScalar::retain(&binary(c.literal[3], w)?, m, w)?);
    }
    Ok(EncodedMeridian {
        component,
        band,
        piece,
        phase,
        radius,
        weights,
        bounds,
    })
}
#[cfg(test)]
mod tests;
