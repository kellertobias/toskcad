//! Opt-in exact-law coupled coefficient candidates. No material or file admission.
use super::*;
use super::{domains::N, surfaces::RowKey};
use num_bigint::{BigInt, Sign};
use num_rational::BigRational;
use num_traits::{One, Signed, Zero};
mod admitted;
mod aliases;
mod number;
mod reduced;
mod rows;
const BITS: usize = 96;
const SCRATCH: usize = 16 * 1024 * 1024;
/// Actual dyadic candidate coefficient; the adjacent-float cache is not its
/// manufacturing authority. Original bytes and the fixed denominator are.
pub struct Coefficient {
    sign: i8,
    numerator: Vec<u8>,
    bounds: [f64; 2],
}
impl Coefficient {
    pub fn signum(&self) -> i8 {
        self.sign
    }
    pub fn numerator_magnitude_le(&self) -> &[u8] {
        &self.numerator
    }
    pub fn fractional_bits(&self) -> usize {
        BITS
    }
    pub fn numeric_bounds_mm(&self) -> [f64; 2] {
        self.bounds
    }
    pub fn decimal(
        &self,
        maximum_work: usize,
        maximum_bytes: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<String, Error> {
        let mut w = Work::new(maximum_work, interrupt)?;
        if maximum_bytes == 0 || maximum_bytes > 64 * 1024 * 1024 {
            return Err(Error::Invalid);
        }
        let mut memory = Memory::new(size_of::<Self>() + self.numerator.capacity(), maximum_bytes)?;
        self.decimal_in(&mut memory, &mut w)
    }
    pub(crate) fn decimal_in(
        &self,
        memory: &mut Memory,
        w: &mut Work<'_>,
    ) -> Result<String, Error> {
        memory.add(64 * 1024)?;
        w.charge(128)?;
        let value = number::Grid(BigInt::from_bytes_le(
            if self.sign < 0 {
                Sign::Minus
            } else if self.sign > 0 {
                Sign::Plus
            } else {
                Sign::NoSign
            },
            &self.numerator,
        ));
        let text = number::dyadic_text(&value, BITS, w)?;
        memory.add(text.capacity())?;
        memory.used -= 64 * 1024;
        Ok(text)
    }
    fn from_grid(
        value: &number::Grid,
        memory: &mut Memory,
        w: &mut Work<'_>,
    ) -> Result<Self, Error> {
        let bytes = encoded(&value.0, memory, w)?;
        let bounds = N(value.actual(0., BITS)).enclosure(w)?;
        Ok(Self {
            sign: sign(&value.0),
            numerator: bytes,
            bounds: [bounds.lo, bounds.hi],
        })
    }
}
/// Exact original-curve parameter retained independently of float interval views.
pub struct Parameter {
    sign: i8,
    numerator: Vec<u8>,
    denominator: Vec<u8>,
}
impl Parameter {
    pub fn signum(&self) -> i8 {
        self.sign
    }
    pub fn numerator_magnitude_le(&self) -> &[u8] {
        &self.numerator
    }
    pub fn denominator_magnitude_le(&self) -> &[u8] {
        &self.denominator
    }
    fn from_number(value: &N, memory: &mut Memory, w: &mut Work<'_>) -> Result<Self, Error> {
        Ok(Self {
            sign: sign(value.0.numer()),
            numerator: encoded(value.0.numer(), memory, w)?,
            denominator: encoded(value.0.denom(), memory, w)?,
        })
    }
}
fn sign(value: &BigInt) -> i8 {
    if value.is_negative() {
        -1
    } else if value.is_zero() {
        0
    } else {
        1
    }
}
fn encoded(value: &BigInt, memory: &mut Memory, w: &mut Work<'_>) -> Result<Vec<u8>, Error> {
    let requested = (value.bits() as usize).div_ceil(8).max(1);
    w.charge(1 + requested)?;
    // Source arithmetic already limits operands to16Kbits. Its caller-owned
    // scratch includes the transient digit export; persistent capacity is actual.
    memory.add(requested.checked_mul(4).ok_or(Error::Limit)?)?;
    let (_, raw) = value.to_bytes_le();
    let mut out = memory.allocate(raw.len(), w)?;
    out.extend_from_slice(&raw);
    memory.used -= requested * 4;
    Ok(out)
}
pub struct CoupledRow {
    key: RowKey,
    owners: [surfaces::VertexOwner; 2],
    parameters: [Parameter; 2],
    coefficients: Vec<[Coefficient; 3]>,
    errors: [f64; 2],
}
impl CoupledRow {
    pub fn key(&self) -> RowKey {
        self.key
    }
    pub fn endpoint_owners(&self) -> [surfaces::VertexOwner; 2] {
        self.owners
    }
    pub fn original_parameters(&self) -> [&Parameter; 2] {
        [&self.parameters[0], &self.parameters[1]]
    }
    pub fn coefficients(&self) -> &[[Coefficient; 3]] {
        &self.coefficients
    }
    pub fn position_error_mm(&self) -> f64 {
        self.errors[0]
    }
    /// Derivative in the ORIGINAL owned source curve parameter, never re-labelled
    /// as the short trim's normalized derivative or an actual-file normal proof.
    pub fn original_curve_derivative_error_mm(&self) -> f64 {
        self.errors[1]
    }
}
/// The SAME immutable authored Boundary is the correspondence authority. This
/// type retains no imported header, inferred Thread, mesh or substituted recipe.
pub struct CoupledPlan<'a> {
    source: &'a Boundary,
    domain: DomainPlan<'a>,
    rows: Vec<CoupledRow>,
    work: usize,
    bytes: usize,
}
impl<'a> CoupledPlan<'a> {
    pub fn source(&self) -> &'a Boundary {
        self.source
    }
    pub fn domain(&self) -> &DomainPlan<'a> {
        &self.domain
    }
    pub fn rows(&self) -> &[CoupledRow] {
        &self.rows
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_bytes(&self) -> usize {
        self.bytes
    }
    pub fn prepare(
        source: &'a Boundary,
        options: CurveOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, Error> {
        options.validate()?;
        let mut w = Work::new(options.work, interrupt)?;
        w.charge(1)?;
        let mut memory = Memory::new(source.retained_bytes(), options.bytes)?;
        memory.add(size_of::<Self>())?;
        let radius = source
            .source()
            .outer()
            .maximum_radius_mm()
            .max(source.source().inner().maximum_radius_mm())
            .max(1.);
        let circle_position =
            options.position_error_mm.min(options.derivative_error_mm) / (8. * radius);
        let circle_derivative = options.derivative_error_mm / (8. * radius);
        let mut circle_options = CurveOptions {
            position_error_mm: circle_position / 2.,
            derivative_error_mm: circle_derivative / 2.,
            ..options
        };
        let (angular, reduced) = loop {
            let angular =
                super::angular::fit_angular([0., 1.], false, circle_options, &mut memory, &mut w)?;
            let reduced = reduced::reduce(&angular, false, &mut memory, &mut w)?;
            if reduced.errors[0] <= circle_position && reduced.errors[1] <= circle_derivative {
                break (angular, reduced);
            }
            memory.release(&angular.knots);
            memory.release(&angular.poles);
            memory.release(&reduced.poles);
            memory.used -= reduced.poles.len() * 22 * 128;
            circle_options.position_error_mm /= 4.;
            circle_options.derivative_error_mm /= 4.;
        };
        let domain = DomainPlan::prepare_with_knots(
            source,
            &angular.knots,
            DomainOptions {
                bytes: options.bytes - (memory.used - source.retained_bytes()),
                work: w.remaining(),
                ..DomainOptions::default()
            },
            interrupt,
        )?;
        w.charge(domain.work_used())?;
        memory.add(domain.retained_bytes() - source.retained_bytes())?;
        memory.add(SCRATCH)?;
        let registry = admitted::Registry::prepare(&domain, &mut memory, &mut w)?;
        let mut cache: Vec<(RowKey, rows::Slice<'_, 'a>)> = memory.allocate(
            domain.stations(Side::Outer).len() + domain.stations(Side::Inner).len(),
            &mut w,
        )?;
        let mut rows: Vec<CoupledRow> = memory.allocate(registry.len(), &mut w)?;
        let mut aliases = aliases::Registry::prepare(source, &mut memory, &mut w)?;
        let mut curve_cache: Option<admitted::Curve<'a>> = None;
        for registered in registry.iter() {
            w.charge(1)?;
            let key = registered.key();
            let side = key.curve.owner();
            let component = source.source_boundary(side);
            if curve_cache.as_ref().is_none_or(|c| c.key() != key.curve) {
                curve_cache = Some(admitted::Curve::load(
                    &registered,
                    &mut memory,
                    &mut w,
                    interrupt,
                )?);
            }
            let request = registered
                .original(curve_cache.as_ref().expect("loaded same owner/key"), &mut w)?;
            let mut cached = None;
            for (i, (old, _)) in cache.iter().enumerate() {
                w.charge(1)?;
                if old.curve.owner() == side && old.interval == key.interval {
                    cached = Some(i);
                    break;
                }
            }
            let index = if let Some(index) = cached {
                index
            } else {
                memory.add(22 * 128)?;
                let angular =
                    rows::angular_slice(registered.angular(&angular, &reduced, &mut w)?, &mut w)?;
                cache.push((key, angular));
                cache.len() - 1
            };
            let mut actual = rows::row(&cache[index].1, &request, options, &mut w)?;
            if component.source().spec().left_hand {
                for point in &mut actual.poles {
                    point[1].0 = -point[1].0.clone();
                }
            }
            let curve = source
                .curve(key.curve)
                .expect("registered actual source slab row");
            let owners = std::array::from_fn(|i| {
                if actual.source_trim[i].0.is_zero() {
                    surfaces::VertexOwner::Original(vertex_key(side, curve.endpoints()[0]))
                } else if actual.source_trim[i].0.is_one() {
                    surfaces::VertexOwner::Original(vertex_key(side, curve.endpoints()[1]))
                } else {
                    surfaces::VertexOwner::Station {
                        curve: key.curve,
                        station: key.interval + i as u32,
                    }
                }
            });
            let parameters = [
                Parameter::from_number(&actual.source_trim[0], &mut memory, &mut w)?,
                Parameter::from_number(&actual.source_trim[1], &mut memory, &mut w)?,
            ];
            let mut coefficients = memory.allocate(12, &mut w)?;
            for point in &actual.poles {
                coefficients.push([
                    Coefficient::from_grid(&point[0], &mut memory, &mut w)?,
                    Coefficient::from_grid(&point[1], &mut memory, &mut w)?,
                    Coefficient::from_grid(&point[2], &mut memory, &mut w)?,
                ]);
            }
            let row = CoupledRow {
                key,
                owners,
                parameters,
                coefficients,
                errors: actual.errors,
            };
            if let Some(previous) = rows.last()
                && previous.key.curve == key.curve
                && previous.key.interval + 1 == key.interval
            {
                same_point(&previous.coefficients[11], &row.coefficients[0], &mut w)?;
            }
            for (i, owner) in owners.iter().enumerate() {
                if let surfaces::VertexOwner::Original(vertex) = owner {
                    aliases.bind(*vertex, &row, rows.len(), i, &rows, &mut w)?;
                }
            }
            rows.push(row);
        }
        aliases.complete(&mut w)?;
        aliases.release(&mut memory);
        drop(aliases);
        registry.release(&mut memory);
        drop(registry);
        memory.used -= cache.len() * 22 * 128;
        memory.release(&cache);
        drop(cache);
        memory.release(&angular.knots);
        memory.release(&angular.poles);
        drop(angular);
        memory.used -= reduced.poles.len() * 22 * 128;
        memory.release(&reduced.poles);
        drop(reduced);
        memory.used -= SCRATCH;
        w.charge(1)?;
        Ok(Self {
            source,
            domain,
            rows,
            work: w.used,
            bytes: memory.used,
        })
    }
}
fn same_point(a: &[Coefficient; 3], b: &[Coefficient; 3], w: &mut Work<'_>) -> Result<(), Error> {
    for axis in 0..3 {
        w.charge(1 + a[axis].numerator.len() + b[axis].numerator.len())?;
        if a[axis].sign != b[axis].sign || a[axis].numerator != b[axis].numerator {
            return Err(Error::Unresolved);
        }
    }
    Ok(())
}
#[cfg(test)]
mod tests;
