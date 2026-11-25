//! Bounded literal replay of sealed independently qualified bytes. Only entity
//! identifiers change; numeric coefficients and distinct topology owners stay.
use super::*;
use crate::step::{Entity, Writer};
use std::fmt::Write;
/// A proper orthonormal placement of a qualified local hanger.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct QualifiedPlacement {
    pub origin: [f64; 3],
    pub u: [f64; 3],
    pub v: [f64; 3],
}
impl QualifiedPlacement {
    /// Whether [`Writer::append_placed_qualified_hanger`] admits this frame.
    pub fn normal_is_proper(self) -> bool {
        self.normal().is_ok()
    }
    /// `u×v`, after proving `u` and `v` are finite unit vectors at a right
    /// angle to within a few rounding units, so the STEP operator's
    /// projected basis is the stated one and right-handed.
    fn normal(self) -> Result<[f64; 3], Error> {
        let Self { origin, u, v } = self;
        if origin.into_iter().chain(u).chain(v).any(|x| !x.is_finite()) {
            return Err(Error::Invalid);
        }
        let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        const EPSILON: f64 = 1e-12;
        if (dot(u, u) - 1.).abs() > EPSILON
            || (dot(v, v) - 1.).abs() > EPSILON
            || dot(u, v).abs() > EPSILON
        {
            return Err(Error::Unresolved);
        }
        Ok([
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ])
    }
}
pub struct QualifiedAppendReceipt {
    work: usize,
    append: usize,
    peak: usize,
    volume: [f64; 2],
    sewing: f64,
    source_position: f64,
    support_derivatives: [f64; 2],
    original_derivative: f64,
}
impl QualifiedAppendReceipt {
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn append_work_used(&self) -> usize {
        self.append
    }
    pub fn peak_accounted_bytes(&self) -> usize {
        self.peak
    }
    pub fn healed_volume_bounds_mm3(&self) -> [f64; 2] {
        self.volume
    }
    pub fn sewing_correspondence_bound_mm(&self) -> f64 {
        self.sewing
    }
    pub fn source_to_raw_position_error_mm(&self) -> f64 {
        self.source_position
    }
    pub fn source_support_derivative_error_mm(&self) -> [f64; 2] {
        self.support_derivatives
    }
    pub fn original_coedge_derivative_error_mm(&self) -> f64 {
        self.original_derivative
    }
    pub fn source_to_healed_position_error_mm(&self) -> f64 {
        (I::point(self.source_position) + I::point(self.sewing)).hi
    }
}
fn bounds(file: &QualifiedHangerStep, l: &mut Ledger<'_>) -> Result<(usize, usize), Error> {
    let text = file.as_str();
    l.charge(text.len())?;
    let start = text.find("\nDATA;\n").ok_or(Error::Invalid)? + 7;
    let mut end = start;
    for (expected, line) in (1..).zip(text[start..].split_inclusive('\n')) {
        l.charge(line.len() + 1)?;
        let (id, body) = line.split_once('=').ok_or(Error::Invalid)?;
        if id.strip_prefix('#').and_then(|v| v.parse::<usize>().ok()) != Some(expected) {
            return Err(Error::Invalid);
        }
        end = end.checked_add(line.len()).ok_or(Error::Limit)?;
        if body.starts_with("MANIFOLD_SOLID_BREP(") {
            return Ok((end, expected));
        }
    }
    Err(Error::Invalid)
}
impl Writer {
    /// Atomic ID rebasing of the exact qualified geometry records. The result
    /// reports healed support-derived mass and sewing correspondence, never raw
    /// STEP occupied mass. A failed/cancelled replay preserves prior solids.
    pub fn append_qualified_hanger(
        &mut self,
        candidate: &QualifiedHangerStep,
        options: QualifiedExportOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<QualifiedAppendReceipt, Error> {
        self.append_qualified_hanger_with(candidate, options, None, interrupt)
    }
    /// The same literal local records, held by one declared STEP map onto a
    /// proper orthonormal placement. No replayed coordinate is rewritten:
    /// the local body keeps the exact qualified bytes and the placement is
    /// one `MAPPED_ITEM` whose operator states `origin`, `u`, `v` and `u×v`.
    pub fn append_placed_qualified_hanger(
        &mut self,
        candidate: &QualifiedHangerStep,
        options: QualifiedExportOptions,
        placement: QualifiedPlacement,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<QualifiedAppendReceipt, Error> {
        self.append_qualified_hanger_with(candidate, options, Some(placement), interrupt)
    }
    fn append_qualified_hanger_with(
        &mut self,
        candidate: &QualifiedHangerStep,
        options: QualifiedExportOptions,
        placement: Option<QualifiedPlacement>,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<QualifiedAppendReceipt, Error> {
        let normal = placement.map(QualifiedPlacement::normal).transpose()?;
        if options.maximum_work == 0
            || options.maximum_work > 1_000_000_000
            || options.maximum_bytes == 0
            || options.maximum_bytes > 64 * 1024 * 1024
            || !options.maximum_healed_volume_width_mm3.is_finite()
            || options.maximum_healed_volume_width_mm3 <= 0.
        {
            return Err(Error::Invalid);
        }
        if (I::point(candidate.volume[1]) - I::point(candidate.volume[0])).hi
            > options.maximum_healed_volume_width_mm3
        {
            return Err(Error::Unresolved);
        }
        let mut l = Ledger {
            used: 0,
            left: options.maximum_work,
            bytes: options.maximum_bytes,
            interrupt,
        };
        l.charge(candidate.work_used())?;
        let inherited = l.used;
        let (end, solid) = bounds(candidate, &mut l)?;
        let start = candidate.as_str().find("\nDATA;\n").ok_or(Error::Invalid)? + 7;
        let text = &candidate.as_str()[start..end];
        let mut peak = 0;
        self.append_typed(|file| {
            let offset = file.next.checked_sub(1).ok_or(Error::Invalid)?;
            let maximum = offset.checked_add(solid).ok_or(Error::Limit)?;
            let digits = maximum.to_string().len();
            l.charge(text.len())?;
            let hashes = text.bytes().filter(|&c| c == b'#').count();
            let growth = text
                .len()
                .checked_add(hashes.checked_mul(digits).ok_or(Error::Limit)?)
                .ok_or(Error::Limit)?;
            let requested = file.data.len().checked_add(growth).ok_or(Error::Limit)?;
            let old = file.data.capacity();
            let held = candidate
                .candidate
                .retained_bytes()
                .checked_add(
                    file.points
                        .capacity()
                        .checked_mul(128)
                        .ok_or(Error::Limit)?,
                )
                .and_then(|v| v.checked_add(64 * 1024))
                .ok_or(Error::Limit)?;
            let conservative = held
                .checked_add(old)
                .and_then(|v| v.checked_add(requested))
                .ok_or(Error::Limit)?;
            if conservative > l.bytes {
                return Err(Error::Limit);
            }
            file.data
                .try_reserve_exact(growth)
                .map_err(|_| Error::Limit)?;
            peak = held + old + file.data.capacity();
            if peak > l.bytes {
                return Err(Error::Limit);
            }
            for line in text.split_inclusive('\n') {
                l.charge(line.len() + digits)?;
                let (id, body) = line.split_once('=').ok_or(Error::Invalid)?;
                let id = id
                    .strip_prefix('#')
                    .and_then(|v| v.parse::<usize>().ok())
                    .ok_or(Error::Invalid)?;
                write!(file.data, "#{}=", offset + id).map_err(|_| Error::Limit)?;
                let bytes = body.as_bytes();
                let mut cursor = 0;
                let mut copied = 0;
                let mut quoted = false;
                while cursor < bytes.len() {
                    if bytes[cursor] == b'\'' {
                        if quoted && bytes.get(cursor + 1) == Some(&b'\'') {
                            cursor += 2;
                            continue;
                        }
                        quoted = !quoted;
                        cursor += 1;
                        continue;
                    }
                    if bytes[cursor] != b'#' || quoted {
                        cursor += 1;
                        continue;
                    }
                    file.data.push_str(&body[copied..cursor]);
                    let a = cursor + 1;
                    cursor = a;
                    while bytes.get(cursor).is_some_and(u8::is_ascii_digit) {
                        cursor += 1;
                    }
                    let reference = body[a..cursor]
                        .parse::<usize>()
                        .map_err(|_| Error::Invalid)?;
                    if reference == 0 || reference > solid {
                        return Err(Error::Invalid);
                    }
                    write!(file.data, "#{}", offset + reference).map_err(|_| Error::Limit)?;
                    copied = cursor;
                }
                if quoted {
                    return Err(Error::Invalid);
                }
                file.data.push_str(&body[copied..]);
            }
            file.next = maximum.checked_add(1).ok_or(Error::Limit)?;
            let item = match (placement, normal) {
                (Some(p), Some(n)) => {
                    l.charge(4096)?;
                    file.mapped_singleton(Entity(maximum), [p.origin, p.u, p.v, n])
                }
                _ => Entity(maximum),
            };
            // Last callback/work charge occurs inside the atomic transaction.
            l.charge(1)?;
            Ok::<_, Error>(item)
        })?;
        Ok(QualifiedAppendReceipt {
            work: l.used,
            append: l.used - inherited,
            peak,
            volume: candidate.volume,
            sewing: candidate.sewing,
            source_position: candidate.source_position,
            support_derivatives: candidate.support_derivatives,
            original_derivative: candidate.original_derivative,
        })
    }
}
