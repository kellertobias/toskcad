//! Atomic literal replay. Only entity numbering changes; numerical coefficients
//! and all distinct vertex/edge/face owners retain their sealed identities.
use super::*;
use std::fmt::Write;
#[derive(Clone, Copy)]
pub struct AppendOptions {
    pub maximum_work: usize,
    pub maximum_bytes: usize,
}
impl Default for AppendOptions {
    fn default() -> Self {
        Self {
            maximum_work: 64_000_000,
            maximum_bytes: 64 * 1024 * 1024,
        }
    }
}
pub struct QualifiedAppendReceipt {
    work: usize,
    append: usize,
    peak: usize,
    solid: u32,
    original_volume: [f64; 2],
    healed_volume: [f64; 2],
    position: f64,
    derivative: f64,
    sewing: f64,
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
    pub fn solid_entity(&self) -> u32 {
        self.solid
    }
    pub fn original_source_volume_bounds_mm3(&self) -> [f64; 2] {
        self.original_volume
    }
    pub fn healed_volume_bounds_mm3(&self) -> [f64; 2] {
        self.healed_volume
    }
    pub fn source_to_healed_position_error_mm(&self) -> f64 {
        self.position
    }
    pub fn source_to_healed_normalized_derivative_error_mm(&self) -> f64 {
        self.derivative
    }
    pub fn raw_support_sewing_error_mm(&self) -> f64 {
        self.sewing
    }
}
fn bounds(
    candidate: &QualifiedOrthogonalSeatsStep,
    w: &mut Ledger<'_>,
) -> Result<(usize, usize), QualificationError> {
    let text = candidate.as_str();
    w.charge(text.len())?;
    let start = text.find("\nDATA;\n").ok_or(QualificationError::Invalid)? + 7;
    let mut end = start;
    for (expected, line) in (1..).zip(text[start..].split_inclusive('\n')) {
        w.charge(line.len() + 1)?;
        let (id, body) = line.split_once('=').ok_or(QualificationError::Invalid)?;
        if id.strip_prefix('#').and_then(|id| id.parse::<usize>().ok()) != Some(expected) {
            return Err(QualificationError::Invalid);
        }
        end = end
            .checked_add(line.len())
            .ok_or(QualificationError::Limit)?;
        if body.starts_with("MANIFOLD_SOLID_BREP(") {
            return Ok((end, expected));
        }
    }
    Err(QualificationError::Invalid)
}
impl Writer {
    /// Append only independently qualified bytes. Budgets include the one
    /// immutable qualification prerequisite; append_work_used reports the new
    /// encoding work separately for a caller's bounded compound-file ledger.
    pub fn append_qualified_orthogonal_seats(
        &mut self,
        candidate: &QualifiedOrthogonalSeatsStep,
        options: AppendOptions,
        interrupt: &dyn Fn() -> bool,
    ) -> Result<QualifiedAppendReceipt, QualificationError> {
        let mut w = Ledger {
            used: 0,
            left: options.maximum_work,
            maximum_bytes: options.maximum_bytes,
            peak: 0,
            interrupt,
        };
        w.charge(0)?;
        if options.maximum_work == 0
            || options.maximum_work > 256_000_000
            || options.maximum_bytes == 0
            || options.maximum_bytes > 64 * 1024 * 1024
        {
            return Err(QualificationError::Invalid);
        }
        w.charge(candidate.work_used())?;
        let inherited = w.used;
        let (end, solid) = bounds(candidate, &mut w)?;
        let start = candidate
            .as_str()
            .find("\nDATA;\n")
            .ok_or(QualificationError::Invalid)?
            + 7;
        let text = &candidate.as_str()[start..end];
        let mut actual_solid = 0u32;
        self.append_typed(|file| {
            let offset = file
                .next
                .checked_sub(1)
                .ok_or(QualificationError::Invalid)?;
            let maximum = offset.checked_add(solid).ok_or(QualificationError::Limit)?;
            if maximum > crate::step_read::MAX_STEP_ENTITIES {
                return Err(QualificationError::Limit);
            }
            actual_solid = u32::try_from(maximum).map_err(|_| QualificationError::Limit)?;
            let digits = maximum.to_string().len();
            w.charge(text.len())?;
            let hashes = text.bytes().filter(|b| *b == b'#').count();
            let growth = text
                .len()
                .checked_add(
                    hashes
                        .checked_mul(digits)
                        .ok_or(QualificationError::Limit)?,
                )
                .ok_or(QualificationError::Limit)?;
            let requested = file
                .data
                .len()
                .checked_add(growth)
                .ok_or(QualificationError::Limit)?;
            let old = file.data.capacity();
            // HashMap capacity bound follows the pinned existing literal replay
            // convention. Fixed64KiB covers <=1024 solid entries and formatting
            // scratch; no reference interning merges the replayed topology.
            let held = candidate
                .text
                .capacity()
                .checked_add(candidate.source.capacity() * std::mem::size_of::<u64>())
                .and_then(|v| v.checked_add(file.points.capacity() * 128))
                .and_then(|v| v.checked_add(64 * 1024))
                .ok_or(QualificationError::Limit)?;
            w.bytes(
                held.checked_add(old)
                    .and_then(|v| v.checked_add(requested))
                    .ok_or(QualificationError::Limit)?,
            )?;
            file.data
                .try_reserve_exact(growth)
                .map_err(|_| QualificationError::Limit)?;
            w.bytes(
                held.checked_add(old)
                    .and_then(|v| v.checked_add(file.data.capacity()))
                    .ok_or(QualificationError::Limit)?,
            )?;
            for line in text.split_inclusive('\n') {
                w.charge(line.len() + digits)?;
                let (id, body) = line.split_once('=').ok_or(QualificationError::Invalid)?;
                let id = id
                    .strip_prefix('#')
                    .and_then(|id| id.parse::<usize>().ok())
                    .ok_or(QualificationError::Invalid)?;
                write!(file.data, "#{}=", offset + id).map_err(|_| QualificationError::Limit)?;
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
                        cursor += 1
                    }
                    let reference = body[a..cursor]
                        .parse::<usize>()
                        .map_err(|_| QualificationError::Invalid)?;
                    if reference == 0 || reference > solid {
                        return Err(QualificationError::Invalid);
                    }
                    write!(file.data, "#{}", offset + reference)
                        .map_err(|_| QualificationError::Limit)?;
                    copied = cursor;
                }
                if quoted {
                    return Err(QualificationError::Invalid);
                }
                file.data.push_str(&body[copied..]);
            }
            file.next = maximum.checked_add(1).ok_or(QualificationError::Limit)?;
            w.charge(1)?;
            Ok::<_, QualificationError>(Entity(maximum))
        })?;
        Ok(QualifiedAppendReceipt {
            work: w.used,
            append: w.used - inherited,
            peak: w.peak,
            solid: actual_solid,
            original_volume: candidate.original_volume,
            healed_volume: candidate.healed_volume,
            position: candidate.composed_position,
            derivative: candidate.composed_derivative,
            sewing: candidate
                .healed
                .maximum_sewing_error_mm()
                .max(candidate.original.maximum_sewing_error_mm()),
        })
    }
}
