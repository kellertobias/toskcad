//! Source-owned exact parameter leaves. These are not foreign geometry headers.
use super::*;
#[cfg(test)]
mod tests;
/// A canonical reduced rational emitted by the source's audited Q encoder.
/// Magnitudes use little-endian bytes; the denominator is strictly positive.
#[derive(Clone, Copy, Debug)]
pub struct ExactBoundaryRational<'a> {
    sign: i8,
    numerator: &'a [u8],
    denominator: &'a [u8],
}
impl<'a> ExactBoundaryRational<'a> {
    pub fn signum(self) -> i8 {
        self.sign
    }
    pub fn numerator_magnitude_le(self) -> &'a [u8] {
        self.numerator
    }
    pub fn denominator_magnitude_le(self) -> &'a [u8] {
        self.denominator
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ExactBoundaryDomain {
    Cell(usize),
    Curve(u32),
}
/// Compact exact values tied to one immutable original graph. Cell rows are
/// original [q,v] corners; curve rows are original affine [r,z,v] base/delta.
pub struct ExactBoundaryData<'a, S> {
    boundary: &'a Boundary<S>,
    domain: ExactBoundaryDomain,
    rows: Vec<Box<[u8]>>,
    columns: usize,
    work: usize,
}
impl<'a, S> ExactBoundaryData<'a, S> {
    pub fn boundary(&self) -> &'a Boundary<S> {
        self.boundary
    }
    pub fn domain(&self) -> ExactBoundaryDomain {
        self.domain
    }
    pub fn row_count(&self) -> usize {
        self.rows.len()
    }
    pub fn column_count(&self) -> usize {
        self.columns
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    /// Only this receipt's owned payload/header; EXCLUDES borrowed graph and
    /// caller-reserved arenas, which callers must count exactly once.
    pub fn owned_bytes(&self) -> usize {
        size_of::<Self>()
            + self.rows.capacity() * size_of::<Box<[u8]>>()
            + self.rows.iter().map(|r| r.len()).sum::<usize>()
    }
    pub fn coordinate(&self, row: usize, column: usize) -> Option<ExactBoundaryRational<'_>> {
        if column >= self.columns {
            return None;
        }
        let bytes = self.rows.get(row)?;
        let mut offset = 0;
        for _ in 0..column {
            let n = u16::from_le_bytes([bytes[offset + 1], bytes[offset + 2]]) as usize;
            let d = u16::from_le_bytes([bytes[offset + 3], bytes[offset + 4]]) as usize;
            offset += 5 + n + d;
        }
        let n = u16::from_le_bytes([bytes[offset + 1], bytes[offset + 2]]) as usize;
        let d = u16::from_le_bytes([bytes[offset + 3], bytes[offset + 4]]) as usize;
        Some(ExactBoundaryRational {
            sign: bytes[offset] as i8 - 1,
            numerator: &bytes[offset + 5..offset + 5 + n],
            denominator: &bytes[offset + 5 + n..offset + 5 + n + d],
        })
    }
}
// Q has at most8192bits in each operand. No BigInt survives the query: fresh
// boxed Q encoding has <=2053bytes/value, plus its actual Vec descriptor arena.
const LEAF_MAX: usize = 5 + 2 * (8192 / 8);
impl<S> Boundary<S> {
    fn exact_query_work<'a>(
        &self,
        maximum_work: usize,
        interrupt: Option<&'a dyn Fn() -> bool>,
    ) -> Result<Work<'a>, Error> {
        let mut work = Work {
            remaining: maximum_work,
            used: 0,
            interrupt,
        };
        work.step()?;
        if maximum_work > 16_000_000 {
            return Err(Error::Invalid);
        }
        Ok(work)
    }
    fn exact_query_rows(
        &self,
        rows: usize,
        columns: usize,
        reserved_bytes: usize,
        w: &mut Work<'_>,
    ) -> Result<Vec<Box<[u8]>>, Error> {
        // reserved_bytes EXCLUDES this graph. Graph::check adds its original
        // arena. query_headers adds Boundary's extension and live Vec headers.
        let reserved = reserved_bytes
            .checked_add(self.query_headers())
            .and_then(|n| n.checked_add(size_of::<ExactBoundaryData<'_, S>>()))
            .and_then(|n| n.checked_add(EXACT_SCRATCH))
            .and_then(|n| n.checked_add(rows * columns * LEAF_MAX))
            .ok_or(Error::Limit)?;
        self.graph.buffer(rows, reserved, w)
    }
    pub fn exact_cell_parameters(
        &self,
        index: usize,
        maximum_work: usize,
        reserved_bytes: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<ExactBoundaryData<'_, S>, Error> {
        let mut work = self.exact_query_work(maximum_work, interrupt)?;
        let cell = self.graph.cells.get(index).ok_or(Error::Invalid)?;
        if !matches!(cell.kind, Kind::Wall) {
            return Err(Error::Unresolved);
        }
        let mut rows = self.exact_query_rows(cell.count, 2, reserved_bytes, &mut work)?;
        let pitch = Q::float(self.graph.source.thread.spec().pitch_mm, &mut work)?;
        let lead = pitch.mul(
            &Q::integer(i64::from(self.graph.source.thread.spec().starts), &mut work)?,
            &mut work,
        )?;
        for corner in 0..cell.count {
            work.step()?;
            let point = &self.graph.vertices[cell.vertices[corner] as usize];
            let mut q = close::scalar(&point.exact, 0, &mut work)?;
            let mut v = close::scalar(&point.exact, 1, &mut work)?;
            // Q.enclosure returns the exact singleton for the producer's
            // representable endpoint1. Only that marker underwent q+=lead,
            // v1->0 in Patch::owned; no rounded near-end value is snapped.
            if cell.parameters[corner][1].lo == 1. && cell.parameters[corner][1].hi == 1. {
                q = q.sub(&lead, &mut work)?;
                v = Q::integer(1, &mut work)?;
            }
            rows.push(Q::encode(&[&q, &v], &mut work)?);
        }
        work.step()?;
        Ok(ExactBoundaryData {
            boundary: self,
            domain: ExactBoundaryDomain::Cell(index),
            rows,
            columns: 2,
            work: work.used,
        })
    }
    pub fn exact_curve_coefficients(
        &self,
        index: u32,
        maximum_work: usize,
        reserved_bytes: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<ExactBoundaryData<'_, S>, Error> {
        let mut work = self.exact_query_work(maximum_work, interrupt)?;
        let curve = self
            .graph
            .curves
            .get(index as usize)
            .ok_or(Error::Invalid)?;
        let mut rows = self.exact_query_rows(2, 3, reserved_bytes, &mut work)?;
        for row in 0..2 {
            let a = close::scalar(&curve.exact, row * 3, &mut work)?;
            let b = close::scalar(&curve.exact, row * 3 + 1, &mut work)?;
            let c = close::scalar(&curve.exact, row * 3 + 2, &mut work)?;
            rows.push(Q::encode(&[&a, &b, &c], &mut work)?);
        }
        work.step()?;
        Ok(ExactBoundaryData {
            boundary: self,
            domain: ExactBoundaryDomain::Curve(index),
            rows,
            columns: 3,
            work: work.used,
        })
    }
}
