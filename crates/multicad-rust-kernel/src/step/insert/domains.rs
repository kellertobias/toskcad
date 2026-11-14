//! Exact angular decomposition of the ORIGINAL source cells; no file surfaces.
use super::*;
use crate::{
    authored_flat_thread::AuthoredThread,
    authored_internal_treatment::{Error as SourceError, ExactBoundaryData, ExactBoundaryRational},
};
use std::cmp::Ordering;
mod exact;
pub(super) use exact::N;
#[cfg(test)]
mod tests;
const SCRATCH: usize = 16 * 1024 * 1024;
const QUERY_OUTPUT: usize = 64 * 1024;
#[derive(Clone, Copy, Debug)]
pub struct DomainOptions {
    pub maximum_stations: usize,
    pub maximum_slabs: usize,
    pub work: usize,
    pub bytes: usize,
}
impl Default for DomainOptions {
    fn default() -> Self {
        Self {
            maximum_stations: 4096,
            maximum_slabs: 32768,
            work: 16_000_000,
            bytes: 64 * 1024 * 1024,
        }
    }
}
/// Exact subdivision of the original angular coordinate; all neighboring cells
/// use its ID. Original corner events and canonical quarter turns stay distinct.
pub struct Station {
    sign: i8,
    numerator: Vec<u8>,
    denominator: Vec<u8>,
}
impl Station {
    pub fn signum(&self) -> i8 {
        self.sign
    }
    pub fn numerator_magnitude_le(&self) -> &[u8] {
        &self.numerator
    }
    pub fn denominator_magnitude_le(&self) -> &[u8] {
        &self.denominator
    }
    pub(super) fn number(&self, w: &mut Work<'_>) -> Result<N, Error> {
        N::parts(self.sign, &self.numerator, &self.denominator, w)
    }
    pub(super) fn cardinal(&self) -> Option<u8> {
        cardinal_parts(self.sign, &self.numerator, &self.denominator)
    }
    fn from_parts(
        sign: i8,
        n: &[u8],
        d: &[u8],
        memory: &mut Memory,
        w: &mut Work<'_>,
    ) -> Result<Self, Error> {
        let mut numerator = memory.allocate(n.len(), w)?;
        numerator.extend_from_slice(n);
        let mut denominator = memory.allocate(d.len(), w)?;
        denominator.extend_from_slice(d);
        Ok(Self {
            sign,
            numerator,
            denominator,
        })
    }
}
#[derive(Clone, Copy, Debug)]
pub struct Slab {
    pub side: Side,
    pub cell: usize,
    pub stations: [u32; 2],
    /// Actual source curves, ordered by q on the WHOLE open slab.
    pub rows: [BoundaryCurveKey; 2],
    /// Exact source endpoint-width zero. Not permission to omit a literal edge:
    /// future carrier construction must prove the ENTIRE row is one owned vertex.
    pub collapsed: [bool; 2],
}
pub struct DomainPlan<'a> {
    boundary: &'a Boundary,
    stations: [Vec<Station>; 2],
    cells: [Vec<ExactBoundaryData<'a, AuthoredThread>>; 2],
    slabs: Vec<Slab>,
    work: usize,
    bytes: usize,
}
impl<'a> DomainPlan<'a> {
    pub fn boundary(&self) -> &'a Boundary {
        self.boundary
    }
    pub fn stations(&self, side: Side) -> &[Station] {
        &self.stations[index(side)]
    }
    pub fn slabs(&self) -> &[Slab] {
        &self.slabs
    }
    pub fn exact_cell(
        &self,
        side: Side,
        cell: usize,
    ) -> Option<&ExactBoundaryData<'a, AuthoredThread>> {
        self.cells[index(side)].get(cell)
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_bytes(&self) -> usize {
        self.bytes
    }
    pub fn prepare(
        boundary: &'a Boundary,
        options: DomainOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, Error> {
        let mut w = Work::new(options.work, interrupt)?;
        w.charge(1)?;
        if options.maximum_stations == 0
            || options.maximum_stations > 16384
            || options.maximum_slabs == 0
            || options.maximum_slabs > 131072
            || options.bytes == 0
            || options.bytes > 64 * 1024 * 1024
        {
            return Err(Error::Invalid);
        }
        let mut memory = Memory::new(boundary.retained_bytes(), options.bytes)?;
        memory.add(size_of::<Self>() + SCRATCH)?;
        let mut stations = [
            memory.allocate(options.maximum_stations, &mut w)?,
            memory.allocate(options.maximum_stations, &mut w)?,
        ];
        let mut cells = [
            memory.allocate(boundary.cells(Side::Outer).len(), &mut w)?,
            memory.allocate(boundary.cells(Side::Inner).len(), &mut w)?,
        ];
        let mut slabs = memory.allocate(options.maximum_slabs, &mut w)?;
        for side in [Side::Outer, Side::Inner] {
            // Exact authored-turn chart subdivisions, not moved source events.
            for (sign, n, d) in [(0, 0, 1), (1, 1, 4), (1, 1, 2), (1, 3, 4), (1, 1, 1)] {
                insert_parts(
                    &mut stations[index(side)],
                    sign,
                    &[n],
                    &[d],
                    options.maximum_stations,
                    &mut memory,
                    &mut w,
                )?;
            }
            let source = boundary.source_boundary(side);
            for cell in 0..source.cells().len() {
                // Scratch is reused: no decoded BigInts are live during the
                // child query. Its8MiB reserve fits our16MiB reserved scratch.
                memory.add(QUERY_OUTPUT)?;
                let reserved = memory.used - source.retained_bytes() - SCRATCH;
                let data = source
                    .exact_cell_parameters(cell, w.remaining(), reserved, interrupt)
                    .map_err(failure)?;
                w.charge(data.work_used())?;
                memory.used -= QUERY_OUTPUT;
                memory.add(data.owned_bytes())?;
                for row in 0..data.row_count() {
                    insert_station(
                        &mut stations[index(side)],
                        data.coordinate(row, 1)
                            .expect("sealed two-column source query"),
                        options.maximum_stations,
                        &mut memory,
                        &mut w,
                    )?;
                }
                cells[index(side)].push(data);
            }
        }
        for side in [Side::Outer, Side::Inner] {
            let source = boundary.source_boundary(side);
            for (cell, data) in cells[index(side)].iter().enumerate() {
                let mut corners: [[N; 2]; 6] =
                    std::array::from_fn(|_| [N::integer(0), N::integer(0)]);
                for (row, p) in corners.iter_mut().enumerate().take(data.row_count()) {
                    *p = [
                        N::view(data.coordinate(row, 0).unwrap(), &mut w)?,
                        N::view(data.coordinate(row, 1).unwrap(), &mut w)?,
                    ];
                }
                let polygon = &corners[..data.row_count()];
                let mut low = &polygon[0][1];
                let mut high = low;
                for p in polygon {
                    if p[1].cmp(low, &mut w)? == Ordering::Less {
                        low = &p[1];
                    }
                    if p[1].cmp(high, &mut w)? == Ordering::Greater {
                        high = &p[1];
                    }
                }
                let mut curves = [curve_key(side, 0); 6];
                // The public occupied coedge traversal reverses the original
                // polygon on RH bare sources. Undo only that documented order.
                for (i, c) in source
                    .cell_coedges(cell)
                    .expect("enumerated original source cell")
                    .enumerate()
                {
                    let original = if source.cell_surface_same_sense() {
                        i
                    } else {
                        polygon.len() - 1 - i
                    };
                    curves[original] = curve_key(side, c.curve as u32);
                }
                for (station, pair) in stations[index(side)].windows(2).enumerate() {
                    w.charge(1)?;
                    let a = pair[0].number(&mut w)?;
                    let b = pair[1].number(&mut w)?;
                    if a.cmp(low, &mut w)? == Ordering::Less
                        || b.cmp(high, &mut w)? == Ordering::Greater
                    {
                        continue;
                    }
                    let (edges, collapsed) = section(polygon, &a, &b, &mut w)?;
                    if slabs.len() == options.maximum_slabs {
                        return Err(Error::Limit);
                    }
                    slabs.push(Slab {
                        side,
                        cell,
                        stations: [station as u32, station as u32 + 1],
                        rows: edges.map(|i| curves[i]),
                        collapsed,
                    });
                }
            }
        }
        w.charge(1)?;
        memory.used -= SCRATCH;
        Ok(Self {
            boundary,
            stations,
            cells,
            slabs,
            work: w.used,
            bytes: memory.used,
        })
    }
}
fn index(side: Side) -> usize {
    match side {
        Side::Outer => 0,
        Side::Inner => 1,
    }
}
pub(super) fn failure(e: SourceError) -> Error {
    match e {
        SourceError::Invalid => Error::Invalid,
        SourceError::Unresolved => Error::Unresolved,
        SourceError::Limit => Error::Limit,
        SourceError::Cancelled => Error::Cancelled,
    }
}
fn insert_station(
    stations: &mut Vec<Station>,
    value: ExactBoundaryRational<'_>,
    maximum: usize,
    memory: &mut Memory,
    w: &mut Work<'_>,
) -> Result<(), Error> {
    insert_parts(
        stations,
        value.signum(),
        value.numerator_magnitude_le(),
        value.denominator_magnitude_le(),
        maximum,
        memory,
        w,
    )
}
fn insert_parts(
    stations: &mut Vec<Station>,
    sign: i8,
    n: &[u8],
    d: &[u8],
    maximum: usize,
    memory: &mut Memory,
    w: &mut Work<'_>,
) -> Result<(), Error> {
    let number = N::parts(sign, n, d, w)?;
    let mut low = 0;
    let mut high = stations.len();
    while low < high {
        let mid = low + (high - low) / 2;
        match stations[mid].number(w)?.cmp(&number, w)? {
            Ordering::Equal => return Ok(()),
            Ordering::Less => low = mid + 1,
            Ordering::Greater => high = mid,
        }
    }
    if stations.len() == maximum {
        return Err(Error::Limit);
    }
    w.charge(1 + stations.len() - low)?;
    stations.insert(low, Station::from_parts(sign, n, d, memory, w)?);
    Ok(())
}

fn at(polygon: &[[N; 2]], edge: usize, v: &N, w: &mut Work<'_>) -> Result<N, Error> {
    let a = &polygon[edge];
    let b = &polygon[(edge + 1) % polygon.len()];
    a[0].add(
        &b[0]
            .sub(&a[0], w)?
            .mul(&v.sub(&a[1], w)?.div(&b[1].sub(&a[1], w)?, w)?, w)?,
        w,
    )
}
fn section(
    polygon: &[[N; 2]],
    a: &N,
    b: &N,
    w: &mut Work<'_>,
) -> Result<([usize; 2], [bool; 2]), Error> {
    let mid = a.add(b, w)?.div(&N::integer(2), w)?;
    // The producer clips a positive q/v rectangle by the two parallel,
    // distinct height planes. Its nonempty cells are convex and positive-area.
    // All their corner-v events occur in the shared station set, so an OPEN
    // slab midpoint crosses exactly two nonhorizontal edges, never a corner.
    let mut first = None;
    let mut last = None;
    for edge in 0..polygon.len() {
        let x = polygon[edge][1].cmp(&mid, w)?;
        let y = polygon[(edge + 1) % polygon.len()][1].cmp(&mid, w)?;
        if (x == Ordering::Less && y == Ordering::Greater)
            || (x == Ordering::Greater && y == Ordering::Less)
        {
            first.get_or_insert(edge);
            last = Some(edge);
        }
    }
    let mut crossings = [
        first.expect("positive source slab first edge"),
        last.expect("positive source slab last edge"),
    ];
    if at(polygon, crossings[0], &mid, w)?.cmp(&at(polygon, crossings[1], &mid, w)?, w)?
        == Ordering::Greater
    {
        crossings.swap(0, 1);
    }
    let mut collapsed = [false; 2];
    for (i, v) in [a, b].into_iter().enumerate() {
        // Exact convex source provenance and the absence of an interior event
        // imply nonnegative width at BOTH ends; two zero ends would contradict
        // the source cell's positive area. This is not foreign polygon repair.
        collapsed[i] = at(polygon, crossings[0], v, w)?
            .cmp(&at(polygon, crossings[1], v, w)?, w)?
            == Ordering::Equal;
    }
    Ok((crossings, collapsed))
}

pub(super) fn cardinal_parts(sign: i8, n: &[u8], d: &[u8]) -> Option<u8> {
    match (sign, n, d) {
        (0, _, _) => Some(0),
        (1, [1], [4]) => Some(1),
        (1, [1], [2]) => Some(2),
        (1, [3], [4]) => Some(3),
        (1, [1], [1]) => Some(4),
        _ => None,
    }
}
