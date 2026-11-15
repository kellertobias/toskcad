//! Exact multi-turn tensor charts, grouping original cells by meridian/strand.
//! This changes only chart partitioning. Trims, file and material admission
//! remain with the closed boundary producer.
use super::*;
pub struct CompactCell {
    surface: usize,
    turn_shift: i64,
}
impl CompactCell {
    pub fn surface(&self) -> usize {
        self.surface
    }
    /// Compact V = original local V + this exact integer.
    pub fn turn_shift(&self) -> i64 {
        self.turn_shift
    }
}
pub struct CompactSurface<'f, 's> {
    field: &'f EncodedField<'s>,
    meridian: usize,
    strand: i32,
    first_turn: i32,
    dimensions: [usize; 2],
    points: Vec<[EncodedScalar; 3]>,
    knots: Vec<f64>,
}
impl CompactSurface<'_, '_> {
    pub fn field(&self) -> &EncodedField<'_> {
        self.field
    }
    pub fn meridian(&self) -> usize {
        self.meridian
    }
    pub fn strand(&self) -> i32 {
        self.strand
    }
    pub fn first_turn(&self) -> i32 {
        self.first_turn
    }
    pub fn dimensions(&self) -> [usize; 2] {
        self.dimensions
    }
    pub fn degrees(&self) -> [usize; 2] {
        [self.dimensions[0] - 1, 10]
    }
    pub fn points(&self) -> &[[EncodedScalar; 3]] {
        &self.points
    }
    /// One positive defining weight per U row, constant along compact V.
    pub fn row_weights(&self) -> &[EncodedScalar] {
        &self.field.meridians[self.meridian].weights
    }
    /// Distinct V knots; degree10 internal multiplicity is10, end multiplicity11.
    pub fn knots(&self) -> &[f64] {
        &self.knots
    }
    pub fn source_bounds(&self) -> SupportBounds {
        self.field.meridians[self.meridian].bounds
    }
}
pub struct CompactSurfaces<'f, 's> {
    field: &'f EncodedField<'s>,
    surfaces: Vec<CompactSurface<'f, 's>>,
    cells: Vec<CompactCell>,
    work: usize,
    bytes: usize,
}
impl<'f, 's> CompactSurfaces<'f, 's> {
    pub fn field(&self) -> &EncodedField<'s> {
        self.field
    }
    pub fn surfaces(&self) -> &[CompactSurface<'f, 's>] {
        &self.surfaces
    }
    /// Entries are indexed by the ORIGINAL immutable field cell ID.
    pub fn cells(&self) -> &[CompactCell] {
        &self.cells
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_bytes(&self) -> usize {
        self.bytes
    }
    pub fn prepare(
        field: &'f EncodedField<'s>,
        options: DomainOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, Error> {
        let (mut w, mut m) =
            crate::step::hanger::domain::ledgers(options, field.retained_bytes(), interrupt)?;
        w.charge(field.work_used())?;
        m.add(SCRATCH)?;
        build(field, &mut m, &mut w)
    }
}
struct Group {
    meridian: usize,
    strand: i32,
    first: i32,
    last: i32,
    reference: usize,
}
fn component(field: &EncodedField<'_>, meridian: usize) -> (f64, i32) {
    let source = field.source.source_boundary().source();
    if field.meridians[meridian].component == Component::Metric {
        let s = source.metric().spec();
        (s.pitch_mm, s.starts as i32)
    } else {
        let s = source.wood().spec();
        (s.pitch_mm, s.starts as i32)
    }
}
pub(super) fn build<'f, 's>(
    field: &'f EncodedField<'s>,
    m: &mut Memory,
    w: &mut Work<'_>,
) -> Result<CompactSurfaces<'f, 's>, Error> {
    let base = m.used;
    m.add(size_of::<CompactSurfaces<'_, '_>>())?;
    let mut groups: Vec<Group> = m.allocate(field.cells.len(), w)?;
    let mut cells = m.allocate(field.cells.len(), w)?;
    for (index, cell) in field.cells.iter().enumerate() {
        let (_, starts) = component(field, cell.meridian);
        let strand = cell.source_key.total_period.rem_euclid(starts);
        let turn = cell.source_key.total_period.div_euclid(starts);
        let mut found = None;
        for (i, g) in groups.iter().enumerate() {
            w.charge(1)?;
            if g.meridian == cell.meridian && g.strand == strand {
                found = Some(i);
                break;
            }
        }
        let group = if let Some(i) = found {
            let g = &mut groups[i];
            if turn < g.first {
                g.first = turn;
                g.reference = index
            }
            g.last = g.last.max(turn);
            i
        } else {
            groups.push(Group {
                meridian: cell.meridian,
                strand,
                first: turn,
                last: turn,
                reference: index,
            });
            groups.len() - 1
        };
        cells.push(CompactCell {
            surface: group,
            turn_shift: i64::from(turn),
        });
    }
    let mut surfaces = m.allocate(groups.len(), w)?;
    for group in &groups {
        let row = &field.meridians[group.meridian];
        let (pitch, starts) = component(field, group.meridian);
        let turns = usize::try_from(i64::from(group.last) - i64::from(group.first) + 1)
            .map_err(|_| Error::Limit)?;
        let pieces = field.angular.coefficients().len();
        let per_turn = pieces.checked_mul(10).ok_or(Error::Limit)?;
        let width = turns
            .checked_mul(per_turn)
            .and_then(|x| x.checked_add(1))
            .ok_or(Error::Limit)?;
        let dimensions = [row.degree() + 1, width];
        let count = dimensions[0].checked_mul(width).ok_or(Error::Limit)?;
        let mut points = m.allocate(count, w)?;
        let knot_count = turns
            .checked_mul(pieces)
            .and_then(|x| x.checked_add(1))
            .ok_or(Error::Limit)?;
        let mut knots = m.allocate(knot_count, w)?;
        for i in 0..knot_count {
            w.charge(1)?;
            let value = if i + 1 == knot_count {
                turns as f64
            } else {
                (i / pieces) as f64 + field.angular.knots()[i % pieces]
            };
            knots.push(value)
        }
        let template_width = per_turn + 1;
        let mut template = m.allocate(dimensions[0] * template_width, w)?;
        for u in 0..dimensions[0] {
            for j in 0..template_width {
                let (piece, v) = if j == per_turn {
                    (pieces - 1, 10)
                } else {
                    (j / 10, j % 10)
                };
                let (_, h) = field.homogeneous_control(group.reference, u, piece, v, w)?;
                template.push([
                    h[0].div(&h[3], w)?,
                    h[1].div(&h[3], w)?,
                    h[2].div(&h[3], w)?,
                ]);
            }
        }
        let lead = binary(pitch, w)?.mul(&N::integer(starts.into()), w)?;
        for u in 0..dimensions[0] {
            for j in 0..width {
                let (turn, local) = if j + 1 == width {
                    (turns - 1, per_turn)
                } else {
                    (j / per_turn, j % per_turn)
                };
                let p = &template[u * template_width + local];
                let z = p[2].add(&lead.mul(&N::integer(turn as i64), w)?, w)?;
                points.push([
                    EncodedScalar::retain(&p[0], m, w)?,
                    EncodedScalar::retain(&p[1], m, w)?,
                    EncodedScalar::retain(&z, m, w)?,
                ]);
            }
        }
        m.release(&template);
        drop(template);
        surfaces.push(CompactSurface {
            field,
            meridian: group.meridian,
            strand: group.strand,
            first_turn: group.first,
            dimensions,
            points,
            knots,
        });
    }
    for cell in &mut cells {
        w.charge(1)?;
        cell.turn_shift -= i64::from(groups[cell.surface].first);
    }
    m.release(&groups);
    drop(groups);
    w.charge(1)?;
    Ok(CompactSurfaces {
        field,
        surfaces,
        cells,
        work: w.used,
        bytes: m.used - base,
    })
}
mod boundary;
#[cfg(test)]
mod tests;
pub use boundary::{
    ChartCoedge, ChartEdge, ChartEdgeKind, ChartFace, ChartFaceKind, CompactBoundary,
};
