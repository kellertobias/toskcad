//! Exact Cartesian/weight encoding of actual tensor isocurves. This does not
//! establish a finite trim, endpoint incidence, closed shell or occupied body.
use super::*;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SideSupportKind {
    Meridian {
        turn_end: End,
    },
    Helix {
        phase_end: End,
        angular_piece: usize,
    },
}
/// Exact numeric owner and its finite decimal representation, prepared together.
pub struct SidePole {
    point: [EncodedScalar; 3],
    weight: EncodedScalar,
    step_point: [FiniteDecimal; 3],
    step_weight: FiniteDecimal,
}
impl SidePole {
    pub fn point(&self) -> &[EncodedScalar; 3] {
        &self.point
    }
    pub fn weight(&self) -> &EncodedScalar {
        &self.weight
    }
    pub fn step_point(&self) -> &[FiniteDecimal; 3] {
        &self.step_point
    }
    pub fn step_weight(&self) -> &FiniteDecimal {
        &self.step_weight
    }
}
pub struct SideSupport<'f, 's> {
    field: &'f EncodedField<'s>,
    cell: usize,
    kind: SideSupportKind,
    parameter: [f64; 2],
    poles: Vec<SidePole>,
    work: usize,
    bytes: usize,
}
impl SideSupport<'_, '_> {
    pub fn field(&self) -> &EncodedField<'_> {
        self.field
    }
    pub fn source_cell(&self) -> usize {
        self.cell
    }
    pub fn kind(&self) -> SideSupportKind {
        self.kind
    }
    /// ORIGINAL surface u for a meridian, physical unsigned turn v for a helix.
    pub fn parameter(&self) -> [f64; 2] {
        self.parameter
    }
    pub fn degree(&self) -> usize {
        self.poles.len() - 1
    }
    pub fn poles(&self) -> &[SidePole] {
        &self.poles
    }
    pub fn position_error_mm(&self) -> f64 {
        self.field.meridians[self.field.cells[self.cell].meridian]
            .bounds
            .position_error_mm
    }
    /// Original surface u or global physical-turn v, matching parameter().
    pub fn derivative_error_mm(&self) -> f64 {
        let axis = usize::from(matches!(self.kind, SideSupportKind::Helix { .. }));
        self.field.meridians[self.field.cells[self.cell].meridian]
            .bounds
            .derivative_error_mm[axis]
    }
    /// Borrowed-field preparation is separate; only this query's work is counted.
    pub fn work_used(&self) -> usize {
        self.work
    }
    /// Conservative owned metadata accounting; borrowed field/scratch are separate.
    pub fn retained_bytes(&self) -> usize {
        self.bytes
    }
}
/// Source-owned exact supports, sharing only mathematical constants. Distinct
/// requests retain distinct source owners even if their numeric poles coincide.
pub struct SideSupportBatch<'f, 's> {
    supports: Vec<SideSupport<'f, 's>>,
    work: usize,
    bytes: usize,
}
impl<'f, 's> SideSupportBatch<'f, 's> {
    pub fn supports(&self) -> &[SideSupport<'f, 's>] {
        &self.supports
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_bytes(&self) -> usize {
        self.bytes
    }
}
impl<'s> EncodedField<'s> {
    /// All requests use one caller work/byte balance and one exact constant
    /// cache. Preparation of the borrowed field is charged separately.
    pub fn side_support_batch(
        &self,
        requests: &[(usize, SideSupportKind)],
        options: DomainOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<SideSupportBatch<'_, 's>, Error> {
        let (mut w, mut m) = crate::step::hanger::domain::ledgers(options, self.bytes, interrupt)?;
        m.add(super::super::SCRATCH)?;
        let base = m.used;
        m.add(size_of::<SideSupportBatch<'_, '_>>())?;
        let mut supports = m.allocate(requests.len(), &mut w)?;
        let mut cache = super::super::decimal::DecimalCache::new(&mut m, &mut w)?;
        for &(cell, kind) in requests {
            supports.push(prepare_cached(
                self, cell, kind, &mut cache, &mut m, &mut w,
            )?);
        }
        cache.release(&mut m);
        w.charge(1)?;
        Ok(SideSupportBatch {
            supports,
            work: w.used,
            bytes: m.used - base,
        })
    }
    pub fn side_support(
        &self,
        cell: usize,
        kind: SideSupportKind,
        options: DomainOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<SideSupport<'_, 's>, Error> {
        let (mut w, mut m) = crate::step::hanger::domain::ledgers(options, self.bytes, interrupt)?;
        m.add(super::super::SCRATCH)?;
        let mut support = prepare(self, cell, kind, &mut m, &mut w)?;
        support.work = w.used;
        Ok(support)
    }
}
// Aggregate callers reserve the SAME shared scratch once, and pass their live
// memory/work ledgers. No nested query resets or inherited-field recharging.
fn prepare<'f, 's>(
    field: &'f EncodedField<'s>,
    cell: usize,
    kind: SideSupportKind,
    m: &mut Memory,
    w: &mut Work<'_>,
) -> Result<SideSupport<'f, 's>, Error> {
    let mut cache = super::super::decimal::DecimalCache::new(m, w)?;
    let output = prepare_cached(field, cell, kind, &mut cache, m, w)?;
    cache.release(m);
    Ok(output)
}
fn prepare_cached<'f, 's>(
    field: &'f EncodedField<'s>,
    cell: usize,
    kind: SideSupportKind,
    cache: &mut super::super::decimal::DecimalCache,
    m: &mut Memory,
    w: &mut Work<'_>,
) -> Result<SideSupport<'f, 's>, Error> {
    let cache_before = cache.retained_bytes();
    let before = m.used;
    let start_work = w.used;
    w.charge(1)?;
    let owner = field.cells.get(cell).ok_or(Error::Invalid)?;
    let row = &field.meridians[owner.meridian];
    let (fixed_u, piece, column, count, parameter) = match kind {
        SideSupportKind::Meridian {
            turn_end: End::Start,
        } => (None, 0, 0, row.degree() + 1, [0., 1.]),
        SideSupportKind::Meridian { turn_end: End::End } => (
            None,
            field.angular.coefficients().len() - 1,
            10,
            row.degree() + 1,
            [0., 1.],
        ),
        SideSupportKind::Helix {
            phase_end,
            angular_piece,
        } => {
            if angular_piece >= field.angular.coefficients().len() {
                return Err(Error::Invalid);
            }
            let u = if phase_end == End::Start {
                0
            } else {
                row.degree()
            };
            (
                Some(u),
                angular_piece,
                0,
                11,
                [
                    field.angular.knots()[angular_piece],
                    field.angular.knots()[angular_piece + 1],
                ],
            )
        }
    };
    m.add(size_of::<SideSupport<'_, '_>>())?;
    let mut poles = m.allocate(count, w)?;
    for i in 0..count {
        let (u, v) = if let Some(u) = fixed_u {
            (u, i)
        } else {
            (i, column)
        };
        let (_, h) = field.homogeneous_control(cell, u, piece, v, w)?;
        // Sealed field weights are exactly 1 or 1/2. The retained phase,
        // radius and angular scalars therefore have terminating Cartesian
        // quotients here, unlike a general restricted rational cap support.
        let point = [
            h[0].div(&h[3], w)?,
            h[1].div(&h[3], w)?,
            h[2].div(&h[3], w)?,
        ];
        let encoded = [
            EncodedScalar::retain(&point[0], m, w)?,
            EncodedScalar::retain(&point[1], m, w)?,
            EncodedScalar::retain(&point[2], m, w)?,
        ];
        let weight = EncodedScalar::retain(&h[3], m, w)?;
        let step_point = [
            cache.format(&point[0], m, w)?,
            cache.format(&point[1], m, w)?,
            cache.format(&point[2], m, w)?,
        ];
        let step_weight = cache.format(&h[3], m, w)?;
        poles.push(SidePole {
            point: encoded,
            weight,
            step_point,
            step_weight,
        });
    }
    w.charge(1)?;
    Ok(SideSupport {
        field,
        cell,
        kind,
        parameter,
        poles,
        work: w.used - start_work,
        bytes: m.used - before - (cache.retained_bytes() - cache_before),
    })
}
mod literal;
#[cfg(test)]
mod tests;
pub use literal::LiteralSideSupport;
