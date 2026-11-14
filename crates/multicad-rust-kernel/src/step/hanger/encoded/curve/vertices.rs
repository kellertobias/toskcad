//! One chosen finite point per exact closed-graph vertex owner. Point/support
//! incidence is bounded, never inferred from shared IDs or a proximity weld.
use super::literal::{Precision, quantize_at};
use super::*;

pub struct LiteralVertex {
    owner: ClosedVertex,
    point: [EncodedScalar; 3],
    field_error: f64,
    source_error: f64,
}
impl LiteralVertex {
    pub fn owner(&self) -> ClosedVertex {
        self.owner
    }
    pub fn coordinates(&self) -> &[EncodedScalar; 3] {
        &self.point
    }
    /// Euclidean error from the actual encoded field at EVERY incident source
    /// footpoint. Initial numeric binding uses L1; peer-owner transfer uses
    /// the two earned Euclidean authored-source receipts.
    pub fn field_incidence_error_mm(&self) -> f64 {
        self.field_error
    }
    /// Includes the original authored source-to-encoded-field receipt.
    pub fn source_incidence_error_mm(&self) -> f64 {
        self.source_error
    }
}
pub struct LiteralVertices<'c, 'd, 'f, 's> {
    source: &'c ClosedField<'d, 'f, 's>,
    vertices: Vec<LiteralVertex>,
    work: usize,
    bytes: usize,
}
impl<'c, 'd, 'f, 's> LiteralVertices<'c, 'd, 'f, 's> {
    pub fn source(&self) -> &'c ClosedField<'d, 'f, 's> {
        self.source
    }
    pub fn vertices(&self) -> &[LiteralVertex] {
        &self.vertices
    }
    /// Includes the inherited source/domain/closed-incidence work ONCE.
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_bytes(&self) -> usize {
        self.bytes
    }
}
fn evaluate<const D: usize>(input: &[[I; D]], t: I, w: &mut Work<'_>) -> Result<[I; D], Error> {
    w.charge(32 * D * input.len() * (input.len() - 1) / 2 + input.len() * D)?;
    // Sealed rows have <=11 controls. This tiny Vec and its interval values
    // live inside the caller's existing exact-arithmetic scratch reserve.
    let mut values = input.to_vec();
    let one = p(1.) - t;
    for n in (1..values.len()).rev() {
        for i in 0..n {
            let (lo, hi) = values.split_at_mut(i + 1);
            for (a, b) in lo[i].iter_mut().zip(&hi[0]) {
                *a = *a * one + *b * t;
            }
        }
    }
    Ok(values.remove(0))
}
fn footpoint(
    f: &EncodedField<'_>,
    cell: &FiniteCell,
    corner: &FiniteCorner,
    w: &mut Work<'_>,
) -> Result<([I; 2], N, f64), Error> {
    let owner = f.cells[cell.source_cell];
    let row = &f.meridians[owner.meridian];
    let source = f.source.source_boundary().source();
    let (pitch, starts, left) = if row.component == Component::Metric {
        let s = source.metric().spec();
        (s.pitch_mm, s.starts, s.left_hand)
    } else {
        let s = source.wood().spec();
        (s.pitch_mm, s.starts, s.left_hand)
    };
    let pitch = binary(pitch, w)?;
    let q = corner.phase.number(w)?;
    let v = corner.turn.number(w)?;
    let u = parameter(row, &pitch, &q, w)?;
    let mut meridian = Vec::with_capacity(row.degree() + 1);
    for i in 0..=row.degree() {
        meridian.push([
            row.radius[i].number(w)?.enclosure(w)?,
            row.weights[i].number(w)?.enclosure(w)?,
        ]);
    }
    let h = evaluate(&meridian, u, w)?;
    // Actual u is in[0,1]; the SAME sealed field weights are[1,1] or
    // [1,1/2,1]. Bernstein convexity proves W in[1/2,1] independently
    // of interval dependency. Intersect with this earned positive range.
    let weight = I::new(h[1].lo.max(0.5), h[1].hi.min(1.));
    let radius = h[0] / weight;
    // Exact knot comparisons choose the actual angular owner. An endpoint
    // shared by two pieces has the same exact coefficient owner on both.
    let mut piece = 0;
    while piece + 1 < f.angular.coefficients().len()
        && v.cmp(&binary(f.angular.knots()[piece + 1], w)?, w)? != Ordering::Less
    {
        piece += 1;
    }
    let a = binary(f.angular.knots()[piece], w)?;
    let b = binary(f.angular.knots()[piece + 1], w)?;
    let t = v.sub(&a, w)?.div(&b.sub(&a, w)?, w)?.enclosure(w)?;
    let index = if left {
        f.angular.coefficients().len() - 1 - piece
    } else {
        piece
    };
    let mut angular = Vec::with_capacity(11);
    for j in 0..=10 {
        let c = &f.angular.coefficients()[index][if left { 10 - j } else { j }];
        angular.push([N(c[0].clone()).enclosure(w)?, N(c[1].clone()).enclosure(w)?]);
    }
    let cs = evaluate(&angular, t, w)?;
    // q(u)=q is the earned inverse relation of the SAME exact phase/weight
    // scalar owner. Preserve its exact axial value instead of independently
    // rounding an evaluation of q(u). In particular every cap Z stays exact.
    let mut z = q
        .add(
            &pitch.mul(&N::integer(i64::from(owner.source_key.total_period)), w)?,
            w,
        )?
        .add(
            &pitch.mul(&N::integer(i64::from(starts)), w)?.mul(&v, w)?,
            w,
        )?;
    if row.component == Component::Metric {
        z = z
            .sub(&binary(source.metric().full_length_mm(), w)?, w)?
            .sub(&binary(source.spec().plain_length_mm, w)?, w)?;
    }
    Ok((
        [radius * cs[0], radius * cs[1]],
        z,
        row.bounds.position_error_mm,
    ))
}
fn proposal(actual: I, w: &mut Work<'_>) -> Result<N, Error> {
    // The enclosure is authority; its midpoint is only a finite proposal.
    if !actual.valid() {
        return Err(Error::Unresolved);
    }
    quantize_at(&binary(actual.midpoint(), w)?, Precision::Bits64, w)
}
fn bind(point: &[EncodedScalar; 3], xy: [I; 2], z: &N, w: &mut Work<'_>) -> Result<f64, Error> {
    let mut error = p(0.);
    for (chosen, actual) in point[..2].iter().zip(xy) {
        error = error + p((chosen.number(w)?.enclosure(w)? - actual).magnitude());
    }
    // Both are exact chosen scalar owners, not an interval-zero equality.
    if point[2].number(w)?.cmp(z, w)? != Ordering::Equal {
        return Err(Error::Unresolved);
    }
    if !error.valid() {
        return Err(Error::Unresolved);
    }
    Ok(error.hi)
}
impl<'d, 'f, 's> ClosedField<'d, 'f, 's> {
    pub fn literal_vertices<'c>(
        &'c self,
        options: CurveOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<LiteralVertices<'c, 'd, 'f, 's>, Error> {
        options.validate()?;
        let mut w = Work::new(options.work, interrupt)?;
        w.charge(self.work_used())?;
        let mut m = Memory::new(
            self.retained_bytes() + size_of::<LiteralVertices>(),
            options.bytes,
        )?;
        m.add(SCRATCH)?;
        let domains = self.domains();
        let f = domains.field();
        let source = f.source.source_boundary().source();
        let mut vertices = m.allocate(self.vertices().len(), &mut w)?;
        // FiniteDomains creates every Thread vertex while appending a genuine
        // positive cell corner. Its sealed arena contains no orphan vertices;
        // ClosedField then proves every retained vertex link is a cycle.
        // Thus every temporary slot below is overwritten by an actual corner
        // before any lookup. This exact arena map replaces repeated searches.
        let mut references = m.allocate(domains.vertices().len(), &mut w)?;
        references.resize(domains.vertices().len(), (0usize, 0usize));
        for (cell_index, cell) in domains.cells().iter().enumerate() {
            for (corner_index, corner) in cell.corners().iter().enumerate() {
                w.charge(1)?;
                references[corner.vertex] = (cell_index, corner_index);
            }
        }
        for &owner in self.vertices() {
            w.charge(1)?;
            let (point, field_box) = match owner {
                ClosedVertex::Thread(vertex) => {
                    let (cell_index, corner_index) = references[vertex];
                    let cell = &domains.cells()[cell_index];
                    let corner = &cell.corners()[corner_index];
                    let (xy, z, source_error) = footpoint(f, cell, corner, &mut w)?;
                    let mut chosen = std::array::from_fn(|_| N::integer(0));
                    for (out, actual) in chosen[..2].iter_mut().zip(xy) {
                        *out = proposal(actual, &mut w)?;
                    }
                    chosen[2] = z.clone();
                    (chosen, Some((xy, z, source_error)))
                }
                ClosedVertex::Plain { end, quarter } => {
                    // The sealed graph creates exactly four cardinal owners.
                    let cs = [[1, 0], [0, 1], [-1, 0], [0, -1]][usize::from(quarter)];
                    let radius = binary(source.spec().plain_radius_mm, &mut w)?;
                    (
                        [
                            radius.mul(&N::integer(cs[0]), &mut w)?,
                            radius.mul(&N::integer(cs[1]), &mut w)?,
                            if end == End::Start {
                                N::integer(0)
                                    .sub(&binary(source.spec().plain_length_mm, &mut w)?, &mut w)?
                            } else {
                                N::integer(0)
                            },
                        ],
                        None,
                    )
                }
            };
            let point = [
                EncodedScalar::retain(&point[0], &mut m, &mut w)?,
                EncodedScalar::retain(&point[1], &mut m, &mut w)?,
                EncodedScalar::retain(&point[2], &mut m, &mut w)?,
            ];
            let (field_error, source_error) = if let Some((xy, z, source_error)) = field_box {
                let error = bind(&point, xy, &z, &mut w)?;
                (error, (p(error) + p(source_error)).hi)
            } else {
                (0., 0.)
            };
            vertices.push(LiteralVertex {
                owner,
                point,
                field_error,
                source_error,
            });
        }
        m.release(&references);
        drop(references);
        // Each finite corner has EXACTLY the source identity created in
        // finite::vertex: q=p transfers one period; v=1 transfers starts
        // periods. Both are mathematical identities of the same true-period
        // authored source, including the crest seam. We do not assert that
        // separately approximated neighboring tensor supports are identical.
        // For every actual incident field P_i and the shared authored point S:
        // |chosen-P_i| <= |chosen-S| + |S-P_i|. Both bounds are earned receipts.
        for cell in domains.cells() {
            let row = &f.meridians[f.cells[cell.source_cell].meridian];
            for corner in cell.corners() {
                w.charge(1)?;
                let point = &mut vertices[corner.vertex];
                point.field_error = point
                    .field_error
                    .max((p(point.source_error) + p(row.bounds.position_error_mm)).hi);
                if point.source_error > options.position_error_mm
                    || point.field_error > options.position_error_mm
                {
                    return Err(Error::Unresolved);
                }
            }
        }
        w.charge(1)?;
        Ok(LiteralVertices {
            source: self,
            vertices,
            work: w.used,
            bytes: m.used - SCRATCH,
        })
    }
}
#[cfg(test)]
mod tests;
