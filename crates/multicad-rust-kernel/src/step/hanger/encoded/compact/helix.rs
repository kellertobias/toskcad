//! A long helix is an ACTUAL boundary row of a compact tensor surface.
//! Finite end parameters come from the original shared graph stations.
use super::*;
use num_traits::One;
use std::cmp::Ordering;
pub struct HelixPcurve {
    face: usize,
    u: f64,
    v_shift: EncodedScalar,
}
impl HelixPcurve {
    pub fn face(&self) -> usize {
        self.face
    }
    pub fn u(&self) -> f64 {
        self.u
    }
    pub fn v_shift(&self) -> &EncodedScalar {
        &self.v_shift
    }
}
pub struct CompactHelix<'b, 'g, 'c, 'd, 'f, 's> {
    vertices: &'g LiteralVertices<'c, 'd, 'f, 's>,
    source: &'b CompactBoundary<'c, 'd, 'f, 's>,
    edge: usize,
    surface: usize,
    row: usize,
    range: [EncodedScalar; 2],
    pcurves: [HelixPcurve; 2],
    residual: [f64; 2],
    same_sense: bool,
    work: usize,
    bytes: usize,
}
impl CompactHelix<'_, '_, '_, '_, '_, '_> {
    pub fn source(&self) -> &CompactBoundary<'_, '_, '_, '_> {
        self.source
    }
    pub fn vertices_source(&self) -> &LiteralVertices<'_, '_, '_, '_> {
        self.vertices
    }
    pub fn edge(&self) -> usize {
        self.edge
    }
    pub fn points(&self) -> &[[EncodedScalar; 3]] {
        let s = &self.source.surfaces().surfaces()[self.surface];
        let width = s.dimensions()[1];
        &s.points()[self.row * width..(self.row + 1) * width]
    }
    pub fn weight(&self) -> &EncodedScalar {
        &self.source.surfaces().surfaces()[self.surface].row_weights()[self.row]
    }
    pub fn knots(&self) -> &[f64] {
        self.source.surfaces().surfaces()[self.surface].knots()
    }
    pub fn range(&self) -> &[EncodedScalar; 2] {
        &self.range
    }
    pub fn pcurves(&self) -> &[HelixPcurve; 2] {
        &self.pcurves
    }
    pub fn vertex_residual_mm(&self) -> [f64; 2] {
        self.residual
    }
    pub fn same_sense(&self) -> bool {
        self.same_sense
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_bytes(&self) -> usize {
        self.bytes
    }
}
fn boundary_row(
    source: &ClosedField<'_, '_, '_>,
    surfaces: &CompactSurfaces<'_, '_>,
    edge: usize,
    old_face: usize,
    w: &mut Work<'_>,
) -> Result<(usize, usize), Error> {
    let ClosedCurve::Thread(domain_edge) = source.edges()[edge].curve else {
        return Err(Error::Invalid);
    };
    let ClosedFaceKind::ThreadCell(domain_cell) = source.faces()[old_face].kind else {
        return Err(Error::Invalid);
    };
    let cell = &source.domains().cells()[domain_cell];
    let field = source.domains().field();
    let row = &field.meridians[field.cells[cell.source_cell].meridian];
    let corner = cell
        .corners()
        .iter()
        .find(|c| c.vertex == source.domains().edges()[domain_edge].vertices[0])
        .ok_or(Error::Unresolved)?;
    let q = corner.phase.number(w)?;
    let lower = row.phase[0].number(w)?.div(&row.weights[0].number(w)?, w)?;
    let u = if q.cmp(&lower, w)? == Ordering::Equal {
        0
    } else {
        let upper = row.phase[row.degree()]
            .number(w)?
            .div(&row.weights[row.degree()].number(w)?, w)?;
        if q.cmp(&upper, w)? != Ordering::Equal {
            return Err(Error::Unresolved);
        }
        row.degree()
    };
    Ok((surfaces.cells()[cell.source_cell].surface(), u))
}
impl<'c, 'd, 'f, 's> CompactBoundary<'c, 'd, 'f, 's> {
    pub fn helix<'b, 'g>(
        &'b self,
        edge: usize,
        vertices: &'g LiteralVertices<'c, 'd, 'f, 's>,
        options: DomainOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Option<CompactHelix<'b, 'g, 'c, 'd, 'f, 's>>, Error> {
        if !std::ptr::eq(vertices.source(), self.source()) {
            return Err(Error::Invalid);
        }
        let owner = self.edges().get(edge).ok_or(Error::Invalid)?;
        let ChartEdgeKind::Helix { source_edge } = owner.kind else {
            return Ok(None);
        };
        let field = self.surfaces().field();
        let inherited =
            self.retained_bytes() + self.surfaces().retained_bytes() + vertices.retained_bytes();
        let (mut w, mut m) = crate::step::hanger::domain::ledgers(options, inherited, interrupt)?;
        w.charge(self.work_used() + vertices.work_used() - self.source().work_used())?;
        m.add(SCRATCH)?;
        let base = m.used;
        m.add(size_of::<CompactHelix<'_, '_, '_, '_, '_, '_>>())?;
        let mut faces = Vec::new();
        for (i, face) in self.source().faces().iter().enumerate() {
            for c in face.loops().iter().flatten() {
                w.charge(1)?;
                if c.edge == source_edge {
                    faces.push(i)
                }
            }
        }
        if faces.len() != 2 {
            return Err(Error::Unresolved);
        }
        let (surface, row) = boundary_row(
            self.source(),
            self.surfaces(),
            source_edge,
            faces[0],
            &mut w,
        )?;
        let (peer, peer_row) = boundary_row(
            self.source(),
            self.surfaces(),
            source_edge,
            faces[1],
            &mut w,
        )?;
        let s = &self.surfaces().surfaces()[surface];
        let peer_surface = &self.surfaces().surfaces()[peer];
        let width = s.dimensions()[1];
        let pwidth = peer_surface.dimensions()[1];
        let z0 = s.points()[row * width][2].number(&mut w)?;
        let peer_z0 = peer_surface.points()[peer_row * pwidth][2].number(&mut w)?;
        let (pitch, starts) = super::component(field, s.meridian());
        let lead = binary(pitch, &mut w)?.mul(&N::integer(starts.into()), &mut w)?;
        let shift = z0.sub(&peer_z0, &mut w)?.div(&lead, &mut w)?;
        if !shift.0.denom().is_one() {
            return Err(Error::Unresolved);
        }
        let mut ranges: [N; 2] = std::array::from_fn(|_| N::integer(0));
        let mut residual = [0.; 2];
        for (e, &v) in owner.vertices.iter().enumerate() {
            let point = vertices.vertices()[v].coordinates();
            ranges[e] = point[2]
                .number(&mut w)?
                .sub(&z0, &mut w)?
                .div(&lead, &mut w)?;
            let t = ranges[e].enclosure(&mut w)?;
            let limit = *s.knots().last().ok_or(Error::Unresolved)?;
            if ranges[e].cmp(&N::integer(0), &mut w)? == Ordering::Less
                || ranges[e].cmp(&binary(limit, &mut w)?, &mut w)? == Ordering::Greater
            {
                return Err(Error::Unresolved);
            }
            let approximate = t.midpoint();
            let mut span = s
                .knots()
                .partition_point(|&k| k < approximate)
                .saturating_sub(1)
                .min(s.knots().len() - 2);
            while ranges[e].cmp(&binary(s.knots()[span], &mut w)?, &mut w)? == Ordering::Less {
                w.charge(1)?;
                span -= 1;
            }
            while ranges[e].cmp(&binary(s.knots()[span + 1], &mut w)?, &mut w)? == Ordering::Greater
            {
                w.charge(1)?;
                span += 1;
            }
            let start = span * 10;
            let a = binary(s.knots()[span], &mut w)?;
            let b = binary(s.knots()[span + 1], &mut w)?;
            let normalized = ranges[e]
                .sub(&a, &mut w)?
                .div(&b.sub(&a, &mut w)?, &mut w)?;
            let controls = s.points()[row * width + start..row * width + start + 11]
                .iter()
                .map(|p| {
                    Ok([
                        p[0].number(&mut w)?,
                        p[1].number(&mut w)?,
                        p[2].number(&mut w)?,
                    ])
                })
                .collect::<Result<Vec<_>, Error>>()?;
            let actual = super::super::curve::eval_rational_curve(&controls, &normalized, &mut w)?;
            let mut error = p(0.);
            for axis in 0..3 {
                error = error
                    + p(actual[axis]
                        .sub(&point[axis].number(&mut w)?, &mut w)?
                        .enclosure(&mut w)?
                        .magnitude())
            }
            if !error.valid() || error.hi > 1e-7 {
                return Err(Error::Unresolved);
            }
            residual[e] = error.hi;
        }
        let range = [
            EncodedScalar::retain(&ranges[0], &mut m, &mut w)?,
            EncodedScalar::retain(&ranges[1], &mut m, &mut w)?,
        ];
        let same_sense = ranges[0].cmp(&ranges[1], &mut w)? == Ordering::Less;
        let pcurves = [
            HelixPcurve {
                face: self.original_faces()[faces[0]],
                u: if row == 0 { 0. } else { 1. },
                v_shift: EncodedScalar::retain(&N::integer(0), &mut m, &mut w)?,
            },
            HelixPcurve {
                face: self.original_faces()[faces[1]],
                u: if peer_row == 0 { 0. } else { 1. },
                v_shift: EncodedScalar::retain(&shift, &mut m, &mut w)?,
            },
        ];
        w.charge(1)?;
        Ok(Some(CompactHelix {
            vertices,
            source: self,
            edge,
            surface,
            row,
            range,
            pcurves,
            residual,
            same_sense,
            work: w.used,
            bytes: m.used - base,
        }))
    }
}
#[cfg(test)]
mod tests;
