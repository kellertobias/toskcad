//! Exact affine-row cap support, preserving the original source and chart.
//! This support has no finite literal trim, shared-vertex or file admission.
use super::*;
pub(super) struct AffineCap {
    pub(super) controls: Vec<[N; 4]>,
    pub(super) trim: [N; 2],
    pub(super) uv: [[N; 2]; 2],
    // Original meridian-u and angular-v derivatives in basis t.
    pub(super) derivative: [N; 2],
    pub(super) reversed: bool,
}
pub(super) fn relation(h: &[[N; 3]], w: &mut Work<'_>) -> Result<Option<[[N; 2]; 2]>, Error> {
    if h.len() != 2
        || h[0][2].cmp(&h[1][2], w)? != Ordering::Equal
        || h[0][2].cmp(&N::integer(0), w)? != Ordering::Greater
    {
        return Ok(None);
    }
    Ok(Some([
        [h[0][0].div(&h[0][2], w)?, h[0][1].div(&h[0][2], w)?],
        [h[1][0].div(&h[1][2], w)?, h[1][1].div(&h[1][2], w)?],
    ]))
}
pub(super) fn build(
    field: &EncodedField<'_>,
    index: usize,
    end: End,
    piece: usize,
    w: &mut Work<'_>,
) -> Result<Option<AffineCap>, Error> {
    w.charge(1)?;
    let cell = field.cells[index];
    let row = &field.meridians[cell.meridian];
    // ONE actual coefficient gate: a real rounded-root row reaches the
    // same fallback as any non-affine/projective row. Do not preclassify by
    // a degree label and leave its retained-coefficient check unreachable.
    let input = (0..=row.degree())
        .map(|i| {
            Ok([
                row.radius[i].number(w)?,
                row.phase[i].number(w)?,
                row.weights[i].number(w)?,
            ])
        })
        .collect::<Result<Vec<_>, Error>>()?;
    let Some([[r0, q0], [r1, q1]]) = relation(&input, w)? else {
        return Ok(None);
    };
    let source = field.source.source_boundary().source();
    let (pitch, starts, left, length) = if row.component == Component::Metric {
        let s = source.metric().spec();
        (s.pitch_mm, s.starts, s.left_hand, s.length_mm)
    } else {
        let s = source.wood().spec();
        (s.pitch_mm, s.starts, s.left_hand, s.length_mm)
    };
    let pitch = binary(pitch, w)?;
    let lead = pitch.mul(&N::integer(i64::from(starts)), w)?;
    let station = if end == End::Start {
        N::integer(0)
    } else {
        binary(length, w)?
    };
    let shift = station.sub(
        &pitch.mul(&N::integer(i64::from(cell.source_key.total_period)), w)?,
        w,
    )?;
    let a = binary(field.angular.knots()[piece], w)?;
    let b = binary(field.angular.knots()[piece + 1], w)?;
    let from = shift.sub(&q1, w)?.div(&lead, w)?;
    let to = shift.sub(&q0, w)?.div(&lead, w)?;
    let lo = if from.cmp(&a, w)? == Ordering::Greater {
        from
    } else {
        a.clone()
    };
    let hi = if to.cmp(&b, w)? == Ordering::Less {
        to
    } else {
        b.clone()
    };
    if lo.cmp(&hi, w)? != Ordering::Less {
        return Ok(None);
    }
    let width = b.sub(&a, w)?;
    let qspan = q1.sub(&q0, w)?;
    let ua = shift
        .sub(&lead.mul(&a, w)?, w)?
        .sub(&q0, w)?
        .div(&qspan, w)?;
    let ub = shift
        .sub(&lead.mul(&b, w)?, w)?
        .sub(&q0, w)?
        .div(&qspan, w)?;
    let ra = r0.add(&r1.sub(&r0, w)?.mul(&ua, w)?, w)?;
    let rb = r0.add(&r1.sub(&r0, w)?.mul(&ub, w)?, w)?;
    let ap = if left {
        field.angular.coefficients().len() - 1 - piece
    } else {
        piece
    };
    let coefficients = &field.angular.coefficients()[ap];
    let mut z = station;
    if row.component == Component::Metric {
        z = z
            .sub(&binary(source.metric().full_length_mm(), w)?, w)?
            .sub(&binary(source.spec().plain_length_mm, w)?, w)?;
    }
    let mut controls = Vec::with_capacity(12);
    for i in 0usize..=11 {
        let mut xy = [N::integer(0), N::integer(0)];
        for (j, r, scale) in [(i, &ra, 11 - i), (i.wrapping_sub(1), &rb, i)] {
            if scale == 0 {
                continue;
            }
            let j = if left { 10 - j } else { j };
            for (out, c) in xy.iter_mut().zip(&coefficients[j]) {
                w.charge(512)?;
                *out = out.add(
                    &r.mul(&N(c.clone()), w)?.mul(&N::integer(scale as i64), w)?,
                    w,
                )?;
            }
        }
        controls.push([
            xy[0].clone(),
            xy[1].clone(),
            z.mul(&N::integer(11), w)?,
            N::integer(11),
        ]);
    }
    let trim = [
        lo.sub(&a, w)?.div(&width, w)?,
        hi.sub(&a, w)?.div(&width, w)?,
    ];
    let du = ub.sub(&ua, w)?;
    Ok(Some(AffineCap {
        controls,
        trim,
        uv: [[ua, a], [ub, b]],
        derivative: [du, width],
        reversed: true,
    }))
}

pub struct AffineCapPiece<'f, 's> {
    field: &'f EncodedField<'s>,
    cell: usize,
    end: End,
    piece: usize,
    controls: Vec<[EncodedScalar; 4]>,
    trim: [EncodedScalar; 2],
    uv: [[EncodedScalar; 2]; 2],
    reversed: bool,
    derivative: [EncodedScalar; 2],
    work: usize,
    bytes: usize,
}
impl AffineCapPiece<'_, '_> {
    pub fn field(&self) -> &EncodedField<'_> {
        self.field
    }
    pub fn source_cell(&self) -> usize {
        self.cell
    }
    pub fn end(&self) -> End {
        self.end
    }
    pub fn angular_piece(&self) -> usize {
        self.piece
    }
    pub fn controls(&self) -> &[[EncodedScalar; 4]] {
        &self.controls
    }
    pub fn trim(&self) -> &[EncodedScalar; 2] {
        &self.trim
    }
    pub fn surface_endpoints(&self) -> &[[EncodedScalar; 2]; 2] {
        &self.uv
    }
    pub fn chart_derivative(&self) -> &[EncodedScalar; 2] {
        &self.derivative
    }
    pub fn original_phase_reversed(&self) -> bool {
        self.reversed
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_bytes(&self) -> usize {
        self.bytes
    }
}
impl<'s> EncodedField<'s> {
    pub fn affine_cap_piece(
        &self,
        cell: usize,
        end: End,
        piece: usize,
        options: DomainOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Option<AffineCapPiece<'_, 's>>, Error> {
        let (mut w, mut m) = crate::step::hanger::domain::ledgers(options, self.bytes, interrupt)?;
        w.charge(1)?;
        if cell >= self.cells.len() || piece >= self.angular.coefficients().len() {
            return Err(Error::Invalid);
        }
        // Conservative live exact-arithmetic scratch for this sealed degree11
        // construction. It remains charged while retained coefficients are built.
        // Scalar formula depth in build is bounded: source differences/divisions,
        // one affine radius interpolation, then one degree10 angular product and
        // multiplication by integer11. 64*B+2048 bits bounds numerator/denominator
        // intermediates; <=512 live BigInts covers48 output rationals, locals and
        // arithmetic scratch. Eightfold limb capacity bounds the pinned bigint
        // allocation strategy. Binary64 source scalars need at most1076 bits.
        let mut bits = 1076_usize;
        for row in &self.meridians {
            for scalar in row.radius.iter().chain(&row.phase).chain(&row.weights) {
                w.charge(1 + scalar.n.len() + scalar.d.len())?;
                bits = bits.max(8 * (scalar.n.len() + scalar.d.len()));
            }
        }
        for piece in self.angular.coefficients() {
            for row in piece {
                for scalar in row {
                    w.charge(
                        1 + (scalar.numer().bits() + scalar.denom().bits()).div_ceil(8) as usize,
                    )?;
                    bits = bits.max((scalar.numer().bits() + scalar.denom().bits()) as usize);
                }
            }
        }
        let scratch = 512 * 128 + 512 * 8 * (64 * bits + 2048).div_ceil(8);
        m.add(scratch)?;
        let Some(raw) = build(self, cell, end, piece, &mut w)? else {
            return Ok(None);
        };
        let before = m.used;
        m.add(size_of::<AffineCapPiece<'_, '_>>())?;
        let mut controls = m.allocate(raw.controls.len(), &mut w)?;
        for row in &raw.controls {
            controls.push([
                EncodedScalar::retain(&row[0], &mut m, &mut w)?,
                EncodedScalar::retain(&row[1], &mut m, &mut w)?,
                EncodedScalar::retain(&row[2], &mut m, &mut w)?,
                EncodedScalar::retain(&row[3], &mut m, &mut w)?,
            ]);
        }
        let mut retain = |x: &N| EncodedScalar::retain(x, &mut m, &mut w);
        let trim = [retain(&raw.trim[0])?, retain(&raw.trim[1])?];
        let uv = [
            [retain(&raw.uv[0][0])?, retain(&raw.uv[0][1])?],
            [retain(&raw.uv[1][0])?, retain(&raw.uv[1][1])?],
        ];
        let derivative = [retain(&raw.derivative[0])?, retain(&raw.derivative[1])?];
        w.charge(1)?;
        Ok(Some(AffineCapPiece {
            field: self,
            cell,
            end,
            piece,
            controls,
            trim,
            uv,
            reversed: raw.reversed,
            derivative,
            work: w.used,
            bytes: m.used - before,
        }))
    }
}
