use super::restrict::amplification;
use super::*;
pub(super) fn eval_interval<const D: usize>(
    input: &[[I; D]],
    t: I,
    w: &mut Work<'_>,
) -> Result<[I; D], Error> {
    w.charge(32 * D * input.len() * (input.len() - 1) / 2 + input.len() * D)?;
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
fn norm(c: &[N], w: &mut Work<'_>) -> Result<I, Error> {
    let mut s = p(0.);
    for x in c {
        s = s + p(x.enclosure(w)?.magnitude());
    }
    Ok(s)
}
pub(super) fn field_bounds(
    f: &EncodedField<'_>,
    index: usize,
    end: End,
    angular_piece: usize,
    vrange: I,
    w: &mut Work<'_>,
) -> Result<[f64; 2], Error> {
    let cell = f.cells[index];
    let r = &f.meridians[cell.meridian];
    let source = f.source.source_boundary().source();
    let (pitch, starts, left, length) = if r.component == Component::Metric {
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
    let shift = pitch
        .mul(&N::integer(i64::from(cell.source_key.total_period)), w)?
        .sub(&station, w)?;
    let a = binary(f.angular.knots()[angular_piece], w)?;
    let b = binary(f.angular.knots()[angular_piece + 1], w)?;
    let width = b.sub(&a, w)?;
    let t = (vrange - a.enclosure(w)?) / width.enclosure(w)?;
    let factors = [
        amplification(t, 10, w)?,
        amplification(t, 9, w)?,
        amplification(t, 8, w)?,
    ];
    let ap = if left {
        f.angular.coefficients().len() - 1 - angular_piece
    } else {
        angular_piece
    };
    let mut net = Vec::new();
    let mut weights = Vec::new();
    for i in 0..=r.degree() {
        let radius = r.radius[i].number(w)?;
        let q = r.phase[i].number(w)?;
        let weight = r.weights[i].number(w)?;
        let mut row = Vec::new();
        for j in 0..=10 {
            let c = &f.angular.coefficients()[ap][if left { 10 - j } else { j }];
            let v = a.add(
                &width
                    .mul(&N::integer(j as i64), w)?
                    .div(&N::integer(10), w)?,
                w,
            )?;
            w.charge(512)?;
            row.push([
                radius.mul(&N(c[0].clone()), w)?,
                radius.mul(&N(c[1].clone()), w)?,
                q.add(&weight.mul(&shift.add(&lead.mul(&v, w)?, w)?, w)?, w)?,
            ]);
        }
        net.push(row);
        weights.push(weight);
    }
    let du = r.degree();
    let mut raw = [0_f64; 6];
    for row in &net {
        for v in row {
            raw[0] = raw[0].max(norm(v, w)?.hi);
        }
    }
    for rows in net.windows(2) {
        for (upper, lower) in rows[1].iter().zip(&rows[0]) {
            let d = std::array::from_fn::<_, 3, _>(|_| N::integer(0));
            let mut d = d;
            for (axis, out) in d.iter_mut().enumerate() {
                *out = upper[axis]
                    .sub(&lower[axis], w)?
                    .mul(&N::integer(du as i64), w)?;
            }
            raw[1] = raw[1].max(norm(&d, w)?.hi);
        }
    }
    for row in &net {
        for p in row.windows(2) {
            let mut d: [N; 3] = std::array::from_fn(|_| N::integer(0));
            for axis in 0..3 {
                d[axis] = p[1][axis]
                    .sub(&p[0][axis], w)?
                    .mul(&N::integer(10), w)?
                    .div(&width, w)?;
            }
            raw[2] = raw[2].max(norm(&d, w)?.hi);
        }
    }
    for rows in net.windows(3) {
        for ((upper, middle), lower) in rows[2].iter().zip(&rows[1]).zip(&rows[0]) {
            let mut d: [N; 3] = std::array::from_fn(|_| N::integer(0));
            for axis in 0..3 {
                d[axis] = upper[axis]
                    .sub(&middle[axis].mul(&N::integer(2), w)?, w)?
                    .add(&lower[axis], w)?
                    .mul(&N::integer((du * (du - 1)) as i64), w)?;
            }
            raw[3] = raw[3].max(norm(&d, w)?.hi);
        }
    }
    for rows in net.windows(2) {
        for j in 0..10 {
            let mut d: [N; 3] = std::array::from_fn(|_| N::integer(0));
            for axis in 0..3 {
                d[axis] = rows[1][j + 1][axis]
                    .sub(&rows[1][j][axis], w)?
                    .sub(&rows[0][j + 1][axis].sub(&rows[0][j][axis], w)?, w)?
                    .mul(&N::integer((du * 10) as i64), w)?
                    .div(&width, w)?;
            }
            raw[4] = raw[4].max(norm(&d, w)?.hi);
        }
    }
    let square = width.mul(&width, w)?;
    for row in &net {
        for p in row.windows(3) {
            let mut d: [N; 3] = std::array::from_fn(|_| N::integer(0));
            for axis in 0..3 {
                d[axis] = p[2][axis]
                    .sub(&p[1][axis].mul(&N::integer(2), w)?, w)?
                    .add(&p[0][axis], w)?
                    .mul(&N::integer(90), w)?
                    .div(&square, w)?;
            }
            raw[5] = raw[5].max(norm(&d, w)?.hi);
        }
    }
    let mut minimum = f64::INFINITY;
    for x in &weights {
        minimum = minimum.min(x.enclosure(w)?.lo);
    }
    let mut wu = 0_f64;
    for p in weights.windows(2) {
        wu = wu.max(
            p[1].sub(&p[0], w)?
                .mul(&N::integer(du as i64), w)?
                .enclosure(w)?
                .magnitude(),
        );
    }
    let mut wuu = 0_f64;
    for p in weights.windows(3) {
        wuu = wuu.max(
            p[2].sub(&p[1].mul(&N::integer(2), w)?, w)?
                .add(&p[0], w)?
                .mul(&N::integer((du * (du - 1)) as i64), w)?
                .enclosure(w)?
                .magnitude(),
        );
    }
    w.charge(1024)?;
    let a = p(raw[0]) * factors[0];
    let au = p(raw[1]) * factors[0];
    let av = p(raw[2]) * factors[1];
    let auu = p(raw[3]) * factors[0];
    let auv = p(raw[4]) * factors[1];
    let avv = p(raw[5]) * factors[2];
    let lo = p(minimum);
    let fu = au / lo + a * p(wu) / lo.square();
    let fv = av / lo;
    let fuu = auu / lo
        + p(2.) * au * p(wu) / lo.square()
        + a * p(wuu) / lo.square()
        + p(2.) * a * p(wu).square() / (lo.square() * lo);
    let fuv = auv / lo + av * p(wu) / lo.square();
    let fvv = avv / lo;
    Ok([fu.hi.max(fv.hi), fuu.hi.max(fuv.hi).max(fvv.hi)])
}
