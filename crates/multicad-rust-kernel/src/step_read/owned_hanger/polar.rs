//! Exact polar factorization of decoded supports. Finite occupied material
//! still requires the independently checked cap/side domain transfer.
use super::*;
use num_bigint::{BigInt, Sign};
use num_rational::BigRational;
use num_traits::{One, ToPrimitive, Zero};
#[derive(Clone, Copy)]
pub struct PolarOptions {
    pub maximum_work: usize,
    pub maximum_bytes: usize,
}
impl Default for PolarOptions {
    fn default() -> Self {
        Self {
            maximum_work: 64_000_000,
            maximum_bytes: 64 * 1024 * 1024,
        }
    }
}
pub struct PolarStrip {
    pub(super) tensor: usize,
    pub(super) radius: Vec<I>,
    pub(super) phase: Vec<I>,
    pub(super) weights: Vec<I>,
    pub(super) phase_start: I,
    pub(super) phase_width: I,
    pub(super) lower: usize,
    pub(super) left: bool,
    pub(super) r_integer: Vec<BigInt>,
    pub(super) q_integer: Vec<BigInt>,
}
impl PolarStrip {
    pub fn tensor_index(&self) -> usize {
        self.tensor
    }
    pub fn radius_controls_mm(&self) -> Vec<[f64; 2]> {
        self.radius.iter().map(|v| [v.lo, v.hi]).collect()
    }
    pub fn phase_controls_mm(&self) -> Vec<[f64; 2]> {
        self.phase.iter().map(|v| [v.lo, v.hi]).collect()
    }
    pub fn phase_start_mm(&self) -> [f64; 2] {
        [self.phase_start.lo, self.phase_start.hi]
    }
    pub fn phase_width_mm(&self) -> [f64; 2] {
        [self.phase_width.lo, self.phase_width.hi]
    }
    pub fn lower_plane_index(&self) -> usize {
        self.lower
    }
    pub fn left_hand(&self) -> bool {
        self.left
    }
    pub fn weight_bounds(&self) -> Vec<[f64; 2]> {
        self.weights.iter().map(|v| [v.lo, v.hi]).collect()
    }
}
pub(super) struct AngularPiece {
    pub(super) knots: [BigRational; 2],
    pub(super) points: [[BigInt; 2]; 11],
    pub(super) bounds: [[I; 2]; 11],
}
pub struct PolarSupports<'g, 'd, 'i> {
    pub(super) grid: &'g ExactPoleGrid<'d, 'i>,
    pub(super) planes: Vec<(u32, BigInt)>,
    pub(super) angular: Vec<AngularPiece>,
    pub(super) plain_radius: BigInt,
    pub(super) strips: Vec<PolarStrip>,
    pub(super) leads: [BigInt; 2],
    pub(super) origins: [BigInt; 2],
    work: usize,
    bytes: usize,
}
fn p(x: f64) -> I {
    I::point(x)
}
fn integer(value: &PoleInteger, b: &mut Budget<'_>) -> Result<BigInt, Error> {
    b.charge(value.numerator_magnitude_le().len() + 1)?;
    Ok(BigInt::from_bytes_le(
        match value.signum() {
            -1 => Sign::Minus,
            0 => Sign::NoSign,
            _ => Sign::Plus,
        },
        value.numerator_magnitude_le(),
    ))
}
fn xyz(grid: &ExactPoleGrid<'_, '_>, id: u32, b: &mut Budget<'_>) -> Result<[BigInt; 3], Error> {
    b.charge(64 * (1 + (grid.points().len() + 1).ilog2() as usize))?;
    let q = grid.points().get(&id).ok_or(Error::Invalid)?;
    Ok([integer(&q[0], b)?, integer(&q[1], b)?, integer(&q[2], b)?])
}
fn multiply(a: &BigInt, c: &BigInt, b: &mut Budget<'_>) -> Result<BigInt, Error> {
    let n = (a.bits() as usize).div_ceil(64);
    let m = (c.bits() as usize).div_ceil(64);
    b.charge(16 * (1 + n * m))?;
    Ok(a * c)
}
fn interval(n: &BigInt, den: &BigInt, b: &mut Budget<'_>) -> Result<I, Error> {
    let words = (n.bits().max(den.bits()) as usize).div_ceil(64);
    b.charge(128 * (1 + words * words))?;
    let q = BigRational::new(n.clone(), den.clone());
    let value = q.to_f64().ok_or(Error::Unsupported)?;
    if !value.is_finite() {
        return Err(Error::Unsupported);
    }
    let cached = BigRational::from_float(value).ok_or(Error::Unsupported)?;
    Ok(match q.cmp(&cached) {
        std::cmp::Ordering::Equal => p(value),
        std::cmp::Ordering::Less => I::new(value.next_down(), value),
        std::cmp::Ordering::Greater => I::new(value, value.next_up()),
    })
}
fn face_for(d: &Decoded<'_>, surface: u32, b: &mut Budget<'_>) -> Result<usize, Error> {
    let mut found = None;
    for (i, f) in d.faces.iter().enumerate() {
        b.charge(1)?;
        if f.surface == surface {
            if found.is_some() {
                return Err(Error::Unsupported);
            }
            found = Some(i);
        }
    }
    found.ok_or(Error::Unsupported)
}
fn cap_pair(
    d: &Decoded<'_>,
    surface: u32,
    planes: &[(u32, BigInt)],
    b: &mut Budget<'_>,
) -> Result<usize, Error> {
    let mut owners = BTreeSet::new();
    for edge in d.edges.values() {
        b.charge(1)?;
        if edge.pcurves.iter().any(|pc| pc.surface == surface) {
            for pc in edge.pcurves {
                if let Some(i) = planes.iter().position(|(id, _)| *id == pc.surface) {
                    owners.insert(i);
                }
            }
        }
    }
    if owners.len() != 2 {
        return Err(Error::Unsupported);
    }
    let mut iter = owners.into_iter();
    let lower = iter.next().ok_or(Error::Invalid)?;
    let upper = iter.next().ok_or(Error::Invalid)?;
    if upper != lower + 1 {
        return Err(Error::Unsupported);
    }
    Ok(lower)
}
fn storage_bound(grid: &ExactPoleGrid<'_, '_>, b: &mut Budget<'_>) -> Result<usize, Error> {
    let d = grid.source();
    let mut bytes = grid.denominator().numerator_magnitude_le().len();
    for point in grid.points().values() {
        b.charge(16)?;
        for coordinate in point {
            bytes = bytes.max(coordinate.numerator_magnitude_le().len());
        }
    }
    let mut kb = 1_usize;
    for t in &d.tensors {
        for k in &t.knots_v {
            b.charge(k.text.len() + 1)?;
            let exponent = k
                .text
                .split_once('E')
                .map_or("0", |(_, e)| e)
                .parse::<i32>()
                .map_err(|_| Error::Unsupported)?;
            let bits = k
                .text
                .len()
                .checked_add(exponent.unsigned_abs() as usize)
                .and_then(|v| v.checked_mul(4))
                .ok_or(Error::Limit)?;
            kb = kb.max(bits);
        }
    }
    let coefficient = 8_usize.checked_mul(32 + bytes).ok_or(Error::Limit)?;
    let derived_bits = bytes
        .checked_mul(8)
        .and_then(|v| v.checked_add(kb.checked_mul(2)?))
        .and_then(|v| v.checked_add(64))
        .ok_or(Error::Limit)?;
    let derived = 8_usize
        .checked_mul(32 + derived_bits.div_ceil(8))
        .ok_or(Error::Limit)?;
    let mut bound = coefficient.checked_mul(32).ok_or(Error::Limit)?;
    for t in &d.tensors {
        let field = 1024_usize
            .checked_add(
                coefficient
                    .checked_mul(t.dimensions[0] * 2)
                    .ok_or(Error::Limit)?,
            )
            .and_then(|v| v.checked_add(derived.checked_mul(6)?))
            .ok_or(Error::Limit)?;
        bound = bound.checked_add(field).ok_or(Error::Limit)?;
        let a = grid.points().get(&t.points[0]).ok_or(Error::Invalid)?;
        let z = grid.points().get(&t.points[10]).ok_or(Error::Invalid)?;
        b.charge(a[2].numerator_magnitude_le().len() + z[2].numerator_magnitude_le().len() + 1)?;
        if a[2].signum() == z[2].signum()
            && a[2].numerator_magnitude_le() == z[2].numerator_magnitude_le()
        {
            let piece = size_of::<AngularPiece>() + 22 * coefficient + 4 * derived;
            bound = bound
                .checked_add(
                    (t.knots_v.len() - 1)
                        .checked_mul(piece)
                        .ok_or(Error::Limit)?,
                )
                .ok_or(Error::Limit)?;
        }
    }
    Ok(bound)
}
impl<'g, 'd, 'i> PolarSupports<'g, 'd, 'i> {
    pub fn prepare(
        grid: &'g ExactPoleGrid<'d, 'i>,
        options: PolarOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, Error> {
        if options.maximum_work == 0
            || options.maximum_work > 128_000_000
            || options.maximum_bytes == 0
            || options.maximum_bytes > 64 * 1024 * 1024
        {
            return Err(Error::Invalid);
        }
        let d = grid.source();
        let mut b = Budget {
            work: 0,
            remaining: options.maximum_work,
            bytes: grid.retained_bytes(),
            limit: options.maximum_bytes,
            interrupt,
        };
        b.charge(1)?;
        let retained = storage_bound(grid, &mut b)?;
        b.retain(8 * 1024 * 1024 + retained)?;
        if d.planes.len() != 4 {
            return Err(Error::Unsupported);
        }
        let unit = integer(grid.denominator(), &mut b)?;
        let mut planes = Vec::new();
        for plane in &d.planes {
            let q = xyz(grid, plane.origin, &mut b)?;
            if !q[0].is_zero() || !q[1].is_zero() {
                return Err(Error::Unsupported);
            }
            planes.push((plane.entity, q[2].clone()));
        }
        for i in 1..planes.len() {
            let mut j = i;
            while j > 0 {
                b.charge(32)?;
                if planes[j - 1].1 < planes[j].1 {
                    break;
                }
                planes.swap(j - 1, j);
                j -= 1;
            }
        }
        if planes.windows(2).any(|p| p[0].1 >= p[1].1) {
            return Err(Error::Unsupported);
        }
        let mut raw = Vec::new();
        let mut quarters: BTreeMap<usize, (usize, BigInt)> = BTreeMap::new();
        let mut plain_radius = None;
        for (index, t) in d.tensors.iter().enumerate() {
            let face = face_for(d, t.entity, &mut b)?;
            if d.faces[face].loops.len() != 1 {
                return Err(Error::Unsupported);
            }
            let lower = cap_pair(d, t.entity, &planes, &mut b)?;
            for row in 0..t.dimensions[0] {
                for j in 1..t.dimensions[1] {
                    if numbers::compare(
                        t.weights[row * t.dimensions[1]],
                        t.weights[row * t.dimensions[1] + j],
                        &mut b,
                    )? != std::cmp::Ordering::Equal
                    {
                        return Err(Error::Unsupported);
                    }
                }
            }
            let first = xyz(grid, t.points[0], &mut b)?;
            let end = xyz(grid, t.points[10], &mut b)?;
            let ka = numbers::rational(t.knots_v[0], &mut b)?;
            let kb = numbers::rational(t.knots_v[1], &mut b)?;
            let width = &kb - &ka;
            let lead = BigRational::from_integer(&end[2] - &first[2]) / width;
            if !lead.is_integer() {
                return Err(Error::Unsupported);
            }
            let lead = lead.to_integer();
            if lead.is_zero() {
                if lower != 1 || t.degrees[0] != 1 || d.faces[face].same_sense {
                    return Err(Error::Unsupported);
                }
                let quarter = &ka * BigRational::from_integer(4.into());
                if !quarter.is_integer() {
                    return Err(Error::Unsupported);
                }
                let quarter = quarter.to_integer().to_usize().ok_or(Error::Unsupported)?;
                if quarter >= 4
                    || numbers::rational(*t.knots_v.last().ok_or(Error::Invalid)?, &mut b)? - &ka
                        != BigRational::new(1.into(), 4.into())
                {
                    return Err(Error::Unsupported);
                }
                let radius = match quarter {
                    0 => first[0].clone(),
                    1 => first[1].clone(),
                    2 => -&first[0],
                    _ => -&first[1],
                };
                if radius <= BigInt::zero() {
                    return Err(Error::Unsupported);
                }
                if plain_radius.as_ref().is_some_and(|r| *r != radius) {
                    return Err(Error::Unsupported);
                }
                plain_radius = Some(radius.clone());
                for j in 0..t.dimensions[1] {
                    let a = xyz(grid, t.points[j], &mut b)?;
                    let z = xyz(grid, t.points[t.dimensions[1] + j], &mut b)?;
                    if a[0] != z[0] || a[1] != z[1] || a[2] != planes[1].1 || z[2] != planes[2].1 {
                        return Err(Error::Unsupported);
                    }
                }
                if quarters.insert(quarter, (index, radius)).is_some() {
                    return Err(Error::Unsupported);
                }
            } else {
                if lead < BigInt::zero() || lower == 1 || !ka.is_zero() {
                    return Err(Error::Unsupported);
                }
                let mut radius = Vec::new();
                let mut phase = Vec::new();
                for row in 0..t.dimensions[0] {
                    let first = xyz(grid, t.points[row * t.dimensions[1]], &mut b)?;
                    if first[0] <= BigInt::zero() || !first[1].is_zero() {
                        return Err(Error::Unsupported);
                    }
                    radius.push(first[0].clone());
                    phase.push(first[2].clone());
                }
                if phase.windows(2).any(|q| q[0] >= q[1]) {
                    return Err(Error::Unsupported);
                }
                raw.push((index, lower, lead, radius, phase));
            }
        }
        if quarters.len() != 4 {
            return Err(Error::Unsupported);
        }
        let plain_radius = plain_radius.ok_or(Error::Unsupported)?;
        let mut angular = Vec::new();
        let cardinals = [[1, 0], [0, 1], [-1, 0], [0, -1], [1, 0]];
        for quarter in 0..4 {
            let (index, _) = &quarters[&quarter];
            let t = &d.tensors[*index];
            for (pole, cardinal) in [
                (0, cardinals[quarter]),
                (t.dimensions[1] - 1, cardinals[quarter + 1]),
            ] {
                let q = xyz(grid, t.points[pole], &mut b)?;
                for axis in 0..2 {
                    if q[axis] != &plain_radius * BigInt::from(cardinal[axis]) {
                        return Err(Error::Unsupported);
                    }
                }
            }
            for span in 0..t.knots_v.len() - 1 {
                let mut points =
                    std::array::from_fn::<_, 11, _>(|_| std::array::from_fn(|_| BigInt::zero()));
                let mut bounds = [[p(0.); 2]; 11];
                for j in 0..=10 {
                    let q = xyz(grid, t.points[span * 10 + j], &mut b)?;
                    for axis in 0..2 {
                        bounds[j][axis] = interval(&q[axis], &plain_radius, &mut b)?;
                        points[j][axis] = q[axis].clone();
                    }
                }
                angular.push(AngularPiece {
                    knots: [
                        numbers::rational(t.knots_v[span], &mut b)?,
                        numbers::rational(t.knots_v[span + 1], &mut b)?,
                    ],
                    points,
                    bounds,
                });
            }
        }
        let mut strips = Vec::new();
        let mut leads: [Option<BigInt>; 2] = [None, None];
        let mut origins: [Option<BigInt>; 2] = [None, None];
        let mut hands = [None, None];
        for (index, lower, lead, radius, phase) in raw {
            let group = if lower == 0 { 0 } else { 1 };
            if leads[group].as_ref().is_some_and(|v| *v != lead) {
                return Err(Error::Unsupported);
            }
            leads[group] = Some(lead.clone());
            if origins[group].as_ref().is_none_or(|v| phase[0] < *v) {
                origins[group] = Some(phase[0].clone());
            }
            let t = &d.tensors[index];
            let first_y = xyz(grid, t.points[1], &mut b)?[1].clone();
            if first_y.is_zero() {
                return Err(Error::Unsupported);
            }
            let left = first_y < BigInt::zero();
            if hands[group].is_some_and(|h| h != left) {
                return Err(Error::Unsupported);
            }
            hands[group] = Some(left);
            let face = face_for(d, t.entity, &mut b)?;
            if d.faces[face].same_sense != left {
                return Err(Error::Unsupported);
            }
            for (span, knots) in t.knots_v.windows(2).enumerate() {
                let a = numbers::rational(knots[0], &mut b)?;
                let z = numbers::rational(knots[1], &mut b)?;
                let period = a.to_integer();
                let fa = &a - BigRational::from_integer(period.clone());
                let fz = &z - BigRational::from_integer(period);
                let pair = if left {
                    [BigRational::one() - &fz, BigRational::one() - &fa]
                } else {
                    [fa.clone(), fz.clone()]
                };
                let piece = angular
                    .iter()
                    .find(|p| p.knots == pair)
                    .ok_or(Error::Unsupported)?;
                for row in 0..t.dimensions[0] {
                    for j in 0..=10 {
                        let q = xyz(
                            grid,
                            t.points[row * t.dimensions[1] + span * 10 + j],
                            &mut b,
                        )?;
                        let angular_point = &piece.points[if left { 10 - j } else { j }];
                        for axis in 0..2 {
                            if multiply(&q[axis], &plain_radius, &mut b)?
                                != multiply(&angular_point[axis], &radius[row], &mut b)?
                            {
                                return Err(Error::Unsupported);
                            }
                        }
                        let parameter =
                            &a + (&z - &a) * BigRational::new((j as i64).into(), 10.into());
                        if multiply(&(&q[2] - &phase[row]), parameter.denom(), &mut b)?
                            != multiply(&lead, parameter.numer(), &mut b)?
                        {
                            return Err(Error::Unsupported);
                        }
                    }
                }
            }
            let last = numbers::rational(*t.knots_v.last().ok_or(Error::Invalid)?, &mut b)?;
            if !last.is_integer()
                || phase.last().ok_or(Error::Invalid)? > &planes[lower].1
                || &phase[0] + &lead * last.to_integer() < planes[lower + 1].1
            {
                return Err(Error::Unsupported);
            }
            let r = radius
                .iter()
                .map(|n| interval(n, &unit, &mut b))
                .collect::<Result<Vec<_>, _>>()?;
            let q = phase
                .iter()
                .map(|n| interval(&(n - &phase[0]), &unit, &mut b))
                .collect::<Result<Vec<_>, _>>()?;
            let max = t.weights.iter().map(|w| w.interval().hi).fold(1., f64::max);
            let weights = (0..t.dimensions[0])
                .map(|row| t.weights[row * t.dimensions[1]].interval() / p(max))
                .collect();
            strips.push(PolarStrip {
                tensor: index,
                radius: r,
                phase: q,
                weights,
                phase_start: p(0.),
                phase_width: p(0.),
                lower,
                left,
                r_integer: radius,
                q_integer: phase,
            });
        }
        let leads = [
            leads[0].take().ok_or(Error::Unsupported)?,
            leads[1].take().ok_or(Error::Unsupported)?,
        ];
        let origins = [
            origins[0].take().ok_or(Error::Unsupported)?,
            origins[1].take().ok_or(Error::Unsupported)?,
        ];
        for group in 0..2 {
            let mut ordered = Vec::new();
            for (i, strip) in strips.iter().enumerate() {
                if strip.lower == if group == 0 { 0 } else { 2 } {
                    let start = (&strip.q_integer[0] - &origins[group]) % &leads[group];
                    let width = strip.q_integer.last().unwrap() - &strip.q_integer[0];
                    ordered.push((i, start, width));
                }
            }
            for i in 1..ordered.len() {
                let mut j = i;
                while j > 0 {
                    b.charge(32)?;
                    if ordered[j - 1].1 < ordered[j].1 {
                        break;
                    }
                    ordered.swap(j - 1, j);
                    j -= 1;
                }
            }
            let mut next = BigInt::zero();
            for k in 0..ordered.len() {
                let (i, start, width) = &ordered[k];
                if *start != next || *width <= BigInt::zero() {
                    return Err(Error::Unsupported);
                }
                next = start + width;
                let following = ordered[(k + 1) % ordered.len()].0;
                if strips[*i].r_integer.last() != strips[following].r_integer.first() {
                    return Err(Error::Unsupported);
                }
                strips[*i].phase_start = interval(start, &unit, &mut b)?;
                strips[*i].phase_width = interval(width, &unit, &mut b)?;
            }
            if next != leads[group] {
                return Err(Error::Unsupported);
            }
        }
        b.bytes -= 8 * 1024 * 1024;
        Ok(Self {
            grid,
            planes,
            angular,
            plain_radius,
            strips,
            leads,
            origins,
            work: b.work,
            bytes: b.bytes,
        })
    }
    pub fn grid(&self) -> &ExactPoleGrid<'d, 'i> {
        self.grid
    }
    pub fn strips(&self) -> &[PolarStrip] {
        &self.strips
    }
    pub fn preparation_work_used(&self) -> usize {
        self.work
    }
    pub fn aggregate_known_work_used(&self) -> usize {
        self.grid.aggregate_known_work_used() + self.work
    }
    pub fn retained_bytes(&self) -> usize {
        self.bytes
    }
    pub fn plane_entities(&self) -> Vec<u32> {
        self.planes.iter().map(|p| p.0).collect()
    }
    pub fn plane_heights_mm(&self) -> Result<Vec<[f64; 2]>, Error> {
        let mut b = Budget {
            work: 0,
            remaining: 16_000_000,
            bytes: 0,
            limit: 64 * 1024 * 1024,
            interrupt: None,
        };
        let unit = integer(self.grid.denominator(), &mut b)?;
        self.planes
            .iter()
            .map(|(_, z)| interval(z, &unit, &mut b).map(|v| [v.lo, v.hi]))
            .collect()
    }
    pub fn angular_coefficients(&self) -> Vec<[[[f64; 2]; 2]; 11]> {
        self.angular
            .iter()
            .map(|p| p.bounds.map(|v| v.map(|x| [x.lo, x.hi])))
            .collect()
    }
    pub fn angular_piece_count(&self) -> usize {
        self.angular.len()
    }
    pub fn plain_radius_mm(&self) -> Result<[f64; 2], Error> {
        let mut b = Budget {
            work: 0,
            remaining: 16_000_000,
            bytes: 0,
            limit: 64 * 1024 * 1024,
            interrupt: None,
        };
        let unit = integer(self.grid.denominator(), &mut b)?;
        interval(&self.plain_radius, &unit, &mut b).map(|v| [v.lo, v.hi])
    }
    pub fn thread_leads_mm(&self) -> Result<[[f64; 2]; 2], Error> {
        let mut b = Budget {
            work: 0,
            remaining: 16_000_000,
            bytes: 0,
            limit: 64 * 1024 * 1024,
            interrupt: None,
        };
        let unit = integer(self.grid.denominator(), &mut b)?;
        let a = interval(&self.leads[0], &unit, &mut b)?;
        let z = interval(&self.leads[1], &unit, &mut b)?;
        Ok([[a.lo, a.hi], [z.lo, z.hi]])
    }
    pub fn phase_origins_mm(&self) -> Result<[[f64; 2]; 2], Error> {
        let mut b = Budget {
            work: 0,
            remaining: 16_000_000,
            bytes: 0,
            limit: 64 * 1024 * 1024,
            interrupt: None,
        };
        let unit = integer(self.grid.denominator(), &mut b)?;
        let a = interval(&self.origins[0], &unit, &mut b)?;
        let z = interval(&self.origins[1], &unit, &mut b)?;
        Ok([[a.lo, a.hi], [z.lo, z.hi]])
    }
}
