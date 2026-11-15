use super::*;
pub(super) struct StockEnvelope {
    columns: [[Q; 3]; 3],
    points: Vec<[Q; 3]>,
    caches: Vec<[I; 3]>,
    spans: [Q; 3],
    jacobian: Q,
}
fn cache(value: &Q) -> R<I> {
    let v = value.to_f64().ok_or(Error::Unsupported)?;
    if q(v)? == *value {
        Ok(ip(v))
    } else {
        enclosure(value)
    }
}
fn identity() -> Frame {
    Frame {
        origin_mm: [0.; 3],
        axes: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
    }
}
fn linear(columns: &[[Q; 3]; 3], p: &[Q; 3], work: &mut Work<'_>) -> R<[Q; 3]> {
    work.charge(24)?;
    let mut result = std::array::from_fn(|_| Q::zero());
    for i in 0..3 {
        result[i] = bounded((0..3).map(|j| &columns[j][i] * &p[j]).sum())?;
    }
    Ok(result)
}
fn apply(origin: &[Q; 3], columns: &[[Q; 3]; 3], p: &[Q; 3], work: &mut Work<'_>) -> R<[Q; 3]> {
    let mut result = linear(columns, p, work)?;
    work.charge(3)?;
    for i in 0..3 {
        result[i] = bounded(&result[i] + &origin[i])?;
    }
    Ok(result)
}
fn determinant(frame: Frame) -> R<Q> {
    let a = vector(frame.axes[0])?;
    let b = vector(frame.axes[1])?;
    let c = vector(frame.axes[2])?;
    dot(&a, &cross(&b, &c)?)
}
impl StockEnvelope {
    pub fn new(
        nut: ClockedLine8Placement<'_>,
        receiver: ReceivingSweep<'_>,
        work: &mut Work<'_>,
    ) -> R<Self> {
        work.charge(2048)?;
        scale(nut.stock_frame)?;
        scale(nut.outer_frame)?;
        scale(receiver.frame)?;
        scale(receiver.occurrence)?;
        // Reuse the established literal-frame inverse. Every original operator
        // remains separate: receiver base^-1 * occurrence^-1 * nut outer * stock.
        // No rounded replacement origin/axes are manufactured.
        let (a, ac) = affine(receiver.occurrence, nut.outer_frame, work.interrupt)?;
        let (b, bc) = affine(receiver.frame, identity(), work.interrupt)?;
        let p = apply(&a, &ac, &vector(nut.stock_frame.origin_mm)?, work)?;
        let origin = apply(&b, &bc, &p, work)?;
        let mut columns = std::array::from_fn(|_| std::array::from_fn(|_| Q::zero()));
        for (j, column) in columns.iter_mut().enumerate() {
            let p = linear(&ac, &vector(nut.stock_frame.axes[j])?, work)?;
            *column = linear(&bc, &p, work)?;
        }
        let stock = nut.source.stock();
        let ring = &stock.base()[0];
        let mut points = Vec::with_capacity(ring.len() * 2);
        for x in [0., stock.height_mm()] {
            for yz in ring {
                points.push(apply(
                    &origin,
                    &columns,
                    &[q(x)?, q(yz[0])?, q(yz[1])?],
                    work,
                )?);
            }
        }
        let mut lo = [f64::INFINITY; 2];
        let mut hi = [f64::NEG_INFINITY; 2];
        for yz in ring {
            for i in 0..2 {
                lo[i] = lo[i].min(yz[i]);
                hi[i] = hi[i].max(yz[i]);
            }
        }
        let spans = [
            q(stock.height_mm())?,
            bounded(q(hi[0])? - q(lo[0])?)?,
            bounded(q(hi[1])? - q(lo[1])?)?,
        ];
        let jacobian = bounded(determinant(nut.stock_frame)? * determinant(nut.outer_frame)?)?;
        if jacobian <= Q::zero() {
            return Err(Error::Unsupported);
        }
        let mut caches = Vec::with_capacity(points.len());
        for p in &points {
            work.charge(48)?;
            caches.push([cache(&p[0])?, cache(&p[1])?, cache(&p[2])?]);
        }
        Ok(Self {
            columns,
            caches,
            points,
            spans,
            jacobian,
        })
    }
    pub fn source_box(
        frame: Frame,
        bounds: [[Q; 2]; 3],
        receiver: ReceivingSweep<'_>,
        work: &mut Work<'_>,
    ) -> R<Self> {
        work.charge(2048)?;
        scale(frame)?;
        scale(receiver.frame)?;
        scale(receiver.occurrence)?;
        let (a, ac) = affine(receiver.occurrence, frame, work.interrupt)?;
        let (b, bc) = affine(receiver.frame, identity(), work.interrupt)?;
        let origin = apply(&b, &bc, &a, work)?;
        let mut columns = std::array::from_fn(|_| std::array::from_fn(|_| Q::zero()));
        for j in 0..3 {
            columns[j] = linear(&bc, &ac[j], work)?;
        }
        let mut points = Vec::with_capacity(8);
        let mut caches = Vec::with_capacity(8);
        for bits in 0..8 {
            let p = std::array::from_fn(|i| bounds[i][(bits >> i) & 1].clone());
            let p = apply(&origin, &columns, &p, work)?;
            work.charge(48)?;
            caches.push([cache(&p[0])?, cache(&p[1])?, cache(&p[2])?]);
            points.push(p);
        }
        let spans = [
            bounded(&bounds[0][1] - &bounds[0][0])?,
            bounded(&bounds[1][1] - &bounds[1][0])?,
            bounded(&bounds[2][1] - &bounds[2][0])?,
        ];
        let jacobian = determinant(frame)?;
        if spans.iter().any(|x| x < &Q::zero()) || jacobian <= Q::zero() {
            return Err(Error::Unsupported);
        }
        Ok(Self {
            columns,
            points,
            caches,
            spans,
            jacobian,
        })
    }
    pub fn projected_ranges(&self, work: &mut Work<'_>) -> R<[[Q; 2]; 2]> {
        work.charge(self.points.len() * 4)?;
        Ok(std::array::from_fn(|i| {
            [
                self.points.iter().map(|p| p[i].clone()).min().unwrap(),
                self.points.iter().map(|p| p[i].clone()).max().unwrap(),
            ]
        }))
    }
    fn volume_box(&self) -> R<Q> {
        bounded(&self.spans[0] * &self.spans[1] * &self.spans[2] * &self.jacobian)
    }
    /// Exact upper volume of original stock in plane(x)>0. Projected values
    /// have no metric conversion: the source-coordinate slab Jacobian is exact.
    fn slab(&self, plane: &[Q; 4], work: &mut Work<'_>) -> R<Q> {
        let mut maximum = Q::zero();
        work.charge(64)?;
        let coefficients = [
            cache(&plane[0])?,
            cache(&plane[1])?,
            cache(&plane[2])?,
            cache(&plane[3])?,
        ];
        for (p, cache) in self.points.iter().zip(&self.caches) {
            work.charge(16)?;
            // Only an earned outward nonpositive enclosure skips exact work.
            // Near contact the original rational values still own the result.
            let upper = (0..3).fold(coefficients[3], |sum, i| sum + coefficients[i] * cache[i]);
            if finite(upper)?.hi <= 0. {
                continue;
            }
            work.charge(12)?;
            maximum = maximum.max(bounded(
                (0..3).map(|i| &plane[i] * &p[i]).sum::<Q>() + &plane[3],
            )?);
        }
        if maximum.is_zero() {
            return Ok(Q::zero());
        }
        let mut axis = 0;
        let mut coefficient = Q::zero();
        for j in 0..3 {
            work.charge(12)?;
            let c = bounded((0..3).map(|i| &plane[i] * &self.columns[j][i]).sum::<Q>())?.abs();
            if c > coefficient {
                coefficient = c;
                axis = j;
            }
        }
        let full = self.volume_box()?;
        if coefficient.is_zero() {
            return Ok(full);
        }
        work.charge(16)?;
        let area = bounded(&self.spans[(axis + 1) % 3] * &self.spans[(axis + 2) % 3])?;
        Ok(bounded(maximum / coefficient * area * &self.jacobian)?.min(full))
    }
    pub fn outer_cover_bound(&self, receiver: ReceivingSweep<'_>, work: &mut Work<'_>) -> R<Q> {
        let mut bounds: Option<[[Q; 2]; 2]> = None;
        for edge in &receiver.source.rings()[0] {
            let next = cells::edge_bounds(*edge, work)?;
            if let Some(old) = &mut bounds {
                for i in 0..2 {
                    old[i][0] = old[i][0].clone().min(next[i][0].clone());
                    old[i][1] = old[i][1].clone().max(next[i][1].clone());
                }
            } else {
                bounds = Some(next);
            }
        }
        let xy = bounds.ok_or(Error::Unsupported)?;
        let h = q(receiver.source.height_mm())?;
        let z = if receiver.reverse {
            [-h, Q::zero()]
        } else {
            [Q::zero(), h]
        };
        let ranges = [xy[0].clone(), xy[1].clone(), z];
        let mut best = self.volume_box()?;
        for axis in 0..3 {
            for upper in [false, true] {
                let mut plane = std::array::from_fn(|_| Q::zero());
                plane[axis] = if upper {
                    -Q::from_integer(1.into())
                } else {
                    Q::from_integer(1.into())
                };
                plane[3] = if upper {
                    ranges[axis][1].clone()
                } else {
                    -&ranges[axis][0]
                };
                best = best.min(self.slab(&plane, work)?);
            }
        }
        Ok(best)
    }
    pub fn cell_bound(&self, cell: &Cell, work: &mut Work<'_>) -> R<Q> {
        let mut bound = Q::zero();
        for p in &cell.planes {
            bound = bounded(
                bound + self.slab(&[p[0].clone(), p[1].clone(), Q::zero(), p[2].clone()], work)?,
            )?;
        }
        Ok(bound.min(self.volume_box()?))
    }
}
