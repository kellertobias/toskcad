//! Whole finite Curve/PCurve correspondence derived only from decoded leaves.
use super::bernstein::{self as poly, Algebra, Polynomial};
use super::*;
use std::cmp::Ordering;
#[derive(Clone, Copy)]
pub struct LiftOptions {
    pub position_limit_mm: f64,
    pub derivative_limit_mm: f64,
    pub maximum_work: usize,
    pub maximum_edge_work: usize,
    pub maximum_bytes: usize,
    pub maximum_subdivisions: usize,
}
impl Default for LiftOptions {
    fn default() -> Self {
        Self {
            position_limit_mm: 1e-7,
            derivative_limit_mm: 1e-6,
            maximum_work: 160_000_000,
            maximum_edge_work: 16_000_000,
            maximum_bytes: 64 * 1024 * 1024,
            maximum_subdivisions: 4096,
        }
    }
}
pub struct EdgeLift {
    edge: u32,
    bounds: [f64; 2],
    vertices: [f64; 2],
    range: [I; 2],
    work: usize,
}
impl EdgeLift {
    pub fn edge_entity(&self) -> u32 {
        self.edge
    }
    pub fn field_bounds(&self) -> [f64; 2] {
        self.bounds
    }
    /// End residuals in increasing curve-parameter order.
    pub fn vertex_errors_mm(&self) -> [f64; 2] {
        self.vertices
    }
    pub fn parameter_bounds(&self) -> [[f64; 2]; 2] {
        self.range.map(|v| [v.lo, v.hi])
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
}
pub struct SurfaceLifts<'d, 'i> {
    source: &'d Decoded<'i>,
    edges: Vec<EdgeLift>,
    work: usize,
    bytes: usize,
}
impl SurfaceLifts<'_, '_> {
    pub fn source(&self) -> &Decoded<'_> {
        self.source
    }
    pub fn edges(&self) -> &[EdgeLift] {
        &self.edges
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_bytes(&self) -> usize {
        self.bytes
    }
}
fn p(x: f64) -> I {
    I::point(x)
}
fn point<'a, 'i>(d: &'a Decoded<'i>, id: u32) -> Result<&'a [Scalar<'i>], Error> {
    Ok(d.points.get(&id).ok_or(Error::Invalid)?.coordinates())
}
fn curve<'a, 'i>(d: &'a Decoded<'i>, id: u32) -> Result<&'a Curve<'i>, Error> {
    d.curves.get(&id).ok_or(Error::Invalid)
}
fn weight<'i>(curve: &Curve<'i>, i: usize) -> Scalar<'i> {
    curve
        .weights
        .get(i)
        .copied()
        .unwrap_or(Scalar::new("1.", 1.))
}
fn same(a: Scalar<'_>, c: Scalar<'_>, b: &mut Budget<'_>) -> Result<bool, Error> {
    Ok(numbers::compare(a, c, b)? == Ordering::Equal)
}
fn cartesian(
    d: &Decoded<'_>,
    c: &Curve<'_>,
    dim: usize,
    b: &mut Budget<'_>,
) -> Result<Vec<[I; 4]>, Error> {
    let scale = c.weights.iter().map(|w| w.interval().hi).fold(1., f64::max);
    let mut out = Vec::with_capacity(c.points.len());
    for (i, &id) in c.points.iter().enumerate() {
        b.charge(96)?;
        let q = point(d, id)?;
        if q.len() != dim {
            return Err(Error::Invalid);
        }
        let w = weight(c, i).interval() / p(scale);
        let mut row = [p(0.); 4];
        for axis in 0..dim {
            row[axis] = q[axis].interval() * w;
        }
        row[3] = w;
        out.push(row);
    }
    Ok(out)
}
fn columns(controls: &[[I; 4]]) -> [Polynomial; 4] {
    std::array::from_fn(|axis| controls.iter().map(|v| v[axis]).collect())
}
fn scalar_range<'i>(
    d: &Decoded<'i>,
    edge: &Edge,
    b: &mut Budget<'_>,
) -> Result<[Scalar<'i>; 2], Error> {
    let c = curve(d, edge.curve)?;
    let mut lo = c.knots[0];
    let mut hi = *c.knots.last().ok_or(Error::Invalid)?;
    for pc in edge.pcurves {
        let pc = curve(d, pc.curve)?;
        if numbers::compare(pc.knots[0], lo, b)? == Ordering::Greater {
            lo = pc.knots[0]
        }
        if numbers::compare(*pc.knots.last().ok_or(Error::Invalid)?, hi, b)? == Ordering::Less {
            hi = *pc.knots.last().ok_or(Error::Invalid)?
        }
    }
    if lo.interval().hi >= hi.interval().lo {
        return Err(Error::Unsupported);
    }
    Ok([lo, hi])
}
fn endpoint(
    d: &Decoded<'_>,
    c: &Curve<'_>,
    parameter: I,
    target: u32,
    b: &mut Budget<'_>,
) -> Result<f64, Error> {
    let point = point(d, *d.vertices.get(&target).ok_or(Error::Invalid)?)?;
    let all = cartesian(d, c, 3, b)?;
    let mut found = false;
    let mut maximum = 0_f64;
    for (knot, pair) in c.knots.windows(2).enumerate() {
        let a = numbers::exact_interval(pair[0], b)?;
        let z = numbers::exact_interval(pair[1], b)?;
        if parameter.hi < a.lo || parameter.lo > z.hi {
            continue;
        }
        let intersection = I::new(parameter.lo.max(a.lo), parameter.hi.min(z.hi));
        let t = (intersection - a) / (z - a);
        let t = I::new(t.lo.max(0.), t.hi.min(1.));
        if !t.valid() {
            return Err(Error::Unsupported);
        }
        let start = knot * c.degree;
        let cols = columns(&all[start..start + c.degree + 1]);
        let den = poly::evaluate(&cols[3], t, b)?;
        if den.lo <= 0. {
            return Err(Error::Unsupported);
        }
        let mut error = p(0.);
        for axis in 0..3 {
            error = error
                + p(
                    (poly::evaluate(&cols[axis], t, b)? / den - point[axis].interval()).magnitude(),
                );
        }
        if !error.valid() {
            return Err(Error::Unsupported);
        }
        maximum = maximum.max(error.hi);
        found = true;
    }
    if !found {
        return Err(Error::Unsupported);
    }
    Ok(maximum)
}
fn plane_lift(
    d: &Decoded<'_>,
    c: &Curve<'_>,
    pc: &Curve<'_>,
    plane: &Plane,
    b: &mut Budget<'_>,
) -> Result<[f64; 2], Error> {
    if c.degree != pc.degree || c.points.len() != pc.points.len() || c.knots.len() != pc.knots.len()
    {
        return Err(Error::Unsupported);
    }
    for (&a, &v) in c.knots.iter().zip(&pc.knots) {
        if !same(a, v, b)? {
            return Err(Error::Unsupported);
        }
    }
    let origin = point(d, plane.origin)?;
    for i in 0..c.points.len() {
        if !same(weight(c, i), weight(pc, i), b)? {
            return Err(Error::Unsupported);
        }
        let q = point(d, c.points[i])?;
        let uv = point(d, pc.points[i])?;
        if !same(q[0], uv[0], b)? || !same(q[1], uv[1], b)? || !same(q[2], origin[2], b)? {
            return Err(Error::Unsupported);
        }
    }
    Ok([0., 0.])
}
fn row_lift(
    d: &Decoded<'_>,
    c: &Curve<'_>,
    pc: &Curve<'_>,
    surface: &Tensor<'_>,
    range: [I; 2],
    b: &mut Budget<'_>,
) -> Result<Option<[f64; 2]>, Error> {
    if pc.degree != 1 || pc.points.len() != 2 || c.degree != 10 {
        return Ok(None);
    }
    let first = point(d, pc.points[0])?;
    let last = point(d, pc.points[1])?;
    if !same(first[0], last[0], b)? || !(first[0].is_zero() || first[0].is_one()) {
        return Ok(None);
    }
    if !same(weight(pc, 0), weight(pc, 1), b)? {
        return Ok(None);
    }
    for i in 1..c.points.len() {
        if !same(weight(c, 0), weight(c, i), b)? {
            return Ok(None);
        }
    }
    let ka = numbers::rational(pc.knots[0], b)?;
    let kb = numbers::rational(*pc.knots.last().ok_or(Error::Invalid)?, b)?;
    let va = numbers::rational(first[1], b)?;
    let vb = numbers::rational(last[1], b)?;
    if vb - va.clone() != kb - ka.clone() {
        return Ok(None);
    }
    let shift = va - ka;
    let row = if first[0].is_zero() {
        0
    } else {
        surface.dimensions[0] - 1
    };
    let width = surface.dimensions[1];
    for i in 1..width {
        if !same(
            surface.weights[row * width],
            surface.weights[row * width + i],
            b,
        )? {
            return Ok(None);
        }
    }
    for (span, knots) in c.knots.windows(2).enumerate() {
        let before = numbers::exact_interval(knots[0], b)?;
        let after = numbers::exact_interval(knots[1], b)?;
        if after.hi <= range[0].lo || before.lo >= range[1].hi {
            continue;
        }
        let a = numbers::rational(knots[0], b)? + &shift;
        let z = numbers::rational(knots[1], b)? + &shift;
        use num_traits::ToPrimitive;
        let estimate = a.to_f64().ok_or(Error::Unsupported)?;
        let candidate = surface
            .knots_v
            .partition_point(|k| k.value <= estimate)
            .saturating_sub(1)
            .min(surface.knots_v.len() - 2);
        let mut owner = None;
        for j in candidate.saturating_sub(1)..=(candidate + 1).min(surface.knots_v.len() - 2) {
            b.charge(1)?;
            let pair = &surface.knots_v[j..=j + 1];
            if numbers::rational(pair[0], b)? == a && numbers::rational(pair[1], b)? == z {
                owner = Some(j);
                break;
            }
        }
        let owner = owner.ok_or(Error::Unsupported)?;
        for j in 0..=10 {
            b.charge(1)?;
            let curve_point = c.points[span * 10 + j];
            let surface_point = surface.points[row * width + owner * 10 + j];
            if curve_point != surface_point {
                let a = point(d, curve_point)?;
                let z = point(d, surface_point)?;
                for axis in 0..3 {
                    if !same(a[axis], z[axis], b)? {
                        return Err(Error::Unsupported);
                    }
                }
            }
        }
    }
    Ok(Some([0., 0.]))
}
fn tensor_lift(
    d: &Decoded<'_>,
    c: &Curve<'_>,
    pc: &Curve<'_>,
    surface: &Tensor<'_>,
    algebra: &Algebra,
    options: LiftOptions,
    b: &mut Budget<'_>,
) -> Result<[f64; 2], Error> {
    if c.knots.len() != 2 || pc.knots.len() != 2 || pc.degree > 3 {
        return Err(Error::Unsupported);
    }
    let controls = cartesian(d, c, 3, b)?;
    let original = columns(&controls);
    let a = numbers::unit(pc.knots[0], c.knots[0], c.knots[1], b)?;
    let z = numbers::unit(pc.knots[1], c.knots[0], c.knots[1], b)?;
    let restricted = original
        .iter()
        .map(|v| poly::restrict(v, a, z, b))
        .collect::<Result<Vec<_>, _>>()?;
    let scale = pc
        .weights
        .iter()
        .map(|w| w.interval().hi)
        .fold(1., f64::max);
    let mut u = Vec::new();
    let mut uc = Vec::new();
    let mut v = Vec::new();
    let mut vc = Vec::new();
    let mut den = Vec::new();
    let mut chosen = None;
    for span in 0..surface.knots_v.len() - 1 {
        let mut valid = true;
        for &id in &pc.points {
            let q = point(d, id)?;
            if numbers::compare(q[1], surface.knots_v[span], b)? == Ordering::Less
                || numbers::compare(q[1], surface.knots_v[span + 1], b)? == Ordering::Greater
            {
                valid = false;
                break;
            }
        }
        if valid {
            chosen = Some(span);
            break;
        }
    }
    #[cfg(test)]
    if chosen.is_none() {
        eprintln!(
            "no source span pc={} V={:?} surface={} knots={:?}",
            pc.entity,
            pc.points
                .iter()
                .map(|id| point(d, *id).unwrap()[1].text)
                .collect::<Vec<_>>(),
            surface.entity,
            surface.knots_v.iter().map(|k| k.value).collect::<Vec<_>>()
        );
    }
    let span = chosen.ok_or(Error::Unsupported)?;
    for (i, &id) in pc.points.iter().enumerate() {
        let q = point(d, id)?;
        let uu = numbers::unit_coordinate(q[0], b)?;
        let vv = numbers::unit(q[1], surface.knots_v[span], surface.knots_v[span + 1], b)?;
        let w = weight(pc, i).interval() / p(scale);
        u.push(uu * w);
        uc.push(I::new((p(1.) - uu).lo.max(0.), (p(1.) - uu).hi.min(1.)) * w);
        v.push(vv * w);
        vc.push(I::new((p(1.) - vv).lo.max(0.), (p(1.) - vv).hi.min(1.)) * w);
        den.push(w);
    }
    let [du, dv] = surface.degrees;
    let degree = (du + dv) * pc.degree;
    let mut actual = std::array::from_fn::<_, 4, _>(|_| vec![p(0.); degree + 1]);
    let scale = surface
        .weights
        .iter()
        .map(|w| w.interval().hi)
        .fold(1., f64::max);
    for i in 0..=du {
        let bu = algebra.basis(&u, &uc, du, i, b)?;
        for j in 0..=dv {
            let bv = algebra.basis(&v, &vc, dv, j, b)?;
            let basis = algebra.multiply(&bu, &bv, b)?;
            let id = i * surface.dimensions[1] + span * dv + j;
            let q = point(d, surface.points[id])?;
            let w = surface.weights[id].interval() / p(scale);
            for axis in 0..4 {
                let mut contribution = basis.clone();
                poly::scale(
                    &mut contribution,
                    if axis == 3 { w } else { q[axis].interval() * w },
                    b,
                )?;
                poly::add(&mut actual[axis], &contribution, b)?;
            }
        }
    }
    let denominator = algebra.multiply(&restricted[3], &actual[3], b)?;
    let mut error = Vec::new();
    let mut derivative = Vec::new();
    let dd = poly::derivative(&denominator, b)?;
    for axis in 0..3 {
        let a = algebra.multiply(&restricted[axis], &actual[3], b)?;
        let z = algebra.multiply(&actual[axis], &restricted[3], b)?;
        let e = poly::subtract(&a, &z, b)?;
        let ed = poly::derivative(&e, b)?;
        let a = algebra.multiply(&ed, &denominator, b)?;
        let z = algebra.multiply(&e, &dd, b)?;
        derivative.push(poly::subtract(&a, &z, b)?);
        error.push(e);
    }
    let denominator_squared = algebra.multiply(&denominator, &denominator, b)?;
    let mut pending = vec![(error, derivative, denominator, denominator_squared, 0_usize)];
    let mut visited = 0;
    let mut out = [0_f64; 2];
    while let Some((error, derivative, denominator, denominator_squared, depth)) = pending.pop() {
        visited += 1;
        if visited > options.maximum_subdivisions {
            return Err(Error::Limit);
        }
        b.charge(1)?;
        let lo = poly::lower(&denominator);
        let dlo = poly::lower(&denominator_squared);
        let c0 = (error
            .iter()
            .fold(p(0.), |sum, e| sum + p(poly::magnitude(e)))
            / p(lo))
        .hi;
        let d1 = (derivative
            .iter()
            .fold(p(0.), |sum, e| sum + p(poly::magnitude(e)))
            / p(dlo))
        .hi;
        if lo > 0.
            && dlo > 0.
            && c0.is_finite()
            && d1.is_finite()
            && c0 <= options.position_limit_mm
            && d1 <= options.derivative_limit_mm
        {
            out[0] = out[0].max(c0);
            out[1] = out[1].max(d1);
            continue;
        }
        if depth >= 12 {
            #[cfg(test)]
            eprintln!(
                "lift residual refusal c0={c0} d1={d1} denominator={lo} derivative_den={dlo}"
            );
            return Err(Error::Unsupported);
        }
        let es = error
            .iter()
            .map(|e| poly::split(e, p(0.5), b))
            .collect::<Result<Vec<_>, _>>()?;
        let ds = derivative
            .iter()
            .map(|e| poly::split(e, p(0.5), b))
            .collect::<Result<Vec<_>, _>>()?;
        let ns = poly::split(&denominator, p(0.5), b)?;
        let qs = poly::split(&denominator_squared, p(0.5), b)?;
        for i in 0..2 {
            pending.push((
                es.iter().map(|v| v[i].clone()).collect(),
                ds.iter().map(|v| v[i].clone()).collect(),
                ns[i].clone(),
                qs[i].clone(),
                depth + 1,
            ));
        }
    }
    Ok(out)
}
fn monotone_parameter(
    d: &Decoded<'_>,
    c: &Curve<'_>,
    vertex: u32,
    limit: f64,
    b: &mut Budget<'_>,
) -> Result<I, Error> {
    let target = point(d, *d.vertices.get(&vertex).ok_or(Error::Invalid)?)?;
    let controls = cartesian(d, c, 3, b)?;
    // These long rows have one constant positive rational weight, so Cartesian
    // Bezier controls and their strictly positive Z differences prove monotonicity.
    for i in 1..c.points.len() {
        if !same(weight(c, 0), weight(c, i), b)? {
            return Err(Error::Unsupported);
        }
    }
    // Prove the complete long curve is monotone, not just the selected span.
    for span in 0..c.knots.len() - 1 {
        for j in 0..c.degree {
            b.charge(32)?;
            let a = point(d, c.points[span * c.degree + j])?[2].interval();
            let z = point(d, c.points[span * c.degree + j + 1])?[2].interval();
            if z.lo <= a.hi {
                return Err(Error::Unsupported);
            }
        }
    }
    for (span, knots) in c.knots.windows(2).enumerate() {
        let first = point(d, c.points[span * c.degree])?;
        let last = point(d, c.points[(span + 1) * c.degree])?;
        if numbers::compare(target[2], first[2], b)? == Ordering::Less
            || numbers::compare(target[2], last[2], b)? == Ordering::Greater
        {
            continue;
        }
        if same(target[2], first[2], b)? {
            return numbers::exact_interval(knots[0], b);
        }
        if same(target[2], last[2], b)? {
            return numbers::exact_interval(knots[1], b);
        }
        let xyz = std::array::from_fn::<_, 3, _>(|axis| {
            controls[span * c.degree..=(span + 1) * c.degree]
                .iter()
                .map(|p| p[axis] / p[3])
                .collect::<Polynomial>()
        });
        let dz = poly::derivative(&xyz[2], b)?;
        let min = poly::lower(&dz);
        let max = dz.iter().map(|v| v.hi).fold(0., f64::max);
        if min <= 0. || !max.is_finite() {
            return Err(Error::Unsupported);
        }
        let mut interval = I::new(0., 1.);
        for _ in 0..64 {
            let mut residual = p(0.);
            for axis in 0..3 {
                residual = residual
                    + p(
                        (poly::evaluate(&xyz[axis], interval, b)? - target[axis].interval())
                            .magnitude(),
                    );
            }
            if residual.hi <= limit / 8. {
                let a = numbers::exact_interval(knots[0], b)?;
                let z = numbers::exact_interval(knots[1], b)?;
                return Ok(a + (z - a) * interval);
            }
            let mid = interval.midpoint();
            let value = poly::evaluate(&xyz[2], p(mid), b)?;
            if value.hi < target[2].interval().lo {
                interval = I::new(mid, interval.hi)
            } else if value.lo > target[2].interval().hi {
                interval = I::new(interval.lo, mid)
            } else {
                let enclosed = p(mid) + (target[2].interval() - value) / I::new(min, max);
                interval = I::new(interval.lo.max(enclosed.lo), interval.hi.min(enclosed.hi));
                if !interval.valid() {
                    return Err(Error::Unsupported);
                }
            }
        }
        return Err(Error::Unsupported);
    }
    Err(Error::Unsupported)
}
impl<'i> Decoded<'i> {
    pub fn surface_lifts<'d>(
        &'d self,
        options: LiftOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<SurfaceLifts<'d, 'i>, Error> {
        if !options.position_limit_mm.is_finite()
            || options.position_limit_mm <= 0.
            || options.position_limit_mm > 1e-7
            || !options.derivative_limit_mm.is_finite()
            || options.derivative_limit_mm <= 0.
            || options.derivative_limit_mm > 1e-6
            || options.maximum_work == 0
            || options.maximum_work > 256_000_000
            || options.maximum_edge_work == 0
            || options.maximum_edge_work > 16_000_000
            || options.maximum_bytes == 0
            || options.maximum_bytes > 64 * 1024 * 1024
            || options.maximum_subdivisions == 0
            || options.maximum_subdivisions > 8192
        {
            return Err(Error::Invalid);
        }
        let mut b = Budget {
            work: 0,
            remaining: options.maximum_work,
            bytes: self.bytes,
            limit: options.maximum_bytes,
            interrupt,
        };
        b.charge(self.work)?;
        b.retain(
            4 * 1024 * 1024
                + size_of::<SurfaceLifts<'_, '_>>()
                + self.edges.len() * size_of::<EdgeLift>(),
        )?;
        let algebra = Algebra::new(&mut b)?;
        let mut output = Vec::with_capacity(self.edges.len());
        for edge in self.edges.values() {
            let mut eb = Budget {
                work: 0,
                remaining: b.remaining.min(options.maximum_edge_work),
                bytes: b.bytes,
                limit: b.limit,
                interrupt,
            };
            eb.charge(1)?;
            let c = curve(self, edge.curve)?;
            let complex = c.knots.len() > 2
                && edge
                    .pcurves
                    .iter()
                    .all(|pc| self.tensors.iter().any(|s| s.entity == pc.surface));
            let (range, scalar) = if complex {
                let a = monotone_parameter(
                    self,
                    c,
                    edge.vertices[0],
                    options.position_limit_mm,
                    &mut eb,
                )?;
                let z = monotone_parameter(
                    self,
                    c,
                    edge.vertices[1],
                    options.position_limit_mm,
                    &mut eb,
                )?;
                if edge.same_sense {
                    ([a, z], None)
                } else {
                    ([z, a], None)
                }
            } else {
                let scalar = scalar_range(self, edge, &mut eb)?;
                (
                    [
                        numbers::exact_interval(scalar[0], &mut eb)?,
                        numbers::exact_interval(scalar[1], &mut eb)?,
                    ],
                    Some(scalar),
                )
            };
            if range[0].hi >= range[1].lo {
                return Err(Error::Unsupported);
            }
            let mut vertex_error = [0.; 2];
            for (i, out) in vertex_error.iter_mut().enumerate() {
                let vertex = edge.vertices[if edge.same_sense { i } else { 1 - i }];
                *out = endpoint(self, c, range[i], vertex, &mut eb)?;
                if *out > options.position_limit_mm {
                    #[cfg(test)]
                    eprintln!(
                        "vertex refusal edge={} value={} range={:?}",
                        edge.entity, *out, range
                    );
                    return Err(Error::Unsupported);
                }
            }
            let mut bounds = [0_f64; 2];
            for pc in edge.pcurves {
                let pc_curve = curve(self, pc.curve)?;
                let bound = if let Some(plane) = self.planes.iter().find(|s| s.entity == pc.surface)
                {
                    plane_lift(self, c, pc_curve, plane, &mut eb)?
                } else {
                    let tensor = self
                        .tensors
                        .iter()
                        .find(|s| s.entity == pc.surface)
                        .ok_or(Error::Invalid)?;
                    if let Some(bound) = row_lift(self, c, pc_curve, tensor, range, &mut eb)? {
                        bound
                    } else {
                        if scalar.is_none() {
                            return Err(Error::Unsupported);
                        }
                        tensor_lift(self, c, pc_curve, tensor, &algebra, options, &mut eb)?
                    }
                };
                for axis in 0..2 {
                    bounds[axis] = bounds[axis].max(bound[axis]);
                }
            }
            b.charge(eb.work)?;
            output.push(EdgeLift {
                edge: edge.entity,
                bounds,
                vertices: vertex_error,
                range,
                work: eb.work,
            });
        }
        b.bytes -= 4 * 1024 * 1024;
        Ok(SurfaceLifts {
            source: self,
            edges: output,
            work: b.work,
            bytes: b.bytes,
        })
    }
}
