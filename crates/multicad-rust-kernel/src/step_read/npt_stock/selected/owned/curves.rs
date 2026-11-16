//! Finite original C11 edge geometry; display facets never define topology.
use super::*;
const POINT_LIMIT: usize = 131_072;
const DEPTH: u8 = 48;
#[derive(Clone)]
pub struct LiteralVertexSample {
    source: OwnedMaterial,
    vertex: usize,
    point: [f64; 3],
    bounds: [[f64; 2]; 3],
    work: usize,
    peak: usize,
}
impl LiteralVertexSample {
    pub fn source(&self) -> &OwnedMaterial {
        &self.source
    }
    pub fn vertex_index(&self) -> usize {
        self.vertex
    }
    pub fn vertex_entity(&self) -> u32 {
        self.source.vertices()[self.vertex].entity()
    }
    pub fn point_entity(&self) -> u32 {
        self.source.0.points[self.source.vertices()[self.vertex].point].id
    }
    pub fn position_mm(&self) -> [f64; 3] {
        self.point
    }
    pub fn position_bounds_mm(&self) -> [[f64; 2]; 3] {
        self.bounds
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn accounted_peak_bytes(&self) -> usize {
        self.peak
    }
}
pub struct LiteralEdgePolyline {
    source: OwnedMaterial,
    edge: usize,
    points: Vec<[f64; 3]>,
    parameters: Vec<f64>,
    error: f64,
    length: I,
    work: usize,
    remaining: usize,
    peak: usize,
}
impl LiteralEdgePolyline {
    pub fn source(&self) -> &OwnedMaterial {
        &self.source
    }
    pub fn edge_index(&self) -> usize {
        self.edge
    }
    pub fn edge_entity(&self) -> u32 {
        self.source.edges()[self.edge].entity()
    }
    pub fn curve_entity(&self) -> u32 {
        self.source.curves()[self.source.edges()[self.edge].curve].entity()
    }
    pub fn vertex_entities(&self) -> [u32; 2] {
        self.source.edges()[self.edge]
            .vertices
            .map(|v| self.source.vertices()[v].entity())
    }
    pub fn points_mm(&self) -> &[[f64; 3]] {
        &self.points
    }
    /// Original uniform knot parameter, increasing from the actual edge's
    /// first vertex to its second. Every original knot is retained.
    pub fn parameters(&self) -> &[f64] {
        &self.parameters
    }
    pub fn position_error_mm(&self) -> f64 {
        self.error
    }
    /// Bounds the original curve's arclength, using each accepted source
    /// leaf's endpoint chord and Bernstein control-polygon inequality.
    pub fn length_bounds_mm(&self) -> [f64; 2] {
        bounds(self.length)
    }
    /// Actual query only; immutable source construction is reported separately.
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn remaining_work(&self) -> usize {
        self.remaining
    }
    pub fn accounted_peak_bytes(&self) -> usize {
        self.peak
    }
}
#[derive(Clone, Copy)]
struct Cell {
    controls: [[I; 3]; 12],
    parameter: [f64; 2],
    depth: u8,
}
fn local(source: &OwnedMaterial, index: usize) -> [I; 3] {
    source.0.points[index]
        .coordinates
        .each_ref()
        .map(|s| interval(s.bounds))
}
fn world(source: &OwnedMaterial, point: [I; 3]) -> Result<[I; 3]> {
    source
        .map_box_mm(point.map(bounds))
        .map(|p| p.map(interval))
}
fn vector(source: &OwnedMaterial, point: [I; 3]) -> [I; 3] {
    std::array::from_fn(|a| {
        (0..3).fold(I::point(0.), |sum, j| {
            sum + source.0.linear[j][a] * point[j]
        })
    })
}
fn l1(v: [I; 3]) -> f64 {
    v.into_iter()
        .fold(I::point(0.), |sum, x| sum + I::point(x.magnitude()))
        .hi
}
fn length(v: [I; 3]) -> I {
    v.into_iter()
        .fold(I::point(0.), |sum, x| sum + x.square())
        .sqrt()
}
fn subtract(a: [I; 3], b: [I; 3]) -> [I; 3] {
    std::array::from_fn(|i| a[i] - b[i])
}
fn begin<'a>(
    source: &OwnedMaterial,
    options: Options,
    interrupt: Option<&'a dyn Fn() -> bool>,
) -> Result<Budget<'a>> {
    let mut budget = Budget::new(options, interrupt)?;
    hold(&mut budget, source.retained_budget_bytes())?;
    Ok(budget)
}
fn push(
    points: &mut Vec<[f64; 3]>,
    parameters: &mut Vec<f64>,
    value: [f64; 3],
    parameter: f64,
    limit: usize,
    budget: &mut Budget<'_>,
) -> Result<()> {
    budget.step(8)?;
    if points.len() == limit {
        return Err(Error::Limit);
    }
    if points.len() == points.capacity() {
        let next = points.capacity().max(8).saturating_mul(2).min(limit);
        let mut new_points = budget.allocate(next)?;
        // Parameter storage must cover the ACTUAL point allocation capacity,
        // even on an allocator that grants more than the requested capacity.
        let mut new_parameters = budget.allocate(new_points.capacity())?;
        budget.step(points.len() * 4 + 1)?;
        new_points.extend_from_slice(points);
        new_parameters.extend_from_slice(parameters);
        budget.release(points);
        budget.release(parameters);
        *points = new_points;
        *parameters = new_parameters;
    }
    points.push(value);
    parameters.push(parameter);
    Ok(())
}
impl OwnedMaterial {
    pub fn sample_vertex(
        &self,
        vertex: usize,
        options: Options,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<LiteralVertexSample> {
        let owner = self.vertices().get(vertex).ok_or(Error::Invalid)?;
        let mut budget = begin(self, options, interrupt)?;
        budget.step(128)?;
        let point = world(self, local(self, owner.point))?;
        budget.step(1)?;
        Ok(LiteralVertexSample {
            source: self.clone(),
            vertex,
            point: point.map(I::midpoint),
            bounds: point.map(bounds),
            work: budget.used,
            peak: budget.peak,
        })
    }
    /// Query the original degree11 curve, including every original knot and
    /// exact topological endpoint owner. Original decimal controls and STEP
    /// mapping bound every COMPLETE source span; no nominal reconstruction.
    pub fn sample_edge(
        &self,
        edge: usize,
        deflection_mm: f64,
        point_limit: usize,
        options: Options,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<LiteralEdgePolyline> {
        if !deflection_mm.is_finite()
            || deflection_mm <= 0.
            || !(2..=POINT_LIMIT).contains(&point_limit)
        {
            return Err(Error::Invalid);
        }
        let owner = self.edges().get(edge).ok_or(Error::Invalid)?;
        let curve = &self.curves()[owner.curve];
        let mut budget = begin(self, options, interrupt)?;
        budget.step(32)?;
        if point_limit < curve.spans + 1 {
            return Err(Error::Limit);
        }
        let mut points = Vec::new();
        let mut parameters = Vec::new();
        let first = world(self, local(self, curve.poles[0]))?.map(I::midpoint);
        push(
            &mut points,
            &mut parameters,
            first,
            0.,
            point_limit,
            &mut budget,
        )?;
        let mut pending = budget.allocate(DEPTH as usize + 1)?;
        let mut error = 0f64;
        let mut arclength = I::point(0.);
        for span in 0..curve.spans {
            budget.step(128)?;
            pending.push(Cell {
                controls: std::array::from_fn(|i| local(self, curve.poles[11 * span + i])),
                parameter: [
                    span as f64 / curve.spans as f64,
                    (span + 1) as f64 / curve.spans as f64,
                ],
                depth: 0,
            });
            while let Some(cell) = pending.pop() {
                budget.step(4096)?;
                let a = *points.last().expect("original first edge endpoint");
                let first = world(self, cell.controls[0])?;
                let last = world(self, cell.controls[11])?;
                let b = last.map(I::midpoint);
                let rounding =
                    l1(subtract(first, a.map(I::point))).max(l1(subtract(last, b.map(I::point))));
                let delta = subtract(cell.controls[11], cell.controls[0]);
                let mut chord_error = 0f64;
                for (i, pole) in cell.controls.into_iter().enumerate() {
                    let weight = I::point(i as f64) / I::point(11.);
                    let residual =
                        std::array::from_fn(|k| pole[k] - cell.controls[0][k] - weight * delta[k]);
                    chord_error = chord_error.max(l1(vector(self, residual)));
                }
                // Bernstein partition of unity and sum(i/11)*B_i=t
                // make these residuals the COMPLETE source-minus-endpoint
                // chord polynomial. The actual target endpoints contribute
                // their independently enclosed rounding error by convexity.
                let enclosure = (I::point(chord_error) + I::point(rounding)).hi;
                if enclosure <= deflection_mm {
                    let lower = length(vector(self, delta));
                    let mut upper = I::point(0.);
                    for pair in cell.controls.windows(2) {
                        upper = upper + length(vector(self, subtract(pair[1], pair[0])));
                    }
                    let bound = I::new(lower.lo, upper.hi);
                    if !bound.valid() || !enclosure.is_finite() {
                        return Err(Error::Unsupported);
                    }
                    arclength = arclength + bound;
                    error = error.max(enclosure);
                    push(
                        &mut points,
                        &mut parameters,
                        b,
                        cell.parameter[1],
                        point_limit,
                        &mut budget,
                    )?;
                } else {
                    if cell.depth == DEPTH {
                        return Err(Error::Unsupported);
                    }
                    let mid = cell.parameter[0].midpoint(cell.parameter[1]);
                    if mid <= cell.parameter[0] || mid >= cell.parameter[1] {
                        return Err(Error::Unsupported);
                    }
                    budget.step(2048)?;
                    let [left, right] = crate::polar_loop11::split_half_controls(cell.controls);
                    pending.push(Cell {
                        controls: right,
                        parameter: [mid, cell.parameter[1]],
                        depth: cell.depth + 1,
                    });
                    pending.push(Cell {
                        controls: left,
                        parameter: [cell.parameter[0], mid],
                        depth: cell.depth + 1,
                    });
                }
            }
        }
        budget.release(&pending);
        budget.step(1)?;
        Ok(LiteralEdgePolyline {
            source: self.clone(),
            edge,
            points,
            parameters,
            error,
            length: arclength,
            work: budget.used,
            remaining: budget.remaining,
            peak: budget.peak,
        })
    }
}
#[cfg(test)]
mod tests;
