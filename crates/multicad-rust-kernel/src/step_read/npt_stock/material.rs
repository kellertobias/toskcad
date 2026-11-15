//! Projected literal geometry. This module never rebuilds a thread recipe.
use super::*;
use crate::{
    polar_loop11,
    tapered_thread::{TensorCell11, calculation_cell11, interval::I},
};
const MATERIAL_SCRATCH: usize = 512 * 1024;
#[derive(Clone, Copy, Debug)]
pub struct MaterialOptions {
    pub maximum_work: usize,
    pub maximum_bytes: usize,
    pub maximum_nodes: usize,
    pub maximum_depth: u8,
}
impl Default for MaterialOptions {
    fn default() -> Self {
        Self {
            maximum_work: 64_000_000,
            maximum_bytes: 64 * 1024 * 1024,
            maximum_nodes: 65_536,
            maximum_depth: 32,
        }
    }
}
struct State<'a> {
    budget: Budget<'a>,
    nodes: usize,
    leaves: usize,
    maximum_nodes: usize,
    maximum_depth: u8,
}
pub struct Projection<'parsed, 'document, 'input> {
    parsed: &'parsed Parsed<'document, 'input>,
    caps: [usize; 2],
    cap_heights: [[f64; 2]; 2],
    radial_squared: [f64; 2],
    work: usize,
    bytes: usize,
    nodes: usize,
    leaves: usize,
}
impl<'parsed, 'document, 'input> Projection<'parsed, 'document, 'input> {
    pub fn parsed(&self) -> &'parsed Parsed<'document, 'input> {
        self.parsed
    }
    pub fn cap_face_indices(&self) -> [usize; 2] {
        self.caps
    }
    pub fn cap_height_bounds_mm(&self) -> [[f64; 2]; 2] {
        self.cap_heights
    }
    pub fn lateral_radius_squared_bounds(&self) -> [f64; 2] {
        self.radial_squared
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn query_budget_bytes(&self) -> usize {
        self.bytes
    }
    pub fn nodes(&self) -> usize {
        self.nodes
    }
    pub fn leaves(&self) -> usize {
        self.leaves
    }
}
// The parser checks every initial interval bound is finite and within±1e7.
// Fixed axis0/1 subdivision at at most32 exact dyadic levels stays below the
// numeric cell's64-level limit; rounded convex combinations and their first
// differences are far from overflow. These internal numeric operations have
// neither a work budget nor an interrupt callback. Their errors would signal
// an internal invariant break; public budgets and interrupts remain fallible.
fn loop_error(error: polar_loop11::Error) -> Error {
    match error {
        polar_loop11::Error::Limit => Error::Limit,
        polar_loop11::Error::Cancelled => Error::Cancelled,
        _ => Error::Unsupported,
    }
}
fn scalar(bounds: [f64; 2]) -> I {
    I::new(bounds[0], bounds[1])
}
fn controls(parsed: &Parsed<'_, '_>, tensor: &Tensor, u: usize, v: usize) -> [[[I; 3]; 12]; 12] {
    let height = 11 * tensor.spans[1] + 1;
    std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            parsed.points[tensor.poles[(11 * u + i) * height + 11 * v + j]]
                .coordinates
                .map(|s| scalar(s.bounds()))
        })
    })
}
impl State<'_> {
    fn cell(&mut self, root: TensorCell11) -> Result<f64> {
        let mut radial_lower = f64::INFINITY;
        let mut pending = Vec::new();
        pending.try_reserve_exact(33).map_err(|_| Error::Limit)?;
        pending.push((root, 0u8));
        while let Some((cell, depth)) = pending.pop() {
            self.budget.step(4096)?;
            self.nodes += 1;
            if self.nodes > self.maximum_nodes {
                return Err(Error::Limit);
            }
            let first = cell
                .first_derivative_bounds()
                .expect("sealed finite literal cell");
            let p = first.position.map(scalar);
            let du = first.du.map(scalar);
            let dv = first.dv.map(scalar);
            let angular_u = p[0] * du[1] - p[1] * du[0];
            let angular_v = p[0] * dv[1] - p[1] * dv[0];
            let jacobian = angular_u * dv[2] - angular_v * du[2];
            // J!=0 also excludes rho=0. Retain the actual mixed z_u term;
            // neither a source profile nor an approximately affine z is used
            // for the geometric projection proof.
            if jacobian.hi < 0. && dv[2].lo > 0. {
                self.budget.step(32)?;
                let rho = p[0].square() + p[1].square();
                radial_lower = radial_lower.min(rho.lo.max(0.));
                self.leaves += 1;
                continue;
            }
            if depth == self.maximum_depth {
                return Err(Error::Unsupported);
            }
            self.budget.step(12_288)?;
            let [a, b] = cell
                .subdivide(usize::from(depth % 2 != 0))
                .expect("sealed finite literal cell");
            pending.push((b, depth + 1));
            pending.push((a, depth + 1));
        }
        Ok(radial_lower)
    }
}
impl<'document, 'input> Parsed<'document, 'input> {
    /// Prove the actual interval-valued tensor charts are local polar/axial
    /// homeomorphisms and both actual cap rims wind once. Occupied material
    /// and independently integrated volume are a subsequent typed gate.
    pub fn projection<'parsed>(
        &'parsed self,
        options: MaterialOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Projection<'parsed, 'document, 'input>> {
        if options.maximum_work == 0
            || options.maximum_work > 64_000_000
            || options.maximum_bytes == 0
            || options.maximum_bytes > 64 * 1024 * 1024
            || options.maximum_nodes > 65_536
            || options.maximum_depth > 32
        {
            return Err(Error::Invalid);
        }
        let bytes = self.bytes + MATERIAL_SCRATCH;
        if bytes > options.maximum_bytes {
            return Err(Error::Limit);
        }
        let mut state = State {
            budget: Budget {
                used: 0,
                remaining: options.maximum_work,
                bytes,
                peak: bytes,
                limit: options.maximum_bytes,
                interrupt,
            },
            nodes: 0,
            leaves: 0,
            maximum_nodes: options.maximum_nodes,
            maximum_depth: options.maximum_depth,
        };
        state.budget.step(self.work)?;
        // parse::body seals ten faces with exactly eight unique Tensor supports;
        // the remaining supports are Plane. Capture their height indices while
        // matching the type, rather than inventing errors for a forged Parsed.
        let mut caps = [usize::MAX; 2];
        let mut heights = [usize::MAX; 2];
        let mut count = 0;
        for (index, face) in self.faces.iter().enumerate() {
            if let Support::Plane { height } = face.support {
                assert!(count < 2, "sealed graph has exactly two plane faces");
                caps[count] = index;
                heights[count] = height;
                count += 1;
            }
        }
        assert_eq!(count, 2, "sealed graph has exactly two plane faces");
        if self.points[heights[0]].coordinates[2].value()
            > self.points[heights[1]].coordinates[2].value()
        {
            caps.swap(0, 1);
            heights.swap(0, 1);
        }
        let lower = self.points[heights[0]].coordinates[2];
        let upper = self.points[heights[1]].coordinates[2];
        // parse::body already seals false lateral face senses, matching
        // the outward -S_u x S_v used by J<0 and z_v>0. Cap senses remain
        // structurally admissible in either direction and are checked here.
        if lower.bounds()[1] >= upper.bounds()[0]
            || self.faces[caps[0]].forward
            || !self.faces[caps[1]].forward
        {
            return Err(Error::Unsupported);
        }
        // Bernstein convex hulls earn the GLOBAL closed z range. At a cap,
        // compare actual decimal leaves for equality rather than independent
        // overlapping caches. The disjoint-bound route covers interior poles.
        for tensor in &self.tensors {
            for index in &tensor.poles {
                state.budget.step(8)?;
                let z = self.points[*index].coordinates[2];
                if (z.bounds()[0] < lower.bounds()[1]
                    && !numeric::equal(z.spelling(), lower.spelling(), &mut state.budget)?)
                    || (z.bounds()[1] > upper.bounds()[0]
                        && !numeric::equal(z.spelling(), upper.spelling(), &mut state.budget)?)
                {
                    return Err(Error::Unsupported);
                }
            }
        }
        let mut points = state.budget.allocate(self.points.len())?;
        for point in &self.points {
            points.push([point.coordinates[0].bounds(), point.coordinates[1].bounds()]);
        }
        for index in caps {
            let face = &self.faces[index];
            let mut indices: [Vec<usize>; 4] = std::array::from_fn(|_| Vec::new());
            for (index, coedge) in face.coedges.iter().enumerate() {
                let curve = &self.curves[self.edges[coedge.edge].curve];
                indices[index] = state.budget.allocate(curve.poles.len())?;
                if coedge.forward {
                    indices[index].extend_from_slice(&curve.poles);
                } else {
                    indices[index].extend(curve.poles.iter().rev().copied());
                }
            }
            let curves: [polar_loop11::IntervalCurve<'_>; 4] = std::array::from_fn(|i| {
                let edge = &self.edges[face.coedges[i].edge];
                polar_loop11::IntervalCurve {
                    curve_owner: edge.id as usize,
                    spans: self.curves[edge.curve].spans,
                    point_indices: &indices[i],
                }
            });
            let proof = polar_loop11::certify_intervals(
                polar_loop11::IntervalInput {
                    points: &points,
                    curves: &curves,
                    counter_clockwise: face.forward,
                },
                polar_loop11::Options {
                    maximum_work: state.budget.remaining,
                    maximum_nodes: options.maximum_nodes - state.nodes,
                    maximum_bytes: 16 * 1024 * 1024,
                    maximum_depth: options.maximum_depth,
                },
                interrupt,
            )
            .map_err(loop_error)?;
            state.budget.step(proof.work_used())?;
            state.nodes += proof.nodes();
            state.leaves += proof.leaves();
            for values in &indices {
                state.budget.release(values);
            }
            drop(indices);
        }
        state.budget.release(&points);
        drop(points);
        let mut radial_squared = [f64::INFINITY, f64::NEG_INFINITY];
        for tensor in &self.tensors {
            for u in 0..tensor.spans[0] {
                for v in 0..tensor.spans[1] {
                    state.budget.step(1024)?;
                    let coefficients = controls(self, tensor, u, v);
                    // Squared XY norm is convex: every Bernstein image lies
                    // below the largest actual interval control-point norm.
                    // This avoids the loose corner of a componentwise hull.
                    state.budget.step(2048)?;
                    for point in coefficients.iter().flatten() {
                        let rho = point[0].square() + point[1].square();
                        radial_squared[1] = radial_squared[1].max(rho.hi);
                    }
                    let cell =
                        calculation_cell11(std::array::from_fn(|i| coefficients[i / 12][i % 12]))
                            .expect("sealed finite literal cell");
                    radial_squared[0] = radial_squared[0].min(state.cell(cell)?);
                }
            }
        }
        state.budget.step(1)?;
        Ok(Projection {
            parsed: self,
            caps,
            cap_heights: [lower.bounds(), upper.bounds()],
            radial_squared,
            work: state.budget.used,
            bytes: state.budget.peak,
            nodes: state.nodes,
            leaves: state.leaves,
        })
    }
}
#[cfg(test)]
pub(super) mod tests;

pub struct Occupied<'parsed, 'document, 'input> {
    projection: Projection<'parsed, 'document, 'input>,
    volume: [f64; 2],
    areas: [[f64; 2]; 2],
    work: usize,
    bytes: usize,
}
impl<'parsed, 'document, 'input> Occupied<'parsed, 'document, 'input> {
    pub fn projection(&self) -> &Projection<'parsed, 'document, 'input> {
        &self.projection
    }
    pub fn volume_mm3(&self) -> [f64; 2] {
        self.volume
    }
    pub fn cap_areas_mm2(&self) -> [[f64; 2]; 2] {
        self.areas
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn query_budget_bytes(&self) -> usize {
        self.bytes
    }
}
fn maximum_derivative_error(error: &[[I; 12]; 12], axis: usize) -> f64 {
    let mut maximum = 0_f64;
    for i in 0..12 {
        for j in 0..12 {
            if (if axis == 0 { i } else { j }) < 11 {
                let next = if axis == 0 {
                    error[i + 1][j]
                } else {
                    error[i][j + 1]
                };
                maximum = maximum.max((I::point(11.) * (next - error[i][j])).magnitude());
            }
        }
    }
    maximum
}
// Integrate the actual literal polynomial. The fixed affine z polynomial is
// only an integration proxy; BOTH derivatives of its full residual are bounded.
fn cell_flux(
    poles: &[[[I; 3]; 12]; 12],
    base: f64,
    slope: f64,
    matrix: &super::moments::Moments,
    budget: &mut Budget<'_>,
) -> Result<I> {
    let cell = calculation_cell11(std::array::from_fn(|i| poles[i / 12][i % 12]))
        .expect("sealed finite literal cell");
    let first = cell
        .first_derivative_bounds()
        .expect("sealed finite literal cell");
    let p = first.position.map(scalar);
    let du = first.du.map(scalar);
    let dv = first.dv.map(scalar);
    let slope = I::point(slope);
    let base = I::point(base);
    let error = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            poles[i][j][2] - (base + slope * I::point(j as f64) / I::point(11.))
        })
    });
    let eu = maximum_derivative_error(&error, 0);
    let ev = maximum_derivative_error(&error, 1);
    let cu = p[0] * du[1] - p[1] * du[0];
    let cv = p[0] * dv[1] - p[1] * dv[0];
    let residual =
        (I::point(cu.magnitude()) * I::point(ev) + I::point(cv.magnitude()) * I::point(eu)).hi;
    let angular = matrix.tensor_angular(poles, budget)?;
    // div(x/2,y/2,0)=1. Horizontal cap flux is exactly zero;
    // outward side flux is -J/2, including the mixed z_u term.
    Ok(-I::point(0.5) * (slope * angular + I::new(-residual, residual)))
}
impl<'parsed, 'document, 'input> Projection<'parsed, 'document, 'input> {
    /// The parsed connected quad sphere has cyclic3/4-valent vertex links.
    /// Removing its two disjoint planar Jordan disks leaves a cylinder. The
    /// literal J sign proves a local polar/axial homeomorphism on every cell;
    /// exact shared edge polynomials glue the orientations. This also applies
    /// to internal C0 Bernstein knot seams: adjacent cells reuse the identical
    /// boundary coefficient row, and their closed-cell J signs agree. Internal
    /// knot crossings have four cyclic incident quadrants, just as the parsed
    /// face vertices have three or four. At each interior vertex every
    /// nondegenerate projected corner is a wedge strictly smaller than pi,
    /// so the closed link has degree1 (positive and strictly below2 turns).
    /// The proper cylinder map has degree1 on either cap and is therefore
    /// one-to-one globally. The projection gate also bounds EVERY literal
    /// tensor z coefficient between the actual cap heights, so boundary
    /// neighborhoods cannot wrap outside the target cylinder. This earns the
    /// actual radial occupied solid,
    /// without an authored source model or an embedding claim from a mesh.
    pub fn occupied(
        self,
        options: MaterialOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Occupied<'parsed, 'document, 'input>> {
        if options.maximum_work == 0
            || options.maximum_work > 64_000_000
            || options.maximum_bytes == 0
            || options.maximum_bytes > 64 * 1024 * 1024
            || options.maximum_nodes > 65_536
            || options.maximum_depth > 32
        {
            return Err(Error::Invalid);
        }
        if self.bytes > options.maximum_bytes {
            return Err(Error::Limit);
        }
        let mut budget = Budget {
            used: 0,
            remaining: options.maximum_work,
            bytes: self.bytes,
            peak: self.bytes,
            limit: options.maximum_bytes,
            interrupt,
        };
        budget.step(self.work)?;
        let matrix = super::moments::Moments::new(&mut budget)?;
        let mut volume = I::point(0.);
        for tensor in &self.parsed.tensors {
            for u in 0..tensor.spans[0] {
                for v in 0..tensor.spans[1] {
                    budget.step(8192)?;
                    let poles = controls(self.parsed, tensor, u, v);
                    let height = 11 * tensor.spans[1] + 1;
                    let a = self.parsed.points[tensor.poles[(11 * u) * height + 11 * v]]
                        .coordinates[2]
                        .value();
                    let b = self.parsed.points[tensor.poles[(11 * u) * height + 11 * v + 11]]
                        .coordinates[2]
                        .value();
                    volume = volume + cell_flux(&poles, a, b - a, &matrix, &mut budget)?;
                }
            }
        }
        let mut areas = [[0.; 2]; 2];
        for (index, face_index) in self.caps.into_iter().enumerate() {
            let face = &self.parsed.faces[face_index];
            let mut integral = I::point(0.);
            for coedge in face.coedges {
                let curve = &self.parsed.curves[self.parsed.edges[coedge.edge].curve];
                for span in 0..curve.spans {
                    let poles = std::array::from_fn(|i| {
                        self.parsed.points[curve.poles[11 * span + i]]
                            .coordinates
                            .map(|s| scalar(s.bounds()))
                    });
                    let value = matrix.curve_angular(&poles, &mut budget)?;
                    integral = if coedge.forward {
                        integral + value
                    } else {
                        integral - value
                    };
                }
            }
            let area = if face.forward {
                I::point(0.5) * integral
            } else {
                -I::point(0.5) * integral
            };
            if !area.valid() || area.lo <= 0. {
                return Err(Error::Unsupported);
            }
            areas[index] = [area.lo, area.hi];
        }
        if !volume.valid() || volume.lo <= 0. {
            return Err(Error::Unsupported);
        }
        budget.step(1)?;
        Ok(Occupied {
            projection: self,
            volume: [volume.lo, volume.hi],
            areas,
            work: budget.used,
            bytes: budget.peak,
        })
    }
}
