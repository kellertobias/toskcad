//! Whole-box queries of the original literal radial material, never a recipe.
use super::*;
use crate::tapered_thread::TensorCell11;
// At most 65 pending C11 nets (144 interval points plus row ownership),
// the current net and subdivision/evaluation temporaries fit within this
// conservative charged peak. The immutable literal arena is charged separately.
const QUERY_SCRATCH: usize = 1024 * 1024;
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Location {
    Inside,
    Outside,
    Unresolved,
}
#[derive(Clone, Copy, Debug)]
pub struct QueryOptions {
    pub maximum_work: usize,
    pub maximum_bytes: usize,
    pub maximum_nodes: usize,
    pub maximum_depth: u8,
}
impl Default for QueryOptions {
    fn default() -> Self {
        Self {
            maximum_work: 16_000_000,
            maximum_bytes: 64 * 1024 * 1024,
            maximum_nodes: 65_536,
            maximum_depth: 48,
        }
    }
}
pub struct MaterialQuery<'source> {
    source: &'source OwnedMaterial,
    local_box: [[f64; 2]; 3],
    location: Location,
    work: usize,
    bytes: usize,
    nodes: usize,
}
impl<'source> MaterialQuery<'source> {
    pub fn source(&self) -> &'source OwnedMaterial {
        self.source
    }
    pub fn local_bounds_mm(&self) -> [[f64; 2]; 3] {
        self.local_box
    }
    pub fn location(&self) -> Location {
        self.location
    }
    /// Actual query work only. Immutable construction was already certified;
    /// it is available separately from source().work_used(), never repeated.
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn query_budget_bytes(&self) -> usize {
        self.bytes
    }
    pub fn nodes(&self) -> usize {
        self.nodes
    }
}
struct State<'a> {
    budget: Budget<'a>,
    point: [I; 3],
    radius: I,
    nodes: usize,
    options: QueryOptions,
    inside: bool,
    outside: bool,
    leaves: usize,
    strict_z: bool,
}
impl State<'_> {
    fn skip(&self, position: [I; 3]) -> bool {
        if position[2].hi < self.point[2].lo || position[2].lo > self.point[2].hi {
            return true;
        }
        let cross = self.point[0] * position[1] - self.point[1] * position[0];
        let dot = self.point[0] * position[0] + self.point[1] * position[1];
        cross.lo > 0. || cross.hi < 0. || dot.hi <= 0.
    }
    fn certify_root(
        &mut self,
        cell: &TensorCell11,
        position: [I; 3],
        location: Location,
    ) -> Result<Option<Location>> {
        if location == Location::Inside && !self.strict_z {
            return Ok(None);
        }
        let dot = self.point[0] * position[0] + self.point[1] * position[1];
        if dot.lo <= 0. {
            return Ok(None);
        }
        self.budget.step(65_536)?;
        let first = cell
            .first_derivative_bounds()
            .expect("finite literal query derivatives");
        let du = first.du.map(interval);
        let dv = first.dv.map(interval);
        let j = [
            [
                self.point[0] * du[1] - self.point[1] * du[0],
                self.point[0] * dv[1] - self.point[1] * dv[0],
            ],
            [du[2], dv[2]],
        ];
        let m = j.map(|row| row.map(I::midpoint));
        let determinant = m[0][0] * m[1][1] - m[0][1] * m[1][0];
        if !determinant.is_finite() || determinant == 0. {
            return Ok(None);
        }
        let c = [
            [m[1][1] / determinant, -m[0][1] / determinant],
            [-m[1][0] / determinant, m[0][0] / determinant],
        ]
        .map(|row| row.map(I::point));
        let det = c[0][0] * c[1][1] - c[0][1] * c[1][0];
        if !det.valid() || (det.lo <= 0. && det.hi >= 0.) {
            return Ok(None);
        }
        let middle = cell
            .point_bounds([[0.5, 0.5]; 2])
            .expect("literal unit midpoint")
            .map(interval);
        let f = [
            self.point[0] * middle[1] - self.point[1] * middle[0],
            middle[2] - self.point[2],
        ];
        let remainder = I::new(-0.5, 0.5);
        for (i, row) in c.iter().enumerate() {
            let mut image = I::point(0.5) - row[0] * f[0] - row[1] * f[1];
            for (axis, (a, b)) in j[0].iter().zip(j[1]).enumerate() {
                let e = I::point(f64::from(i == axis)) - row[0] * *a - row[1] * b;
                image = image + e * remainder;
            }
            if !image.valid() || image.lo <= 0. || image.hi >= 1. {
                return Ok(None);
            }
        }
        // For EVERY query point, x-C*F(x) maps the closed unit square strictly
        // inside itself. Brouwer gives a root; the directed determinant earns
        // nonsingular defining C. Positive dot selects the forward ray. The
        // retained global radial material proof makes that intersection unique.
        Ok(Some(location))
    }
    fn cell(&mut self, root: TensorCell11) -> Result<Option<Location>> {
        let mut pending = Vec::new();
        pending.try_reserve_exact(65).map_err(|_| Error::Limit)?;
        pending.push((root, 0u8));
        while let Some((cell, depth)) = pending.pop() {
            self.budget.step(4096)?;
            self.nodes += 1;
            if self.nodes > self.options.maximum_nodes {
                return Err(Error::Limit);
            }
            let position = cell
                .bounds()
                .expect("finite retained C11 hull")
                .map(interval);
            if self.skip(position) {
                continue;
            }
            let radial = position[0].square() + position[1].square();
            if self.radius.hi < radial.lo {
                if let Some(location) = self.certify_root(&cell, position, Location::Inside)? {
                    return Ok(Some(location));
                }
                self.outside = false;
                self.leaves += 1;
            } else if self.radius.lo > radial.hi {
                if let Some(location) = self.certify_root(&cell, position, Location::Outside)? {
                    return Ok(Some(location));
                }
                self.inside = false;
                self.leaves += 1;
            } else {
                if depth == self.options.maximum_depth {
                    self.inside = false;
                    self.outside = false;
                    continue;
                }
                self.budget.step(12_288)?;
                let [a, b] = cell
                    .subdivide(usize::from(depth % 2 != 0))
                    .expect("bounded dyadic literal query subdivision");
                pending.push((b, depth + 1));
                pending.push((a, depth + 1));
            }
        }
        Ok(None)
    }
}
impl OwnedMaterial {
    /// Classify EVERY point of a closed local-mm box against the original
    /// decoded material. Inside/Outside are strict whole-box conclusions.
    /// Boundary/precision ambiguity returns Unresolved; exhausting a caller
    /// work/node/storage allowance returns Limit without changing the source.
    pub fn classify_local_box(
        &self,
        input: [[f64; 2]; 3],
        options: QueryOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<MaterialQuery<'_>> {
        self.query(input, options, interrupt, false)
    }
    /// Enclose the world box through the retained original STEP inverse, then
    /// classify that local enclosure. A wider inverse enclosure can reduce
    /// precision but cannot turn an ambiguous world box into an admission.
    pub fn classify_box(
        &self,
        input: [[f64; 2]; 3],
        options: QueryOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<MaterialQuery<'_>> {
        self.query(input, options, interrupt, true)
    }
    fn query(
        &self,
        input: [[f64; 2]; 3],
        options: QueryOptions,
        interrupt: Option<&dyn Fn() -> bool>,
        world: bool,
    ) -> Result<MaterialQuery<'_>> {
        if options.maximum_work == 0
            || options.maximum_work > 16_000_000
            || options.maximum_bytes == 0
            || options.maximum_bytes > 64 * 1024 * 1024
            || options.maximum_nodes > 65_536
            || options.maximum_depth > 64
            || input.map(interval).iter().any(|i| !i.valid())
        {
            return Err(Error::Invalid);
        }
        let bytes = self.retained_budget_bytes() + QUERY_SCRATCH;
        if bytes > options.maximum_bytes {
            return Err(Error::Limit);
        }
        let mut budget = Budget {
            used: 0,
            remaining: options.maximum_work,
            bytes,
            peak: bytes,
            limit: options.maximum_bytes,
            interrupt,
        };
        budget.step(256)?;
        let input = if world {
            self.inverse_box_mm(input)?
        } else {
            input
        };
        let point = input.map(interval);
        let mut state = State {
            budget,
            point,
            radius: I::point(0.),
            nodes: 0,
            options,
            inside: true,
            outside: true,
            leaves: 0,
            strict_z: false,
        };
        let location = self.locate(&mut state)?;
        state.budget.step(1)?;
        Ok(MaterialQuery {
            source: self,
            local_box: input,
            location,
            work: state.budget.used,
            bytes,
            nodes: state.nodes,
        })
    }
    fn locate(&self, state: &mut State<'_>) -> Result<Location> {
        if state
            .point
            .iter()
            .zip(self.0.local_bounds)
            .any(|(q, b)| q.hi < b[0] || q.lo > b[1])
        {
            return Ok(Location::Outside);
        }
        let [lower, upper] = self.0.cap_heights.map(interval);
        // The sealed tensor hull reaches both caps and every original pole
        // lies between them; the preceding local-box test already covers
        // all strictly exterior cap heights.
        let strict_z = state.point[2].lo > lower.hi && state.point[2].hi < upper.lo;
        state.strict_z = strict_z;
        // Every strict interior axial point belongs to the earned radial
        // solid. Resolve this common case before squared-norm arithmetic,
        // which can underflow for genuinely tiny decoded decimal geometry.
        if strict_z && state.point[..2].iter().all(|q| q.lo == 0. && q.hi == 0.) {
            return Ok(Location::Inside);
        }
        state.radius = state.point[0].square() + state.point[1].square();
        if !state.radius.valid() {
            return Ok(Location::Unresolved);
        }
        if state.radius.lo > self.0.radial_squared[1] {
            return Ok(Location::Outside);
        }
        if strict_z && state.radius.hi < self.0.radial_squared[0] {
            return Ok(Location::Inside);
        }
        for (tensor_index, tensor) in self.0.tensors.iter().enumerate() {
            for u in 0..tensor.spans[0] {
                for v in 0..tensor.spans[1] {
                    state.budget.step(4096)?;
                    let height = 11 * tensor.spans[1] + 1;
                    let mut hull = [I::new(f64::INFINITY, f64::NEG_INFINITY); 3];
                    for i in 0..12 {
                        for j in 0..12 {
                            let point =
                                &self.0.points[tensor.poles[(11 * u + i) * height + 11 * v + j]];
                            for (axis, bound) in hull.iter_mut().enumerate() {
                                bound.lo = bound.lo.min(point.coordinates[axis].bounds[0]);
                                bound.hi = bound.hi.max(point.coordinates[axis].bounds[1]);
                            }
                        }
                    }
                    if state.skip(hull) {
                        continue;
                    }
                    if let Some(location) = state.cell(self.cell(tensor_index, [u, v])?)? {
                        return Ok(location);
                    }
                }
            }
        }
        // Proper radial projection guarantees a boundary intersection for
        // every nonzero ray at every strict interior height. We classify only
        // when ALL retained possible intersections give the same strict sign.
        if state.leaves == 0 {
            return Ok(Location::Unresolved);
        }
        Ok(if state.inside && strict_z {
            Location::Inside
        } else if state.outside {
            Location::Outside
        } else {
            Location::Unresolved
        })
    }
}
#[cfg(test)]
mod tests;
