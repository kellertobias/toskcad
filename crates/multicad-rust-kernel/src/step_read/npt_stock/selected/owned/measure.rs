//! Integrals of retained literal supports, independent of display facets.
use super::*;
#[derive(Clone, Copy)]
pub struct FaceMeasureOptions {
    pub area_width_mm2: f64,
    pub centroid_width_mm: f64,
    pub maximum_work: usize,
    pub maximum_bytes: usize,
}
impl Default for FaceMeasureOptions {
    fn default() -> Self {
        Self {
            area_width_mm2: 1e-6,
            centroid_width_mm: 1e-7,
            maximum_work: 512_000_000,
            maximum_bytes: 64 * 1024 * 1024,
        }
    }
}
pub struct FaceMeasure {
    source: OwnedMaterial,
    face: usize,
    area: I,
    centroid: [I; 3],
    centered_moments: [I; 3],
    anchor: [f64; 3],
    work: usize,
    bytes: usize,
}
impl FaceMeasure {
    pub fn source(&self) -> &OwnedMaterial {
        &self.source
    }
    pub fn face_index(&self) -> usize {
        self.face
    }
    pub fn area_mm2(&self) -> [f64; 2] {
        bounds(self.area)
    }
    pub fn centroid_mm(&self) -> [[f64; 2]; 3] {
        self.centroid.map(bounds)
    }
    pub fn moment_anchor_mm(&self) -> [f64; 3] {
        self.anchor
    }
    pub fn first_moments_about_anchor_mm3(&self) -> [[f64; 2]; 3] {
        self.centered_moments.map(bounds)
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.bytes
    }
}
struct Work<'a> {
    left: usize,
    used: usize,
    bytes: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl Work<'_> {
    fn step(&mut self, n: usize) -> Result<()> {
        if self.interrupt.is_some_and(|f| f()) {
            return Err(Error::Cancelled);
        }
        self.left = self.left.checked_sub(n).ok_or(Error::Limit)?;
        self.used += n;
        Ok(())
    }
}
fn p(x: f64) -> I {
    I::point(x)
}
fn binomials<const N: usize>(work: &mut Work<'_>) -> Result<[f64; N]> {
    let mut c = [1.; N];
    for i in 1..N {
        work.step(4)?;
        c[i] = c[i - 1] * (N - i) as f64 / i as f64;
    }
    Ok(c)
}
struct CapWeights {
    skew: [[I; 12]; 12],
    triple: [[[I; 12]; 12]; 12],
}
impl CapWeights {
    fn new(work: &mut Work<'_>) -> Result<Self> {
        let a = binomials::<12>(work)?;
        let b = binomials::<11>(work)?;
        let c = binomials::<22>(work)?;
        let d = binomials::<33>(work)?;
        // All binomial values and products here are integers below2^53.
        // Only division and subsequent interval algebra require rounding.
        let pair = |i: usize, j: usize| p(a[i] * b[j]) / p(22. * c[i + j]);
        let derivative = |i: usize, j: usize| {
            p(11.)
                * ((if j == 0 { p(0.) } else { pair(i, j - 1) })
                    - (if j == 11 { p(0.) } else { pair(i, j) }))
        };
        let mut skew = [[p(0.); 12]; 12];
        let mut triple = [[[p(0.); 12]; 12]; 12];
        for i in 0..12 {
            for j in 0..12 {
                work.step(32)?;
                skew[i][j] = derivative(i, j) - derivative(j, i);
                for (k, value) in triple[i][j].iter_mut().enumerate() {
                    work.step(32)?;
                    let term = |l: usize| p(a[i] * a[j] * b[l]) / p(33. * d[i + j + l]);
                    *value = p(11.)
                        * ((if k == 0 { p(0.) } else { term(k - 1) })
                            - (if k == 11 { p(0.) } else { term(k) }));
                }
            }
        }
        Ok(Self { skew, triple })
    }
    fn curve(&self, poles: [[I; 3]; 12], work: &mut Work<'_>) -> Result<[I; 3]> {
        let mut value = [p(0.); 3];
        for i in 0..12 {
            for j in 0..12 {
                work.step(8)?;
                value[0] = value[0] + p(0.5) * poles[i][0] * self.skew[i][j] * poles[j][1];
                for k in 0..12 {
                    work.step(16)?;
                    let w = p(0.5) * self.triple[i][j][k];
                    value[1] = value[1] + poles[i][0] * poles[j][0] * poles[k][1] * w;
                    value[2] = value[2] - poles[i][1] * poles[j][1] * poles[k][0] * w;
                }
            }
        }
        Ok(value)
    }
}
impl OwnedMaterial {
    /// Actual planar cap area and centroid from the original finite degree11
    /// boundary curves. Tensor faces require the separate surface integrator.
    pub fn cap_face_measure(
        &self,
        face: usize,
        options: FaceMeasureOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<FaceMeasure> {
        if !options.area_width_mm2.is_finite()
            || options.area_width_mm2 <= 0.
            || !options.centroid_width_mm.is_finite()
            || options.centroid_width_mm <= 0.
            || options.maximum_work == 0
            || options.maximum_work > 512_000_000
            || options.maximum_bytes == 0
            || options.maximum_bytes > 64 * 1024 * 1024
        {
            return Err(Error::Invalid);
        }
        // Includes all fixed weight tables, curve controls, intervals and
        // result storage, while retaining the immutable source arena once.
        let bytes = self
            .retained_budget_bytes()
            .checked_add(64 * 1024)
            .ok_or(Error::Limit)?;
        if bytes > options.maximum_bytes {
            return Err(Error::Limit);
        }
        let mut work = Work {
            left: options.maximum_work,
            used: 0,
            bytes,
            interrupt,
        };
        work.step(1)?;
        let owner = self.faces().get(face).ok_or(Error::Invalid)?;
        let Support::Plane { height } = owner.support() else {
            return Err(Error::Unsupported);
        };
        let weights = CapWeights::new(&mut work)?;
        let mut local = [p(0.); 3];
        for coedge in owner.coedges() {
            let curve = &self.0.curves[self.0.edges[coedge.edge_index()].curve_index()];
            for span in 0..curve.spans() {
                let poles = std::array::from_fn(|i| {
                    self.0.points[curve.point_indices()[11 * span + i]]
                        .coordinates
                        .each_ref()
                        .map(|s| interval(s.bounds))
                });
                let value = weights.curve(poles, &mut work)?;
                for i in 0..3 {
                    local[i] = if coedge.forward() {
                        local[i] + value[i]
                    } else {
                        local[i] - value[i]
                    };
                }
            }
        }
        if !owner.forward() {
            local = local.map(|v| -v);
        }
        if !local[0].valid() || local[0].lo <= 0. {
            return Err(Error::Unsupported);
        }
        let original_area = interval(self.0.local_areas[usize::from(owner.forward())]);
        if local[0].lo > original_area.hi || local[0].hi < original_area.lo {
            return Err(Error::Unsupported);
        }
        let z = interval(self.0.points[height].coordinates[2].bounds);
        let moments = [local[1], local[2], z * local[0]];
        let [u, v, _] = self.0.linear;
        let normal = [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ];
        let scale = (normal[0].square() + normal[1].square() + normal[2].square()).sqrt();
        let area = local[0] * scale;
        let anchor = self.world_enclosure_mm().map(|b| b[0].midpoint(b[1]));
        let centered_moments: [I; 3] = std::array::from_fn(|i| {
            scale
                * (self
                    .0
                    .linear
                    .iter()
                    .zip(moments)
                    .fold(p(0.), |sum, (column, m)| sum + column[i] * m)
                    + (self.0.origin[i] - p(anchor[i])) * local[0])
        });
        // A positive surface-density centroid is a convex combination of the
        // original support points, so it also lies in the proved source box.
        // Intersect that independent enclosure after interval division.
        let source_box = self.world_enclosure_mm();
        let centroid: [I; 3] = std::array::from_fn(|i| {
            let value = p(anchor[i]) + centered_moments[i] / area;
            I::new(
                value.lo.max(source_box[i][0]),
                value.hi.min(source_box[i][1]),
            )
        });
        if !area.valid()
            || area.lo <= 0.
            || area.width() > options.area_width_mm2
            || centroid
                .iter()
                .any(|v| !v.valid() || v.width() > options.centroid_width_mm)
        {
            return Err(Error::Unsupported);
        }
        work.step(1)?;
        Ok(FaceMeasure {
            source: self.clone(),
            face,
            area,
            centroid,
            centered_moments,
            anchor,
            work: work.used,
            bytes: work.bytes,
        })
    }
}
#[cfg(test)]
mod tests;
