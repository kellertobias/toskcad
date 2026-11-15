//! Distinct exact-periodic radial source. Cached Thread slopes/rate are metadata,
//! not the defining geometry. The four positive radial bands share their REAL
//! phase and radius nodes. The final phase is the expression first+pitch, so
//! theta=hand*2pi*(z-phase)/pitch closes by identity, not a floating comparison.
//! At fixed z, phase is strictly increasing through exactly one pitch and r>0:
//! every radial ray meets exactly once. Axial coordinates are literal z. This
//! proves a radial Jordan side and two disk caps for this thread source alone;
//! no head/cavity union or raw-file exact-material claim follows.
use super::*;
mod correspondence;
pub use correspondence::PeriodicThreadCorrespondence;
#[derive(Clone, Copy, Debug)]
pub struct PeriodicThreadSource {
    phase: [f64; 4],
    major: f64,
    root: f64,
    pitch: f64,
    length: f64,
    hand: f64,
}
#[derive(Clone, Copy, Debug)]
pub struct PeriodicThreadBounds {
    pub position_mm: [[f64; 2]; 3],
    /// Original band parameter u, then physical axial coordinate z.
    pub derivative_mm: [[[f64; 2]; 3]; 2],
}
impl PeriodicThreadSource {
    pub(super) fn new(thread: crate::helical_thread::Thread, work: &mut Work<'_>) -> Result<Self> {
        work.step(64)?;
        let [_, a, b, c, d, _] = thread.profile_breaks();
        let s = thread.spec();
        let source = Self {
            phase: [a, b, c, d],
            major: s.major_radius_mm,
            root: thread.root_radius_mm(),
            pitch: s.pitch_mm,
            length: s.length_mm,
            hand: if s.left_hand { -1. } else { 1. },
        };
        let p = source.phases();
        require(p.windows(2).all(|q| (q[1] - q[0]).lo > 0.) && source.root > 0.)?;
        Ok(source)
    }
    pub fn phase_node_expressions_mm(&self) -> [[f64; 2]; 5] {
        [
            [self.phase[0], 0.],
            [self.phase[1], 0.],
            [self.phase[2], 0.],
            [self.phase[3], 0.],
            [self.phase[0], self.pitch],
        ]
    }
    pub fn radius_nodes_mm(&self) -> [f64; 5] {
        [self.major, self.root, self.root, self.major, self.major]
    }
    pub fn pitch_mm(&self) -> f64 {
        self.pitch
    }
    pub fn length_mm(&self) -> f64 {
        self.length
    }
    pub fn left_hand(&self) -> bool {
        self.hand < 0.
    }
    pub(super) fn phases(&self) -> [I; 5] {
        self.phase_node_expressions_mm()
            .map(|p| I::point(p[0]) + I::point(p[1]))
    }
    pub(super) fn radii(&self, band: usize, u: I) -> I {
        let r = self.radius_nodes_mm();
        I::point(r[band]) + u * (I::point(r[band + 1]) - I::point(r[band]))
    }
    pub fn bounds(
        &self,
        band: usize,
        u: [f64; 2],
        z_mm: [f64; 2],
        maximum_work: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> std::result::Result<PeriodicThreadBounds, PhillipsReadError> {
        let mut work = Work {
            remaining: maximum_work,
            interrupt,
            cancelled: Cell::new(false),
        };
        let result = (|| {
            work.step(256)?;
            let u = I::new(u[0], u[1]);
            let z = I::new(z_mm[0], z_mm[1]);
            if band >= 4
                || !u.valid()
                || u.lo < 0.
                || u.hi > 1.
                || !z.valid()
                || z.lo < 0.
                || z.hi > self.length
                || maximum_work > WORK
            {
                return Err(ReadError::Invalid);
            }
            let rotation = rotation::Rotation::new(&mut work)?;
            let p = self.phases();
            let width = p[band + 1] - p[band];
            let phase = p[band] + u * width;
            let rate =
                I::point(self.hand) * I::point(2.) * rotation.pi(&mut work)? / I::point(self.pitch);
            let angle = rate * (z - phase);
            let [c, s] = if angle.width() > 1. {
                [I::new(-1., 1.); 2]
            } else {
                rotation.sin_cos(angle, &mut work)?
            };
            let r = self.radii(band, u);
            let nodes = self.radius_nodes_mm();
            let dr = I::point(nodes[band + 1]) - I::point(nodes[band]);
            let du = -rate * width;
            let position = [r * c, r * s, z];
            let derivative = [
                [dr * c - r * s * du, dr * s + r * c * du, I::point(0.)],
                [-r * s * rate, r * c * rate, I::point(1.)],
            ];
            Ok(PeriodicThreadBounds {
                position_mm: position.map(|x| [x.lo, x.hi]),
                derivative_mm: derivative.map(|d| d.map(|x| [x.lo, x.hi])),
            })
        })();
        result.map_err(|e| bound_error(e, &work))
    }
}
/// Same original Body/finite-field ownership, freshly compared with the exact
/// periodic shared-node source. This type cannot be created by the cached-rate
/// field-only producer and does not itself qualify a whole screw material.
pub struct PeriodicThreadRepair<'a, 'input> {
    field: AnalyticThreadRepair<'a, 'input>,
}
impl<'a, 'input> PeriodicThreadRepair<'a, 'input> {
    pub fn source(&self) -> &PeriodicThreadSource {
        self.field
            .periodic
            .as_ref()
            .expect("periodic producer source")
    }
    pub fn body(&self) -> &'a Body<'input> {
        self.field.body()
    }
    pub fn frame(&self) -> [Point; 3] {
        self.field.frame()
    }
    pub fn transfers(&self) -> &[AnalyticThreadTransfer] {
        self.field.transfers()
    }
    pub fn new_work_used(&self) -> usize {
        self.field.new_work_used()
    }
    pub fn inherited_head_work_used(&self) -> usize {
        self.field.inherited_head_work_used()
    }
    pub fn total_work_used(&self) -> usize {
        self.field.total_work_used()
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.field.retained_budget_bytes()
    }
}
impl<'h, 'q, 'r, 'a, 'input> ClosedRepairedHead<'h, 'q, 'r, 'a, 'input> {
    pub fn periodic_thread_repair(
        &self,
        options: AnalyticThreadOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> std::result::Result<PeriodicThreadRepair<'a, 'input>, PhillipsReadError> {
        self.build_thread_repair(options, interrupt, true, None)
            .map(|field| PeriodicThreadRepair { field })
    }
}
#[cfg(test)]
mod tests;
