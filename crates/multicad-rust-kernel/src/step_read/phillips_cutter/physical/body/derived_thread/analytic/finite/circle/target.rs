//! Constant-radius arc tied to the same exact periodic phase-node expressions.
//! A shoulder complement selects0→3; a side crest selects3→4. This is geometry,
//! not an inferred authoring recipe, raw trimming or owning-graph admission.
use super::*;
#[derive(Clone, Copy)]
pub(in super::super::super) struct ArcTarget {
    source: PeriodicThreadSource,
    nodes: [usize; 2],
    radius: f64,
    z: f64,
}
impl ArcTarget {
    pub(in super::super::super) fn new(
        source: PeriodicThreadSource,
        nodes: [usize; 2],
        radius: f64,
        z: f64,
        work: &mut Work<'_>,
    ) -> FResult<Self> {
        work.step(64)?;
        if nodes[0] >= nodes[1]
            || nodes[1] > 4
            || !radius.is_finite()
            || radius <= 0.
            || radius > 1e7
            || !z.is_finite()
            || z < 0.
            || z > source.length_mm()
        {
            return Err(PhillipsReadError::Unsupported);
        }
        Ok(Self {
            source,
            nodes,
            radius,
            z,
        })
    }
    pub(in super::super::super) fn radius(&self) -> f64 {
        self.radius
    }
    pub(in super::super::super) fn z(&self) -> f64 {
        self.z
    }
    pub(in super::super::super) fn left_hand(&self) -> bool {
        self.source.left_hand()
    }
    pub(in super::super::super) fn angular_span(&self, pi: I) -> I {
        let p = self.source.phases();
        I::point(2.) * pi * (p[self.nodes[1]] - p[self.nodes[0]]) / I::point(self.source.pitch_mm())
    }
    pub(in super::super::super) fn endpoint(
        &self,
        end: bool,
        work: &mut Work<'_>,
    ) -> FResult<[I; 3]> {
        work.step(256)?;
        let rotation = rotation::Rotation::new(work)?;
        let pi = rotation.pi(work)?;
        let phase = self.source.phases()[self.nodes[usize::from(end)]];
        let hand = I::point(if self.source.left_hand() { -1. } else { 1. });
        let angle = hand * I::point(2.) * pi * (I::point(self.z) - phase)
            / I::point(self.source.pitch_mm());
        let [c, s] = rotation.sin_cos(angle, work)?;
        Ok([
            I::point(self.radius) * c,
            I::point(self.radius) * s,
            I::point(self.z),
        ])
    }
}
