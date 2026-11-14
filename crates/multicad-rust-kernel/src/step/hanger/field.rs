//! Canonical actual-source ownership for periodic profile seam points.
//! These are approximation points, not a positive feature/sewn-body proof.
use super::*;
use crate::authored_hanger_bolt::{Boundary, Component, Knot};
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum PhaseNode {
    Knot(Knot),
    /// Computational split inside the constant-radius crest; not a new
    /// physical source face/edge or a replacement for the original knots.
    CrestSplit,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
/// Identity inside ONE SourceField owner. A compound/occurrence registry must
/// additionally retain that actual source and placement namespace.
pub struct StationKey {
    pub component: Component,
    pub total_period: i32,
    pub node: PhaseNode,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SeamPoint {
    pub key: StationKey,
    pub point_mm: [f64; 3],
    pub source_bounds_mm: [[f64; 2]; 3],
    pub source_error_mm: f64,
    pub work_used: usize,
}
#[derive(Clone, Copy)]
struct Cover {
    starts: u32,
    last_cell: i32,
    last_turn: i32,
}
pub struct SourceField<'a> {
    boundary: &'a Boundary,
    covers: [Cover; 2],
    work: usize,
    bytes: usize,
}
fn index(c: Component) -> usize {
    match c {
        Component::Metric => 0,
        Component::Wood => 1,
    }
}
// Both sealed component constructors require .001<=p<=L<=1e7. Their
// integer_decode exponent is in[-62,-29], so this exact common dyadic scale
// fits below2^86. Multiplying by <=8 starts fits below2^89; ceil numerator
// stays below2^87. No float quotient is used to decide a genuine last cell.
fn scaled_dimension(value: f64) -> u128 {
    let (m, e, _) = num_traits::Float::integer_decode(value);
    u128::from(m) << ((e + 62) as u32)
}
fn cover(pitch: f64, length: f64, starts: u32) -> Cover {
    let p = scaled_dimension(pitch);
    let l = scaled_dimension(length);
    let lead = p * u128::from(starts);
    // Source admission bounds the rounded quotient L/p+2N<=128; exact
    // quotient is consequently <127 and these exact ceilings <=127.
    Cover {
        starts,
        last_cell: l.div_ceil(p) as i32,
        last_turn: l.div_ceil(lead) as i32,
    }
}
impl<'a> SourceField<'a> {
    pub fn prepare(
        boundary: &'a Boundary,
        options: DomainOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, Error> {
        let bytes = boundary.retained_bytes() + size_of::<Self>();
        let mut w = super::domain::work(options, bytes, interrupt)?;
        w.charge(512)?;
        let source = boundary.source();
        let m = source.metric().spec();
        let wood = source.wood().spec();
        let covers = [
            cover(m.pitch_mm, m.length_mm, m.starts),
            cover(wood.pitch_mm, wood.length_mm, wood.starts),
        ];
        w.charge(1)?;
        Ok(Self {
            boundary,
            covers,
            work: w.used,
            bytes,
        })
    }
    pub fn source_boundary(&self) -> &'a Boundary {
        self.boundary
    }
    /// One physical-v support cover for EVERY source band of this component,
    /// including both genuine crest halves. Later finite face trims remove
    /// the support outside original z0..L; it is not extra material.
    pub fn physical_turns(&self, c: Component) -> [i32; 2] {
        [-1, self.covers[index(c)].last_turn]
    }
    pub fn cell_periods(&self, c: Component) -> [i32; 2] {
        let x = self.covers[index(c)];
        [-(x.starts as i32), x.last_cell]
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_bytes(&self) -> usize {
        self.bytes
    }
    /// Source cell(j,v=1) and cell(j+N,v=0) use the SAME total-period
    /// expression/key/point, not merely overlapping numeric intervals.
    pub fn seam_point(
        &self,
        c: Component,
        cell: i32,
        upper_turn: bool,
        node: PhaseNode,
        options: DomainOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<SeamPoint, Error> {
        let mut w = super::domain::work(options, self.bytes, interrupt)?;
        let cover = self.covers[index(c)];
        if cell < -(cover.starts as i32) || cell > cover.last_cell {
            return Err(Error::Invalid);
        }
        let total = cell + if upper_turn { cover.starts as i32 } else { 0 };
        w.charge(512)?;
        let s = self.boundary.source();
        let (pitch, phase, radius) = match c {
            Component::Metric => {
                let m = s.metric();
                let spec = m.spec();
                let (q, r) = match node {
                    PhaseNode::CrestSplit => (p(0.), p(spec.major_radius_mm)),
                    PhaseNode::Knot(k) => {
                        use crate::authored_metric_thread::Knot as K;
                        let k = match k {
                            Knot::CrestEnd => K::CrestEnd,
                            Knot::UpperRoot => K::UpperRoot,
                            Knot::LowerRoot => K::LowerRoot,
                            Knot::CrestStart => K::CrestStart,
                        };
                        let v = m.knot(k);
                        (hull(v.phase_mm), hull(v.radius_mm))
                    }
                };
                (spec.pitch_mm, q, r)
            }
            Component::Wood => {
                let a = s.wood();
                let (q, r) = match node {
                    PhaseNode::CrestSplit => (p(0.), p(a.spec().major_radius_mm)),
                    PhaseNode::Knot(k) => {
                        let j = match k {
                            Knot::CrestEnd => 0,
                            Knot::UpperRoot => 1,
                            Knot::LowerRoot => 2,
                            Knot::CrestStart => 3,
                        };
                        let b = a.bands()[j];
                        (p(b.phase_mm[1]), p(b.radius_mm[1]))
                    }
                };
                (a.spec().pitch_mm, q, r)
            }
        };
        let station = phase + p(pitch) * p(f64::from(total));
        let z = match c {
            Component::Metric => {
                station - p(s.metric().full_length_mm()) - p(s.spec().plain_length_mm)
            }
            Component::Wood => station,
        };
        // The source closes by TRUE2pi and PeriodicAngular owns exactly the
        // same literal xy endpoint. Never call a rounded sin(2pi*j) here.
        let point = [radius.midpoint(), 0., z.midpoint()];
        let error = (p((radius - p(point[0])).magnitude()) + p((z - p(point[2])).magnitude())).hi;
        w.charge(1)?;
        Ok(SeamPoint {
            key: StationKey {
                component: c,
                total_period: total,
                node,
            },
            point_mm: point,
            source_bounds_mm: [[radius.lo, radius.hi], [0., 0.], [z.lo, z.hi]],
            source_error_mm: error,
            work_used: w.used,
        })
    }
}
#[cfg(test)]
mod tests;
