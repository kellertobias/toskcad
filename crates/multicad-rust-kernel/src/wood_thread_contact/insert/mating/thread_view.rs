//! Distinct actual thread fields. A periodic repair never becomes a legacy
//! Thread, cached evaluator or original-source receipt.
use super::*;
use crate::step_read::PhillipsOwnedPeriodicMaterial as Periodic;
#[derive(Clone, Copy)]
pub(in crate::wood_thread_contact) enum MatingThread<'a> {
    Literal(&'a Thread),
    Periodic(&'a Periodic),
}
#[derive(Clone, Copy)]
pub(in crate::wood_thread_contact) struct Definition {
    pub pitch_mm: f64,
    pub starts: u32,
    pub left_hand: bool,
    pub length_mm: f64,
}
impl<'a> From<&'a Thread> for MatingThread<'a> {
    fn from(v: &'a Thread) -> Self {
        Self::Literal(v)
    }
}
pub(in crate::wood_thread_contact) struct Parts<'a> {
    pub thread: MatingThread<'a>,
    pub visible: [Q; 2],
    pub prefix: Q,
}
impl<'a> Parts<'a> {
    pub(super) fn literal(thread: &'a Thread, visible: [f64; 2], prefix: f64) -> R<Self> {
        Ok(Self {
            thread: thread.into(),
            visible: [q(visible[0])?, q(visible[1])?],
            prefix: q(prefix)?,
        })
    }
}
impl MatingThread<'_> {
    pub(in crate::wood_thread_contact) fn spec(self) -> Definition {
        match self {
            Self::Literal(t) => {
                let s = t.spec();
                Definition {
                    pitch_mm: s.pitch_mm,
                    starts: s.starts,
                    left_hand: s.left_hand,
                    length_mm: s.length_mm,
                }
            }
            Self::Periodic(t) => {
                let s = t.thread();
                Definition {
                    pitch_mm: s.pitch_mm(),
                    starts: 1,
                    left_hand: s.left_hand(),
                    length_mm: s.length_mm(),
                }
            }
        }
    }
    pub(in crate::wood_thread_contact) fn phase_rate(self) -> I {
        match self {
            Self::Literal(t) => ip(t.helical_rate()),
            Self::Periodic(t) => {
                let tau = ip(2.)
                    * I::new(
                        std::f64::consts::PI.next_down(),
                        std::f64::consts::PI.next_up(),
                    );
                ip(if t.thread().left_hand() { -1. } else { 1. }) * ip(t.thread().pitch_mm()) / tau
            }
        }
    }
    pub(in crate::wood_thread_contact) fn phase_joins(self) -> Vec<I> {
        match self {
            Self::Literal(t) => t.profile_breaks().map(ip).to_vec(),
            Self::Periodic(t) => {
                let nodes = t.thread().phase_node_expressions_mm();
                let mut p = vec![ip(0.)];
                p.extend(nodes[..4].iter().map(|p| ip(p[0]) + ip(p[1])));
                p.push(ip(t.thread().pitch_mm()));
                p
            }
        }
    }
    pub(in crate::wood_thread_contact) fn source_coordinate_bounds_mm(self) -> [[f64; 2]; 3] {
        match self {
            Self::Literal(t) => t.source_coordinate_bounds_mm(),
            Self::Periodic(t) => {
                let r = t.thread().radius_nodes_mm()[0];
                let cap = t.cap_expressions_mm().map(|p| ip(p[0]) + ip(p[1]));
                [[-r, r], [-r, r], [cap[0].lo, cap[1].hi]]
            }
        }
    }
    pub(in crate::wood_thread_contact) fn tip_bounds_mm(self) -> [[f64; 2]; 3] {
        let mut b = self.source_coordinate_bounds_mm();
        b[2] = match self {
            Self::Literal(t) => [t.spec().length_mm; 2],
            Self::Periodic(t) => {
                let p = t.cap_expressions_mm()[1];
                let z = ip(p[0]) + ip(p[1]);
                [z.lo, z.hi]
            }
        };
        b
    }
    pub(in crate::wood_thread_contact) fn radius_at_phase_interval(self, phase: I) -> Option<I> {
        match self {
            Self::Literal(t) => t.radius_at_phase_interval(phase),
            Self::Periodic(t) => period_radius(t, phase).ok(),
        }
    }
}
fn exact_expression(p: [f64; 2]) -> R<Q> {
    bounded(q(p[0])? + q(p[1])?)
}
pub(super) fn periodic_parts(source: &Periodic) -> R<Parts<'_>> {
    let chart = source.chart();
    // Normal authored material uses this exact map. Decoded composed maps
    // require a separate actual-map contact proof; no rounded composition is
    // substituted into the caller's literal frame.
    if chart.columns != [[0.; 3], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]
        || chart.outer_top_z != 0.
    {
        return Err(Error::Unsupported);
    }
    let caps = source.cap_expressions_mm();
    let visible = [exact_expression(caps[0])?, exact_expression(caps[1])?];
    let head = source.head().view();
    let end = head
        .station(head.segment_count())
        .expect("qualified head retains final endpoint");
    // Every repaired line/arc earned nonnegative axial progress and exact
    // shared endpoints. Thus the final endpoint is the WHOLE head/plain/neck
    // prefix's axial maximum, including zero derivative at an arc endpoint.
    let prefix = exact_expression([end[1], -source.head_z_shift_mm()])?;
    Ok(Parts {
        thread: MatingThread::Periodic(source),
        visible,
        prefix,
    })
}
pub(in crate::wood_thread_contact) fn periodic_prefix_bounds(
    source: &Periodic,
) -> R<[[f64; 2]; 3]> {
    let prefix = periodic_parts(source)?.prefix;
    let mut bounds = source.local_enclosure_bounds();
    bounds[2][1] = exact_interval(&prefix)?.hi;
    Ok(bounds)
}
fn period_radius(source: &Periodic, phase: I) -> R<I> {
    let knots = MatingThread::Periodic(source).phase_joins();
    let r = source.thread().radius_nodes_mm();
    // The crest crosses the phase origin. Split its SAME constant source law
    // at0/pitch; all other pieces use original shared phase/radius nodes.
    let radii = [r[0], r[0], r[1], r[2], r[3], r[4]];
    if !phase.valid() || phase.lo < 0. || phase.hi > source.thread().pitch_mm() {
        return Err(Error::Unsupported);
    }
    let mut result = I::new(f64::INFINITY, f64::NEG_INFINITY);
    for band in 0..5 {
        let p = I::new(
            phase.lo.max(knots[band].lo),
            phase.hi.min(knots[band + 1].hi),
        );
        if p.lo > p.hi {
            continue;
        }
        let denominator = positive(knots[band + 1] - knots[band])?;
        let u = (p - knots[band]) / denominator;
        let radius = finite(ip(radii[band]) + (ip(radii[band + 1]) - ip(radii[band])) * u)?;
        // An affine band stays between its SAME defining endpoint radii.
        let lo = radius.lo.max(radii[band].min(radii[band + 1]));
        let hi = radius.hi.min(radii[band].max(radii[band + 1]));
        result.lo = result.lo.min(lo);
        result.hi = result.hi.max(hi);
    }
    finite(result)
}
pub(super) fn periodic_bound(source: &Periodic, phase: I) -> (I, bool) {
    match periodic_radius(phase, source.thread().pitch_mm(), &mut |p| {
        period_radius(source, p)
    }) {
        Ok(v) => (v, false),
        Err(_) => (I::new(0., source.thread().radius_nodes_mm()[0]), true),
    }
}
