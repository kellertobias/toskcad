//! Ordered source events for the qualified three-dimensional cutter boundary.
//! The shared caches enclose immutable authored expressions, not independently
//! rounded or imported support intersections. Manufacturing remains separate.
use super::*;
use crate::tapered_thread::interval::J;
mod carrier;
mod evaluate;
pub(super) mod patches;
pub use carrier::{
    CubicTrimSpan, TRIM_SPAN_COMPARE_WORK, TRIM_SPAN_SOURCE_WORK, TrimErrorBound,
    TrimSpanComparison, TrimSpanCorrespondence,
};
pub use patches::{
    BoundaryMesh, FaceProperties, FacePropertyOptions, PatchAtlas, PatchBounds,
    PatchCorrespondence, PatchCorrespondenceOptions, PatchKey, QuadErrorBound, SurfacePatch,
};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(usize)]
pub enum EventRole {
    MouthEnvelope,
    EntryJoin,
    FloorJoin,
    WallJoin,
    WallContact,
    FloorApex,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum WallKind {
    Inner,
    Outer,
}
#[derive(Clone, Copy, Debug)]
pub struct SourceEvent {
    point: [I; 3],
    envelope: EnvelopeRole,
}
impl SourceEvent {
    pub fn point_mm_bounds(&self) -> [Bounds; 3] {
        self.point.map(bounds)
    }
    pub fn envelope_role(&self) -> EnvelopeRole {
        self.envelope
    }
}
#[derive(Clone, Copy, Debug)]
pub struct FiniteTrim {
    ends: [EventRole; 2],
    envelope: EnvelopeRole,
    wall: WallKind,
    parameters: [I; 2],
}
impl FiniteTrim {
    pub fn end_events(&self) -> [EventRole; 2] {
        self.ends
    }
    pub fn envelope_role(&self) -> EnvelopeRole {
        self.envelope
    }
    pub fn wall_kind(&self) -> WallKind {
        self.wall
    }
    /// Actual source parameter endpoints, possibly decreasing. Torus uses its
    /// meridian angle; cone intersections use positive transverse coordinate y.
    pub fn parameter_bounds(&self) -> [Bounds; 2] {
        self.parameters.map(bounds)
    }
}
/// A canonical positive-corner orbit. Both incident trim pieces reference the
/// same event object. Exact quarter rotations/reflections generate its copies.
#[derive(Clone, Debug)]
pub struct BoundaryEvents<'a> {
    source: &'a Cutter,
    mouth: SectionBoundary<'a>,
    events: [SourceEvent; 6],
    trims: [FiniteTrim; 4],
    pub(super) beta_tangent: I,
    pub(super) base: I,
    pub(super) cone_slopes: [I; 2],
    pub(super) cone_intercepts: [I; 2],
}
impl<'a> BoundaryEvents<'a> {
    pub(crate) fn retained_heap_bytes(&self) -> usize {
        self.mouth.retained_heap_bytes()
    }

    pub fn source(&self) -> &'a Cutter {
        self.source
    }
    pub fn mouth(&self) -> &SectionBoundary<'a> {
        &self.mouth
    }
    pub fn event(&self, role: EventRole) -> &SourceEvent {
        &self.events[role as usize]
    }
    pub fn trims(&self) -> &[FiniteTrim; 4] {
        &self.trims
    }
    pub fn event_orbit_point_mm_bounds(
        &self,
        role: EventRole,
        quarter: u8,
        lower: bool,
    ) -> Result<[Bounds; 3], Error> {
        if quarter > 3 {
            return Err(Error::Invalid);
        }
        let [x, y, z] = self.event(role).point;
        let [x, y] = rotate([x, if lower { -y } else { y }], quarter);
        Ok([bounds(x), bounds(y), bounds(z)])
    }
    /// Bounded point on the actual source trim, normalized from its first
    /// shared event to its last. Endpoint queries reuse the shared event cache.
    pub fn trim_point_mm_bounds(&self, index: usize, t: f64) -> Result<[Bounds; 3], Error> {
        if index >= 4 || !t.is_finite() || !(0. ..=1.).contains(&t) {
            return Err(Error::Invalid);
        }
        let trim = self.trims[index];
        if t == 0. || t == 1. {
            return Ok(self.event(trim.ends[t as usize]).point_mm_bounds());
        }
        evaluate::trim_state::<I>(self, index, I::point(t)).map(|p| p.map(bounds))
    }
}
fn checkpoint(cancel: &dyn Fn() -> bool) -> Result<(), Error> {
    if cancel() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn j(x: I) -> J {
    J {
        v: x,
        ..J::constant(0.)
    }
}
fn b([lo, hi]: Bounds) -> I {
    I::new(lo, hi)
}
fn cone_point(base: I, beta: I, slope: I, intercept: I, y: I, offset: I) -> Result<[I; 3], Error> {
    evaluate::cone_point(base, beta, slope, intercept, y, offset)
}

fn xy(source: &Cutter, radius: I, rb: I, wall: WallKind) -> Result<[I; 2], Error> {
    let k = &source.corners[0];
    let result = match wall {
        WallKind::Outer => {
            let slope = k.cos_theta / k.sin_theta / k.cos_beta;
            let offset = rb + (k.join_s - k.join_y * k.cos_theta / k.sin_theta) / k.cos_beta;
            volume::line_intersection(j(radius), j(offset), j(slope))
        }
        WallKind::Inner => match k.inner {
            SourceInner::ShallowV { sine, cosine } => {
                volume::line_intersection(j(radius), j(rb), j(cosine / sine / k.cos_beta))
            }
            SourceInner::Rounded { radius_mm } => {
                volume::rounded_intersection(j(radius), j(rb), radius_mm, k)
            }
        },
    }
    .ok_or(Error::Unresolved)?
    .map(|v| v.v);
    if result.into_iter().all(|v| v.valid() && v.lo > 0.) {
        Ok(result)
    } else {
        Err(Error::Unresolved)
    }
}
fn trim(
    ends: [EventRole; 2],
    envelope: EnvelopeRole,
    wall: WallKind,
    events: &[SourceEvent; 6],
    torus_end: I,
) -> FiniteTrim {
    let parameters = if envelope == EnvelopeRole::EntryRoundover {
        [I::point(0.), torus_end]
    } else {
        ends.map(|role| events[role as usize].point[1])
    };
    FiniteTrim {
        ends,
        envelope,
        wall,
        parameters,
    }
}
pub(super) fn construct<'a>(
    source: &'a Cutter,
    cancel: &dyn Fn() -> bool,
) -> Result<BoundaryEvents<'a>, Error> {
    checkpoint(cancel)?;
    volume::Sections::new(source)?;
    let e = &source.envelope;
    let k = &source.corners[0];
    let beta = k.sin_beta / k.cos_beta;
    let base = k.b_half + k.b_depth * beta;
    let slopes = [e.sin_a / e.cos_a, e.cos_c / e.sin_c];
    let intercepts = [
        I::point(source.spec.mouth_diameter_mm) / I::point(2.),
        e.tip * slopes[1],
    ];
    // These prove increasing envelope/wall residual with depth, since x<=R:
    // dF/dd=2*R*|R'|-2*x*tan(beta)>0. The torus has |R'|>=tan(a).
    for index in 0..2 {
        let q = beta / slopes[index];
        if !q.valid() || q.lo <= 0. || q.hi >= 1. || (base - q * intercepts[index]).lo <= 0. {
            return Err(Error::Unresolved);
        }
    }
    let rb_entry = base - e.join * beta;
    let rb_floor = base - e.flank * beta;
    let entry_radius = e.entry - I::point(source.spec.mouth_radius_mm) * e.cos_a;
    let floor_radius = I::point(source.spec.g_diameter_mm) / I::point(2.);
    let join_at_mouth = base + k.join_s / k.cos_beta;
    // The entire entry roundover meets the OUTER plane, including its endpoint:
    // the smallest envelope circle exceeds even the largest join front.
    if (entry_radius.square() - join_at_mouth.square() - k.join_y.square()).lo <= 0.
        || floor_radius.lo <= rb_floor.hi
    {
        return Err(Error::Unresolved);
    }
    let floor_join_test =
        floor_radius.square() - (rb_floor + k.join_s / k.cos_beta).square() - k.join_y.square();
    let floor_wall = if floor_join_test.lo > 0. {
        WallKind::Outer
    } else if floor_join_test.hi < 0. {
        WallKind::Inner
    } else {
        return Err(Error::Unresolved);
    };
    checkpoint(cancel)?;
    // Authored coupling: floorR(tb)=b/2 and rb(tb)=b/2 exactly. Retain
    // that common expression instead of solving the same event twice.
    let contact = [k.b_half, I::point(0.), k.b_depth];
    let join_envelope = if floor_wall == WallKind::Inner {
        EnvelopeRole::OuterFlank
    } else {
        EnvelopeRole::Floor
    };
    let join_index = if join_envelope == EnvelopeRole::OuterFlank {
        0
    } else {
        1
    };
    let join = cone_point(
        base,
        beta,
        slopes[join_index],
        intercepts[join_index],
        k.join_y,
        k.join_s / k.cos_beta,
    )?;
    let mouth = source.finite_section(0.)?;
    let mouth_vertex = mouth
        .vertices()
        .iter()
        .find(|v| v.source_identity() == (0, VertexRole::UpperEnvelope))
        .ok_or(Error::Unresolved)?;
    let [mx, my] = mouth_vertex.point_mm_bounds().map(b);
    let [ex, ey] = xy(source, entry_radius, rb_entry, WallKind::Outer)?;
    let [fx, fy] = xy(source, floor_radius, rb_floor, floor_wall)?;
    checkpoint(cancel)?;
    // Strict source order is also a representability proof for every finite trim.
    if contact[0].lo <= 0.
        || contact[2].lo <= e.flank.hi
        || contact[2].hi >= e.tip.lo
        || join[2].lo <= e.join.hi
        || join[2].hi >= contact[2].lo
        || my.lo <= ey.hi
        || ey.lo <= fy.hi.max(k.join_y.hi)
        || (floor_wall == WallKind::Outer && (fy.lo <= k.join_y.hi || join[2].lo <= e.flank.hi))
        || (floor_wall == WallKind::Inner && (fy.hi >= k.join_y.lo || join[2].hi >= e.flank.lo))
    {
        return Err(Error::Unresolved);
    }
    let events = [
        SourceEvent {
            point: [mx, my, I::point(0.)],
            envelope: EnvelopeRole::EntryRoundover,
        },
        SourceEvent {
            point: [ex, ey, e.join],
            envelope: EnvelopeRole::EntryRoundover,
        },
        SourceEvent {
            point: [fx, fy, e.flank],
            envelope: EnvelopeRole::OuterFlank,
        },
        SourceEvent {
            point: join,
            envelope: join_envelope,
        },
        SourceEvent {
            point: contact,
            envelope: EnvelopeRole::Floor,
        },
        SourceEvent {
            point: [I::point(0.), I::point(0.), e.tip],
            envelope: EnvelopeRole::Floor,
        },
    ];
    let torus_end = (I::point(90.) - I::point(source.spec.outer_flank_angle_deg))
        * I::new(PI.next_down(), PI.next_up())
        / I::point(180.);
    let middle = if floor_wall == WallKind::Inner {
        [
            (
                [EventRole::EntryJoin, EventRole::WallJoin],
                EnvelopeRole::OuterFlank,
                WallKind::Outer,
            ),
            (
                [EventRole::WallJoin, EventRole::FloorJoin],
                EnvelopeRole::OuterFlank,
                WallKind::Inner,
            ),
            (
                [EventRole::FloorJoin, EventRole::WallContact],
                EnvelopeRole::Floor,
                WallKind::Inner,
            ),
        ]
    } else {
        [
            (
                [EventRole::EntryJoin, EventRole::FloorJoin],
                EnvelopeRole::OuterFlank,
                WallKind::Outer,
            ),
            (
                [EventRole::FloorJoin, EventRole::WallJoin],
                EnvelopeRole::Floor,
                WallKind::Outer,
            ),
            (
                [EventRole::WallJoin, EventRole::WallContact],
                EnvelopeRole::Floor,
                WallKind::Inner,
            ),
        ]
    };
    let trims = [
        trim(
            [EventRole::MouthEnvelope, EventRole::EntryJoin],
            EnvelopeRole::EntryRoundover,
            WallKind::Outer,
            &events,
            torus_end,
        ),
        trim(middle[0].0, middle[0].1, middle[0].2, &events, torus_end),
        trim(middle[1].0, middle[1].1, middle[1].2, &events, torus_end),
        trim(middle[2].0, middle[2].1, middle[2].2, &events, torus_end),
    ];
    checkpoint(cancel)?;
    Ok(BoundaryEvents {
        source,
        mouth,
        events,
        trims,
        beta_tangent: beta,
        base,
        cone_slopes: slopes,
        cone_intercepts: intercepts,
    })
}
#[cfg(test)]
mod tests;

pub use patches::{
    MeshCurveKey, MeshPointKey, MeshRow, MouthFootpoint, MouthTraceOptions, TracedMouthMesh,
};
