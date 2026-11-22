//! Radial clearance between an authored external tapered pipe thread and a
//! pipe sleeve's authored cavity on opposing coaxial axes.
//!
//! Both laws are radial graphs `r = a(q) + k z` of the phase `q = z - theta
//! p / 2pi`. With opposing axes the cavity phase is `c0 - q` and its taper
//! term cancels the male's, so the radial gap depends on the phase alone.
//! It is enclosed from below on a complete partition of one pitch. Actual
//! boundaries deviate from their laws by recorded Euclidean bounds; a
//! Lipschitz radial transfer of those bounds, the axis misalignment and a
//! fixed tolerance reserve are subtracted before any clearance is admitted.
use super::*;
use multicad_rust_kernel::axial_pipe_thread::authored::AuthoredSection;

/// One authored law placed in the world, with its actual-boundary bound.
#[derive(Clone, Copy)]
struct Law {
    thread: AuthoredThread,
    frame: PlaneFrame,
    error_mm: f64,
}
const CELLS: usize = 4096;

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|i| a[i] * b[i]).sum()
}
fn compose(outer: PlaneFrame, inner: PlaneFrame) -> PlaneFrame {
    PlaneFrame {
        origin: position(outer, inner.origin),
        u: rotate(outer, inner.u),
        v: rotate(outer, inner.v),
    }
}
/// The external thread's authored law, without any stock crest limit: a
/// clip only removes male material, so the uncapped law bounds it above.
fn male(solid: &Solid) -> Option<Law> {
    let (thread, frame, error_mm) = if let Some(body) = solid.clipped_npt_body() {
        // The clipped body's material is its authored plan predicate itself.
        let (thread, pose) = body.authored_law();
        let r = pose.rotation;
        let frame = PlaneFrame {
            origin: pose.translation_mm,
            u: [r[0][0], r[1][0], r[2][0]],
            v: [r[0][1], r[1][1], r[2][1]],
        };
        (thread, frame, 0.)
    } else {
        let (source, occurrence) = solid.rust_authored_axial_pipe_thread_source()?;
        let (body, at) = solid.radial_boundary_body()?;
        let boundary = body.placed_boundary(at).ok()?;
        (
            source.transport,
            compose(occurrence, source.construction_frame),
            (source.local_source_error_mm + boundary.boundary_error_mm()).next_up(),
        )
    };
    let spec = thread.spec();
    if spec.radial_taper <= 0. {
        return None;
    }
    let uncapped = AuthoredThread::new(
        AuthoredSpec {
            crest_radius_limit_mm: None,
            ..spec
        },
        thread.tolerance_mm(),
        thread.minimum_feature_mm(),
    )
    .ok()?;
    Some(Law {
        thread: uncapped,
        frame,
        error_mm,
    })
}
/// Upper bound of `|grad r|` of an authored radial graph in Euclidean
/// cylindrical coordinates, over its whole length.
fn slope(thread: &AuthoredThread) -> Option<f64> {
    let spec = thread.spec();
    let k = spec.radial_taper.abs();
    let flank = match spec.section {
        AuthoredSection::NptFlat { .. } => 3_f64.sqrt(),
        // cot(27.5 deg) bounds every rounded-band slope of this section.
        AuthoredSection::BsptRounded { .. } => 1.921,
    };
    let depth = match spec.section {
        AuthoredSection::NptFlat { depth_mm, .. } | AuthoredSection::BsptRounded { depth_mm } => {
            depth_mm
        }
    };
    let minor = spec.major_radius_mm - depth - k * spec.length_mm;
    if minor.is_nan() || minor <= 0. {
        return None;
    }
    let along = flank + 2. * k;
    let around = (flank + k) * spec.pitch_mm / (2. * PI * minor);
    let g = along.portable_hypot(around) * (1. + 1e-12);
    g.is_finite().then_some(g)
}
/// Lower bound of the cavity minus male radius over the complete phase
/// domain, or `None` when an enclosure is unresolved.
fn gap(male: &AuthoredThread, cavity: &AuthoredThread, c0: I, station: f64) -> Option<f64> {
    let p = male.spec().pitch_mm;
    let mut lowest = f64::INFINITY;
    for j in 0..CELLS {
        let q0 = p * j as f64 / CELLS as f64;
        let q1 = if j + 1 == CELLS {
            p
        } else {
            p * (j + 1) as f64 / CELLS as f64
        };
        let outer = male.radial_bounds([q0, q1], [0., 0.]).ok()?[1];
        // Cavity phase c0 - q over this cell, wrapped into [0, p].
        let lo = (c0.lo - q1).next_down();
        let hi = (c0.hi - q0).next_up();
        let start = lo.rem_euclid(p);
        let width = hi - lo;
        if width.is_nan() || width < 0. || width >= p {
            return None;
        }
        let mut inner = f64::INFINITY;
        let mut cover = |a: f64, b: f64| -> Option<()> {
            let a = a.clamp(0., p);
            let b = b.clamp(0., p);
            inner = inner.min(cavity.radial_bounds([a, b], [0., 0.]).ok()?[0]);
            Some(())
        };
        if start + width <= p {
            cover(start, start + width)?;
        } else {
            cover(start, p)?;
            cover(0., start + width - p)?;
        }
        let value = (I::point(inner) - I::point(outer)
            + I::point(cavity.spec().radial_taper) * I::point(station))
        .lo;
        lowest = lowest.min(value);
    }
    lowest.is_finite().then_some(lowest)
}

/// Proven positive radial clearance of `male` inside the sleeve's cavity.
fn prove(male: Law, cavity: Law, height: f64, tolerance: f64) -> Option<f64> {
    let (m, c) = (male.thread.spec(), cavity.thread.spec());
    if m.pitch_mm != c.pitch_mm
        || c.radial_taper != -m.radial_taper
        || c.crest_radius_limit_mm.is_some()
        || std::mem::discriminant(&m.section) != std::mem::discriminant(&c.section)
    {
        return None;
    }
    let p = m.pitch_mm;
    let nm = male.frame.normal();
    let nc = cavity.frame.normal();
    let antiparallel = dot(nm, nc);
    if antiparallel.is_nan() || antiparallel > -1. + 1e-12 {
        return None;
    }
    // Small tilt angle of the opposing axes, bounded above.
    let tilt = (2. * (1. + antiparallel)).max(0.).sqrt() * (1. + 1e-9) + 1e-15;
    let offset: [f64; 3] = std::array::from_fn(|i| cavity.frame.origin[i] - male.frame.origin[i]);
    let station = dot(offset, nm);
    let radial: [f64; 3] = std::array::from_fn(|i| offset[i] - nm[i] * station);
    let lateral = dot(radial, radial).sqrt();
    let phi = dot(cavity.frame.u, male.frame.v).portable_atan2(dot(cavity.frame.u, male.frame.u));
    // Exact pose errors are carried by the misalignment reserve below; the
    // floating phase constant gets its own directed enclosure.
    let c0 = station - phi * p / (2. * PI);
    let rounding = 64. * f64::EPSILON * (station.abs() + p + 1.);
    let c0 = I::new(c0 - rounding, c0 + rounding);
    let reach = m.major_radius_mm.max(c.major_radius_mm)
        + (m.radial_taper.abs() * (m.length_mm + height))
        + station.abs()
        + m.length_mm
        + height;
    let misalignment = lateral + tilt * reach;
    let g = slope(&male.thread)?.max(slope(&cavity.thread)?);
    let transfer = 1. + 2. * g;
    let margin = transfer * (male.error_mm + cavity.error_mm + misalignment) + 64. * tolerance;
    let clearance = gap(&male.thread, &cavity.thread, c0, station)?;
    (margin.is_finite() && clearance > margin).then_some(clearance - margin)
}

impl Solid {
    /// An external authored tapered pipe thread clear of a sleeve's cavity
    /// shares no material with the sleeve: their intersection is empty and
    /// each difference is the unchanged left operand.
    pub(in crate::geometry) fn pipe_sleeve_clearance_boolean(
        &self,
        other: &Self,
        operation: BooleanOperation,
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        if policy.boolean_fuzzy_mm != 0.
            || self.policy != policy
            || other.policy != policy
            || operation == BooleanOperation::Union
        {
            return None;
        }
        let ((sleeve, frame), thread) = match (self.pipe_sleeve_body(), other.pipe_sleeve_body()) {
            (Some(sleeve), None) => (sleeve, other),
            (None, Some(sleeve)) => (sleeve, self),
            _ => return None,
        };
        let law = sleeve.law()?;
        let cavity = Law {
            thread: law.transport,
            frame,
            error_mm: law.error_mm,
        };
        prove(
            male(thread)?,
            cavity,
            sleeve.height_mm(),
            policy.linear_tolerance_mm,
        )?;
        Some(match operation {
            BooleanOperation::Subtract => self.share().map(Some),
            _ => Ok(None),
        })
    }
}
