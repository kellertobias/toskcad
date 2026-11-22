//! Shared-ring tessellation of the authoritative clipped radial solid.
use super::{
    interval::I,
    phase::{Knot, Segment},
    *,
};
use crate::display::{DisplayMesh, MeshError, Tessellation};
mod builder;
mod cell;
mod geometry;
use geometry::Bounds;
#[derive(Clone, Copy)]
struct Sample {
    phase: I,
    nominal: f64,
    threshold: bool,
    glued: bool,
}
#[derive(Clone, Copy)]
struct Leaf {
    segment: Segment,
    ends: [Sample; 2],
    dq: I,
    thread: Bounds,
    cone: Option<Bounds>,
    trim: [f64; 2],
    normals: [geometry::Normal; 2],
}
struct Plan {
    leaves: Vec<Leaf>,
    angular_step: f64,
    angular_origin: f64,
    angular_count: usize,
}
impl Segment {
    fn sample(self, u: f64) -> Sample {
        if u == 0. {
            return self.ends[0].into();
        }
        if u == 1. {
            return self.ends[1].into();
        }
        Sample {
            phase: (I::point(1.) - I::point(u)) * self.ends[0].phase
                + I::point(u) * self.ends[1].phase,
            nominal: (self.ends[1].nominal - self.ends[0].nominal).mul_add(u, self.ends[0].nominal),
            threshold: false,
            glued: false,
        }
    }
}
impl From<Knot> for Sample {
    fn from(k: Knot) -> Self {
        Self {
            phase: k.phase,
            nominal: k.nominal,
            threshold: k.threshold,
            glued: k.glued,
        }
    }
}
/// Indexed manufacturing boundary. All vertices and normals are checked
/// against the retained analytic material, including shared trim vertices.
pub fn tessellate(body: &Chamfer, options: Tessellation) -> Result<DisplayMesh, MeshError> {
    tessellate_with_interrupt(body, options, None).map_err(|e| match e {
        Error::Invalid => MeshError::InvalidOptions,
        Error::Limit | Error::Cancelled => MeshError::BudgetExceeded,
        Error::Unresolved => MeshError::UnresolvedBoundary,
    })
}
pub fn tessellate_with_interrupt(
    body: &Chamfer,
    options: Tessellation,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<DisplayMesh, Error> {
    options.check().map_err(|_| Error::Invalid)?;
    let mut work = Work {
        remaining: 16_000_000,
        interrupt,
    };
    work.step()?;
    let floor = (I::point(body.tool.spec().major_radius_mm) * I::point(1.)
        + I::point(body.full_length_mm()) * I::point(1.))
    .width();
    if options.deflection * 0.125 <= floor {
        return Err(Error::Unresolved);
    }
    let segments = mesh_segments(body, options)?;
    for segment in &segments {
        let error = body.glue_residual(*segment)?;
        if error[0] > options.deflection * 0.125 || error[1] > options.angle * 0.1 {
            return Err(Error::Unresolved);
        }
    }
    let plan = make_plan(body, &segments, options, 1., &mut work)?;
    build_adaptive(body, plan, options, &mut work, 8)
}
fn build_adaptive(
    body: &Chamfer,
    mut plan: Plan,
    options: Tessellation,
    work: &mut Work<'_>,
    passes: usize,
) -> Result<DisplayMesh, Error> {
    for _ in 0..passes {
        match builder::build(body, &plan, options, work) {
            Ok(mesh) => return Ok(mesh),
            Err(builder::Failure::Geometry(cause)) => {
                return Err(cause);
            }
            Err(builder::Failure::Refine(failed)) => {
                let mut refined = Vec::with_capacity(plan.leaves.len() + failed.len());
                for (index, leaf) in plan.leaves.iter().enumerate() {
                    work.step()?;
                    let ends = leaf.ends.map(|sample| Knot {
                        phase: sample.phase,
                        nominal: sample.nominal,
                        threshold: sample.threshold,
                        glued: sample.glued,
                    });
                    let source = Segment {
                        ends,
                        ..leaf.segment
                    };
                    if failed.binary_search(&index).is_ok() {
                        let middle = source.sample(0.5);
                        if middle.nominal <= ends[0].nominal || middle.nominal >= ends[1].nominal {
                            return Err(Error::Unresolved);
                        }
                        let knot = Knot {
                            phase: middle.phase,
                            nominal: middle.nominal,
                            threshold: false,
                            glued: false,
                        };
                        refined.push(Segment {
                            ends: [ends[0], knot],
                            ..source
                        });
                        refined.push(Segment {
                            ends: [knot, ends[1]],
                            ..source
                        });
                    } else {
                        refined.push(source);
                    }
                }
                plan = make_plan(body, &refined, options, 1., work)?;
            }
        }
    }
    Err(Error::Unresolved)
}
fn mesh_segments(body: &Chamfer, options: Tessellation) -> Result<Vec<Segment>, Error> {
    let mut segments = body.phase_segments()?;
    if let Some(root) = body.tool.root_circle() {
        let k = body.cone_slope()?;
        let critical = (I::point(root.phase)
            + I::point(root.radius) * k / (I::point(1.) + k.square()).sqrt())
        .checked()?;
        if !segments.iter().any(|segment| {
            segment.band == 2
                && segment.active
                && critical.lo > segment.ends[0].phase.hi
                && critical.hi < segment.ends[1].phase.lo
        }) {
            return Ok(segments);
        }
        let q = critical.lo.midpoint(critical.hi);
        let resolution = (I::point(options.deflection) * I::point(0.125)
            / (I::point(1.)
                + I::point(body.tool.spec().major_radius_mm) / body.rate_interval()?.abs()
                + I::point(body.profile(2, I::point(q))?.first.magnitude())))
        .checked()?
        .lo;
        let mut split = Vec::new();
        for segment in segments {
            let left = (critical - segment.ends[0].phase).lo;
            let right = (segment.ends[1].phase - critical).lo;
            if segment.band == 2 && segment.active && left > resolution && right > resolution {
                let knot = Knot {
                    phase: I::point(q),
                    nominal: q,
                    threshold: false,
                    glued: false,
                };
                split.push(Segment {
                    ends: [segment.ends[0], knot],
                    ..segment
                });
                split.push(Segment {
                    ends: [knot, segment.ends[1]],
                    ..segment
                });
            } else {
                split.push(segment);
            }
        }
        segments = split;
    }
    Ok(segments)
}
fn make_plan(
    body: &Chamfer,
    segments: &[Segment],
    options: Tessellation,
    scale: f64,
    work: &mut Work<'_>,
) -> Result<Plan, Error> {
    let mut leaves = Vec::new();
    let mut stack = Vec::new();
    let phase_budget = options.deflection * 0.25 * scale;
    let normal_budget = options.angle * 0.4 * scale;
    let rate = body.rate_interval()?.abs();
    for &segment in segments {
        stack.push([0., 1.]);
        while let Some(u) = stack.pop() {
            work.step()?;
            let ends = u.map(|u| segment.sample(u));
            if ends[0].nominal >= ends[1].nominal {
                return Err(Error::Unresolved);
            }
            let q = I::new(ends[0].phase.lo, ends[1].phase.hi).checked()?;
            let dq = (ends[1].phase - ends[0].phase).positive().checked()?;
            let thread = body.patch_bounds(segment, q, false)?;
            let cone = if segment.active {
                Some(body.patch_bounds(segment, q, true)?)
            } else {
                None
            };
            let profile = body.profile(segment.band, q)?;
            let cap = I::point(profile.second.magnitude())
                + I::point(2. * profile.first.magnitude()) / rate
                + I::point(body.tool.spec().major_radius_mm) / rate.square();
            let q_curvature = thread.hessian[0].max(cap.hi);
            let q_normal = thread.angular_q.max((I::point(1.) / rate).hi);
            let chord = I::point(q_curvature) * dq.square() / I::point(8.);
            let angle = I::point(q_normal) * dq;
            let trim = body.trim_error(segment, q, dq)?;
            if chord.hi > phase_budget
                || angle.hi > normal_budget
                || trim[0] > options.deflection * 0.125
                || trim[1] > options.angle * 0.1
            {
                let mid = u[0].midpoint(u[1]);
                let center = segment.sample(mid);
                if center.nominal <= ends[0].nominal || center.nominal >= ends[1].nominal {
                    return Err(Error::Unresolved);
                }
                if 4 * (leaves.len() + stack.len() + 2) as u64 * u64::from(body.tool.spec().starts)
                    > options.triangle_limit as u64
                {
                    return Err(Error::Limit);
                }
                stack.push([mid, u[1]]);
                stack.push([u[0], mid]);
                continue;
            }
            let qmid = (ends[0].phase + ends[1].phase) / I::point(2.);
            leaves.push(Leaf {
                segment,
                ends,
                dq,
                thread,
                cone,
                trim,
                normals: [
                    geometry::normal_field(body, segment, qmid, false)?,
                    geometry::normal_field(body, segment, qmid, true)?,
                ],
            });
        }
    }
    let curvature = I::point(body.tool.spec().major_radius_mm)
        + I::point(2.) * body.cone_slope()?.abs() * body.rate_interval()?.abs();
    let step = (I::point(options.deflection) / (I::point(2.) * curvature))
        .sqrt()
        .checked()?
        .lo
        .min(options.angle * 0.2);
    if step <= 0. {
        return Err(Error::Unresolved);
    }
    let cycle = body.tool.spec().pitch_mm * f64::from(body.tool.spec().starts)
        / body.tool.helical_rate().abs();
    let count = (cycle / step).ceil();
    if !count.is_finite() || count > options.triangle_limit as f64 {
        return Err(Error::Limit);
    }
    let angular_step = cycle / count;
    // The angular grid is a sampling proposal, not a source boundary. Select
    // its phase to avoid near-coincident artificial cuts at real source knots.
    // Every resulting facet still passes the same full source certificates.
    let mut angular_origin = 0.;
    let mut best = -1.;
    for fraction in [
        0.123456789,
        0.3819660112501051,
        0.6180339887498949,
        0.876543211,
    ] {
        let mut clearance = 0.5_f64;
        for period in 0..body.tool.spec().starts {
            for leaf in &leaves {
                work.step()?;
                let q = leaf.ends[0].nominal + f64::from(period) * body.tool.spec().pitch_mm;
                let seam = body.seam_station(leaf.ends[0].nominal);
                for z in [0., seam, body.material_length_mm()] {
                    let row = (z - q) / body.tool.helical_rate().abs() / angular_step - fraction;
                    clearance = clearance.min((row - row.round()).abs());
                }
            }
        }
        if clearance > best {
            best = clearance;
            angular_origin = fraction * angular_step;
        }
    }
    Ok(Plan {
        leaves,
        angular_step,
        angular_origin,
        angular_count: count as usize,
    })
}

#[cfg(test)]
mod tests;
