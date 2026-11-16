//! Source-qualified surface measures transferred from an independently
//! integrated proposal. Cross-product and first-moment errors are integrated
//! over every original parameter domain; no proposal scalar is relabelled.
use super::*;
use crate::phillips_drive::{FacePropertyOptions, ScrewFace};
#[derive(Clone, Copy, Debug)]
pub struct PeriodicFaceProperties {
    pub source_face: usize,
    pub area_mm2: [f64; 2],
    pub first_moment_mm3: [[f64; 2]; 3],
    pub centroid_mm: [[f64; 2]; 3],
    pub area_transfer_error_mm2: f64,
    pub first_moment_transfer_error_mm3: f64,
}
fn bounds(i: I) -> [f64; 2] {
    [i.lo, i.hi]
}
fn widen(i: I, e: f64) -> I {
    i + I::new(-e, e)
}
fn magnitude(a: [I; 3]) -> I {
    a.into_iter()
        .fold(I::point(0.), |s, v| s + I::point(v.magnitude()))
}
impl PeriodicAuthoredApproximation<'_> {
    /// The budget reserves the entire delegated integration cell allowance,
    /// including unsuccessful cells, before invoking the existing integrator.
    pub fn source_face_properties(
        &self,
        options: FacePropertyOptions,
        maximum_work: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> FResult<Vec<PeriodicFaceProperties>> {
        let mut work = Work {
            remaining: maximum_work,
            interrupt,
            cancelled: Cell::new(false),
        };
        let result = (|| {
            work.step(1)?;
            if maximum_work > WORK || options.maximum_cells == 0 || options.maximum_cells > 65536 {
                return Err(PhillipsReadError::Invalid);
            }
            let owner = self.owner;
            let options_values = [
                options.area_tolerance_mm2,
                options.first_moment_tolerance_mm3,
                options.centroid_tolerance_mm,
            ];
            if options_values
                .into_iter()
                .any(|v| !v.is_finite() || v <= 0.)
            {
                return Err(PhillipsReadError::Invalid);
            }
            work.step(options.maximum_cells)?;
            // The delegated integrator already allocates its own strict
            // per-domain tolerances. Tightening those a second time can reach
            // its genuine cached-phase precision floor. Acceptance below is
            // always against the original final source options after transfer.
            let inner = options;
            let proposal = self
                .proposal
                .face_properties_with_cancel(inner, &|| interrupt.is_some_and(|f| f()))
                .map_err(|e| match e {
                    crate::phillips_drive::Error::Invalid => PhillipsReadError::Invalid,
                    crate::phillips_drive::Error::Limit => PhillipsReadError::Limit,
                    crate::phillips_drive::Error::Cancelled => PhillipsReadError::Cancelled,
                    crate::phillips_drive::Error::Unresolved => PhillipsReadError::Unsupported,
                })?;
            let rotation = rotation::Rotation::new(&mut work)?;
            let tau = I::point(2.) * rotation.pi(&mut work)?;
            let omega = tau / I::point(owner.thread.pitch_mm());
            let old_omega = I::point(std::f64::consts::TAU) / I::point(owner.thread.pitch_mm());
            let angular_gap = (tau - I::point(std::f64::consts::TAU)).magnitude();
            let phases = owner.thread.phases();
            let radii = owner.thread.radius_nodes_mm();
            let caps = owner.caps.map(expression);
            let old_caps = self.proposal.visible_thread_interval_mm().map(I::point);
            let station_gap = I::point((caps[0] - old_caps[0]).magnitude())
                + I::point((caps[1] - old_caps[1]).magnitude());
            let height = caps[1] - caps[0];
            let extent = owner
                .enclosure
                .into_iter()
                .map(interval)
                .map(I::magnitude)
                .fold(1_f64, f64::max);
            let mut side_errors = [(0., 0.); 4];
            let mut cap_area_error = I::point(0.);
            let mut cap_position = 0_f64;
            for band in 0..4 {
                work.step(2048)?;
                let dp = phases[band + 1] - phases[band];
                let dr = I::point(radii[band + 1]) - I::point(radii[band]);
                let old_slope = I::point(
                    self.proposal
                        .thread()
                        .slope_at_phase(phases[band].midpoint(), band + 1),
                );
                let anchor = if band == 2 { phases[3] } else { phases[band] };
                let old_radius = if band == 1 { radii[1] } else { radii[0] };
                let old = [phases[band], phases[band + 1]]
                    .map(|p| I::point(old_radius) + old_slope * (p - anchor));
                let radial_gap = (old[0] - I::point(radii[band]))
                    .magnitude()
                    .max((old[1] - I::point(radii[band + 1])).magnitude());
                let radius = I::point(radii[band].max(radii[band + 1]));
                let old_radius = I::point(old[0].magnitude().max(old[1].magnitude()));
                // Include the one wrap of the original split crest and BOTH
                // finite endpoint ranges. Real periodic closure uses identity;
                // the old integrator's cached-angle sewing reserve is retained.
                let range =
                    I::point(
                        caps.into_iter()
                            .chain(old_caps)
                            .map(I::magnitude)
                            .fold(0_f64, f64::max),
                    ) + I::point(phases.into_iter().map(I::magnitude).fold(0_f64, f64::max));
                let angle_error =
                    I::point((omega - old_omega).magnitude()) * range + I::point(angular_gap);
                let old_dr = old_slope * dp;
                let du_error = I::point((dr - old_dr).magnitude())
                    + I::point(old_dr.magnitude()) * angle_error
                    + I::point(radial_gap) * omega * dp
                    + old_radius
                        * (I::point((omega - old_omega).magnitude()) * dp
                            + old_omega * dp * angle_error);
                let dz_error = I::point(radial_gap) * omega
                    + old_radius
                        * (I::point((omega - old_omega).magnitude()) + old_omega * angle_error);
                let position = I::point(radial_gap) + old_radius * angle_error;
                let old_du = I::point(old_dr.magnitude()) + old_radius * old_omega * dp;
                let new_dz = radius * omega + I::point(1.);
                let new_du = I::point(dr.magnitude()) + radius * omega * dp;
                let jacobian_error = du_error * new_dz + old_du * dz_error;
                let area_error = jacobian_error * height
                    + (new_du * new_dz + old_du * (old_radius * old_omega + I::point(1.)))
                        * station_gap;
                side_errors[band] = (area_error.hi, position.hi);
                // Radial disk fan X(u,v)=v*r(u)*e(theta(u)), v in[0,1].
                // Du difference <=v*du_error; Dv difference <=position.
                cap_area_error =
                    cap_area_error + (du_error * radius + old_du * position) / I::point(2.);
                cap_position = cap_position.max(position.hi);
            }
            let head = owner.head.view();
            let head_transfer = head
                .proposal_transfer(work.remaining.min(WORK), interrupt)
                .map_err(
                    crate::step_read::phillips_cutter::physical::body::closed_head::repaired_error,
                )?;
            work.step(head_transfer.work_used)?;
            let mut result = Vec::with_capacity(owner.face_count());
            for index in 0..owner.face_count() {
                work.step(4096)?;
                let face = owner.face_source(index).expect("retained face table");
                let (role, area_error, position_error, moment_weight) = match face {
                    MaterialFaceSource::ThreadSide { face } => {
                        let (e, p) = side_errors[face];
                        (
                            ScrewFace::ThreadWall(
                                self.proposal.thread().face_index(0, face + 1) - 2,
                            ),
                            e,
                            p,
                            None,
                        )
                    }
                    MaterialFaceSource::ThreadTip => (
                        ScrewFace::ThreadTip,
                        cap_area_error.hi,
                        (I::point(cap_position) + station_gap).hi,
                        None,
                    ),
                    MaterialFaceSource::JoinAnnulus => {
                        let radius = I::point(owner.bottom[0]);
                        let disk = tau * radius.square() / I::point(2.);
                        // Subtraction can greatly reduce the annulus area.
                        // Its moment transfer therefore uses disk+inner cap,
                        // never the small resulting annulus as a weight.
                        (
                            ScrewFace::ThreadShoulder,
                            cap_area_error.hi,
                            (I::point(cap_position) + station_gap).hi,
                            Some((disk * I::point(2.)).hi),
                        )
                    }
                    MaterialFaceSource::TopAnnulus => (
                        ScrewFace::HeadTop,
                        0.,
                        (expression(owner.mouth) - I::point(self.proposal.head_top_z()))
                            .magnitude(),
                        None,
                    ),
                    MaterialFaceSource::ReversedCavity { face } => (
                        ScrewFace::Cavity(face),
                        0.,
                        (expression(owner.mouth) - I::point(self.proposal.head_top_z()))
                            .magnitude(),
                        None,
                    ),
                    MaterialFaceSource::HeadMeridian { segment } => {
                        let p=head.evaluate_with_interrupt(segment,[0.,1.],head.work_used()+257,BYTES,interrupt).map_err(crate::step_read::phillips_cutter::physical::body::closed_head::repaired_error)?;
                        let r = I::point(p.position_mm[0][0].abs().max(p.position_mm[0][1].abs()));
                        let speed = magnitude([
                            interval(p.derivative_mm[0]),
                            interval(p.derivative_mm[1]),
                            I::point(0.),
                        ]);
                        let e = I::point(head_transfer.derivative_mm) * tau * r
                            + (speed + I::point(head_transfer.derivative_mm))
                                * tau
                                * I::point(head_transfer.position_mm);
                        (
                            ScrewFace::ProfileWall(segment as u32),
                            e.hi,
                            head_transfer.position_mm,
                            None,
                        )
                    }
                };
                let old = proposal
                    .iter()
                    .find(|p| p.face == role)
                    .expect("same authored piece ownership");
                let area = widen(interval(old.area_mm2), area_error);
                let weight = moment_weight.unwrap_or(area.hi);
                let moment_error = I::point(position_error) * I::point(weight)
                    + (I::point(extent) + I::point(position_error)) * I::point(area_error);
                let moments = old
                    .first_moment_mm3
                    .map(|m| widen(interval(m), moment_error.hi));
                // Retain the old integral's SAME-area centroid correlation.
                // Cnew-Cold=(deltaM-Cold*deltaA)/Anew. Intersect this
                // enclosure with the independent moment/area quotient.
                let centroid = std::array::from_fn(|axis| {
                    let old_c = interval(old.centroid_mm[axis]);
                    let error =
                        (moment_error + I::point(old_c.magnitude()) * I::point(area_error)) / area;
                    let correlated = widen(old_c, error.hi);
                    let quotient = moments[axis] / area;
                    I::new(
                        correlated.lo.max(quotient.lo),
                        correlated.hi.min(quotient.hi),
                    )
                });
                require(
                    area.valid()
                        && area.lo > 0.
                        && area.width() <= options.area_tolerance_mm2
                        && moments
                            .into_iter()
                            .all(|m| m.valid() && m.width() <= options.first_moment_tolerance_mm3)
                        && centroid
                            .into_iter()
                            .all(|m| m.valid() && m.width() <= options.centroid_tolerance_mm),
                )?;
                result.push(PeriodicFaceProperties {
                    source_face: index,
                    area_mm2: bounds(area),
                    first_moment_mm3: moments.map(bounds),
                    centroid_mm: centroid.map(bounds),
                    area_transfer_error_mm2: area_error,
                    first_moment_transfer_error_mm3: moment_error.hi,
                });
            }
            Ok(result)
        })();
        if work.cancelled.get() {
            Err(PhillipsReadError::Cancelled)
        } else {
            result
        }
    }
}
