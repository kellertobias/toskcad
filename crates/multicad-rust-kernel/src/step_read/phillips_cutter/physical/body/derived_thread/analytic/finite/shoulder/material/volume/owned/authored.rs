//! Direct typed authored authority. No STEP serialization, reader, decoded IDs
//! or old cached-rate thread is used to define this periodic geometry.
use super::*;
use crate::{
    head_profile::{Profile, RepairedHeadOptions},
    helical_thread::Thread,
    phillips_drive::{
        ClosedBoundaryCertificate, ClosedBoundaryOptions, PhillipsScrew,
        RepairedHeadContainmentOptions,
    },
    profiled_screw::SourceSpec,
};
#[cfg(test)]
mod tests;
// All typed inputs, fixed repaired caches,88-domain cavity certificate, fixed
// owned source and publication temporaries together fit this reserve. Each
// numerical child frontier is sequential and is separately bounded below.
const FIXED: usize = 256 * 1024;
const CHILD_READONLY: usize = 64 * 1024;
#[derive(Clone, Copy)]
pub struct AuthoredPeriodicOptions {
    pub minimum_clearance_mm: f64,
    pub maximum_width_mm3: f64,
    pub maximum_containment_cells: usize,
    pub maximum_head_cells: usize,
    pub maximum_cavity_evaluations: usize,
    pub maximum_work: usize,
    pub maximum_bytes: usize,
}
impl Default for AuthoredPeriodicOptions {
    fn default() -> Self {
        Self {
            minimum_clearance_mm: 1e-7,
            maximum_width_mm3: 1e-6,
            maximum_containment_cells: 65536,
            maximum_head_cells: 65536,
            maximum_cavity_evaluations: 1_000_000,
            maximum_work: TOTAL,
            maximum_bytes: BYTES,
        }
    }
}
impl OwnedPeriodicMaterial {
    /// Route eligibility only, never a geometry/material certificate. Other
    /// families keep their existing producer. Once eligible, an actual factory
    /// failure must be reported rather than silently falling back to that path.
    pub fn supports_authored_family(thread: &Thread, profile: &Profile, spec: SourceSpec) -> bool {
        thread.spec().starts == 1
            && thread.spec().root_radius_mm.is_none()
            && matches!(spec.head, Some(crate::head_profile::Head::Crowned { .. }))
            && profile.end()[0] == thread.spec().major_radius_mm
    }
    pub fn from_authored(
        source: &PhillipsScrew,
        o: AuthoredPeriodicOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> FResult<Self> {
        Self::from_authored_parts(
            source.thread(),
            source.profile(),
            source.cavity().source(),
            source.spec(),
            o,
            interrupt,
        )
    }
    pub fn from_authored_parts(
        thread: &Thread,
        profile: &Profile,
        cutter: &Cutter,
        spec: SourceSpec,
        o: AuthoredPeriodicOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> FResult<Self> {
        let mut work = Work {
            remaining: o.maximum_work,
            interrupt,
            cancelled: Cell::new(false),
        };
        let result = (|| {
            work.step(1)?;
            if !o.minimum_clearance_mm.is_finite()
                || o.minimum_clearance_mm <= 0.
                || !o.maximum_width_mm3.is_finite()
                || o.maximum_width_mm3 <= 0.
                || o.maximum_containment_cells == 0
                || o.maximum_containment_cells > 65536
                || o.maximum_head_cells == 0
                || o.maximum_head_cells > 65536
                || o.maximum_cavity_evaluations == 0
                || o.maximum_cavity_evaluations > 1_000_000
                || o.maximum_work > TOTAL
                || o.maximum_bytes > BYTES
                || !spec.minimum_wall_mm.is_finite()
                || spec.minimum_wall_mm <= 0.
                || !spec.shank_length_mm.is_finite()
                || spec.shank_length_mm < 0.
            {
                return Err(PhillipsReadError::Invalid);
            }
            work.step(32768)?;
            if !Self::supports_authored_family(thread, profile, spec) {
                return Err(PhillipsReadError::Unsupported);
            }
            crate::profiled_screw::qualify_outer_source(
                thread,
                profile,
                spec.shank_length_mm,
                true,
            )
            .map_err(|_| PhillipsReadError::Unsupported)?;
            let available = o
                .maximum_bytes
                .checked_sub(FIXED)
                .filter(|n| *n >= ClosedBoundaryCertificate::scratch_budget_bytes())
                .ok_or(PhillipsReadError::Limit)?;
            let child_bytes = available.saturating_add(CHILD_READONLY);
            let head = profile
                .closed_repaired_meridian(
                    RepairedHeadOptions {
                        maximum_transfer_mm: 1e-7,
                        maximum_work: work.remaining.min(WORK),
                        maximum_bytes: child_bytes,
                    },
                    interrupt,
                )
                .map_err(
                    crate::step_read::phillips_cutter::physical::body::closed_head::repaired_error,
                )?;
            work.step(head.work_used())?;
            let local = cutter
                .closed_boundary_certificate(
                    ClosedBoundaryOptions {
                        maximum_work: work.remaining.min(WORK),
                        maximum_bytes: child_bytes,
                    },
                    interrupt,
                )
                .map_err(cavity_error)?;
            work.step(local.work_used())?;
            let periodic = PeriodicThreadSource::new(*thread, &mut work)?;
            let bottom = head.station(head.segment_count()).unwrap();
            let caps = [
                [bottom[1], -spec.shank_length_mm],
                [thread.spec().length_mm, 0.],
            ];
            let lower = expression(caps[0]);
            // The same two-term expression is the translated head bottom and thread
            // entry. This is exact sewing by identity, including nonrepresentable sums;
            // the old rounded visible_start value is never substituted for that seam.
            if bottom[0] != periodic.radius_nodes_mm()[0]
                || periodic.radius_nodes_mm()[1] >= bottom[0]
                || !lower.valid()
                || lower.lo < 0.
                || (I::point(thread.spec().length_mm) - lower).lo <= 0.
            {
                return Err(PhillipsReadError::Unsupported);
            }
            let child_prefix = local.work_used() + head.work_used();
            let containment_limit = fitting_limit(
                o.maximum_containment_cells,
                available,
                RepairedHeadContainment::scratch_budget_bytes,
            );
            let containment = local
                .certify_repaired_head(
                    &head,
                    RepairedHeadContainmentOptions {
                        minimum_clearance_mm: o.minimum_clearance_mm.max(spec.minimum_wall_mm),
                        maximum_cells: containment_limit,
                        maximum_work: work.remaining.saturating_add(child_prefix).min(WORK),
                        maximum_bytes: child_bytes,
                    },
                    interrupt,
                )
                .map_err(cavity_error)?;
            work.step(
                containment
                    .work_used()
                    .checked_sub(child_prefix)
                    .expect("authored containment includes both already-paid source prefixes"),
            )?;
            let volume_options = PeriodicVolumeOptions {
                maximum_width_mm3: o.maximum_width_mm3,
                maximum_head_cells: o.maximum_head_cells,
                maximum_cavity_evaluations: o.maximum_cavity_evaluations,
                maximum_work: WORK,
                maximum_total_work: TOTAL,
                maximum_bytes: child_bytes,
            };
            let mut mass_work = Work {
                remaining: work.remaining.min(WORK),
                interrupt,
                cancelled: Cell::new(false),
            };
            let mass_allowance = mass_work.remaining;
            mass_work.step(META).map_err(|e| {
                if mass_work.cancelled.get() {
                    PhillipsReadError::Cancelled
                } else {
                    e.into()
                }
            })?;
            let measured = measure::measure(
                measure::MeasureInput {
                    head: &head,
                    cavity: cutter,
                    thread: &periodic,
                    stations: caps,
                    determinant: [1., 1.],
                },
                volume_options,
                available,
                &mut mass_work,
            )
            .map_err(|e| {
                if mass_work.cancelled.get() {
                    PhillipsReadError::Cancelled
                } else {
                    e
                }
            })?;
            work.step(mass_allowance - mass_work.remaining)?;
            let mut enclosure = [
                [-bottom[0], bottom[0]],
                [-bottom[0], bottom[0]],
                [lower.lo, thread.spec().length_mm],
            ];
            for i in 0..head.segment_count() {
                work.step(257)?;
                let b = head.certified_interval_bounds(i, [0., 1.]);
                let r = b.position_mm[0][1];
                for axis in &mut enclosure[..2] {
                    axis[0] = axis[0].min(-r);
                    axis[1] = axis[1].max(r);
                }
                let z = interval(b.position_mm[1]) - I::point(spec.shank_length_mm);
                enclosure[2][0] = enclosure[2][0].min(z.lo);
                enclosure[2][1] = enclosure[2][1].max(z.hi);
            }
            let mut cavity_faces = [CavityFace::Mouth; 88];
            let mut count = 0;
            for domain in local.patch_domains() {
                work.step(1)?;
                let face = domain.source_face;
                if face != CavityFace::Mouth && !cavity_faces[..count].contains(&face) {
                    cavity_faces[count] = face;
                    count += 1;
                }
            }
            let owned = OwnedPeriodicMaterial {
                provenance: PeriodicMaterialProvenance::Authored {
                    thread: *thread,
                    profile: *profile,
                    source_spec: spec,
                },
                head: head.owned_snapshot(),
                thread: periodic,
                cavity: cutter.clone(),
                chart: CommonMouthChart {
                    columns: [[0.; 3], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
                    outer_top_z: 0.,
                },
                determinant: [1., 1.],
                caps,
                head_shift: spec.shank_length_mm,
                mouth: [head.station(0).unwrap()[1], -spec.shank_length_mm],
                bottom,
                containment: [
                    containment.minimum_radial_clearance_mm(),
                    containment.axial_floor_clearance_mm(),
                ],
                cavity_faces,
                cavity_count: count,
                local_volume: measured.local,
                components: measured.components,
                enclosure,
                inherited: 0,
                new_work: 0,
            };
            // All retained fields are fixed inline caches in the preflighted reserve.
            assert!(
                owned.retained_budget_bytes() <= FIXED,
                "fixed authored periodic source publication"
            );
            work.step(1)?;
            Ok(OwnedPeriodicMaterial {
                new_work: o.maximum_work - work.remaining,
                ..owned
            })
        })();
        if work.cancelled.get() {
            Err(PhillipsReadError::Cancelled)
        } else {
            result
        }
    }
}
