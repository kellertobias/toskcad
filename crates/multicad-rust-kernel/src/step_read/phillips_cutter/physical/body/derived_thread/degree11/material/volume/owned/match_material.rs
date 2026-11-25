//! Whole-boundary correspondence for unchanged LOCAL repaired D11 sources.
//! Exact binary64 defining coefficients are compared, never decimal tokens or
//! sampled positions. Each operand independently owns its original-file repair,
//! occupied material and volume proofs. This does not certify exact raw shells.
use super::*;
mod world;
const EXTRA: usize = 2 * 1024 * 1024;
#[derive(Clone, Copy)]
pub struct Degree11ReplayMatchOptions {
    pub maximum_distance_mm: f64,
    pub maximum_derivative_mm: f64,
    pub maximum_volume_discrepancy_mm3: f64,
    /// Both independently earned owning prefixes, plus fresh comparison work.
    pub maximum_work: usize,
    /// Both live owners, the receipt and complete atlas/numerical workspace.
    pub maximum_bytes: usize,
}
impl Default for Degree11ReplayMatchOptions {
    fn default() -> Self {
        Self {
            maximum_distance_mm: 1e-7,
            maximum_derivative_mm: 1e-7,
            maximum_volume_discrepancy_mm3: 1e-5,
            maximum_work: WORK_LIMIT,
            maximum_bytes: MAX_BYTES,
        }
    }
}
pub struct Degree11ReplayMatch<'s, 't> {
    source: &'s OwnedDegree11Material,
    target: &'t OwnedDegree11Material,
    distance: f64,
    derivative: f64,
    mass: f64,
    inherited: usize,
    fresh: usize,
    peak: usize,
}
impl Degree11ReplayMatch<'_, '_> {
    pub fn source(&self) -> &OwnedDegree11Material {
        self.source
    }
    pub fn target(&self) -> &OwnedDegree11Material {
        self.target
    }
    pub fn maximum_world_distance_mm(&self) -> f64 {
        self.distance
    }
    /// Complete first derivatives in each ORIGINAL thread root-cell unit
    /// square and each source-atlas normalized patch chart. These are not
    /// derivatives after rescaling a multi-root thread face to one unit square,
    /// and are not unit-normal angle bounds.
    pub fn maximum_world_chart_derivative_mm(&self) -> f64 {
        self.derivative
    }
    pub fn independent_volume_discrepancy_mm3(&self) -> f64 {
        self.mass
    }
    pub fn inherited_work_used(&self) -> usize {
        self.inherited
    }
    pub fn new_work_used(&self) -> usize {
        self.fresh
    }
    pub fn total_work_used(&self) -> usize {
        self.inherited + self.fresh
    }
    pub fn accounted_peak_bytes(&self) -> usize {
        self.peak
    }
}
fn compatible(value: bool) -> std::result::Result<(), PhillipsReadError> {
    if value {
        Ok(())
    } else {
        Err(PhillipsReadError::Unsupported)
    }
}
fn charge(work: &mut Work<'_>, n: usize) -> std::result::Result<(), PhillipsReadError> {
    work.step(n).map_err(|e| bound_error(e, work))
}
fn magnitude(values: impl IntoIterator<Item = I>) -> I {
    values
        .into_iter()
        .fold(I::point(0.), |s, x| s + I::point(x.magnitude()))
}
fn maximum(a: I, b: I) -> I {
    I::point(a.hi.max(b.hi))
}
fn role(a: Degree11MaterialFace, b: Degree11MaterialFace) -> bool {
    use Degree11MaterialFace::*;
    match (a, b) {
        (HeadMeridian { segment: a, .. }, HeadMeridian { segment: b, .. }) => a == b,
        (ThreadSide { index: a, .. }, ThreadSide { index: b, .. }) => a == b,
        (ThreadTip { .. }, ThreadTip { .. })
        | (JoinAnnulus { .. }, JoinAnnulus { .. })
        | (TopAnnulus { .. }, TopAnnulus { .. }) => true,
        (ReversedCavity { face: a, .. }, ReversedCavity { face: b, .. }) => a == b,
        _ => false,
    }
}
impl OwnedDegree11Material {
    /// Narrow retained-replay contract: local source charts must be unchanged.
    /// Entity IDs may be rebased, but complete original rail incidence, actual
    /// repaired P11 coefficients, cap stations, rim uses and face roles must
    /// agree. Arbitrary merely close/edited D11 materials are refused.
    pub fn match_replayed_material<'s, 't>(
        &'s self,
        target: &'t Self,
        occurrence: [[f64; 3]; 4],
        options: Degree11ReplayMatchOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> std::result::Result<Degree11ReplayMatch<'s, 't>, PhillipsReadError> {
        let identity =
            [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]].map(|p| p.map(|v| [v, v]));
        self.match_replayed_material_through(target, occurrence, identity, options, interrupt)
    }
    /// The same contract for a target decoded from a mapped item's LOCAL
    /// records: `target_operator` encloses that item's exact file map
    /// (origin, then the images of the local axes). Local coefficients must
    /// still be unchanged; every world bound transports through the map.
    pub fn match_replayed_material_through<'s, 't>(
        &'s self,
        target: &'t Self,
        occurrence: [[f64; 3]; 4],
        target_operator: [[[f64; 2]; 3]; 4],
        options: Degree11ReplayMatchOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> std::result::Result<Degree11ReplayMatch<'s, 't>, PhillipsReadError> {
        let target_operator = target_operator.map(|p| p.map(|[lo, hi]| I::new(lo, hi)));
        if [options.maximum_distance_mm, options.maximum_derivative_mm]
            .into_iter()
            .any(|x| !x.is_finite() || x <= 0. || x > 1e-7)
            || !options.maximum_volume_discrepancy_mm3.is_finite()
            || options.maximum_volume_discrepancy_mm3 <= 0.
            || options.maximum_volume_discrepancy_mm3 > 1e-5
            || options.maximum_work > WORK_LIMIT
            || options.maximum_bytes > MAX_BYTES
        {
            return Err(PhillipsReadError::Invalid);
        }
        let same = std::ptr::eq(self, target);
        let inherited = self
            .work_used()
            .checked_add(if same { 0 } else { target.work_used() })
            .ok_or(PhillipsReadError::Limit)?;
        let allowance = options
            .maximum_work
            .checked_sub(inherited)
            .ok_or(PhillipsReadError::Limit)?;
        let mut work = Work {
            remaining: allowance,
            interrupt,
            cancelled: Cell::new(false),
        };
        charge(&mut work, 1)?;
        let peak = self
            .retained_budget_bytes()
            .checked_add(if same {
                0
            } else {
                target.retained_budget_bytes()
            })
            .and_then(|x| x.checked_add(EXTRA + std::mem::size_of::<Degree11ReplayMatch<'_, '_>>()))
            .filter(|x| *x <= options.maximum_bytes)
            .ok_or(PhillipsReadError::Limit)?;
        charge(&mut work, 8192)?;
        let world = world::compare_through(self.chart, target.chart, occurrence, target_operator)?;
        compatible(
            self.inventory.len() == target.inventory.len()
                && self.thread.len() == target.thread.len()
                && self.rim.len() == target.rim.len()
                && self.stations == target.stations
                && self.cavity == target.cavity,
        )?;
        for (a, b) in self.inventory.iter().zip(&target.inventory) {
            charge(&mut work, 32)?;
            compatible(role(*a, *b))?;
        }
        for (a, b) in self.rim.iter().zip(&target.rim) {
            charge(&mut work, 32)?;
            compatible(a == b)?;
        }
        let head = self.head.view();
        let other = target.head.view();
        compatible(head.segment_count() == other.segment_count())?;
        for i in 0..=head.segment_count() {
            charge(&mut work, 64)?;
            compatible(head.station(i) == other.station(i))?;
        }
        let h = head
            .translated_correspondence(&other, 0., 0., work.remaining.min(16_000_000), interrupt)
            .map_err(
                crate::step_read::phillips_cutter::physical::body::closed_head::repaired_error,
            )?;
        charge(&mut work, h.work_used)?;
        let mut speed = I::point(0.);
        for (a, b) in self.thread.iter().zip(&target.thread) {
            charge(&mut work, 256)?;
            let at = a.transfer;
            let bt = b.transfer;
            compatible(
                a.forward == b.forward
                    && a.source.spans() == b.source.spans()
                    && at.source_face_index == bt.source_face_index
                    && at.original_v_start_cell == bt.original_v_start_cell
                    && at.original_cells == bt.original_cells
                    && at.source_cells == bt.source_cells,
            )?;
            for (p, q) in a.source.poles().iter().zip(b.source.poles()) {
                charge(&mut work, 32)?;
                compatible(p == q)?;
            }
            let [nu, nv] = a.source.spans();
            for u in 0..nu {
                for v in 0..nv {
                    // Original cell construction+first derivative grids/hulls,
                    // all copies/validation included; <16384 action units.
                    charge(&mut work, 16384)?;
                    let cell = a
                        .source
                        .cell([u, v])
                        .map_err(|_| PhillipsReadError::Unsupported)?;
                    let d = cell
                        .first_derivative_bounds()
                        .map_err(|_| PhillipsReadError::Unsupported)?;
                    for x in [d.du, d.dv] {
                        speed = maximum(speed, magnitude(x.map(interval)));
                    }
                }
            }
        }
        // Preserve the complete equality pattern of original rail owners,
        // not coincident endpoints or a guessed four-band convention.
        for a in 0..self.thread.len() {
            for b in 0..self.thread.len() {
                for i in 0..2 {
                    for j in 0..2 {
                        charge(&mut work, 32)?;
                        compatible(
                            (self.thread[a].transfer.original_rail_entities[i]
                                == self.thread[b].transfer.original_rail_entities[j])
                                == (target.thread[a].transfer.original_rail_entities[i]
                                    == target.thread[b].transfer.original_rail_entities[j]),
                        )?;
                    }
                }
            }
        }
        // Same validated Cutter spec+domain define identical analytic source
        // charts, including rounded walls and trims (Cutter::PartialEq).
        // Build the ACTUAL atlas to earn complete source derivative speeds;
        // no normal/mesh or raw-decimal equality supplies those bounds.
        charge(&mut work, 1_000_000)?;
        let atlas = self
            .cavity
            .surface_patches_with_cancel(&|| interrupt.is_some_and(|f| f()))
            .map_err(cavity_error)?;
        if atlas.retained_budget_bytes() > 1024 * 1024 {
            return Err(PhillipsReadError::Limit);
        }
        let mut cavity_speed = I::point(0.);
        let mut mouth_radius = f64::INFINITY;
        let mut mouth_count = 0usize;
        for i in 0..atlas.patches().len() {
            charge(&mut work, 65536)?;
            let d = atlas
                .patch_bounds(i, [[0., 1.], [0., 1.]])
                .map_err(cavity_error)?;
            for x in [d.du_mm, d.dv_mm] {
                cavity_speed = maximum(cavity_speed, magnitude(x.map(interval)));
            }
            if matches!(
                atlas.patches()[i].source_identity(),
                crate::phillips_drive::PatchKey::MouthTrace { .. }
                    | crate::phillips_drive::PatchKey::MouthGap { .. }
            ) {
                // V=1 is the COMPLETE actual mouth curve of these source
                // charts. The atlas radial-margin reserve is an intercept,
                // not assumed to be a lower bound on this geometric norm.
                charge(&mut work, 65536)?;
                let mouth = atlas
                    .patch_bounds(i, [[0., 1.], [1., 1.]])
                    .map_err(cavity_error)?;
                let r = (interval(mouth.point_mm[0]).square()
                    + interval(mouth.point_mm[1]).square())
                .sqrt();
                compatible(r.valid() && r.lo > 0. && r.hi.is_finite())?;
                mouth_radius = mouth_radius.min(r.lo);
                mouth_count += 1;
            }
        }
        charge(&mut work, 16384)?;
        let rho = I::point(self.core.min(mouth_radius));
        compatible(mouth_count > 0 && mouth_radius.is_finite())?;
        compatible(rho.valid() && rho.lo > 0.)?;
        let mut radius = I::point(0.);
        let mut meridian_speed = I::point(0.);
        for i in 0..head.segment_count() {
            charge(&mut work, 4096)?;
            let q = head.certified_interval_bounds(i, [0., 1.]);
            radius = maximum(radius, I::point(q.position_mm[0][1]));
            meridian_speed = maximum(meridian_speed, magnitude(q.derivative_mm.map(interval)));
        }
        // 2pi<7. Exact source circle/trapezoid charts: full-rim derivatives
        // normalized per original U cell; ray derivative norm<=speed/rho.
        // This covers tip disks, actual join annuli and the shaped mouth/top.
        let extent = magnitude(self.enclosure.map(interval));
        let source_speed = [
            meridian_speed,
            I::point(7.) * radius,
            speed,
            cavity_speed,
            speed + radius * speed / rho,
            cavity_speed + radius * cavity_speed / rho,
            I::point(2.) * extent,
        ]
        .into_iter()
        .fold(I::point(0.), maximum);
        let local_distance = I::point(h.position_mm);
        let local_derivative = maximum(I::point(h.derivative_mm), I::point(7.) * local_distance);
        let distance =
            world.translation + world.linear * extent + world.target_norm * local_distance;
        let derivative = world.linear * source_speed + world.target_norm * local_derivative;
        let mass = I::point(
            (interval(self.world_volume_mm3()) * world.occurrence_determinant
                - interval(target.world_volume_mm3()) * world.target_determinant)
                .magnitude(),
        );
        compatible(
            [distance, derivative, mass]
                .into_iter()
                .all(|x| x.valid() && x.hi.is_finite())
                && distance.hi <= options.maximum_distance_mm
                && derivative.hi <= options.maximum_derivative_mm
                && mass.hi <= options.maximum_volume_discrepancy_mm3,
        )?;
        charge(&mut work, 1)?;
        Ok(Degree11ReplayMatch {
            source: self,
            target,
            distance: distance.hi,
            derivative: derivative.hi,
            mass: mass.hi,
            inherited,
            fresh: allowance - work.remaining,
            peak,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::step_read::phillips_cutter::physical::body::{
        derived_thread::actual_m3, read_degree11_phillips_material,
    };
    const IDENTITY: [[f64; 3]; 4] = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
    fn read(input: &str) -> OwnedDegree11Material {
        read_degree11_phillips_material(input, Default::default(), None)
            .unwrap()
            .into_material()
    }
    #[test]
    fn replay_match_compares_independently_decoded_actual_rounded_material_and_mass() {
        for left in [false, true] {
            let input = actual_m3(true, left);
            let a = read(&input);
            let b = read(&input);
            assert_ne!(
                a.thread[0].source.poles().as_ptr(),
                b.thread[0].source.poles().as_ptr()
            );
            assert!(a.original_input() == b.original_input());
            let proof = a
                .match_replayed_material(&b, IDENTITY, Default::default(), None)
                .unwrap();
            eprintln!(
                "D11_REPLAY_MATCH left={left} work={} fresh={} bytes={} c0={} d1={} mass={}",
                proof.total_work_used(),
                proof.new_work_used(),
                proof.accounted_peak_bytes(),
                proof.maximum_world_distance_mm(),
                proof.maximum_world_chart_derivative_mm(),
                proof.independent_volume_discrepancy_mm3()
            );
            assert!(proof.maximum_world_distance_mm() <= 1e-7);
            assert!(proof.maximum_world_chart_derivative_mm() <= 1e-7);
            assert!(proof.independent_volume_discrepancy_mm3() <= 1e-5);
            assert_eq!(proof.inherited_work_used(), a.work_used() + b.work_used());
            assert_eq!(
                proof.total_work_used(),
                proof.inherited_work_used() + proof.new_work_used()
            );
        }
    }
    #[test]
    fn replay_match_refuses_actual_other_sources_frames_limits_and_cancellation_atomically() {
        let input = actual_m3(true, false);
        let a = read(&input);
        let b = read(&input);
        let options = Degree11ReplayMatchOptions::default();
        let calls = Cell::new(0usize);
        let count = || {
            calls.set(calls.get() + 1);
            false
        };
        let earned = a
            .match_replayed_material(&b, IDENTITY, options, Some(&count))
            .unwrap();
        let total = earned.total_work_used();
        let bytes = earned.accounted_peak_bytes();
        let distance = earned.maximum_world_distance_mm();
        for bad in [
            Degree11ReplayMatchOptions {
                maximum_work: total - 1,
                ..options
            },
            Degree11ReplayMatchOptions {
                maximum_bytes: bytes - 1,
                ..options
            },
        ] {
            assert!(matches!(
                a.match_replayed_material(&b, IDENTITY, bad, None),
                Err(PhillipsReadError::Limit)
            ));
        }
        let exact = a
            .match_replayed_material(
                &b,
                IDENTITY,
                Degree11ReplayMatchOptions {
                    maximum_work: total,
                    maximum_bytes: bytes,
                    ..options
                },
                None,
            )
            .unwrap();
        assert_eq!(exact.total_work_used(), total);
        assert_eq!(
            exact.maximum_world_distance_mm().to_bits(),
            distance.to_bits()
        );
        assert!(matches!(
            a.match_replayed_material(
                &b,
                IDENTITY,
                Degree11ReplayMatchOptions {
                    maximum_distance_mm: f64::NAN,
                    ..options
                },
                None
            ),
            Err(PhillipsReadError::Invalid)
        ));
        assert!(matches!(
            a.match_replayed_material(
                &b,
                IDENTITY,
                Degree11ReplayMatchOptions {
                    maximum_distance_mm: distance * 0.5,
                    ..options
                },
                None
            ),
            Err(PhillipsReadError::Unsupported)
        ));
        let mut wrong = IDENTITY;
        wrong[0][0] = 1.;
        assert!(matches!(
            a.match_replayed_material(&b, wrong, options, None),
            Err(PhillipsReadError::Unsupported)
        ));
        let mut invalid = IDENTITY;
        invalid[1][0] = 2.;
        assert!(matches!(
            a.match_replayed_material(&b, invalid, options, None),
            Err(PhillipsReadError::Invalid)
        ));
        let other = read(&actual_m3(false, false));
        assert!(matches!(
            a.match_replayed_material(&other, IDENTITY, options, None),
            Err(PhillipsReadError::Unsupported)
        ));
        let other = read(&actual_m3(true, true));
        assert!(matches!(
            a.match_replayed_material(&other, IDENTITY, options, None),
            Err(PhillipsReadError::Unsupported)
        ));
        assert!(matches!(
            a.match_replayed_material(&b, IDENTITY, options, Some(&|| true)),
            Err(PhillipsReadError::Cancelled)
        ));
        let late_calls = Cell::new(0usize);
        let cancel = || {
            late_calls.set(late_calls.get() + 1);
            late_calls.get() >= calls.get()
        };
        assert!(matches!(
            a.match_replayed_material(&b, IDENTITY, options, Some(&cancel)),
            Err(PhillipsReadError::Cancelled)
        ));
        let retry = a
            .match_replayed_material(&b, IDENTITY, options, None)
            .unwrap();
        assert_eq!(retry.total_work_used(), total);
        assert_eq!(
            retry.maximum_world_distance_mm().to_bits(),
            distance.to_bits()
        );
        assert_eq!(a.original_input(), input);
    }
}
