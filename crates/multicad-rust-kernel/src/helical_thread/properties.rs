//! Directed integrals of the retained phase coefficients and visible interval.
//! Full-angle domain and finite radial sewing allowances remain explicit;
//! these measures do not independently admit a closed topology or mesh.
use super::*;
use crate::tapered_thread::{Scalar, Series, interval::I};
type S = Series<5>;
type Values = [I; 4];
#[derive(Clone, Copy, Debug)]
pub struct FacePropertyOptions {
    pub area_tolerance_mm2: f64,
    pub first_moment_tolerance_mm3: f64,
    pub centroid_tolerance_mm: f64,
    /// Cumulative attempted cell proofs over the complete query.
    pub maximum_cells: usize,
}
impl Default for FacePropertyOptions {
    fn default() -> Self {
        Self {
            area_tolerance_mm2: 1e-7,
            first_moment_tolerance_mm3: 1e-7,
            centroid_tolerance_mm: 1e-7,
            maximum_cells: 65536,
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FacePropertyError {
    Invalid,
    Unresolved,
    Limit,
    Cancelled,
}
#[derive(Clone, Copy, Debug)]
pub struct FaceProperties {
    /// Original Thread face ordinal: caps0/1, walls as face_index().
    pub face: u32,
    pub area_mm2: [f64; 2],
    pub first_moment_mm3: [[f64; 2]; 3],
    pub centroid_mm: [[f64; 2]; 3],
    /// Already included in area_mm2; includes literal radial sewing reserves.
    pub phase_domain_area_error_mm2: f64,
    /// Already included in first_moment_mm3, including indexed angular drift.
    pub phase_domain_first_moment_error_mm3: [f64; 3],
}
fn n(x: f64) -> S {
    S::constant(x)
}
fn c(x: I) -> S {
    S::interval(x)
}
fn bounds(x: I) -> [f64; 2] {
    [x.lo, x.hi]
}
fn cancelled(f: &dyn Fn() -> bool) -> Result<(), FacePropertyError> {
    if f() {
        Err(FacePropertyError::Cancelled)
    } else {
        Ok(())
    }
}
struct PhasePrecision {
    radius: I,
    angular: I,
    cap_area: I,
    wall_area: I,
    axial: [f64; 2],
    volume: I,
}
impl PhasePrecision {
    fn new(thread: &Thread, axial: [f64; 2]) -> Result<Self, FacePropertyError> {
        let radius = I::point(thread.source_coordinate_bounds_mm()[0][1]);
        let rate = I::point(thread.helical_rate().abs());
        let lead = I::point(f64::from(thread.spec.starts)) * I::point(thread.spec.pitch_mm);
        let tau = I::point(2.) * I::new(PI.next_down(), PI.next_up());
        let delta = I::point((rate * tau - lead).magnitude());
        let mut indexed = 0_f64;
        // Include global clipped stations and the negative/neighbor periods
        // used by finite support and tessellation queries, not only start IDs.
        let margin = 2 * i64::from(thread.spec.starts) + 4;
        let last = (thread.spec.length_mm / thread.spec.pitch_mm).ceil() as i64 + margin;
        for period in -margin..=last {
            for phase in thread.breaks {
                let exact =
                    I::point(period as f64) * I::point(thread.spec.pitch_mm) + I::point(phase);
                for literal in [
                    (period as f64).mul_add(thread.spec.pitch_mm, phase),
                    (period as f64) * thread.spec.pitch_mm + phase,
                ] {
                    indexed = indexed.max((I::point(literal) - exact).magnitude());
                }
            }
        }
        let mut slope = 0_f64;
        let mut jumps = I::point(0.);
        for band in 0..5 {
            let v = profile_jet(
                thread,
                band,
                S::variable(I::new(thread.breaks[band], thread.breaks[band + 1])),
            );
            if !v[0].c[0].valid() || !v[1].c[0].valid() {
                return Err(FacePropertyError::Unresolved);
            }
            slope = slope.max(v[1].c[0].magnitude());
            if band > 0 {
                let phase = n(thread.breaks[band]);
                let left = profile_jet(thread, band - 1, phase)[0].c[0];
                let right = profile_jet(thread, band, phase)[0].c[0];
                jumps = jumps + I::point((left - right).magnitude());
            }
        }
        let m = I::point(slope);
        let height = I::point(axial[1]) - I::point(axial[0]);
        let g = (radius.square() + m.square() * (radius.square() + rate.square())).sqrt();
        let cap_density = radius.square() / (I::point(2.) * rate);
        // True full-angle domain versus one retained lead; additionally retain
        // all rounded indexed boundary shifts rather than treating them exact.
        let cap_area = (delta + I::point(10. * f64::from(thread.spec.starts)) * I::point(indexed))
            * cap_density;
        let wall_strip = (delta + I::point(4.) * I::point(indexed)) * height * g / rate;
        // Any finite radial sewing pieces have positive area. Bound both the
        // angular seam and literal profile-rail residuals, without dropping
        // them or asserting a zero-area discontinuous closure. Allocating this
        // whole reserve to each adjacent semantic face is conservative.
        let angular_seam = (m * delta + jumps) * height;
        let rail_seams = I::point(f64::from(thread.spec.starts))
            * jumps
            * height
            * (I::point(1.) + (radius / rate).square()).sqrt();
        let wall_area = wall_strip + angular_seam + rail_seams;
        let turns = (I::point(axial[0].abs().max(axial[1].abs())) / lead)
            .hi
            .ceil()
            + 2.;
        let angular = (I::point(turns) * delta / rate + I::point(2. * indexed) / rate)
            .hi
            .min(2.);
        let volume = I::point(2.) * delta * cap_density * height;
        Self {
            radius,
            angular: I::point(angular),
            cap_area,
            wall_area,
            axial,
            volume,
        }
        .checked()
    }
    fn checked(self) -> Result<Self, FacePropertyError> {
        if [
            self.radius,
            self.angular,
            self.cap_area,
            self.wall_area,
            self.volume,
        ]
        .iter()
        .any(|v| !v.valid())
        {
            Err(FacePropertyError::Unresolved)
        } else {
            Ok(self)
        }
    }
    fn errors(&self, face: u32, area: f64) -> [f64; 4] {
        let error = if face < 2 {
            self.cap_area
        } else {
            self.wall_area
        };
        let xy = (self.radius * (error + I::point(area) * self.angular)).hi;
        [
            error.hi,
            xy,
            xy,
            (error * I::point(self.axial[0].abs().max(self.axial[1].abs()))).hi,
        ]
    }
}
fn profile_jet(thread: &Thread, band: usize, phase: S) -> [S; 2] {
    let (r, m) = match band {
        1 => (
            n(thread.spec.major_radius_mm) + (phase - n(thread.breaks[1])) * n(thread.flanks[0]),
            n(thread.flanks[0]),
        ),
        3 => (
            n(thread.spec.major_radius_mm) + (phase - n(thread.breaks[4])) * n(thread.flanks[1]),
            n(thread.flanks[1]),
        ),
        2 => {
            if let Some(root) = thread.root {
                let x = phase - n(root.phase);
                let h = (n(root.radius) * n(root.radius) - x * x).sqrt();
                (n(root.center_radius) - h, x / h)
            } else {
                (n(thread.root_radius_mm()), n(0.))
            }
        }
        _ => (n(thread.spec.major_radius_mm), n(0.)),
    };
    [r, m]
}

struct Part<'a> {
    thread: &'a Thread,
    period: u32,
    band: usize,
    axial: [f64; 2],
    cap: Option<usize>,
}
impl Part<'_> {
    fn evaluate(&self, t: I) -> Result<[S; 4], FacePropertyError> {
        let a = self.thread.breaks[self.band];
        let span = I::point(self.thread.breaks[self.band + 1]) - I::point(a);
        let phase = n(a) + S::variable(t) * c(span);
        let [r, m] = profile_jet(self.thread, self.band, phase);
        let rate = I::point(self.thread.helical_rate());
        let q = phase + c(I::point(f64::from(self.period)) * I::point(self.thread.spec.pitch_mm));
        let result = if let Some(cap) = self.cap {
            let theta = (n(self.axial[cap]) - q) / c(rate);
            let [sin, cos] = theta.sin_cos();
            let area = r * r * c(span / (I::point(2.) * I::point(rate.magnitude())));
            let moment = r * r * r * c(span / (I::point(3.) * I::point(rate.magnitude())));
            [area, moment * cos, moment * sin, area * n(self.axial[cap])]
        } else {
            let g = (r * r + m * m * (r * r + c(rate.square()))).sqrt();
            let height = I::point(self.axial[1]) - I::point(self.axial[0]);
            let area = g * c(span * height / I::point(rate.magnitude()));
            let [sa, ca] = ((n(self.axial[0]) - q) / c(rate)).sin_cos();
            let [sb, cb] = ((n(self.axial[1]) - q) / c(rate)).sin_cos();
            let moment = g * r * c(span * rate / I::point(rate.magnitude()));
            [
                area,
                moment * (sb - sa),
                moment * (ca - cb),
                area * c((I::point(self.axial[0]) + I::point(self.axial[1])) / I::point(2.)),
            ]
        };
        for v in result {
            v.checked().map_err(|_| FacePropertyError::Unresolved)?;
        }
        Ok(result)
    }
}
fn cell(part: &Part<'_>, a: f64, b: f64) -> Result<Values, FacePropertyError> {
    let w = I::point(b) - I::point(a);
    let mid = (I::point(a) + I::point(b)) / I::point(2.);
    let offset = w / I::point(12.).sqrt();
    let left = part.evaluate(mid - offset)?;
    let right = part.evaluate(mid + offset)?;
    let whole = part.evaluate(I::new(a, b))?;
    let mut result = [I::point(0.); 4];
    for i in 0..4 {
        let e = (I::point(whole[i].c[4].magnitude()) * w.square().square() * w / I::point(180.)).hi;
        result[i] = ((left[i].c[0] + right[i].c[0]) * w / I::point(2.) + I::new(-e, e))
            .checked()
            .map_err(|_| FacePropertyError::Unresolved)?;
    }
    Ok(result)
}
fn receipt(face: u32, v: Values, o: FacePropertyOptions) -> Option<FaceProperties> {
    if v.iter().any(|v| !v.valid())
        || v[0].lo <= 0.
        || v[0].width() > o.area_tolerance_mm2
        || v[1..]
            .iter()
            .any(|v| v.width() > o.first_moment_tolerance_mm3)
    {
        return None;
    }
    let centroid = [v[1], v[2], v[3]].map(|m| m / v[0]);
    if centroid
        .iter()
        .any(|v| !v.valid() || v.width() > o.centroid_tolerance_mm)
    {
        return None;
    }
    Some(FaceProperties {
        face,
        area_mm2: bounds(v[0]),
        first_moment_mm3: [v[1], v[2], v[3]].map(bounds),
        centroid_mm: centroid.map(bounds),
        phase_domain_area_error_mm2: 0.,
        phase_domain_first_moment_error_mm3: [0.; 3],
    })
}
impl Thread {
    /// Difference between a true full-angle material volume and the retained
    /// one-lead phase mean, without changing the original floating rate.
    pub fn phase_domain_volume_error_mm3(&self, axial: [f64; 2]) -> Result<f64, FacePropertyError> {
        if axial.iter().any(|v| !v.is_finite())
            || axial[0] < 0.
            || axial[0] >= axial[1]
            || axial[1] > self.spec.length_mm
        {
            return Err(FacePropertyError::Invalid);
        }
        Ok(PhasePrecision::new(self, axial)?.volume.hi)
    }
    /// Conservative actual-support box, including the retained root circle.
    /// The complete circle overbounds its finite phase trim without fitting it.
    pub fn source_coordinate_bounds_mm(&self) -> [[f64; 2]; 3] {
        let radial = self.root.map_or(I::point(self.spec.major_radius_mm), |r| {
            I::point(r.center_radius.abs()) + I::point(r.radius)
        });
        let radius = radial.hi.max(self.spec.major_radius_mm);
        [
            [-radius, radius],
            [-radius, radius],
            [0., self.spec.length_mm],
        ]
    }
    /// Integrate the original finite helical faces over this explicit visible
    /// axial interval. Neither the thread length nor its phase is changed.
    pub fn face_properties_with_cancel(
        &self,
        axial: [f64; 2],
        options: FacePropertyOptions,
        interrupt: &dyn Fn() -> bool,
    ) -> Result<Vec<FaceProperties>, FacePropertyError> {
        cancelled(interrupt)?;
        if axial.iter().any(|v| !v.is_finite())
            || axial[0] < 0.
            || axial[0] >= axial[1]
            || axial[1] > self.spec.length_mm
            || options.maximum_cells == 0
            || options.maximum_cells > 65536
            || [
                options.area_tolerance_mm2,
                options.first_moment_tolerance_mm3,
                options.centroid_tolerance_mm,
            ]
            .into_iter()
            .any(|v| !v.is_finite() || v <= 0.)
        {
            return Err(FacePropertyError::Invalid);
        }
        let precision = PhasePrecision::new(self, axial)?;
        let mut remaining = options.maximum_cells;
        let mut result = Vec::with_capacity(2 + 4 * self.spec.starts as usize);
        let extent = self.spec.major_radius_mm.max(self.spec.length_mm).max(1.);
        for face in 0..2 + 4 * self.spec.starts {
            let mut parts = Vec::new();
            for period in 0..self.spec.starts {
                for band in 0..5 {
                    if face < 2 || self.face_index(i64::from(period), band) == face {
                        parts.push(Part {
                            thread: self,
                            period,
                            band,
                            axial,
                            cap: (face < 2).then_some(face as usize),
                        });
                    }
                }
            }
            let count = parts.len() as f64;
            let mut values = [I::point(0.); 4];
            // Relative allocation to exact dyadic parameter width bounds all
            // accepted leaf sums without subtracting historical intervals. For
            // each part, disjoint dyadic leaf widths sum to one, so all part
            // errors sum to at most budget/16. The remaining 15/16 is reserved
            // for interval summation and phase/sewing allowances. Final receipt
            // checks still enforce the complete requested widths. The former
            // budget/32 allocation could fall below evaluation roundoff even
            // for a valid canonical M6 root; subdivision cannot reduce that
            // relative floor.
            let area_budget = options
                .area_tolerance_mm2
                .min(options.centroid_tolerance_mm / (extent + 1.));
            let moment_budget = options
                .first_moment_tolerance_mm3
                .min(options.centroid_tolerance_mm);
            for part in parts {
                let mut pending = vec![(0., 1., 0_u8)];
                while let Some((a, b, depth)) = pending.pop() {
                    cancelled(interrupt)?;
                    remaining = remaining.checked_sub(1).ok_or(FacePropertyError::Limit)?;
                    let proof = cell(&part, a, b);
                    let budgets =
                        [area_budget, moment_budget, moment_budget, moment_budget].map(|budget| {
                            (I::point(budget) * I::point(b - a) / I::point(16. * count)).lo
                        });
                    if let Ok(value) = proof
                        && value.iter().zip(budgets).all(|(v, b)| v.width() <= b)
                    {
                        for i in 0..4 {
                            values[i] = values[i] + value[i];
                        }
                        continue;
                    }
                    if depth >= 48 {
                        return Err(FacePropertyError::Unresolved);
                    }
                    let m = a.midpoint(b);
                    pending.push((m, b, depth + 1));
                    pending.push((a, m, depth + 1));
                }
            }
            let errors = precision.errors(face, values[0].hi);
            for i in 0..4 {
                values[i] = values[i] + I::new(-errors[i], errors[i]);
            }
            let mut measured =
                receipt(face, values, options).ok_or(FacePropertyError::Unresolved)?;
            measured.phase_domain_area_error_mm2 = errors[0];
            measured.phase_domain_first_moment_error_mm3 = [errors[1], errors[2], errors[3]];
            result.push(measured);
        }
        cancelled(interrupt)?;
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn numerical_receipts_refuse_area_and_centroid_widths_independently() {
        let t = tool(false, false, 1);
        let mut precision = PhasePrecision::new(&t, [0., 4.]).unwrap();
        precision.wall_area = I::new(0., f64::INFINITY);
        assert!(precision.checked().is_err());
        let options = FacePropertyOptions::default();
        let mut v = [I::point(1.), I::point(0.), I::point(0.), I::point(0.)];
        v[0] = I::new(1., 2.);
        assert!(receipt(0, v, options).is_none());
        v[0] = I::point(1e-10);
        v[1] = I::new(0., 1e-8);
        assert!(receipt(0, v, options).is_none());
        v[0] = I::new(f64::NAN, 1.);
        assert!(receipt(0, v, options).is_none());
    }
    fn tool(round: bool, left: bool, starts: u32) -> Thread {
        let mut spec = Spec {
            major_radius_mm: 1.47,
            pitch_mm: 0.5,
            length_mm: 4.,
            depth_mm: 0.3,
            crest_width_mm: 0.0625,
            lower_flank_deg: 30.,
            upper_flank_deg: 30.,
            starts,
            left_hand: left,
            root_radius_mm: None,
        };
        if round {
            let radius = 0.05;
            let a = 30_f64.to_radians();
            spec.root_radius_mm = Some(radius);
            spec.depth_mm = (spec.pitch_mm - spec.crest_width_mm) / (2. * a.tan())
                - radius * (1. / a.sin() - 1.);
        }
        Thread::new(spec).unwrap()
    }
    #[test]
    fn retained_flat_and_round_faces_keep_phase_hands_starts_and_visible_interval() {
        for round in [false, true] {
            for starts in [1, 2] {
                let a = tool(round, false, starts);
                let domain = a.source_coordinate_bounds_mm();
                assert!(domain[0][1] >= a.spec.major_radius_mm && domain[2] == [0., 4.]);
                let b = tool(round, true, starts);
                let axial = [0.123, 3.75];
                let options = FacePropertyOptions::default();
                let right = a
                    .face_properties_with_cancel(axial, options, &|| false)
                    .unwrap();
                let left = b
                    .face_properties_with_cancel(axial, options, &|| false)
                    .unwrap();
                assert_eq!(right.len(), 2 + 4 * starts as usize);
                for f in &right {
                    assert!(f.area_mm2[0] > 0.);
                    assert!(f.area_mm2[1] - f.area_mm2[0] <= options.area_tolerance_mm2);
                    for v in f.centroid_mm {
                        assert!(v[1] - v[0] <= options.centroid_tolerance_mm);
                    }
                }
                for (r, l) in right.iter().zip(left) {
                    assert!(r.area_mm2[0] <= l.area_mm2[1] && l.area_mm2[0] <= r.area_mm2[1]);
                    assert!(
                        r.centroid_mm[0][0] <= l.centroid_mm[0][1]
                            && l.centroid_mm[0][0] <= r.centroid_mm[0][1]
                    );
                    assert!(
                        r.centroid_mm[1][0] <= -l.centroid_mm[1][0]
                            && -l.centroid_mm[1][1] <= r.centroid_mm[1][1]
                    );
                }
                let area = right[0].area_mm2;
                assert!((a.volume_mm3() / a.spec.length_mm - area[0]).abs() < 1e-7);
                let crest = right
                    .iter()
                    .filter(|f| f.face >= 2 && (f.face - 2) % 4 == 1)
                    .map(|f| I::new(f.area_mm2[0], f.area_mm2[1]))
                    .fold(I::point(0.), |a, b| a + b);
                let expected = a.spec.major_radius_mm * (axial[1] - axial[0])
                    / a.helical_rate().abs()
                    * a.spec.crest_width_mm
                    * f64::from(starts);
                assert!(crest.lo <= expected && expected <= crest.hi);
                assert_eq!(a.spec(), tool(round, false, starts).spec());
            }
        }
    }
    #[test]
    fn actual_fitted_m6_original_root_meets_native_property_accuracy() {
        for left_hand in [false, true] {
            let thread = Thread::new(Spec {
                major_radius_mm: 2.942,
                pitch_mm: 1.,
                length_mm: 12.,
                depth_mm: 0.5966751345948134,
                crest_width_mm: 0.1443522341791902,
                lower_flank_deg: 30.,
                upper_flank_deg: 30.,
                starts: 1,
                left_hand,
                root_radius_mm: Some(0.14433756729740646),
            })
            .unwrap();
            let before = thread;
            // Actual canonical Core M6x1 6g source; these are the SAME stricter
            // per-thread allocations used by the native screw-property query.
            let options = FacePropertyOptions {
                area_tolerance_mm2: 6.009615384615385e-11,
                first_moment_tolerance_mm3: 7.8125e-10,
                centroid_tolerance_mm: 7.8125e-10,
                maximum_cells: 32768,
            };
            let calls = std::cell::Cell::new(0usize);
            let result = thread
                .face_properties_with_cancel([0., 12.], options, &|| {
                    calls.set(calls.get() + 1);
                    false
                })
                .unwrap();
            let actual_cells = calls.get() - 2; // Initial/final cancellation checkpoints.
            assert!(actual_cells <= options.maximum_cells);
            assert!(matches!(
                thread.face_properties_with_cancel(
                    [0., 12.],
                    FacePropertyOptions {
                        maximum_cells: actual_cells - 1,
                        ..options
                    },
                    &|| false
                ),
                Err(FacePropertyError::Limit)
            ));
            assert!(
                thread
                    .face_properties_with_cancel(
                        [0., 12.],
                        FacePropertyOptions {
                            maximum_cells: actual_cells,
                            ..options
                        },
                        &|| false
                    )
                    .is_ok()
            );
            assert!(matches!(
                thread.face_properties_with_cancel(
                    [0., 12.],
                    FacePropertyOptions {
                        area_tolerance_mm2: 1e-30,
                        ..options
                    },
                    &|| false
                ),
                Err(FacePropertyError::Unresolved | FacePropertyError::Limit)
            ));
            let stop = std::cell::Cell::new(0usize);
            assert!(matches!(
                thread.face_properties_with_cancel([0., 12.], options, &|| {
                    stop.set(stop.get() + 1);
                    stop.get() == actual_cells / 2
                }),
                Err(FacePropertyError::Cancelled)
            ));
            assert!(
                thread
                    .face_properties_with_cancel([0., 12.], options, &|| false)
                    .is_ok()
            );
            assert_eq!(result.len(), 6);
            for face in result {
                assert!(face.area_mm2[1] - face.area_mm2[0] <= options.area_tolerance_mm2);
                for moment in face.first_moment_mm3 {
                    assert!(moment[1] - moment[0] <= options.first_moment_tolerance_mm3);
                }
                for center in face.centroid_mm {
                    assert!(center[1] - center[0] <= options.centroid_tolerance_mm);
                }
            }
            assert_eq!(thread, before);
        }
    }
    #[test]
    fn actual_tiny_retained_root_can_refuse_the_directed_precision_floor_without_mutation() {
        let r = 2_f64.powi(-52).next_up();
        let a = 1_f64.to_radians();
        let depth = (0.5 - 0.0625) / (2. * a.tan()) - r * (1. / a.sin() - 1.);
        let t = Thread::new(Spec {
            major_radius_mm: 20.,
            pitch_mm: 0.5,
            length_mm: 4.,
            depth_mm: depth,
            crest_width_mm: 0.0625,
            lower_flank_deg: 1.,
            upper_flank_deg: 1.,
            starts: 1,
            left_hand: false,
            root_radius_mm: Some(r),
        })
        .unwrap();
        let before = t;
        assert!(matches!(
            t.phase_domain_volume_error_mm3([0., 4.]),
            Err(FacePropertyError::Unresolved)
        ));
        assert!(matches!(
            t.face_properties_with_cancel([0., 4.], FacePropertyOptions::default(), &|| false),
            Err(FacePropertyError::Unresolved)
        ));
        assert_eq!(t, before);
        assert!(matches!(
            t.phase_domain_volume_error_mm3([4., 4.]),
            Err(FacePropertyError::Invalid)
        ));
    }
    #[test]
    fn full_angle_flat_cap_oracle_and_named_original_rate_reserves_agree() {
        for left in [false, true] {
            for starts in [1, 2] {
                let t = tool(false, left, starts);
                let axial = [0., 4.];
                let measured = t
                    .face_properties_with_cancel(axial, FacePropertyOptions::default(), &|| false)
                    .unwrap();
                let mut integral = I::point(0.);
                for band in 0..5 {
                    let a = t.breaks[band];
                    let b = t.breaks[band + 1];
                    let ra = profile_jet(&t, band, n(a))[0].c[0];
                    let rb = profile_jet(&t, band, n(b))[0].c[0];
                    integral = integral
                        + (I::point(b) - I::point(a)) * (ra.square() + ra * rb + rb.square())
                            / I::point(3.);
                }
                let rate = I::point(t.helical_rate().abs());
                let lead = I::point(t.spec.pitch_mm) * I::point(f64::from(starts));
                let tau = I::point(2.) * I::new(PI.next_down(), PI.next_up());
                let end_strip = (rate * tau - lead) * I::point(t.spec.major_radius_mm).square()
                    / (I::point(2.) * rate);
                let oracle =
                    I::point(f64::from(starts)) * integral / (I::point(2.) * rate) + end_strip;
                assert!(
                    measured[0].area_mm2[0] <= oracle.lo && oracle.hi <= measured[0].area_mm2[1]
                );
                assert!(measured[0].phase_domain_area_error_mm2 >= end_strip.magnitude());
                assert!(measured.iter().all(|f| {
                    f.phase_domain_area_error_mm2 > 0.
                        && f.phase_domain_first_moment_error_mm3
                            .iter()
                            .all(|v| v.is_finite() && *v >= 0.)
                }));
                assert!(t.phase_domain_volume_error_mm3(axial).unwrap() > 0.);
            }
        }
    }
    #[test]
    fn invalid_cancelled_tight_and_limited_queries_leave_the_retained_source_unchanged() {
        use std::cell::Cell;
        let t = tool(false, false, 1);
        let before = t;
        let options = FacePropertyOptions::default();
        assert!(matches!(
            t.face_properties_with_cancel([0., 4.], options, &|| true),
            Err(FacePropertyError::Cancelled)
        ));
        let calls = Cell::new(0);
        assert!(matches!(
            t.face_properties_with_cancel([0., 4.], options, &|| {
                calls.set(calls.get() + 1);
                calls.get() > 20
            }),
            Err(FacePropertyError::Cancelled)
        ));
        for axial in [[-1., 4.], [0., 5.], [1., 1.], [f64::NAN, 4.]] {
            assert!(matches!(
                t.face_properties_with_cancel(axial, options, &|| false),
                Err(FacePropertyError::Invalid)
            ));
        }
        for options in [
            FacePropertyOptions {
                maximum_cells: 0,
                ..options
            },
            FacePropertyOptions {
                maximum_cells: 65537,
                ..options
            },
            FacePropertyOptions {
                area_tolerance_mm2: 0.,
                ..options
            },
            FacePropertyOptions {
                first_moment_tolerance_mm3: f64::INFINITY,
                ..options
            },
            FacePropertyOptions {
                centroid_tolerance_mm: f64::NAN,
                ..options
            },
        ] {
            assert!(matches!(
                t.face_properties_with_cancel([0., 4.], options, &|| false),
                Err(FacePropertyError::Invalid)
            ));
        }
        assert!(matches!(
            t.face_properties_with_cancel(
                [0., 4.],
                FacePropertyOptions {
                    maximum_cells: 1,
                    ..options
                },
                &|| false
            ),
            Err(FacePropertyError::Limit)
        ));
        assert!(
            t.face_properties_with_cancel(
                [0., 4.],
                FacePropertyOptions {
                    area_tolerance_mm2: 1e-30,
                    ..options
                },
                &|| false
            )
            .is_err()
        );
        assert_eq!(t, before);
        assert!(
            t.face_properties_with_cancel([0., 4.], options, &|| false)
                .is_ok()
        );
    }
}
