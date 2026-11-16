//! Bounded native helices with flat or circular crest/root axial profiles.
//! These are geometric inputs; standard families and fit tables belong to core.
pub use crate::ffi::ThreadFactoryStats;
use crate::{GeometryError, GeometryErrorKind, GeometryPolicy, Solid, ffi};
use multicad_core::document::EntityId;

/// Diagnostic work counts for the current native owner thread. No native
/// handles or cache mutation cross the thread boundary through these values.
pub fn thread_factory_stats() -> Result<ThreadFactoryStats, GeometryError> {
    ffi::thread_factory_stats().map_err(|error| GeometryError {
        operation: "thread factory metrics",
        kind: GeometryErrorKind::Native,
        message: error.to_string(),
    })
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HelicalThread {
    /// Crest radius at the base plane z=0.
    pub major_radius_mm: f64,
    /// Axial spacing between adjacent threads, independent of start count.
    pub pitch_mm: f64,
    pub length_mm: f64,
    pub depth_mm: f64,
    pub crest_width_mm: f64,
    /// Flank angles from the radial direction in the axial section.
    pub lower_flank_deg: f64,
    pub upper_flank_deg: f64,
    /// Change of radius per unit axial length, rather than diameter taper.
    pub radial_taper: f64,
    pub starts: u32,
    pub left_hand: bool,
}
/// Circular crest/root profile in an axial meridian, with true 55° flanks.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RoundedHelicalThread {
    pub major_radius_mm: f64,
    pub pitch_mm: f64,
    pub length_mm: f64,
    pub depth_mm: f64,
    pub rounding_radius_mm: f64,
    /// Half of the included angle, measured from the radial direction.
    pub flank_angle_deg: f64,
    /// Radius change per axial millimetre: zero or signed 1:16 diameter taper.
    pub radial_taper: f64,
    pub starts: u32,
    pub left_hand: bool,
}
/// Actual axial flank angles with cone-parallel flats. This distinct adapter
/// preserves legacy Custom's existing conical profile/transport convention.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AxialFlatHelicalThread {
    pub profile: HelicalThread,
    /// Optional physical stock radius which clips external crests. It never
    /// changes the nominal cone or its pitch diameter.
    pub crest_radius_limit_mm: Option<f64>,
}
/// Parallel 3/30 or 30/30-degree axial flanks with a flat crest and tangent circular
/// valley. The external source family and load direction belong to core.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RootRoundedHelicalThread {
    pub profile: HelicalThread,
    pub root_radius_mm: f64,
}
/// Parallel 30-degree knuckle form with distinct tangent crest/root circles.
/// The material/cutter convention and nominal source diameter belong to core.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct KnuckleHelicalThread {
    pub major_radius_mm: f64,
    pub pitch_mm: f64,
    pub length_mm: f64,
    pub crest_radius_mm: f64,
    pub root_radius_mm: f64,
    pub left_hand: bool,
}
impl KnuckleHelicalThread {
    fn checked(self, policy: GeometryPolicy) -> Result<ffi::ThreadData, GeometryError> {
        let angle = 15_f64.to_radians();
        let basic = self.pitch_mm * (angle.cos() - angle.sin()) / (4. * (1. - angle.sin()));
        let delta = self.pitch_mm * 0.05 * angle.sin() / (1. - angle.sin());
        let matches = |actual: f64, expected: f64| {
            actual.is_finite() && (actual - expected).abs() <= expected.abs() * 1e-12
        };
        let external = matches(self.crest_radius_mm, basic) && matches(self.root_radius_mm, basic);
        let internal = matches(self.crest_radius_mm, basic - delta)
            && matches(self.root_radius_mm, basic + delta);
        if !(external || internal) {
            return Err(failure(
                GeometryErrorKind::Input,
                "Knuckle thread needs the exact 30-degree tangent external or internal circle pair",
            ));
        }
        policy.length(self.crest_radius_mm)?;
        policy.length(self.root_radius_mm)?;
        let mut data = HelicalThread {
            major_radius_mm: self.major_radius_mm,
            pitch_mm: self.pitch_mm,
            length_mm: self.length_mm,
            depth_mm: self.pitch_mm / 2.,
            crest_width_mm: 2. * self.crest_radius_mm * angle.cos(),
            lower_flank_deg: 15.,
            upper_flank_deg: 15.,
            radial_taper: 0.,
            starts: 1,
            left_hand: self.left_hand,
        }
        .checked(policy)?;
        data.rounding = self.crest_radius_mm;
        data.knuckle_root_radius = self.root_radius_mm;
        Ok(data)
    }
}
fn failure(kind: GeometryErrorKind, message: impl Into<String>) -> GeometryError {
    GeometryError {
        operation: "thread",
        kind,
        message: message.into(),
    }
}
impl HelicalThread {
    pub fn lead_mm(self) -> f64 {
        self.pitch_mm * f64::from(self.starts)
    }
    pub(crate) fn checked(self, policy: GeometryPolicy) -> Result<ffi::ThreadData, GeometryError> {
        policy.validate()?;
        for length in [
            self.major_radius_mm,
            self.pitch_mm,
            self.length_mm,
            self.depth_mm,
            self.crest_width_mm,
        ] {
            policy.length(length)?;
        }
        if !(1..=8).contains(&self.starts)
            || ![self.lower_flank_deg, self.upper_flank_deg]
                .into_iter()
                .all(|angle| angle.is_finite() && (1. ..=80.).contains(&angle))
            || !self.radial_taper.is_finite()
            || self.radial_taper.abs() > 0.1
        {
            return Err(failure(
                GeometryErrorKind::Input,
                "Thread needs 1–8 starts, finite 1–80 degree flanks and radius taper within ±0.1",
            ));
        }
        if self.pitch_mm < 0.001
            || self.length_mm < self.pitch_mm
            || self.length_mm / self.pitch_mm + 2. * f64::from(self.starts) > 128.
        {
            return Err(failure(
                GeometryErrorKind::Limit,
                "Thread needs pitch ≥0.001 mm, at least one pitch of length, and at most 128 total extended turns",
            ));
        }
        let lower_slope = self.lower_flank_deg.to_radians().tan();
        let upper_slope = self.upper_flank_deg.to_radians().tan();
        let overlap = (self.pitch_mm * 0.02).max(policy.linear_tolerance_mm * 64.);
        let root = self.major_radius_mm - self.depth_mm;
        let width = self.crest_width_mm + (self.depth_mm + overlap) * (lower_slope + upper_slope);
        if self.crest_width_mm < policy.linear_tolerance_mm * 100.
            || root - overlap < policy.minimum_feature_mm
            || width >= self.pitch_mm - policy.minimum_feature_mm
        {
            return Err(failure(
                GeometryErrorKind::Input,
                "Thread profile needs a positive root, resolvable crest and separated adjacent starts",
            ));
        }
        for z in [-2. * self.lead_mm(), self.length_mm + 2. * self.lead_mm()] {
            let radius = root + self.radial_taper * z - overlap;
            policy.length(radius)?;
            policy.point([self.major_radius_mm + self.radial_taper * z, 0., z])?;
        }
        Ok(ffi::ThreadData {
            radius: self.major_radius_mm,
            pitch: self.pitch_mm,
            height: self.length_mm,
            depth: self.depth_mm,
            crest_width: self.crest_width_mm,
            rounding: 0.,
            lower_slope,
            upper_slope,
            taper: self.radial_taper,
            overlap,
            tolerance: policy.linear_tolerance_mm,
            crest_limit: 0.,
            root_rounding: 0.,
            knuckle_root_radius: 0.,
            starts: self.starts,
            left: self.left_hand,
            axial_flat: false,
        })
    }
}
impl RootRoundedHelicalThread {
    pub(crate) fn checked(self, policy: GeometryPolicy) -> Result<ffi::ThreadData, GeometryError> {
        let p = self.profile;
        let mut data = p.checked(policy)?;
        policy.length(self.root_radius_mm)?;
        if p.starts != 1
            || p.radial_taper != 0.
            || ![(3., 30.), (30., 3.), (30., 30.)].contains(&(p.lower_flank_deg, p.upper_flank_deg))
        {
            return Err(failure(
                GeometryErrorKind::Input,
                "Root-rounded construction requires a parallel single-start 3/30 or 30/30-degree profile",
            ));
        }
        let lower = p.lower_flank_deg.to_radians();
        let upper = p.upper_flank_deg.to_radians();
        let intersection = (p.pitch_mm - p.crest_width_mm) / (lower.tan() + upper.tan());
        let half = (lower + upper) / 2.;
        let tilt = (upper - lower) / 2.;
        let center_z = p.crest_width_mm / 2. + intersection * upper.tan()
            - self.root_radius_mm * tilt.sin() / half.sin();
        let center_r =
            p.major_radius_mm - intersection + self.root_radius_mm * tilt.cos() / half.sin();
        let depth = p.major_radius_mm - center_r + self.root_radius_mm;
        let lo = center_z - self.root_radius_mm * upper.cos();
        let hi = center_z + self.root_radius_mm * lower.cos();
        if (depth - p.depth_mm).abs() > p.depth_mm.abs() * 1e-12
            || lo <= p.crest_width_mm / 2. + policy.minimum_feature_mm
            || hi >= p.pitch_mm - p.crest_width_mm / 2. - policy.minimum_feature_mm
        {
            return Err(failure(
                GeometryErrorKind::Input,
                "Root circle must be tangent to both flanks and preserve their straight sections",
            ));
        }
        data.root_rounding = self.root_radius_mm;
        Ok(data)
    }
}
impl AxialFlatHelicalThread {
    fn checked(self, policy: GeometryPolicy) -> Result<ffi::ThreadData, GeometryError> {
        if self.profile.lower_flank_deg != 30. || self.profile.upper_flank_deg != 30. {
            return Err(failure(
                GeometryErrorKind::Input,
                "Axial flat pipe construction currently requires the sourced 60-degree profile",
            ));
        }
        let p = self.profile;
        let mut data = p.checked(policy)?;
        if p.starts != 1 || p.radial_taper.abs() != 1. / 32. {
            return Err(failure(
                GeometryErrorKind::Input,
                "Axial flat conical threads require one start and signed 1:16 diameter taper",
            ));
        }
        let lower = 1. / p.lower_flank_deg.to_radians().tan() - p.radial_taper;
        let upper = 1. / p.upper_flank_deg.to_radians().tan() + p.radial_taper;
        if lower <= 0. || upper <= 0. {
            return Err(failure(
                GeometryErrorKind::Input,
                "Axial thread flank crosses its cone",
            ));
        }
        data.lower_slope = 1. / lower;
        data.upper_slope = 1. / upper;
        let width =
            p.crest_width_mm + (p.depth_mm + data.overlap) * (data.lower_slope + data.upper_slope);
        if width >= p.pitch_mm - policy.minimum_feature_mm {
            return Err(failure(
                GeometryErrorKind::Input,
                "Axial thread profile overlaps adjacent turns",
            ));
        }
        if let Some(limit) = self.crest_radius_limit_mm {
            policy.length(limit)?;
            let root = p.major_radius_mm - p.depth_mm;
            let maximum_root = [-2. * p.lead_mm(), p.length_mm + 2. * p.lead_mm()]
                .into_iter()
                .map(|z| root + p.radial_taper * z)
                .fold(f64::NEG_INFINITY, f64::max);
            if limit <= maximum_root + policy.minimum_feature_mm * 2. {
                return Err(failure(
                    GeometryErrorKind::Input,
                    "Crest clipping must preserve the extended thread core",
                ));
            }
            data.crest_limit = limit;
        }
        data.axial_flat = true;
        Ok(data)
    }
}
impl RoundedHelicalThread {
    fn checked(self, policy: GeometryPolicy) -> Result<ffi::ThreadData, GeometryError> {
        let angle = 27.5_f64.to_radians();
        let depth = self.pitch_mm / (3. * angle.tan());
        let k = self.radial_taper;
        let rounding = if k == 0. {
            self.pitch_mm * angle.cos() / (12. * (1. - angle.sin()))
        } else {
            let m = 1. / angle.tan();
            let fundamental = self.pitch_mm * (m * m - k * k) / (2. * m);
            (fundamental - depth) / (2. * (1. / angle.sin() - (1. + k * k).sqrt()))
        };
        let matches = |actual: f64, expected: f64| {
            actual.is_finite() && (actual - expected).abs() <= expected.abs() * 1e-12
        };
        if self.starts != 1
            || ![0., 1. / 32., -1. / 32.].contains(&k)
            || !self.flank_angle_deg.is_finite()
            || (self.flank_angle_deg - 27.5).abs() > 1e-10
            || !matches(self.depth_mm, depth)
            || !matches(self.rounding_radius_mm, rounding)
        {
            return Err(failure(
                GeometryErrorKind::Input,
                "Rounded thread needs the exact 55-degree single-start parallel or 1:16 tapered profile",
            ));
        }
        policy.length(self.rounding_radius_mm)?;
        let mut data = HelicalThread {
            major_radius_mm: self.major_radius_mm,
            pitch_mm: self.pitch_mm,
            length_mm: self.length_mm,
            depth_mm: self.depth_mm,
            crest_width_mm: 2. * self.rounding_radius_mm * angle.cos(),
            lower_flank_deg: self.flank_angle_deg,
            upper_flank_deg: self.flank_angle_deg,
            radial_taper: k,
            starts: self.starts,
            left_hand: self.left_hand,
        }
        .checked(policy)?;
        data.rounding = self.rounding_radius_mm;
        Ok(data)
    }
}
impl Solid {
    pub fn knuckle_helical_thread_named(
        spec: KnuckleHelicalThread,
        feature: EntityId,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        if feature.is_nil() {
            return Err(failure(
                GeometryErrorKind::Input,
                "Thread feature identity must not be nil",
            ));
        }
        let native = ffi::helical_thread(&spec.checked(policy)?)
            .map_err(|error| failure(GeometryErrorKind::Native, error.to_string()))?;
        let mut solid = Self::adopt(native, policy, "knuckle thread")?;
        if solid.solid_count()? != 1 {
            return Err(failure(
                GeometryErrorKind::InvalidResult,
                "Thread must remain one connected solid",
            ));
        }
        solid.name_thread(feature, 1)?;
        let _ = ffi::remember_validated_thread(solid.occt()?);
        Ok(solid)
    }
    pub fn root_rounded_helical_thread_named(
        spec: RootRoundedHelicalThread,
        feature: EntityId,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        if feature.is_nil() {
            return Err(failure(
                GeometryErrorKind::Input,
                "Thread feature identity must not be nil",
            ));
        }
        let data = spec.checked(policy)?;
        let cached = ffi::cached_helical_thread(&data)
            .map_err(|error| failure(GeometryErrorKind::Native, error.to_string()))?;
        let native = if !cached.is_null() {
            cached
        } else {
            ffi::helical_thread_after_miss(&data)
                .map_err(|error| failure(GeometryErrorKind::Native, error.to_string()))?
        };
        let mut solid = Self::adopt(native, policy, "root-rounded thread")?;
        if solid.solid_count()? != 1 {
            return Err(failure(
                GeometryErrorKind::InvalidResult,
                "Thread must remain one connected solid",
            ));
        }
        solid.name_thread(feature, 1)?;
        let _ = ffi::remember_validated_thread(solid.occt()?);
        Ok(solid)
    }
    /// Developer derivation path preserving every native root/physical guard
    /// while bypassing compiled assets and owner-thread factory retention.
    pub fn root_rounded_helical_thread_named_fresh(
        spec: RootRoundedHelicalThread,
        feature: EntityId,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        if feature.is_nil() {
            return Err(failure(
                GeometryErrorKind::Input,
                "Thread feature identity must not be nil",
            ));
        }
        let native = ffi::helical_thread_fresh(&spec.checked(policy)?)
            .map_err(|error| failure(GeometryErrorKind::Native, error.to_string()))?;
        let mut solid = Self::adopt(native, policy, "fresh root-rounded thread")?;
        if solid.solid_count()? != 1 {
            return Err(failure(
                GeometryErrorKind::InvalidResult,
                "Thread must remain one connected solid",
            ));
        }
        solid.name_thread(feature, 1)?;
        Ok(solid)
    }
    pub fn axial_flat_helical_thread_named(
        spec: AxialFlatHelicalThread,
        feature: EntityId,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        if feature.is_nil() {
            return Err(failure(
                GeometryErrorKind::Input,
                "Thread feature identity must not be nil",
            ));
        }
        let native = ffi::helical_thread(&spec.checked(policy)?)
            .map_err(|error| failure(GeometryErrorKind::Native, error.to_string()))?;
        let mut solid = Self::adopt(native, policy, "axial flat thread")?;
        if solid.solid_count()? != 1 {
            return Err(failure(
                GeometryErrorKind::InvalidResult,
                "Thread must remain one connected solid",
            ));
        }
        solid.name_thread(feature, spec.profile.starts)?;
        let _ = ffi::remember_validated_thread(solid.occt()?);
        Ok(solid)
    }
    pub fn rounded_helical_thread_named(
        spec: RoundedHelicalThread,
        feature: EntityId,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        if feature.is_nil() {
            return Err(failure(
                GeometryErrorKind::Input,
                "Thread feature identity must not be nil",
            ));
        }
        let native = ffi::helical_thread(&spec.checked(policy)?)
            .map_err(|error| failure(GeometryErrorKind::Native, error.to_string()))?;
        let mut solid = Self::adopt(native, policy, "rounded thread")?;
        if solid.solid_count()? != 1 {
            return Err(failure(
                GeometryErrorKind::InvalidResult,
                "Thread must remain one connected solid",
            ));
        }
        solid.name_thread(feature, spec.starts)?;
        // The optional factory only remembers geometry after every ordinary
        // adoption, policy and owner-specific naming guard succeeded. Cache
        // allocation/serialization failure must not invalidate a valid solid.
        let _ = ffi::remember_validated_thread(solid.occt()?);
        Ok(solid)
    }
    /// Build a real helical B-rep around +Z with flat cut ends. Native face
    /// ancestry names each start's flanks/crest, the root and both end caps.
    pub fn helical_thread_named(
        spec: HelicalThread,
        feature: EntityId,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        if feature.is_nil() {
            return Err(failure(
                GeometryErrorKind::Input,
                "Thread feature identity must not be nil",
            ));
        }
        let data = spec.checked(policy)?;
        // Record one request, preserving the existing owner-thread factory hit
        // before attempting any compiled prototype reads or measurements.
        let cached = ffi::cached_helical_thread(&data)
            .map_err(|error| failure(GeometryErrorKind::Native, error.to_string()))?;
        let native = if !cached.is_null() {
            cached
        } else {
            ffi::helical_thread_after_miss(&data)
                .map_err(|error| failure(GeometryErrorKind::Native, error.to_string()))?
        };
        let mut solid = Self::adopt(native, policy, "thread")?;
        if solid.solid_count()? != 1 {
            return Err(failure(
                GeometryErrorKind::InvalidResult,
                "Thread must remain one connected solid",
            ));
        }
        solid.name_thread(feature, spec.starts)?;
        let _ = ffi::remember_validated_thread(solid.occt()?);
        Ok(solid)
    }
    /// Developer derivation path: rebuild every native thread and bypass all
    /// compiled prototypes and factory retention while keeping normal guards.
    pub fn helical_thread_named_fresh(
        spec: HelicalThread,
        feature: EntityId,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        if feature.is_nil() {
            return Err(failure(
                GeometryErrorKind::Input,
                "Thread feature identity must not be nil",
            ));
        }
        let data = spec.checked(policy)?;
        let native = ffi::helical_thread_fresh(&data)
            .map_err(|error| failure(GeometryErrorKind::Native, error.to_string()))?;
        let mut solid = Self::adopt(native, policy, "thread")?;
        if solid.solid_count()? != 1 {
            return Err(failure(
                GeometryErrorKind::InvalidResult,
                "Thread must remain one connected solid",
            ));
        }
        solid.name_thread(feature, spec.starts)?;
        Ok(solid)
    }
}
