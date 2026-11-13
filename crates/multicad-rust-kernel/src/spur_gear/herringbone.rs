//! A continuous twisted gear formed from one exact transverse tooth outline.
use super::{
    Edge, Error, Hub, HubExtension, Point, Profile, Spec, TAU, ToothEdge, normal_sweep, rotate,
};
use std::f64::consts::PI;

mod mesh;
mod support;
mod topology;

pub use topology::{Band, Curve, Face, Ring, Topology, TopologyEdge};

/// A single solid with two tooth bands and optional bore. Each tooth
/// section is exactly the shared transverse profile rotated about Z. A
/// gapped chevron separates its bands by a central root-core cylinder.
#[derive(Clone, Debug, PartialEq)]
pub struct Herringbone {
    profile: Profile,
    tip_radius_mm: f64,
    twist_rad_per_mm: f64,
    helical: bool,
    gap_mm: f64,
}

impl Herringbone {
    pub fn new(
        spec: Spec,
        twist_rad_per_mm: f64,
        minimum_mm: f64,
        maximum_coordinate_mm: f64,
    ) -> Result<Self, Error> {
        Self::with_sweep(
            spec,
            twist_rad_per_mm,
            minimum_mm,
            maximum_coordinate_mm,
            false,
        )
    }

    /// Two opposite-hand bands of width (face − gap)/2 at both ends, joined
    /// by the root cylinder through the central gap. Each band ends at the
    /// phase reached by its own twist, so both gap shoulders align.
    pub fn new_gapped(
        spec: Spec,
        twist_rad_per_mm: f64,
        gap_mm: f64,
        minimum_mm: f64,
        maximum_coordinate_mm: f64,
    ) -> Result<Self, Error> {
        let band = (spec.face_width_mm - gap_mm) / 2.;
        if !gap_mm.is_finite()
            || gap_mm < minimum_mm
            || band < minimum_mm
            || (twist_rad_per_mm * band).abs() * spec.tip_radius_mm < minimum_mm
        {
            return Err(Error::InvalidDimensions);
        }
        let mut gear = Self::with_sweep(
            spec,
            twist_rad_per_mm,
            minimum_mm,
            maximum_coordinate_mm,
            false,
        )?;
        gear.gap_mm = gap_mm;
        Ok(gear)
    }

    /// A single-hand helical sweep. The middle ring is shared by two
    /// smoothly joined bands so topology remains stable across both laws.
    pub fn new_helical(
        spec: Spec,
        twist_rad_per_mm: f64,
        minimum_mm: f64,
        maximum_coordinate_mm: f64,
    ) -> Result<Self, Error> {
        Self::with_sweep(
            spec,
            twist_rad_per_mm,
            minimum_mm,
            maximum_coordinate_mm,
            true,
        )
    }

    fn with_sweep(
        spec: Spec,
        twist_rad_per_mm: f64,
        minimum_mm: f64,
        maximum_coordinate_mm: f64,
        helical: bool,
    ) -> Result<Self, Error> {
        if !twist_rad_per_mm.is_finite()
            || twist_rad_per_mm == 0.
            || (twist_rad_per_mm * spec.face_width_mm).abs() > PI
            || spec.face_width_mm / 2. < minimum_mm
            || (twist_rad_per_mm * spec.face_width_mm / 2.).abs() * spec.tip_radius_mm < minimum_mm
        {
            return Err(Error::InvalidDimensions);
        }
        let profile = Profile::new(spec, minimum_mm, maximum_coordinate_mm)?;
        Ok(Self {
            profile,
            tip_radius_mm: spec.tip_radius_mm,
            twist_rad_per_mm,
            helical,
            gap_mm: 0.,
        })
    }

    pub fn profile(&self) -> &Profile {
        &self.profile
    }

    /// Attach a coaxial hub contained radially within every tooth section.
    pub fn with_hub(
        mut self,
        hub: Hub,
        minimum_mm: f64,
        maximum_coordinate_mm: f64,
    ) -> Result<Self, Error> {
        self.profile = self
            .profile
            .with_hub(hub, minimum_mm, maximum_coordinate_mm)?;
        Ok(self)
    }

    pub fn axial_bounds(&self) -> [f64; 2] {
        self.profile.axial_bounds()
    }

    /// Extensions append after the two caps, both tooth bands, the bore
    /// and any central gap faces.
    pub fn hub_extensions(&self) -> [Option<HubExtension>; 2] {
        self.profile
            .hub_extensions_from(self.gap_faces_start() + self.gap_face_count())
    }

    /// Central gap width; zero for continuous bands.
    pub fn gap_mm(&self) -> f64 {
        self.gap_mm
    }

    pub fn is_gapped(&self) -> bool {
        self.gap_mm > 0.
    }

    pub fn teeth(&self) -> u32 {
        self.profile
            .outline()
            .last()
            .map_or(0, |edge| edge.tooth + 1)
    }

    /// Axial length of each tooth band.
    pub fn band_length_mm(&self) -> f64 {
        (self.profile.height_mm() - self.gap_mm) / 2.
    }

    pub fn band_interval_mm(&self, band: Band) -> [f64; 2] {
        let height = self.profile.height_mm();
        let length = self.band_length_mm();
        match band {
            Band::Lower => [0., length],
            Band::Upper => [height - length, height],
        }
    }

    /// The lower and upper gap shoulder planes of a gapped chevron.
    pub fn shoulder_heights_mm(&self) -> Option<[f64; 2]> {
        self.is_gapped()
            .then(|| [self.band_length_mm(), self.band_interval_mm(Band::Upper)[0]])
    }

    pub fn root_radius_mm(&self) -> f64 {
        self.profile.root_radius_mm
    }

    /// First face ordinal of the root core and its gap shoulders.
    pub fn gap_faces_start(&self) -> u32 {
        2 + 2 * self.profile.outline().len() as u32 + u32::from(self.profile.bore().is_some())
    }

    /// The root core and one shoulder per tooth at each band end.
    pub fn gap_face_count(&self) -> u32 {
        if self.is_gapped() {
            1 + 2 * self.teeth()
        } else {
            0
        }
    }

    pub fn tip_radius_mm(&self) -> f64 {
        self.tip_radius_mm
    }

    pub fn twist_rad_per_mm(&self) -> f64 {
        self.twist_rad_per_mm
    }

    pub fn is_helical(&self) -> bool {
        self.helical
    }

    pub fn band_twist_rad_per_mm(&self, band: Band) -> f64 {
        if band == Band::Upper && !self.helical {
            -self.twist_rad_per_mm
        } else {
            self.twist_rad_per_mm
        }
    }

    pub fn volume_mm3(&self) -> f64 {
        // Every transverse section has the same area; a rigid rotation
        // through depth leaves the sweep's Jacobian determinant equal to one.
        // The central gap keeps only the root annulus around the bore.
        let bore = self.profile.bore_radius_mm().unwrap_or(0.);
        let root = self.profile.root_radius_mm;
        self.profile.volume_mm3()
            - self.gap_mm * (self.profile.area_mm2() - PI * (root - bore) * (root + bore))
    }

    pub fn phase_at(&self, z_mm: f64) -> f64 {
        let height = self.profile.height_mm();
        self.twist_rad_per_mm
            * if self.helical {
                z_mm
            } else {
                z_mm.min(height - z_mm).min(self.band_length_mm())
            }
    }

    pub fn section_point(&self, point: Point, z_mm: f64) -> [f64; 3] {
        let [x, y] = rotate(point, self.phase_at(z_mm));
        [x, y, z_mm]
    }

    pub fn bounds(&self) -> [[f64; 3]; 2] {
        [-1., 1.].map(|sign| {
            std::array::from_fn(|i| {
                let mut direction = [0.; 3];
                direction[i] = sign;
                sign * self.support(direction)
            })
        })
    }
}

#[cfg(test)]
mod gap_tests;
#[cfg(test)]
mod hub_tests;
#[cfg(test)]
pub(super) mod tests;
