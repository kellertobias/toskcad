//! Hex stock around or above the literal rounded thread. The hexagon is the
//! regular one of `headed_screw::hex_vertices`, with its first corner at 30°.
use super::{Error, Thread};
use crate::headed_screw::{hex_area, hex_radius, hex_vertices};
use multicad_precision::PortableMath as _;

const MAX_DIMENSION_MM: f64 = 1e7;

/// Which side of the thread the hex prism occupies.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum HexForm {
    /// Hex prism over the full thread length minus the rounded cutter.
    Nut,
    /// Fully threaded shaft `0 <= z <= L` joined at its `z = 0` end to a
    /// solid hex head occupying `-height <= z <= 0`.
    Head { height_mm: f64 },
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HexStock {
    thread: Thread,
    across_flats_mm: f64,
    form: HexForm,
}

impl HexStock {
    /// The flats must strictly clear the thread's major radius: a tangent
    /// zero-width wall or seat is a different, non-regular boundary.
    pub fn new(thread: Thread, across_flats_mm: f64, form: HexForm) -> Result<Self, Error> {
        if !across_flats_mm.is_finite()
            || across_flats_mm / 2. <= thread.spec().major_radius_mm
            || matches!(form, HexForm::Head { height_mm } if !height_mm.is_finite() || height_mm <= 0.)
        {
            return Err(Error::Invalid);
        }
        let height = match form {
            HexForm::Nut => 0.,
            HexForm::Head { height_mm } => height_mm,
        };
        if across_flats_mm / 3_f64.sqrt() > MAX_DIMENSION_MM
            || height + thread.spec().length_mm > MAX_DIMENSION_MM
        {
            return Err(Error::Limit);
        }
        Ok(Self {
            thread,
            across_flats_mm,
            form,
        })
    }
    pub fn thread(&self) -> &Thread {
        &self.thread
    }
    pub fn across_flats_mm(&self) -> f64 {
        self.across_flats_mm
    }
    pub fn form(&self) -> HexForm {
        self.form
    }
    pub fn length_mm(&self) -> f64 {
        self.thread.spec().length_mm
    }
    /// Axial interval of the hex prism.
    pub fn hex_interval_mm(&self) -> [f64; 2] {
        match self.form {
            HexForm::Nut => [0., self.length_mm()],
            HexForm::Head { height_mm } => [-height_mm, 0.],
        }
    }
    pub fn hex_height_mm(&self) -> f64 {
        let [a, b] = self.hex_interval_mm();
        b - a
    }
    pub fn hex_vertices(&self) -> [[f64; 2]; 6] {
        hex_vertices(self.across_flats_mm)
    }
    pub fn hex_area_mm2(&self) -> f64 {
        hex_area(self.across_flats_mm)
    }
    pub fn outer_radius(&self, theta: f64) -> f64 {
        hex_radius(self.across_flats_mm, theta)
    }
    pub fn volume_mm3(&self) -> f64 {
        let prism = self.hex_area_mm2() * self.hex_height_mm();
        match self.form {
            HexForm::Nut => prism - self.thread.volume_mm3(),
            HexForm::Head { .. } => prism + self.thread.volume_mm3(),
        }
    }
    /// Logical faces: both end planes and the five thread bands of the thread
    /// numbering (0..=6), then the six hex walls (7..=12) counterclockwise from
    /// the wall between the 30° and 90° corners. A head adds its top (13).
    pub fn face_count(&self) -> usize {
        match self.form {
            HexForm::Nut => 13,
            HexForm::Head { .. } => 14,
        }
    }
    pub fn wall_face(side: usize) -> u32 {
        7 + side as u32
    }
    /// Material membership with the same radial and axial tolerance as the
    /// thread. Boundaries are material.
    pub fn contains(&self, [x, y, z]: [f64; 3], tolerance: f64) -> bool {
        if !tolerance.is_finite() || tolerance < 0. || ![x, y, z].into_iter().all(f64::is_finite) {
            return false;
        }
        let [lo, hi] = self.hex_interval_mm();
        let radius = x.portable_hypot(y);
        let theta = y.portable_atan2(x);
        let in_hex = z >= lo - tolerance
            && z <= hi + tolerance
            && radius <= self.outer_radius(theta) + tolerance;
        match self.form {
            HexForm::Nut => {
                in_hex
                    && radius
                        >= self.thread.radius(theta, z.clamp(0., self.length_mm())) - tolerance
            }
            HexForm::Head { .. } => in_hex || self.thread.contains([x, y, z], tolerance),
        }
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::rounded_thread::Spec;

    pub(crate) fn thread(left_hand: bool) -> Thread {
        // DIN 405 Rd 12 x 1/10": 30° flanks, equal crest and root rounding.
        let p = 25.4 / 10.;
        let a = 15_f64.to_radians();
        let r = 0.238507 * p;
        let depth = p / (2. * a.portable_tan()) - 2. * r * (1. / a.portable_sin() - 1.);
        Thread::new(Spec {
            major_radius_mm: 6.,
            pitch_mm: p,
            length_mm: 3.37 * p,
            depth_mm: depth,
            crest_radius_mm: r,
            root_radius_mm: r,
            flank_angle_deg: 15.,
            starts: 1,
            left_hand,
        })
        .unwrap()
    }

    #[test]
    fn hex_stock_volume_membership_and_refusals() {
        let t = thread(true);
        let nut = HexStock::new(t, 30., HexForm::Nut).unwrap();
        let head = HexStock::new(t, 30., HexForm::Head { height_mm: 3. }).unwrap();
        let prism = 3_f64.sqrt() / 2. * 900.;
        let l = t.spec().length_mm;
        assert!((nut.volume_mm3() - (prism * l - t.volume_mm3())).abs() < 1e-9);
        assert!((head.volume_mm3() - (prism * 3. + t.volume_mm3())).abs() < 1e-9);
        assert_eq!((nut.face_count(), head.face_count()), (13, 14));
        // Axis is a void of the nut, material of the screw; corners are hex.
        assert!(!nut.contains([0., 0., l / 2.], 1e-9));
        assert!(head.contains([0., 0., l / 2.], 1e-9));
        assert!(head.contains([0., 0., -1.5], 1e-9));
        assert!(nut.contains([14.9, 0., l / 2.], 1e-9));
        assert!(!nut.contains([15.1, 0., l / 2.], 1e-9));
        assert!(nut.contains([0., 17.3, 0.], 1e-9));
        assert!(!head.contains([14.9, 0., 0.5], 1e-9));
        for (width, form) in [
            (12., HexForm::Nut),
            (f64::NAN, HexForm::Nut),
            (30., HexForm::Head { height_mm: 0. }),
        ] {
            assert_eq!(HexStock::new(t, width, form), Err(Error::Invalid));
        }
    }
}
