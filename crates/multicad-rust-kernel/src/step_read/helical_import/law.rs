//! Infer the physical helical law from literal tensor parameter domains and
//! radii. This header is not a substitute for verifying retained control nets.
use super::*;
pub(in super::super) struct Law {
    major: f64,
    pitch: f64,
    phase_start: f64,
    depth: f64,
    root_radius: Option<f64>,
    flank_degrees: [f64; 2],
}
impl Law {
    pub(in super::super) fn infer(
        major: f64,
        upper_radius: f64,
        lower_radius: f64,
        radial: &[[f64; 2]; 4],
        rate: f64,
        pitch: f64,
    ) -> Result<Self> {
        let phase_start = radial[0][0] * rate;
        let phase_upper = radial[0][1] * rate;
        let phase_lower = radial[2][0] * rate;
        let phase_crest = radial[2][1] * rate;
        let upper_run = (phase_upper - phase_start) / (major - upper_radius);
        let lower_run = (phase_crest - phase_lower) / (major - lower_radius);
        let (depth, root_radius, flank_degrees) = if radial[1][0] < 0. {
            let [a, b] = radial[1];
            require(
                a < 0.
                    && b > 0.
                    && a > -std::f64::consts::FRAC_PI_2
                    && b < std::f64::consts::FRAC_PI_2,
            )?;
            let f = (phase_lower - phase_upper) / (b.sin() - a.sin());
            // The literal circular-root domain carries its tangent angles;
            // avoid recovering depth by subtracting distant world radii.
            let lower = std::f64::consts::FRAC_PI_2 - b;
            let upper = std::f64::consts::FRAC_PI_2 + a;
            let half = (lower + upper) / 2.;
            let tilt = (upper - lower) / 2.;
            let intersection = (pitch - 2. * phase_start) / (lower.tan() + upper.tan());
            (
                intersection - f * (tilt.cos() / half.sin() - 1.),
                Some(f),
                [lower.to_degrees(), upper.to_degrees()],
            )
        } else {
            (
                major - upper_radius,
                None,
                [lower_run.atan().to_degrees(), upper_run.atan().to_degrees()],
            )
        };
        Ok(Self {
            major,
            pitch,
            phase_start,
            depth,
            root_radius,
            flank_degrees,
        })
    }
    pub(in super::super) fn thread(
        self,
        length: f64,
        starts: u32,
        left_hand: bool,
    ) -> Result<Thread> {
        Thread::new(Spec {
            major_radius_mm: self.major,
            pitch_mm: self.pitch,
            length_mm: length,
            depth_mm: self.depth,
            crest_width_mm: 2. * self.phase_start,
            lower_flank_deg: self.flank_degrees[0],
            upper_flank_deg: self.flank_degrees[1],
            starts,
            left_hand,
            root_radius_mm: self.root_radius,
        })
        .map_err(|e| match e {
            crate::helical_thread::ThreadError::Invalid => ReadError::Invalid,
            crate::helical_thread::ThreadError::Limit => ReadError::Limit,
        })
    }
}
