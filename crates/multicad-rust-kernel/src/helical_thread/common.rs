//! The common material of an external thread and a coaxial threaded cavity.
//!
//! Both sources are periodic in the phase `z - rate * theta`. The external
//! material is `r <= R_e(phase)` and the cavity material is
//! `r >= R_c(phase - shift)`, so at every height their common material in
//! the meridian half plane `(phase, r)` is the set where
//! `R_c(phase - shift) < r < R_e(phase)`. Away from a rounded root both
//! profiles are piecewise affine, so each connected part of that set is a
//! lens bounded by affine pieces between exact crossings of the two laws.
//! A rounded external root inside a lens, a crossing at a profile vertex
//! and coincident flanks are refused rather than approximated.
use super::Thread;

/// One affine boundary piece `phase -> r` of a lens, from its first to its
/// second endpoint.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MeridianPiece {
    pub phase_mm: [f64; 2],
    pub radius_mm: [f64; 2],
}

/// One closed lens of common material in the meridian half plane: the
/// cavity law from left to right, then the external law back.
#[derive(Clone, Debug, PartialEq)]
pub struct Lens {
    pub pieces: Vec<MeridianPiece>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CommonError {
    /// The two threads do not share pitch, starts and hand.
    Invalid,
    /// A vertex crossing, a tangential contact, a lens thinner than the
    /// requested resolution or a lens on the rounded root.
    Unsupported,
}

impl Thread {
    /// Breaks of one profile period, with the band each interval belongs to.
    fn bands(&self) -> [(f64, f64, usize); 5] {
        let b = self.breaks;
        std::array::from_fn(|i| (b[i], b[i + 1], i))
    }

    /// Whether a phase interval of this profile is affine.
    fn affine_band(&self, band: usize) -> bool {
        band != 2 || self.root.is_none()
    }

    /// The lenses of common material over one phase period, each listed
    /// once per thread start so that they tile one lead. `resolution_mm`
    /// is the smallest radial separation of the two laws, away from their
    /// crossings, that the caller can represent.
    pub fn common_lenses(
        &self,
        cavity: &Thread,
        shift_mm: f64,
        resolution_mm: f64,
    ) -> Result<Vec<Lens>, CommonError> {
        let pitch = self.spec.pitch_mm;
        if cavity.spec.pitch_mm != pitch
            || cavity.spec.starts != self.spec.starts
            || cavity.spec.left_hand != self.spec.left_hand
            || cavity.root.is_some()
            || !shift_mm.is_finite()
            || !(resolution_mm.is_finite() && resolution_mm > 0.)
        {
            return Err(CommonError::Invalid);
        }
        let shift = shift_mm.rem_euclid(pitch);
        let external = |phase: f64| self.radius_at_phase(phase);
        let inner = |phase: f64| cavity.radius_at_phase(phase - shift);
        let gap = |phase: f64| external(phase) - inner(phase);

        // Every phase where either law changes its affine piece.
        let mut knots: Vec<f64> = self
            .breaks
            .into_iter()
            .chain(
                cavity
                    .breaks
                    .into_iter()
                    .map(|b| (b + shift).rem_euclid(pitch)),
            )
            .filter(|&b| b < pitch)
            .collect();
        knots.sort_by(f64::total_cmp);
        knots.dedup_by(|a, b| (*a - *b).abs() <= 64. * f64::EPSILON * pitch);
        // A vertex where the laws meet, or meet within the resolution, is
        // a contact the affine lens cannot separate.
        let scale = 64. * f64::EPSILON * (self.spec.major_radius_mm + pitch);
        let margin = resolution_mm.max(scale);
        if knots.iter().any(|&k| gap(k).abs() <= margin) {
            return Err(CommonError::Unsupported);
        }
        // Start the period outside the common material.
        let Some(start) = knots.iter().position(|&k| gap(k) < 0.) else {
            // Overlap at every phase would leave no external crest.
            return Err(CommonError::Unsupported);
        };
        let count = knots.len();
        let at = |i: usize| knots[(start + i) % count] + pitch * ((start + i) / count) as f64;

        let band_of = |profile: &Thread, phase: f64| {
            let local = phase.rem_euclid(pitch);
            profile
                .bands()
                .into_iter()
                .find(|&(a, b, _)| local >= a && local < b)
                .map_or(4, |(_, _, band)| band)
        };
        // Crossings in order, each entering or leaving the common material.
        let mut crossings = Vec::new();
        for i in 0..count {
            let [a, b] = [at(i), at(i + 1)];
            let middle = a.midpoint(b);
            if !self.affine_band(band_of(self, middle)) {
                // The rounded root is convex, so the gap is convex here and
                // its maximum lies at the interval ends.
                if gap(a) >= 0. || gap(b) >= 0. {
                    return Err(CommonError::Unsupported);
                }
                continue;
            }
            let [ga, gb] = [gap(a), gap(b)];
            if ga.signum() != gb.signum() {
                let phase = a + (b - a) * ga / (ga - gb);
                crossings.push((phase, gb > 0.));
            }
        }
        if crossings.is_empty() {
            return Ok(Vec::new());
        }
        let mut lenses = Vec::new();
        for pair in crossings.chunks(2) {
            let [(enter, true), (leave, false)] = pair else {
                return Err(CommonError::Unsupported);
            };
            let lens = self.lens(cavity, shift, [*enter, *leave], margin)?;
            for start in 0..self.spec.starts {
                let offset = pitch * f64::from(start);
                lenses.push(Lens {
                    pieces: lens
                        .iter()
                        .map(|piece| MeridianPiece {
                            phase_mm: piece.phase_mm.map(|p| p + offset),
                            ..*piece
                        })
                        .collect(),
                });
            }
        }
        Ok(lenses)
    }

    fn lens(
        &self,
        cavity: &Thread,
        shift: f64,
        [enter, leave]: [f64; 2],
        margin: f64,
    ) -> Result<Vec<MeridianPiece>, CommonError> {
        let pitch = self.spec.pitch_mm;
        let external = |phase: f64| self.radius_at_phase(phase);
        let inner = |phase: f64| cavity.radius_at_phase(phase - shift);
        // Profile vertices strictly inside the lens, lifted past `enter`.
        let inside = |breaks: [f64; 6], offset: f64| {
            let mut points: Vec<f64> = (-1..=2)
                .flat_map(|turn| {
                    breaks
                        .into_iter()
                        .map(move |b| b + offset + pitch * f64::from(turn))
                })
                .filter(|&p| p > enter && p < leave)
                .collect();
            points.sort_by(f64::total_cmp);
            points.dedup_by(|a, b| (*a - *b).abs() <= 64. * f64::EPSILON * pitch);
            points
        };
        let chain = |law: &dyn Fn(f64) -> f64, vertices: Vec<f64>, [from, to]: [f64; 2]| {
            let mut stations = vec![from];
            stations.extend(vertices);
            stations.push(to);
            stations
                .windows(2)
                .map(|pair| MeridianPiece {
                    phase_mm: [pair[0], pair[1]],
                    radius_mm: [law(pair[0]), law(pair[1])],
                })
                .collect::<Vec<_>>()
        };
        // Both laws meet at the crossings; take one radius for both chains.
        let ends = [enter, leave].map(|p| external(p).midpoint(inner(p)));
        let mut lower = chain(&inner, inside(cavity.breaks, shift), [enter, leave]);
        let mut upper = chain(&external, inside(self.breaks, 0.), [enter, leave]);
        for chain in [&mut lower, &mut upper] {
            chain.first_mut().expect("lens piece").radius_mm[0] = ends[0];
            chain.last_mut().expect("lens piece").radius_mm[1] = ends[1];
        }
        // Each piece must lie on one affine band of its law, and the lens
        // must be thicker than the resolution between its crossings.
        let pieces: Vec<_> = lower
            .into_iter()
            .chain(upper.into_iter().rev().map(|piece| MeridianPiece {
                phase_mm: [piece.phase_mm[1], piece.phase_mm[0]],
                radius_mm: [piece.radius_mm[1], piece.radius_mm[0]],
            }))
            .collect();
        if pieces
            .iter()
            .any(|piece| (piece.phase_mm[1] - piece.phase_mm[0]).abs() <= margin)
        {
            return Err(CommonError::Unsupported);
        }
        let thickest = inside(self.breaks, 0.)
            .into_iter()
            .chain(inside(cavity.breaks, shift))
            .map(|p| external(p) - inner(p))
            .chain([external(enter.midpoint(leave)) - inner(enter.midpoint(leave))])
            .fold(f64::NEG_INFINITY, f64::max);
        if thickest <= margin {
            return Err(CommonError::Unsupported);
        }
        Ok(pieces)
    }
}

#[cfg(test)]
mod tests;
