//! Exact disjoint spline columns added to or removed from a retained stock cap.
pub use crate::spline_sweep::ConstructionError as Error;
use crate::spline_sweep::{Edge, Sweep, validation};
mod boundary;
mod mesh;
pub use crate::normal_relief::Kind;
pub use crate::spline_sweep::ConstructionError;
pub use boundary::{Boundary, BoundaryEdge, BoundaryFace, Curve, EdgeUse, FaceRole, Surface};
#[cfg(test)]
pub(crate) mod tests;
type Point = [f64; 3];
pub const MAX_MOTIFS: usize = 192;
pub const MAX_EDGES: usize = 16_384;
#[derive(Clone, Debug, PartialEq)]
pub struct Relief {
    base: Sweep,
    motifs: Vec<Sweep>,
    kind: Kind,
    levels: Vec<f64>,
    tolerance: f64,
    minimum: f64,
    boundary_error: f64,
    axial_error: f64,
    volume: f64,
    centroid: [f64; 3],
    boundary: Boundary,
}
fn sum_error(a: f64, b: f64) -> (f64, f64) {
    let value = a + b;
    let part = value - a;
    (value, (a - (value - part)) + (b - part))
}
pub(super) fn remaining_budget(requested: f64, spent: f64) -> f64 {
    let (value, residual) = sum_error(requested, -spent);
    if residual < 0. {
        value.next_down()
    } else {
        value
    }
}

fn finish_metrics(volume: f64, moment: Point, origin: Point) -> Result<Point, Error> {
    let centroid = std::array::from_fn(|i| origin[i] + moment[i] / volume);
    if !volume.is_finite() || volume <= 0. || !centroid.into_iter().all(f64::is_finite) {
        return Err(Error::Unresolved);
    }
    Ok(centroid)
}

impl Relief {
    pub fn new(
        base: Sweep,
        motifs: Vec<Sweep>,
        kind: Kind,
        tolerance: f64,
        minimum: f64,
    ) -> Result<Self, Error> {
        Self::new_with_interrupt(base, motifs, kind, tolerance, minimum, None)
    }
    pub fn new_with_interrupt(
        base: Sweep,
        motifs: Vec<Sweep>,
        kind: Kind,
        tolerance: f64,
        minimum: f64,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, Error> {
        let mut budget = validation::Budget {
            remaining: validation::MAX_PAIR_WORK,
            interrupt,
        };
        validation::consume(&mut budget)?;
        if motifs.is_empty()
            || !tolerance.is_finite()
            || tolerance <= 0.
            || !minimum.is_finite()
            || minimum <= tolerance
            || base.height_mm() < minimum
        {
            return Err(Error::Invalid);
        }
        if motifs.len() > MAX_MOTIFS
            || base
                .edges()
                .count()
                .saturating_add(motifs.iter().map(|s| s.edges().count()).sum::<usize>())
                > MAX_EDGES
        {
            return Err(Error::Limit);
        }
        let mut levels = Vec::with_capacity(motifs.len());
        let join_error = std::iter::once(&base)
            .chain(&motifs)
            .map(Sweep::join_adjustment_mm)
            .fold(0_f64, f64::max);
        let mut axial_error = 0_f64;
        let mut outside = base.holes().to_vec();
        for sweep in std::iter::once(&base).chain(&motifs) {
            validation::refined_region_with_budget(
                sweep.outline(),
                sweep.holes(),
                tolerance,
                tolerance,
                &mut budget,
            )?;
        }
        for motif in &motifs {
            validation::consume(&mut budget)?;
            if motif.height_mm() < minimum {
                return Err(Error::Invalid);
            }
            let signed = if kind == Kind::Raised {
                motif.height_mm()
            } else {
                -motif.height_mm()
            };
            let (level, residual) = sum_error(base.height_mm(), signed);
            if kind == Kind::Recessed && (level < minimum || level == minimum && residual < 0.) {
                return Err(Error::Invalid);
            }
            if !level.is_finite() || level > 1e8 || level == 1e8 && residual > 0. {
                return Err(Error::Limit);
            }
            axial_error = axial_error.max(residual.abs());
            levels.push(level);
            outside.push(motif.outline().iter().rev().map(Edge::reversed).collect());
        }
        // A single certificate proves stock-hole and motif-footprint containment,
        // pairwise separation and nonnesting. Motif counters remain independent
        // certified islands, never holes discarded from the source geometry.
        validation::refined_region_with_budget(
            base.outline(),
            &outside,
            tolerance,
            minimum,
            &mut budget,
        )?;
        let boundary_error = crate::normal_relief::precision::add_up(join_error, axial_error);
        let sign = if kind == Kind::Raised { 1. } else { -1. };
        let mut volume = base.volume_mm3();
        let origin = base.centroid_mm();
        let mut moment = [0.; 3];
        for (motif, &level) in motifs.iter().zip(&levels) {
            let v = sign * motif.volume_mm3();
            volume += v;
            let c = motif.centroid_mm();
            let center = [c[0], c[1], base.height_mm().midpoint(level)];
            for i in 0..3 {
                moment[i] += v * (center[i] - origin[i]);
            }
        }
        let centroid = finish_metrics(volume, moment, origin)?;
        let mut result = Self {
            base,
            motifs,
            kind,
            levels,
            tolerance,
            minimum,
            boundary_error,
            axial_error,
            volume,
            centroid,
            boundary: Boundary::default(),
        };
        result.boundary = boundary::build(&result, &mut budget)?;
        Ok(result)
    }
    pub fn base(&self) -> &Sweep {
        &self.base
    }
    pub fn motifs(&self) -> &[Sweep] {
        &self.motifs
    }
    pub fn kind(&self) -> Kind {
        self.kind
    }
    pub fn motif_level_mm(&self, index: usize) -> Option<f64> {
        self.levels.get(index).copied()
    }
    pub fn tolerance_mm(&self) -> f64 {
        self.tolerance
    }
    pub fn minimum_mm(&self) -> f64 {
        self.minimum
    }
    pub fn boundary_error_mm(&self) -> f64 {
        self.boundary_error
    }
    pub fn boundary(&self) -> &Boundary {
        &self.boundary
    }
    pub fn volume_mm3(&self) -> f64 {
        self.volume
    }
    pub fn centroid_mm(&self) -> [f64; 3] {
        self.centroid
    }
    /// Classify a point away from the source precision enclosure. Points near
    /// a curved boundary return Unresolved, without tolerance-filled material.
    pub fn contains(&self, point: [f64; 3], tolerance: f64) -> Result<bool, Error> {
        if !point.into_iter().all(f64::is_finite) || !tolerance.is_finite() || tolerance <= 0. {
            return Err(Error::Invalid);
        }
        let mut budget = validation::Budget {
            remaining: validation::MAX_PAIR_WORK,
            interrupt: None,
        };
        let xy = [point[0], point[1]];
        let inside = |s: &Sweep, budget: &mut validation::Budget<'_>| -> Result<bool, Error> {
            let winding = s
                .rings()
                .map(|ring| validation::winding_profile(ring, xy, tolerance, budget))
                .sum::<Result<f64, Error>>()?;
            Ok(winding.abs() > std::f64::consts::PI)
        };
        if point[2] < 0. || point[2] > self.bounds()[1][2] || !inside(&self.base, &mut budget)? {
            return Ok(false);
        }
        for (m, motif) in self.motifs.iter().enumerate() {
            let active = if self.kind == Kind::Raised {
                point[2] > self.base.height_mm() && point[2] <= self.levels[m]
            } else {
                point[2] > self.levels[m]
            };
            let uncertain =
                self.axial_error > 0. && (point[2] - self.levels[m]).abs() <= self.axial_error;
            if (active || uncertain) && inside(motif, &mut budget)? {
                if uncertain {
                    return Err(Error::Unresolved);
                }
                return Ok(self.kind == Kind::Raised);
            }
        }
        Ok(point[2] <= self.base.height_mm())
    }
    /// Occupied support from retained source curves, never from a display mesh.
    pub fn support(&self, direction: [f64; 3]) -> f64 {
        let base = self.base.support(direction);
        if self.kind == Kind::Recessed {
            return base;
        }
        self.motifs
            .iter()
            .map(|m| m.support(direction) + direction[2] * self.base.height_mm())
            .fold(base, f64::max)
    }
    /// Exact Green properties of a finite horizontal cap (area and centroid).
    pub fn cap_properties(&self, index: usize) -> Option<(f64, [f64; 3])> {
        let face = self.boundary.faces.get(index)?;
        let Surface::Plane { z_mm, upward } = face.surface else {
            return None;
        };
        let first = self.boundary.edges[face.loops[0][0].edge].curve.point(0.);
        let origin = [first[0], first[1]];
        let mut sums = [crate::revolution::Sum::default(); 3];
        for edge in face.loops.iter().flatten() {
            let sign = if edge.forward == upward { 1. } else { -1. };
            for (sum, value) in sums
                .iter_mut()
                .zip(self.boundary.edges[edge.edge].curve.moments(origin))
            {
                sum.add(sign * value);
            }
        }
        let [area, x, y] = sums.map(crate::revolution::Sum::value);
        Some((area, [origin[0] + x / area, origin[1] + y / area, z_mm]))
    }
    pub fn bounds(&self) -> [[f64; 3]; 2] {
        let mut bounds = self.base.bounds();
        bounds[1][2] = self.levels.iter().copied().fold(bounds[1][2], f64::max);
        bounds
    }
}
