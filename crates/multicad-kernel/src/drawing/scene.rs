//! Exact-curve hidden-line scenes.
//!
//! Every body in a scene is a [`Member`]: it supplies the actual boundary and
//! outline curves of its source geometry, the material intervals along a line
//! of sight, and its exact intersection with a plane facing the viewer. The
//! engine arranges those curves in the paper plane, splits them where the
//! projections meet, classifies each piece against the members' material and
//! applies section, depth, window and detail limits uniformly. No member may
//! answer from a display or manufacturing mesh.
//!
//! All member geometry is expressed in view space: `x` and `y` are the paper
//! axes, `z` is the depth toward the viewer, in millimetres from the drawing
//! frame origin.
use super::*;
use multicad_rust_kernel::curved_projection::Curve;
use std::sync::Arc;

pub(in crate::drawing) type P2 = [f64; 2];
pub(in crate::drawing) type P3 = [f64; 3];

/// World to view-space transform of one drawing frame.
#[derive(Clone, Copy, Debug)]
pub(in crate::drawing) struct View {
    pub origin: [f64; 3],
    /// Paper u, paper v and the direction toward the viewer, orthonormal.
    pub axes: [[f64; 3]; 3],
}

impl View {
    pub(in crate::drawing) fn new(frame: PlaneFrame) -> Self {
        Self {
            origin: frame.origin,
            axes: planar::view_axes(frame),
        }
    }

    pub(in crate::drawing) fn vector(&self, world: [f64; 3]) -> P3 {
        self.axes
            .map(|axis| axis.into_iter().zip(world).map(|(a, b)| a * b).sum())
    }

    pub(in crate::drawing) fn point(&self, world: [f64; 3]) -> P3 {
        self.vector(std::array::from_fn(|i| world[i] - self.origin[i]))
    }
}

/// A finite curve of the source geometry, evaluated exactly in view space.
pub(in crate::drawing) trait Path: Send + Sync {
    /// The point at `t` in `[0, 1]`. Must be continuous and injective.
    fn point(&self, t: f64) -> P3;

    /// Ascending parameters from 0 to 1 inclusive whose chords stay within
    /// `tolerance_mm` of the curve in all three view coordinates.
    fn flatten(&self, tolerance_mm: f64) -> Vec<f64>;

    /// The exact paper curve of `[from, to]` when it is a line or a conic
    /// arc. Other curves return `None` and publish as certified chords.
    fn paper(&self, range: [f64; 2]) -> Option<Curve>;
}

/// Why a candidate curve can appear in a drawing.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(in crate::drawing) enum Kind {
    /// A sharp boundary edge between two faces.
    Edge,
    /// A tangent-continuous seam; drawn only when tangent edges are requested.
    Smooth,
    /// A view-dependent outline: the surface normal is perpendicular to the
    /// line of sight along it.
    Outline,
}

#[derive(Clone)]
pub(in crate::drawing) struct Candidate {
    pub kind: Kind,
    pub path: Arc<dyn Path>,
}

/// One oriented boundary loop of a planar section, material on the left when
/// seen by the viewer: outer loops counterclockwise, holes clockwise.
pub(in crate::drawing) type Loop = Vec<Arc<dyn Path>>;

/// One solid of a scene, already expressed in a fixed view.
pub(in crate::drawing) trait Member: Send + Sync {
    /// Axis-aligned view-space bounds enclosing the solid.
    fn bounds(&self) -> [P3; 2];

    /// Every sharp edge, smooth seam and outline of the unclipped solid for
    /// this view. Outlines are complete: wherever the boundary's normal turns
    /// through perpendicular to the line of sight there is a candidate.
    fn candidates(&self, tolerance_mm: f64) -> Result<Vec<Candidate>, GeometryError>;

    /// Depth intervals `[z0, z1]`, ascending and disjoint, in which the line
    /// of sight through the paper point is strictly inside the solid. Roots
    /// are those of the actual source surfaces. A grazing contact contributes
    /// no interval.
    fn sight(&self, paper: P2) -> Result<Vec<[f64; 2]>, GeometryError>;

    /// Boundary loops of the solid's intersection with the plane `z = depth`.
    /// Empty when the plane misses the solid or only touches it.
    fn section(&self, depth: f64, tolerance_mm: f64) -> Result<Vec<Loop>, GeometryError>;
}

// Shared by the family members' refusals.
#[allow(dead_code)]
pub(in crate::drawing) fn unsupported(message: &'static str) -> GeometryError {
    err(GeometryErrorKind::Unsupported, message)
}

mod arrange;
mod bore;
mod engine;
mod families;
mod paths;
mod prism;
mod profile;
mod publish;
pub(in crate::drawing) use engine::{SceneMember, project};
pub(in crate::drawing) use families::project_solids;
