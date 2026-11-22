//! Owned Rust geometry and typed application boundary.
#![deny(unsafe_code)]

mod mesh;
pub use mesh::{ManufacturingMesh, Mesh};
#[allow(dead_code, unused_variables)]
mod provider_data;
use provider_data as ffi;
type OcctHandle = ffi::UniquePtr<ffi::Shape>;

pub struct Shape(OcctHandle);

/// Compatibility constant for callers checking the removed comparison provider.
pub const OCCT_LINKED: bool = false;

pub use multicad_rust_kernel as rust_kernel;
pub type Error = Box<dyn std::error::Error>;
mod distance;
pub use distance::{BRepSelection, DistanceWitness, MinimumDistance};
mod interference;
pub use interference::{SolidPairInspection, SolidPairRelation};
mod geometry;
pub use geometry::*;
mod affine;
pub use affine::AffineTransform;
mod advanced;
pub use advanced::TreatmentFace;
mod advanced_operations;
mod advanced_profiles;
pub mod thread_prototypes;
mod threads;
pub use threads::{
    AxialFlatHelicalThread, HelicalThread, KnuckleHelicalThread, RootRoundedHelicalThread,
    RoundedHelicalThread, ThreadFactoryStats, thread_factory_stats,
};
mod gears;
pub use gears::{GearHub, GearSweep, InvoluteGear};
mod bevel;
pub use bevel::SphericalBevelGear;
mod worm;
pub use worm::{GeneratedZaWormWheel, ZaWorm};
mod rack;
pub use rack::StraightRack;
mod edge_treatment;
pub use edge_treatment::{EdgeTreatment, SidedChamfer, SidedTreatmentEdge, TreatmentEdge};
mod drawing;
pub use drawing::*;
mod worker;
pub use worker::*;
mod curves;
mod elementary_curves;
mod rust_display;
mod rust_prism_topology;
mod sketch_queries;
pub use sketch_queries::NativeRegionQueries;
mod provenance;
mod revolve;
mod split;
mod taper;
mod to_plane;
pub use curves::*;
pub use revolve::RevolutionAxis;
pub use split::{PlaneSplit, SolidSplit};
/// Identity for exports stamped with the compiled geometry provider.
pub const fn source_identity() -> &'static str {
    "multicad-rust-kernel"
}

/// Compatibility hook for protocol callers. Rust geometry has no C++ output.
pub fn route_native_output_to_stderr() {
    ffi::route_native_output_to_stderr();
}

fn step_path(path: &std::path::Path) -> Result<&str, Error> {
    let value = path.to_str().ok_or("STEP path must be UTF-8")?;
    if value.contains('\0') {
        return Err("STEP path contains a NUL byte".into());
    }
    Ok(value)
}

fn positive(values: &[f64]) -> Result<(), Error> {
    if values
        .iter()
        .all(|v| v.is_finite() && *v > 1e-7 && *v < 1e7)
    {
        Ok(())
    } else {
        Err("Dimensions must be finite, above 1e-7 mm and below 1e7 mm".into())
    }
}
impl Shape {
    pub fn extrude(width: f64, depth: f64, height: f64) -> Result<Self, Error> {
        positive(&[width, depth, height])?;
        Ok(Self(ffi::extrude(width, depth, height)?))
    }
    pub fn cut(&self, x: f64, width: f64, depth: f64, height: f64) -> Result<Self, Error> {
        positive(&[width, depth, height])?;
        if !x.is_finite() || x.abs() >= 1e7 {
            return Err("Invalid cut placement".into());
        }
        Ok(Self(ffi::cut(&self.0, x, width, depth, height)?))
    }
    pub fn mesh(&self, deflection: f64) -> Result<Mesh, Error> {
        positive(&[deflection])?;
        Ok(ffi::mesh(&self.0, deflection)?.into())
    }
    pub fn volume(&self) -> Result<f64, Error> {
        Ok(ffi::volume(&self.0)?)
    }
    pub fn valid(&self) -> Result<bool, Error> {
        Ok(ffi::valid(&self.0)?)
    }
    pub fn write_step(&self, path: &std::path::Path) -> Result<(), Error> {
        ffi::write_step(&self.0, step_path(path)?)?;
        Ok(())
    }
    pub fn read_step(path: &std::path::Path) -> Result<Self, Error> {
        Ok(Self(ffi::read_step(step_path(path)?)?))
    }
}

mod planar_contours;
pub use planar_contours::SourcedPlanarContours;
mod transferred_cuts;
pub use transferred_cuts::TransferHoleAxis;
mod spline_edits;
pub use spline_edits::*;

mod sketch_references;
pub use sketch_references::*;
