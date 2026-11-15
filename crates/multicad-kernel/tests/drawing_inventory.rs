//! Rust drawing capability inventory across body families, view kinds and mixed scenes.
//!
//! Every combination is either a named supported regression, which must produce
//! well-formed vectors, or an explicit unresolved finding, which must refuse with
//! a typed error. Nothing may silently return an empty or malformed drawing.
use multicad_core::{
    document::EntityId,
    drawings::{VectorCurve, VectorRole},
    references::Origin,
};
use multicad_kernel::{
    BooleanOperation, DrawingHatch, DrawingProjection, DrawingSolid, GearSweep,
    GeneratedZaWormWheel, GeometryError, GeometryPolicy, HelicalThread, InvoluteGear, PlaneFrame,
    ProjectedVector, RigidTransform, Solid, SphericalBevelGear, ZaWorm, project_drawing,
    rust_kernel::rounded_thread::Spec as RoundedThreadSpec,
};
use std::f64::consts::PI;

fn spur(sweep: GearSweep, twist: f64) -> InvoluteGear {
    let a = 20_f64.to_radians();
    let pitch = 12.;
    InvoluteGear {
        hub: None,
        teeth: 24,
        pitch_radius_mm: pitch,
        base_radius_mm: pitch * a.cos(),
        tip_radius_mm: pitch + 1.,
        root_radius_mm: pitch - 1.25,
        root_fillet_radius_mm: 0.2,
        half_base_angle_rad: (PI / 2. - 0.05) / (2. * pitch) + a.tan() - a,
        face_width_mm: 4.,
        bore_radius_mm: 3.,
        twist_rad_per_mm: twist,
        sweep,
    }
}

fn ring() -> InvoluteGear {
    let a = 20_f64.to_radians();
    let pitch = 24.;
    InvoluteGear {
        hub: None,
        teeth: 48,
        pitch_radius_mm: pitch,
        base_radius_mm: pitch * a.cos(),
        tip_radius_mm: pitch - 1.,
        root_radius_mm: pitch + 1.25,
        root_fillet_radius_mm: 0.2,
        half_base_angle_rad: (PI / 2. - 0.05) / (2. * pitch) - a.tan() + a,
        face_width_mm: 2.,
        bore_radius_mm: 0.,
        twist_rad_per_mm: 0.,
        sweep: GearSweep::InternalSpur {
            outer_radius_mm: pitch + 4.,
        },
    }
}

fn bevel() -> SphericalBevelGear {
    let (module, teeth, mate) = (2., 20_u32, 30_u32);
    let (z, mate) = (f64::from(teeth), f64::from(mate));
    let alpha = 20_f64.to_radians();
    let delta = 1_f64.atan2(mate / z);
    let distance = module * z / (2. * delta.sin());
    let beta = (delta.sin() * alpha.cos()).asin();
    let inv = |g: f64| {
        let sigma = (g.cos() / beta.cos()).acos();
        sigma / beta.sin() - sigma.sin().atan2(beta.sin() * sigma.cos())
    };
    SphericalBevelGear {
        teeth,
        pitch_cone_angle_rad: delta,
        base_cone_angle_rad: beta,
        tip_cone_angle_rad: delta + (module / distance).atan(),
        root_cone_angle_rad: delta - (1.25 * module / distance).atan(),
        outer_cone_distance_mm: distance,
        face_width_mm: 8.,
        bore_radius_mm: 3.,
        half_base_tooth_angle_rad: (PI * module / 2. - 0.05) / (module * z) + inv(delta),
        root_fillet_radius_mm: module * 0.2,
    }
}

/// The measured outcome of an inventoried combination. Every combination
/// currently draws; an explicit refusal returns here as a variant carrying
/// its exact message only together with a reviewed finding.
#[derive(Clone, Copy, Debug)]
enum Expect {
    Drawn,
}
use Expect::Drawn;

/// Expected outcome of every inventoried combination: each must publish
/// well-formed vectors. Regenerate the table from
/// `print_rust_drawing_inventory` only together with a reviewed change.
const EXPECTED: &[(&str, &str, Expect)] = &[
    ("box", "axial", Drawn),
    ("box", "axial reversed", Drawn),
    ("box", "axial hidden", Drawn),
    ("box", "side", Drawn),
    ("box", "side hidden", Drawn),
    ("box", "oblique", Drawn),
    ("box", "oblique hidden", Drawn),
    ("box", "axial section", Drawn),
    ("box", "axial section-only", Drawn),
    ("box", "axial finite depth", Drawn),
    ("box", "side section", Drawn),
    ("box", "side section-only", Drawn),
    ("box", "side finite depth", Drawn),
    ("box", "side lateral window", Drawn),
    ("box", "oblique section", Drawn),
    ("box", "oblique section-only", Drawn),
    ("box", "axial detail", Drawn),
    ("box", "side detail", Drawn),
    ("bored box", "axial", Drawn),
    ("bored box", "axial reversed", Drawn),
    ("bored box", "axial hidden", Drawn),
    ("bored box", "side", Drawn),
    ("bored box", "side hidden", Drawn),
    ("bored box", "oblique", Drawn),
    ("bored box", "oblique hidden", Drawn),
    ("bored box", "axial section", Drawn),
    ("bored box", "axial section-only", Drawn),
    ("bored box", "axial finite depth", Drawn),
    ("bored box", "side section", Drawn),
    ("bored box", "side section-only", Drawn),
    ("bored box", "side finite depth", Drawn),
    ("bored box", "side lateral window", Drawn),
    ("bored box", "oblique section", Drawn),
    ("bored box", "oblique section-only", Drawn),
    ("bored box", "axial detail", Drawn),
    ("bored box", "side detail", Drawn),
    ("oblique bore", "axial", Drawn),
    ("oblique bore", "axial reversed", Drawn),
    ("oblique bore", "axial hidden", Drawn),
    ("oblique bore", "side", Drawn),
    ("oblique bore", "side hidden", Drawn),
    ("oblique bore", "oblique", Drawn),
    ("oblique bore", "oblique hidden", Drawn),
    ("oblique bore", "axial section", Drawn),
    ("oblique bore", "axial section-only", Drawn),
    ("oblique bore", "axial finite depth", Drawn),
    ("oblique bore", "side section", Drawn),
    ("oblique bore", "side section-only", Drawn),
    ("oblique bore", "side finite depth", Drawn),
    ("oblique bore", "side lateral window", Drawn),
    ("oblique bore", "oblique section", Drawn),
    ("oblique bore", "oblique section-only", Drawn),
    ("oblique bore", "axial detail", Drawn),
    ("oblique bore", "side detail", Drawn),
    ("cylinder", "axial", Drawn),
    ("cylinder", "axial reversed", Drawn),
    ("cylinder", "axial hidden", Drawn),
    ("cylinder", "side", Drawn),
    ("cylinder", "side hidden", Drawn),
    ("cylinder", "oblique", Drawn),
    ("cylinder", "oblique hidden", Drawn),
    ("cylinder", "axial section", Drawn),
    ("cylinder", "axial section-only", Drawn),
    ("cylinder", "axial finite depth", Drawn),
    ("cylinder", "side section", Drawn),
    ("cylinder", "side section-only", Drawn),
    ("cylinder", "side finite depth", Drawn),
    ("cylinder", "side lateral window", Drawn),
    ("cylinder", "oblique section", Drawn),
    ("cylinder", "oblique section-only", Drawn),
    ("cylinder", "axial detail", Drawn),
    ("cylinder", "side detail", Drawn),
    ("annulus", "axial", Drawn),
    ("annulus", "axial reversed", Drawn),
    ("annulus", "axial hidden", Drawn),
    ("annulus", "side", Drawn),
    ("annulus", "side hidden", Drawn),
    ("annulus", "oblique", Drawn),
    ("annulus", "oblique hidden", Drawn),
    ("annulus", "axial section", Drawn),
    ("annulus", "axial section-only", Drawn),
    ("annulus", "axial finite depth", Drawn),
    ("annulus", "side section", Drawn),
    ("annulus", "side section-only", Drawn),
    ("annulus", "side finite depth", Drawn),
    ("annulus", "side lateral window", Drawn),
    ("annulus", "oblique section", Drawn),
    ("annulus", "oblique section-only", Drawn),
    ("annulus", "axial detail", Drawn),
    ("annulus", "side detail", Drawn),
    ("sphere", "axial", Drawn),
    ("sphere", "axial reversed", Drawn),
    ("sphere", "axial hidden", Drawn),
    ("sphere", "side", Drawn),
    ("sphere", "side hidden", Drawn),
    ("sphere", "oblique", Drawn),
    ("sphere", "oblique hidden", Drawn),
    ("sphere", "axial section", Drawn),
    ("sphere", "axial section-only", Drawn),
    ("sphere", "axial finite depth", Drawn),
    ("sphere", "side section", Drawn),
    ("sphere", "side section-only", Drawn),
    ("sphere", "side finite depth", Drawn),
    ("sphere", "side lateral window", Drawn),
    ("sphere", "oblique section", Drawn),
    ("sphere", "oblique section-only", Drawn),
    ("sphere", "axial detail", Drawn),
    ("sphere", "side detail", Drawn),
    ("rounded slot", "axial", Drawn),
    ("rounded slot", "axial reversed", Drawn),
    ("rounded slot", "axial hidden", Drawn),
    ("rounded slot", "side", Drawn),
    ("rounded slot", "side hidden", Drawn),
    ("rounded slot", "oblique", Drawn),
    ("rounded slot", "oblique hidden", Drawn),
    ("rounded slot", "axial section", Drawn),
    ("rounded slot", "axial section-only", Drawn),
    ("rounded slot", "axial finite depth", Drawn),
    ("rounded slot", "side section", Drawn),
    ("rounded slot", "side section-only", Drawn),
    ("rounded slot", "side finite depth", Drawn),
    ("rounded slot", "side lateral window", Drawn),
    ("rounded slot", "oblique section", Drawn),
    ("rounded slot", "oblique section-only", Drawn),
    ("rounded slot", "axial detail", Drawn),
    ("rounded slot", "side detail", Drawn),
    ("spur gear", "axial", Drawn),
    ("spur gear", "axial reversed", Drawn),
    ("spur gear", "axial hidden", Drawn),
    ("spur gear", "side", Drawn),
    ("spur gear", "side hidden", Drawn),
    ("spur gear", "oblique", Drawn),
    ("spur gear", "oblique hidden", Drawn),
    ("spur gear", "axial section", Drawn),
    ("spur gear", "axial section-only", Drawn),
    ("spur gear", "axial finite depth", Drawn),
    ("spur gear", "side section", Drawn),
    ("spur gear", "side section-only", Drawn),
    ("spur gear", "side finite depth", Drawn),
    ("spur gear", "side lateral window", Drawn),
    ("spur gear", "oblique section", Drawn),
    ("spur gear", "oblique section-only", Drawn),
    ("spur gear", "axial detail", Drawn),
    ("spur gear", "side detail", Drawn),
    ("internal ring", "axial", Drawn),
    ("internal ring", "axial reversed", Drawn),
    ("internal ring", "axial hidden", Drawn),
    ("internal ring", "side", Drawn),
    ("internal ring", "side hidden", Drawn),
    ("internal ring", "oblique", Drawn),
    ("internal ring", "oblique hidden", Drawn),
    ("internal ring", "axial section", Drawn),
    ("internal ring", "axial section-only", Drawn),
    ("internal ring", "axial finite depth", Drawn),
    ("internal ring", "side section", Drawn),
    ("internal ring", "side section-only", Drawn),
    ("internal ring", "side finite depth", Drawn),
    ("internal ring", "side lateral window", Drawn),
    ("internal ring", "oblique section", Drawn),
    ("internal ring", "oblique section-only", Drawn),
    ("internal ring", "axial detail", Drawn),
    ("internal ring", "side detail", Drawn),
    ("helical gear", "axial", Drawn),
    ("helical gear", "axial reversed", Drawn),
    ("helical gear", "axial hidden", Drawn),
    ("helical gear", "side", Drawn),
    ("helical gear", "side hidden", Drawn),
    ("helical gear", "oblique", Drawn),
    ("helical gear", "oblique hidden", Drawn),
    ("helical gear", "axial section", Drawn),
    ("helical gear", "axial section-only", Drawn),
    ("helical gear", "axial finite depth", Drawn),
    ("helical gear", "side section", Drawn),
    ("helical gear", "side section-only", Drawn),
    ("helical gear", "side finite depth", Drawn),
    ("helical gear", "side lateral window", Drawn),
    ("helical gear", "oblique section", Drawn),
    ("helical gear", "oblique section-only", Drawn),
    ("helical gear", "axial detail", Drawn),
    ("helical gear", "side detail", Drawn),
    ("herringbone gear", "axial", Drawn),
    ("herringbone gear", "axial reversed", Drawn),
    ("herringbone gear", "axial hidden", Drawn),
    ("herringbone gear", "side", Drawn),
    ("herringbone gear", "side hidden", Drawn),
    ("herringbone gear", "oblique", Drawn),
    ("herringbone gear", "oblique hidden", Drawn),
    ("herringbone gear", "axial section", Drawn),
    ("herringbone gear", "axial section-only", Drawn),
    ("herringbone gear", "axial finite depth", Drawn),
    ("herringbone gear", "side section", Drawn),
    ("herringbone gear", "side section-only", Drawn),
    ("herringbone gear", "side finite depth", Drawn),
    ("herringbone gear", "side lateral window", Drawn),
    ("herringbone gear", "oblique section", Drawn),
    ("herringbone gear", "oblique section-only", Drawn),
    ("herringbone gear", "axial detail", Drawn),
    ("herringbone gear", "side detail", Drawn),
    ("double-helical gear", "axial", Drawn),
    ("double-helical gear", "axial reversed", Drawn),
    ("double-helical gear", "axial hidden", Drawn),
    ("double-helical gear", "side", Drawn),
    ("double-helical gear", "side hidden", Drawn),
    ("double-helical gear", "oblique", Drawn),
    ("double-helical gear", "oblique hidden", Drawn),
    ("double-helical gear", "axial section", Drawn),
    ("double-helical gear", "axial section-only", Drawn),
    ("double-helical gear", "axial finite depth", Drawn),
    ("double-helical gear", "side section", Drawn),
    ("double-helical gear", "side section-only", Drawn),
    ("double-helical gear", "side finite depth", Drawn),
    ("double-helical gear", "side lateral window", Drawn),
    ("double-helical gear", "oblique section", Drawn),
    ("double-helical gear", "oblique section-only", Drawn),
    ("double-helical gear", "axial detail", Drawn),
    ("double-helical gear", "side detail", Drawn),
    ("bevel gear", "axial", Drawn),
    ("bevel gear", "axial reversed", Drawn),
    ("bevel gear", "axial hidden", Drawn),
    ("bevel gear", "side", Drawn),
    ("bevel gear", "side hidden", Drawn),
    ("bevel gear", "oblique", Drawn),
    ("bevel gear", "oblique hidden", Drawn),
    ("bevel gear", "axial section", Drawn),
    ("bevel gear", "axial section-only", Drawn),
    ("bevel gear", "axial finite depth", Drawn),
    ("bevel gear", "side section", Drawn),
    ("bevel gear", "side section-only", Drawn),
    ("bevel gear", "side finite depth", Drawn),
    ("bevel gear", "side lateral window", Drawn),
    ("bevel gear", "oblique section", Drawn),
    ("bevel gear", "oblique section-only", Drawn),
    ("bevel gear", "axial detail", Drawn),
    ("bevel gear", "side detail", Drawn),
    ("thread end relief", "axial", Drawn),
    ("thread end relief", "axial reversed", Drawn),
    ("thread end relief", "axial hidden", Drawn),
    ("thread end relief", "side", Drawn),
    ("thread end relief", "side hidden", Drawn),
    ("thread end relief", "oblique", Drawn),
    ("thread end relief", "oblique hidden", Drawn),
    ("thread end relief", "axial section", Drawn),
    ("thread end relief", "axial section-only", Drawn),
    ("thread end relief", "axial finite depth", Drawn),
    ("thread end relief", "side section", Drawn),
    ("thread end relief", "side section-only", Drawn),
    ("thread end relief", "side finite depth", Drawn),
    ("thread end relief", "side lateral window", Drawn),
    ("thread end relief", "oblique section", Drawn),
    ("thread end relief", "oblique section-only", Drawn),
    ("thread end relief", "axial detail", Drawn),
    ("thread end relief", "side detail", Drawn),
    ("rounded thread", "axial", Drawn),
    ("rounded thread", "axial reversed", Drawn),
    ("rounded thread", "axial hidden", Drawn),
    ("rounded thread", "side", Drawn),
    ("rounded thread", "side hidden", Drawn),
    ("rounded thread", "oblique", Drawn),
    ("rounded thread", "oblique hidden", Drawn),
    ("rounded thread", "axial section", Drawn),
    ("rounded thread", "axial section-only", Drawn),
    ("rounded thread", "axial finite depth", Drawn),
    ("rounded thread", "side section", Drawn),
    ("rounded thread", "side section-only", Drawn),
    ("rounded thread", "side finite depth", Drawn),
    ("rounded thread", "side lateral window", Drawn),
    ("rounded thread", "oblique section", Drawn),
    ("rounded thread", "oblique section-only", Drawn),
    ("rounded thread", "axial detail", Drawn),
    ("rounded thread", "side detail", Drawn),
    ("za worm", "axial", Drawn),
    ("za worm", "axial reversed", Drawn),
    ("za worm", "axial hidden", Drawn),
    ("za worm", "side", Drawn),
    ("za worm", "side hidden", Drawn),
    ("za worm", "oblique", Drawn),
    ("za worm", "oblique hidden", Drawn),
    ("za worm", "axial section", Drawn),
    ("za worm", "axial section-only", Drawn),
    ("za worm", "axial finite depth", Drawn),
    ("za worm", "side section", Drawn),
    ("za worm", "side section-only", Drawn),
    ("za worm", "side finite depth", Drawn),
    ("za worm", "side lateral window", Drawn),
    ("za worm", "oblique section", Drawn),
    ("za worm", "oblique section-only", Drawn),
    ("za worm", "axial detail", Drawn),
    ("za worm", "side detail", Drawn),
    ("worm wheel", "axial", Drawn),
    ("worm wheel", "axial reversed", Drawn),
    ("worm wheel", "axial hidden", Drawn),
    ("worm wheel", "side", Drawn),
    ("worm wheel", "side hidden", Drawn),
    ("worm wheel", "oblique", Drawn),
    ("worm wheel", "oblique hidden", Drawn),
    ("worm wheel", "axial section", Drawn),
    ("worm wheel", "axial section-only", Drawn),
    ("worm wheel", "axial finite depth", Drawn),
    ("worm wheel", "side section", Drawn),
    ("worm wheel", "side section-only", Drawn),
    ("worm wheel", "side finite depth", Drawn),
    ("worm wheel", "side lateral window", Drawn),
    ("worm wheel", "oblique section", Drawn),
    ("worm wheel", "oblique section-only", Drawn),
    ("worm wheel", "axial detail", Drawn),
    ("worm wheel", "side detail", Drawn),
    ("box + cylinder overlapping", "axial", Drawn),
    ("box + cylinder overlapping", "axial reversed", Drawn),
    ("box + cylinder overlapping", "axial hidden", Drawn),
    ("box + cylinder overlapping", "side", Drawn),
    ("box + cylinder overlapping", "side hidden", Drawn),
    ("box + cylinder overlapping", "oblique", Drawn),
    ("box + cylinder overlapping", "oblique hidden", Drawn),
    ("box + cylinder overlapping", "axial section", Drawn),
    ("box + cylinder overlapping", "axial section-only", Drawn),
    ("box + cylinder overlapping", "axial finite depth", Drawn),
    ("box + cylinder overlapping", "side section", Drawn),
    ("box + cylinder overlapping", "side section-only", Drawn),
    ("box + cylinder overlapping", "side finite depth", Drawn),
    ("box + cylinder overlapping", "side lateral window", Drawn),
    ("box + cylinder overlapping", "oblique section", Drawn),
    ("box + cylinder overlapping", "oblique section-only", Drawn),
    ("box + cylinder overlapping", "axial detail", Drawn),
    ("box + cylinder overlapping", "side detail", Drawn),
    ("box + cylinder separated", "axial", Drawn),
    ("box + cylinder separated", "axial reversed", Drawn),
    ("box + cylinder separated", "axial hidden", Drawn),
    ("box + cylinder separated", "side", Drawn),
    ("box + cylinder separated", "side hidden", Drawn),
    ("box + cylinder separated", "oblique", Drawn),
    ("box + cylinder separated", "oblique hidden", Drawn),
    ("box + cylinder separated", "axial section", Drawn),
    ("box + cylinder separated", "axial section-only", Drawn),
    ("box + cylinder separated", "axial finite depth", Drawn),
    ("box + cylinder separated", "side section", Drawn),
    ("box + cylinder separated", "side section-only", Drawn),
    ("box + cylinder separated", "side finite depth", Drawn),
    ("box + cylinder separated", "side lateral window", Drawn),
    ("box + cylinder separated", "oblique section", Drawn),
    ("box + cylinder separated", "oblique section-only", Drawn),
    ("box + cylinder separated", "axial detail", Drawn),
    ("box + cylinder separated", "side detail", Drawn),
    ("box + sphere overlapping", "axial", Drawn),
    ("box + sphere overlapping", "axial reversed", Drawn),
    ("box + sphere overlapping", "axial hidden", Drawn),
    ("box + sphere overlapping", "side", Drawn),
    ("box + sphere overlapping", "side hidden", Drawn),
    ("box + sphere overlapping", "oblique", Drawn),
    ("box + sphere overlapping", "oblique hidden", Drawn),
    ("box + sphere overlapping", "axial section", Drawn),
    ("box + sphere overlapping", "axial section-only", Drawn),
    ("box + sphere overlapping", "axial finite depth", Drawn),
    ("box + sphere overlapping", "side section", Drawn),
    ("box + sphere overlapping", "side section-only", Drawn),
    ("box + sphere overlapping", "side finite depth", Drawn),
    ("box + sphere overlapping", "side lateral window", Drawn),
    ("box + sphere overlapping", "oblique section", Drawn),
    ("box + sphere overlapping", "oblique section-only", Drawn),
    ("box + sphere overlapping", "axial detail", Drawn),
    ("box + sphere overlapping", "side detail", Drawn),
    ("box + sphere separated", "axial", Drawn),
    ("box + sphere separated", "axial reversed", Drawn),
    ("box + sphere separated", "axial hidden", Drawn),
    ("box + sphere separated", "side", Drawn),
    ("box + sphere separated", "side hidden", Drawn),
    ("box + sphere separated", "oblique", Drawn),
    ("box + sphere separated", "oblique hidden", Drawn),
    ("box + sphere separated", "axial section", Drawn),
    ("box + sphere separated", "axial section-only", Drawn),
    ("box + sphere separated", "axial finite depth", Drawn),
    ("box + sphere separated", "side section", Drawn),
    ("box + sphere separated", "side section-only", Drawn),
    ("box + sphere separated", "side finite depth", Drawn),
    ("box + sphere separated", "side lateral window", Drawn),
    ("box + sphere separated", "oblique section", Drawn),
    ("box + sphere separated", "oblique section-only", Drawn),
    ("box + sphere separated", "axial detail", Drawn),
    ("box + sphere separated", "side detail", Drawn),
    ("box + spur gear overlapping", "axial", Drawn),
    ("box + spur gear overlapping", "axial reversed", Drawn),
    ("box + spur gear overlapping", "axial hidden", Drawn),
    ("box + spur gear overlapping", "side", Drawn),
    ("box + spur gear overlapping", "side hidden", Drawn),
    ("box + spur gear overlapping", "oblique", Drawn),
    ("box + spur gear overlapping", "oblique hidden", Drawn),
    ("box + spur gear overlapping", "axial section", Drawn),
    ("box + spur gear overlapping", "axial section-only", Drawn),
    ("box + spur gear overlapping", "axial finite depth", Drawn),
    ("box + spur gear overlapping", "side section", Drawn),
    ("box + spur gear overlapping", "side section-only", Drawn),
    ("box + spur gear overlapping", "side finite depth", Drawn),
    ("box + spur gear overlapping", "side lateral window", Drawn),
    ("box + spur gear overlapping", "oblique section", Drawn),
    ("box + spur gear overlapping", "oblique section-only", Drawn),
    ("box + spur gear overlapping", "axial detail", Drawn),
    ("box + spur gear overlapping", "side detail", Drawn),
    ("box + spur gear separated", "axial", Drawn),
    ("box + spur gear separated", "axial reversed", Drawn),
    ("box + spur gear separated", "axial hidden", Drawn),
    ("box + spur gear separated", "side", Drawn),
    ("box + spur gear separated", "side hidden", Drawn),
    ("box + spur gear separated", "oblique", Drawn),
    ("box + spur gear separated", "oblique hidden", Drawn),
    ("box + spur gear separated", "axial section", Drawn),
    ("box + spur gear separated", "axial section-only", Drawn),
    ("box + spur gear separated", "axial finite depth", Drawn),
    ("box + spur gear separated", "side section", Drawn),
    ("box + spur gear separated", "side section-only", Drawn),
    ("box + spur gear separated", "side finite depth", Drawn),
    ("box + spur gear separated", "side lateral window", Drawn),
    ("box + spur gear separated", "oblique section", Drawn),
    ("box + spur gear separated", "oblique section-only", Drawn),
    ("box + spur gear separated", "axial detail", Drawn),
    ("box + spur gear separated", "side detail", Drawn),
    ("cylinder + sphere overlapping", "axial", Drawn),
    ("cylinder + sphere overlapping", "axial reversed", Drawn),
    ("cylinder + sphere overlapping", "axial hidden", Drawn),
    ("cylinder + sphere overlapping", "side", Drawn),
    ("cylinder + sphere overlapping", "side hidden", Drawn),
    ("cylinder + sphere overlapping", "oblique", Drawn),
    ("cylinder + sphere overlapping", "oblique hidden", Drawn),
    ("cylinder + sphere overlapping", "axial section", Drawn),
    ("cylinder + sphere overlapping", "axial section-only", Drawn),
    ("cylinder + sphere overlapping", "axial finite depth", Drawn),
    ("cylinder + sphere overlapping", "side section", Drawn),
    ("cylinder + sphere overlapping", "side section-only", Drawn),
    ("cylinder + sphere overlapping", "side finite depth", Drawn),
    (
        "cylinder + sphere overlapping",
        "side lateral window",
        Drawn,
    ),
    ("cylinder + sphere overlapping", "oblique section", Drawn),
    (
        "cylinder + sphere overlapping",
        "oblique section-only",
        Drawn,
    ),
    ("cylinder + sphere overlapping", "axial detail", Drawn),
    ("cylinder + sphere overlapping", "side detail", Drawn),
    ("cylinder + sphere separated", "axial", Drawn),
    ("cylinder + sphere separated", "axial reversed", Drawn),
    ("cylinder + sphere separated", "axial hidden", Drawn),
    ("cylinder + sphere separated", "side", Drawn),
    ("cylinder + sphere separated", "side hidden", Drawn),
    ("cylinder + sphere separated", "oblique", Drawn),
    ("cylinder + sphere separated", "oblique hidden", Drawn),
    ("cylinder + sphere separated", "axial section", Drawn),
    ("cylinder + sphere separated", "axial section-only", Drawn),
    ("cylinder + sphere separated", "axial finite depth", Drawn),
    ("cylinder + sphere separated", "side section", Drawn),
    ("cylinder + sphere separated", "side section-only", Drawn),
    ("cylinder + sphere separated", "side finite depth", Drawn),
    ("cylinder + sphere separated", "side lateral window", Drawn),
    ("cylinder + sphere separated", "oblique section", Drawn),
    ("cylinder + sphere separated", "oblique section-only", Drawn),
    ("cylinder + sphere separated", "axial detail", Drawn),
    ("cylinder + sphere separated", "side detail", Drawn),
    ("cylinder + annulus overlapping", "axial", Drawn),
    ("cylinder + annulus overlapping", "axial reversed", Drawn),
    ("cylinder + annulus overlapping", "axial hidden", Drawn),
    ("cylinder + annulus overlapping", "side", Drawn),
    ("cylinder + annulus overlapping", "side hidden", Drawn),
    ("cylinder + annulus overlapping", "oblique", Drawn),
    ("cylinder + annulus overlapping", "oblique hidden", Drawn),
    ("cylinder + annulus overlapping", "axial section", Drawn),
    (
        "cylinder + annulus overlapping",
        "axial section-only",
        Drawn,
    ),
    (
        "cylinder + annulus overlapping",
        "axial finite depth",
        Drawn,
    ),
    ("cylinder + annulus overlapping", "side section", Drawn),
    ("cylinder + annulus overlapping", "side section-only", Drawn),
    ("cylinder + annulus overlapping", "side finite depth", Drawn),
    (
        "cylinder + annulus overlapping",
        "side lateral window",
        Drawn,
    ),
    ("cylinder + annulus overlapping", "oblique section", Drawn),
    (
        "cylinder + annulus overlapping",
        "oblique section-only",
        Drawn,
    ),
    ("cylinder + annulus overlapping", "axial detail", Drawn),
    ("cylinder + annulus overlapping", "side detail", Drawn),
    ("cylinder + annulus separated", "axial", Drawn),
    ("cylinder + annulus separated", "axial reversed", Drawn),
    ("cylinder + annulus separated", "axial hidden", Drawn),
    ("cylinder + annulus separated", "side", Drawn),
    ("cylinder + annulus separated", "side hidden", Drawn),
    ("cylinder + annulus separated", "oblique", Drawn),
    ("cylinder + annulus separated", "oblique hidden", Drawn),
    ("cylinder + annulus separated", "axial section", Drawn),
    ("cylinder + annulus separated", "axial section-only", Drawn),
    ("cylinder + annulus separated", "axial finite depth", Drawn),
    ("cylinder + annulus separated", "side section", Drawn),
    ("cylinder + annulus separated", "side section-only", Drawn),
    ("cylinder + annulus separated", "side finite depth", Drawn),
    ("cylinder + annulus separated", "side lateral window", Drawn),
    ("cylinder + annulus separated", "oblique section", Drawn),
    (
        "cylinder + annulus separated",
        "oblique section-only",
        Drawn,
    ),
    ("cylinder + annulus separated", "axial detail", Drawn),
    ("cylinder + annulus separated", "side detail", Drawn),
    ("box + internal ring overlapping", "axial", Drawn),
    ("box + internal ring overlapping", "axial reversed", Drawn),
    ("box + internal ring overlapping", "axial hidden", Drawn),
    ("box + internal ring overlapping", "side", Drawn),
    ("box + internal ring overlapping", "side hidden", Drawn),
    ("box + internal ring overlapping", "oblique", Drawn),
    ("box + internal ring overlapping", "oblique hidden", Drawn),
    ("box + internal ring overlapping", "axial section", Drawn),
    (
        "box + internal ring overlapping",
        "axial section-only",
        Drawn,
    ),
    (
        "box + internal ring overlapping",
        "axial finite depth",
        Drawn,
    ),
    ("box + internal ring overlapping", "side section", Drawn),
    (
        "box + internal ring overlapping",
        "side section-only",
        Drawn,
    ),
    (
        "box + internal ring overlapping",
        "side finite depth",
        Drawn,
    ),
    (
        "box + internal ring overlapping",
        "side lateral window",
        Drawn,
    ),
    ("box + internal ring overlapping", "oblique section", Drawn),
    (
        "box + internal ring overlapping",
        "oblique section-only",
        Drawn,
    ),
    ("box + internal ring overlapping", "axial detail", Drawn),
    ("box + internal ring overlapping", "side detail", Drawn),
    ("box + internal ring separated", "axial", Drawn),
    ("box + internal ring separated", "axial reversed", Drawn),
    ("box + internal ring separated", "axial hidden", Drawn),
    ("box + internal ring separated", "side", Drawn),
    ("box + internal ring separated", "side hidden", Drawn),
    ("box + internal ring separated", "oblique", Drawn),
    ("box + internal ring separated", "oblique hidden", Drawn),
    ("box + internal ring separated", "axial section", Drawn),
    ("box + internal ring separated", "axial section-only", Drawn),
    ("box + internal ring separated", "axial finite depth", Drawn),
    ("box + internal ring separated", "side section", Drawn),
    ("box + internal ring separated", "side section-only", Drawn),
    ("box + internal ring separated", "side finite depth", Drawn),
    (
        "box + internal ring separated",
        "side lateral window",
        Drawn,
    ),
    ("box + internal ring separated", "oblique section", Drawn),
    (
        "box + internal ring separated",
        "oblique section-only",
        Drawn,
    ),
    ("box + internal ring separated", "axial detail", Drawn),
    ("box + internal ring separated", "side detail", Drawn),
    ("spur gear + spur gear overlapping", "axial", Drawn),
    ("spur gear + spur gear overlapping", "axial reversed", Drawn),
    ("spur gear + spur gear overlapping", "axial hidden", Drawn),
    ("spur gear + spur gear overlapping", "side", Drawn),
    ("spur gear + spur gear overlapping", "side hidden", Drawn),
    ("spur gear + spur gear overlapping", "oblique", Drawn),
    ("spur gear + spur gear overlapping", "oblique hidden", Drawn),
    ("spur gear + spur gear overlapping", "axial section", Drawn),
    (
        "spur gear + spur gear overlapping",
        "axial section-only",
        Drawn,
    ),
    (
        "spur gear + spur gear overlapping",
        "axial finite depth",
        Drawn,
    ),
    ("spur gear + spur gear overlapping", "side section", Drawn),
    (
        "spur gear + spur gear overlapping",
        "side section-only",
        Drawn,
    ),
    (
        "spur gear + spur gear overlapping",
        "side finite depth",
        Drawn,
    ),
    (
        "spur gear + spur gear overlapping",
        "side lateral window",
        Drawn,
    ),
    (
        "spur gear + spur gear overlapping",
        "oblique section",
        Drawn,
    ),
    (
        "spur gear + spur gear overlapping",
        "oblique section-only",
        Drawn,
    ),
    ("spur gear + spur gear overlapping", "axial detail", Drawn),
    ("spur gear + spur gear overlapping", "side detail", Drawn),
    ("spur gear + spur gear separated", "axial", Drawn),
    ("spur gear + spur gear separated", "axial reversed", Drawn),
    ("spur gear + spur gear separated", "axial hidden", Drawn),
    ("spur gear + spur gear separated", "side", Drawn),
    ("spur gear + spur gear separated", "side hidden", Drawn),
    ("spur gear + spur gear separated", "oblique", Drawn),
    ("spur gear + spur gear separated", "oblique hidden", Drawn),
    ("spur gear + spur gear separated", "axial section", Drawn),
    (
        "spur gear + spur gear separated",
        "axial section-only",
        Drawn,
    ),
    (
        "spur gear + spur gear separated",
        "axial finite depth",
        Drawn,
    ),
    ("spur gear + spur gear separated", "side section", Drawn),
    (
        "spur gear + spur gear separated",
        "side section-only",
        Drawn,
    ),
    (
        "spur gear + spur gear separated",
        "side finite depth",
        Drawn,
    ),
    (
        "spur gear + spur gear separated",
        "side lateral window",
        Drawn,
    ),
    ("spur gear + spur gear separated", "oblique section", Drawn),
    (
        "spur gear + spur gear separated",
        "oblique section-only",
        Drawn,
    ),
    ("spur gear + spur gear separated", "axial detail", Drawn),
    ("spur gear + spur gear separated", "side detail", Drawn),
    ("box + bevel gear overlapping", "axial", Drawn),
    ("box + bevel gear overlapping", "axial reversed", Drawn),
    ("box + bevel gear overlapping", "axial hidden", Drawn),
    ("box + bevel gear overlapping", "side", Drawn),
    ("box + bevel gear overlapping", "side hidden", Drawn),
    ("box + bevel gear overlapping", "oblique", Drawn),
    ("box + bevel gear overlapping", "oblique hidden", Drawn),
    ("box + bevel gear overlapping", "axial section", Drawn),
    ("box + bevel gear overlapping", "axial section-only", Drawn),
    ("box + bevel gear overlapping", "axial finite depth", Drawn),
    ("box + bevel gear overlapping", "side section", Drawn),
    ("box + bevel gear overlapping", "side section-only", Drawn),
    ("box + bevel gear overlapping", "side finite depth", Drawn),
    ("box + bevel gear overlapping", "side lateral window", Drawn),
    ("box + bevel gear overlapping", "oblique section", Drawn),
    (
        "box + bevel gear overlapping",
        "oblique section-only",
        Drawn,
    ),
    ("box + bevel gear overlapping", "axial detail", Drawn),
    ("box + bevel gear overlapping", "side detail", Drawn),
    ("box + bevel gear separated", "axial", Drawn),
    ("box + bevel gear separated", "axial reversed", Drawn),
    ("box + bevel gear separated", "axial hidden", Drawn),
    ("box + bevel gear separated", "side", Drawn),
    ("box + bevel gear separated", "side hidden", Drawn),
    ("box + bevel gear separated", "oblique", Drawn),
    ("box + bevel gear separated", "oblique hidden", Drawn),
    ("box + bevel gear separated", "axial section", Drawn),
    ("box + bevel gear separated", "axial section-only", Drawn),
    ("box + bevel gear separated", "axial finite depth", Drawn),
    ("box + bevel gear separated", "side section", Drawn),
    ("box + bevel gear separated", "side section-only", Drawn),
    ("box + bevel gear separated", "side finite depth", Drawn),
    ("box + bevel gear separated", "side lateral window", Drawn),
    ("box + bevel gear separated", "oblique section", Drawn),
    ("box + bevel gear separated", "oblique section-only", Drawn),
    ("box + bevel gear separated", "axial detail", Drawn),
    ("box + bevel gear separated", "side detail", Drawn),
    ("box + za worm overlapping", "axial", Drawn),
    ("box + za worm overlapping", "axial reversed", Drawn),
    ("box + za worm overlapping", "axial hidden", Drawn),
    ("box + za worm overlapping", "side", Drawn),
    ("box + za worm overlapping", "side hidden", Drawn),
    ("box + za worm overlapping", "oblique", Drawn),
    ("box + za worm overlapping", "oblique hidden", Drawn),
    ("box + za worm overlapping", "axial section", Drawn),
    ("box + za worm overlapping", "axial section-only", Drawn),
    ("box + za worm overlapping", "axial finite depth", Drawn),
    ("box + za worm overlapping", "side section", Drawn),
    ("box + za worm overlapping", "side section-only", Drawn),
    ("box + za worm overlapping", "side finite depth", Drawn),
    ("box + za worm overlapping", "side lateral window", Drawn),
    ("box + za worm overlapping", "oblique section", Drawn),
    ("box + za worm overlapping", "oblique section-only", Drawn),
    ("box + za worm overlapping", "axial detail", Drawn),
    ("box + za worm overlapping", "side detail", Drawn),
    ("box + za worm separated", "axial", Drawn),
    ("box + za worm separated", "axial reversed", Drawn),
    ("box + za worm separated", "axial hidden", Drawn),
    ("box + za worm separated", "side", Drawn),
    ("box + za worm separated", "side hidden", Drawn),
    ("box + za worm separated", "oblique", Drawn),
    ("box + za worm separated", "oblique hidden", Drawn),
    ("box + za worm separated", "axial section", Drawn),
    ("box + za worm separated", "axial section-only", Drawn),
    ("box + za worm separated", "axial finite depth", Drawn),
    ("box + za worm separated", "side section", Drawn),
    ("box + za worm separated", "side section-only", Drawn),
    ("box + za worm separated", "side finite depth", Drawn),
    ("box + za worm separated", "side lateral window", Drawn),
    ("box + za worm separated", "oblique section", Drawn),
    ("box + za worm separated", "oblique section-only", Drawn),
    ("box + za worm separated", "axial detail", Drawn),
    ("box + za worm separated", "side detail", Drawn),
    ("box + herringbone gear overlapping", "axial", Drawn),
    (
        "box + herringbone gear overlapping",
        "axial reversed",
        Drawn,
    ),
    ("box + herringbone gear overlapping", "axial hidden", Drawn),
    ("box + herringbone gear overlapping", "side", Drawn),
    ("box + herringbone gear overlapping", "side hidden", Drawn),
    ("box + herringbone gear overlapping", "oblique", Drawn),
    (
        "box + herringbone gear overlapping",
        "oblique hidden",
        Drawn,
    ),
    ("box + herringbone gear overlapping", "axial section", Drawn),
    (
        "box + herringbone gear overlapping",
        "axial section-only",
        Drawn,
    ),
    (
        "box + herringbone gear overlapping",
        "axial finite depth",
        Drawn,
    ),
    ("box + herringbone gear overlapping", "side section", Drawn),
    (
        "box + herringbone gear overlapping",
        "side section-only",
        Drawn,
    ),
    (
        "box + herringbone gear overlapping",
        "side finite depth",
        Drawn,
    ),
    (
        "box + herringbone gear overlapping",
        "side lateral window",
        Drawn,
    ),
    (
        "box + herringbone gear overlapping",
        "oblique section",
        Drawn,
    ),
    (
        "box + herringbone gear overlapping",
        "oblique section-only",
        Drawn,
    ),
    ("box + herringbone gear overlapping", "axial detail", Drawn),
    ("box + herringbone gear overlapping", "side detail", Drawn),
    ("box + herringbone gear separated", "axial", Drawn),
    ("box + herringbone gear separated", "axial reversed", Drawn),
    ("box + herringbone gear separated", "axial hidden", Drawn),
    ("box + herringbone gear separated", "side", Drawn),
    ("box + herringbone gear separated", "side hidden", Drawn),
    ("box + herringbone gear separated", "oblique", Drawn),
    ("box + herringbone gear separated", "oblique hidden", Drawn),
    ("box + herringbone gear separated", "axial section", Drawn),
    (
        "box + herringbone gear separated",
        "axial section-only",
        Drawn,
    ),
    (
        "box + herringbone gear separated",
        "axial finite depth",
        Drawn,
    ),
    ("box + herringbone gear separated", "side section", Drawn),
    (
        "box + herringbone gear separated",
        "side section-only",
        Drawn,
    ),
    (
        "box + herringbone gear separated",
        "side finite depth",
        Drawn,
    ),
    (
        "box + herringbone gear separated",
        "side lateral window",
        Drawn,
    ),
    ("box + herringbone gear separated", "oblique section", Drawn),
    (
        "box + herringbone gear separated",
        "oblique section-only",
        Drawn,
    ),
    ("box + herringbone gear separated", "axial detail", Drawn),
    ("box + herringbone gear separated", "side detail", Drawn),
    ("bored box + rounded slot overlapping", "axial", Drawn),
    (
        "bored box + rounded slot overlapping",
        "axial reversed",
        Drawn,
    ),
    (
        "bored box + rounded slot overlapping",
        "axial hidden",
        Drawn,
    ),
    ("bored box + rounded slot overlapping", "side", Drawn),
    ("bored box + rounded slot overlapping", "side hidden", Drawn),
    ("bored box + rounded slot overlapping", "oblique", Drawn),
    (
        "bored box + rounded slot overlapping",
        "oblique hidden",
        Drawn,
    ),
    (
        "bored box + rounded slot overlapping",
        "axial section",
        Drawn,
    ),
    (
        "bored box + rounded slot overlapping",
        "axial section-only",
        Drawn,
    ),
    (
        "bored box + rounded slot overlapping",
        "axial finite depth",
        Drawn,
    ),
    (
        "bored box + rounded slot overlapping",
        "side section",
        Drawn,
    ),
    (
        "bored box + rounded slot overlapping",
        "side section-only",
        Drawn,
    ),
    (
        "bored box + rounded slot overlapping",
        "side finite depth",
        Drawn,
    ),
    (
        "bored box + rounded slot overlapping",
        "side lateral window",
        Drawn,
    ),
    (
        "bored box + rounded slot overlapping",
        "oblique section",
        Drawn,
    ),
    (
        "bored box + rounded slot overlapping",
        "oblique section-only",
        Drawn,
    ),
    (
        "bored box + rounded slot overlapping",
        "axial detail",
        Drawn,
    ),
    ("bored box + rounded slot overlapping", "side detail", Drawn),
    ("bored box + rounded slot separated", "axial", Drawn),
    (
        "bored box + rounded slot separated",
        "axial reversed",
        Drawn,
    ),
    ("bored box + rounded slot separated", "axial hidden", Drawn),
    ("bored box + rounded slot separated", "side", Drawn),
    ("bored box + rounded slot separated", "side hidden", Drawn),
    ("bored box + rounded slot separated", "oblique", Drawn),
    (
        "bored box + rounded slot separated",
        "oblique hidden",
        Drawn,
    ),
    ("bored box + rounded slot separated", "axial section", Drawn),
    (
        "bored box + rounded slot separated",
        "axial section-only",
        Drawn,
    ),
    (
        "bored box + rounded slot separated",
        "axial finite depth",
        Drawn,
    ),
    ("bored box + rounded slot separated", "side section", Drawn),
    (
        "bored box + rounded slot separated",
        "side section-only",
        Drawn,
    ),
    (
        "bored box + rounded slot separated",
        "side finite depth",
        Drawn,
    ),
    (
        "bored box + rounded slot separated",
        "side lateral window",
        Drawn,
    ),
    (
        "bored box + rounded slot separated",
        "oblique section",
        Drawn,
    ),
    (
        "bored box + rounded slot separated",
        "oblique section-only",
        Drawn,
    ),
    ("bored box + rounded slot separated", "axial detail", Drawn),
    ("bored box + rounded slot separated", "side detail", Drawn),
];

const SCENES: &[(&str, &str)] = &[
    ("box", "cylinder"),
    ("box", "sphere"),
    ("box", "spur gear"),
    ("cylinder", "sphere"),
    ("cylinder", "annulus"),
    ("box", "internal ring"),
    ("spur gear", "spur gear"),
    ("box", "bevel gear"),
    ("box", "za worm"),
    ("box", "herringbone gear"),
    ("bored box", "rounded slot"),
];

pub fn families() -> Vec<(&'static str, Result<Solid, GeometryError>)> {
    let policy = GeometryPolicy::default();
    let id = EntityId::new;
    let boxed = Solid::box_solid([20., 12., 8.], policy).unwrap();
    let cylinder = Solid::cylinder(PlaneFrame::default(), 5., 8., policy).unwrap();
    let drill = Solid::cylinder(
        PlaneFrame {
            origin: [10., 6., -1.],
            ..Default::default()
        },
        2.,
        10.,
        policy,
    )
    .unwrap();
    let bored = boxed
        .boolean(&drill, BooleanOperation::Subtract, policy)
        .unwrap()
        .unwrap();
    let owner = id();
    let slot = Solid::rust_extrude_rounded_slot(
        PlaneFrame::default(),
        [[0., 0.], [6., 0.]],
        2.,
        5.,
        owner,
        std::array::from_fn(|edge| Origin {
            source: owner,
            role: format!("profile.{edge}"),
        }),
        policy,
    )
    .unwrap();
    // A 0.5 mm bore along the main diagonal of a 10 mm cube (CAD-354).
    let oblique = {
        let (sx, ty) = (0.5_f64.sqrt(), (1_f64 / 6.).sqrt());
        Solid::cylinder(
            PlaneFrame {
                origin: [1. / 3. - 3_f64.sqrt().recip(); 3],
                u: [sx, -sx, 0.],
                v: [ty, ty, -2. * ty],
            },
            0.5,
            20.,
            policy,
        )
        .and_then(|tool| {
            Solid::box_solid([10.; 3], policy)?.boolean(&tool, BooleanOperation::Subtract, policy)
        })
        .map(|cut| cut.expect("material remains"))
    };
    vec![
        ("box", Ok(boxed)),
        ("bored box", Ok(bored)),
        ("oblique bore", oblique),
        ("cylinder", Ok(cylinder)),
        (
            "annulus",
            Solid::rust_annulus(3., 6., 2.5, id(), PlaneFrame::default(), policy),
        ),
        (
            "sphere",
            Solid::rust_sphere_named(4., PlaneFrame::default(), id(), policy),
        ),
        ("rounded slot", Ok(slot)),
        (
            "spur gear",
            Solid::rust_spur_gear_named(spur(GearSweep::Spur, 0.), id(), policy),
        ),
        (
            "internal ring",
            Solid::rust_internal_ring_named(ring(), id(), policy),
        ),
        (
            "helical gear",
            Solid::rust_helical_gear_named(
                spur(GearSweep::Helical, 15_f64.to_radians().tan() / 12.),
                id(),
                policy,
            ),
        ),
        (
            "herringbone gear",
            Solid::rust_herringbone_named(
                spur(GearSweep::Chevron, 15_f64.to_radians().tan() / 12.),
                id(),
                policy,
            ),
        ),
        (
            "double-helical gear",
            Solid::rust_double_helical_named(
                spur(
                    GearSweep::GappedChevron { gap_mm: 1. },
                    15_f64.to_radians().tan() / 12.,
                ),
                id(),
                policy,
            ),
        ),
        (
            "bevel gear",
            Solid::rust_spherical_bevel_gear_named(bevel(), id(), policy),
        ),
        (
            "thread end relief",
            Solid::rust_helical_thread_end_relief(
                HelicalThread {
                    major_radius_mm: 3.,
                    pitch_mm: 1.,
                    length_mm: 6.6,
                    depth_mm: (1. - 0.125) / (2. * 30_f64.to_radians().tan())
                        - 1. / (4. * 3_f64.sqrt()),
                    crest_width_mm: 0.125,
                    lower_flank_deg: 30.,
                    upper_flank_deg: 30.,
                    radial_taper: 0.,
                    starts: 1,
                    left_hand: false,
                },
                Some(1. / (4. * 3_f64.sqrt())),
                0.6,
                0.1,
                id(),
                PlaneFrame::default(),
                policy,
            ),
        ),
        ("rounded thread", {
            let p = 1.27;
            let a = 27.5_f64.to_radians();
            let r = p * a.cos() / (12. * (1. - a.sin()));
            Solid::rust_rounded_thread(
                RoundedThreadSpec {
                    major_radius_mm: 3.175,
                    pitch_mm: p,
                    length_mm: 3.7 * p,
                    depth_mm: p / (3. * a.tan()),
                    crest_radius_mm: r,
                    root_radius_mm: r,
                    flank_angle_deg: 27.5,
                    starts: 1,
                    left_hand: false,
                },
                id(),
                PlaneFrame::default(),
                policy,
            )
        }),
        (
            "za worm",
            Solid::rust_za_worm_named(
                ZaWorm {
                    axial_module_mm: 1.,
                    threaded_length_mm: 20.,
                    bore_radius_mm: 0.,
                    axial_thinning_mm: 0.02,
                },
                id(),
                policy,
            ),
        ),
        (
            "worm wheel",
            Solid::rust_generated_za_worm_wheel_named(
                GeneratedZaWormWheel {
                    teeth: 40,
                    axial_module_mm: 2.,
                    face_width_mm: 4.,
                    bore_radius_mm: 6.,
                    pitch_thinning_mm: 0.1,
                },
                id(),
                policy,
            ),
        ),
    ]
}

fn frame(u: [f64; 3], v: [f64; 3]) -> PlaneFrame {
    PlaneFrame {
        origin: [0.; 3],
        u,
        v,
    }
}

pub fn views(center: [f64; 3]) -> Vec<(&'static str, DrawingProjection)> {
    let s = 0.5_f64.sqrt();
    let t = (1_f64 / 6.).sqrt();
    let axial = frame([1., 0., 0.], [0., 1., 0.]);
    let reversed = frame([1., 0., 0.], [0., -1., 0.]);
    let side = frame([1., 0., 0.], [0., 0., 1.]);
    let oblique = frame([s, -s, 0.], [t, t, -2. * t]);
    let base = DrawingProjection::default();
    let at = |frame: PlaneFrame| PlaneFrame {
        origin: center,
        ..frame
    };
    let section = |frame: PlaneFrame| DrawingProjection {
        frame: at(frame),
        section_offset_mm: Some(0.),
        ..base
    };
    vec![
        (
            "axial",
            DrawingProjection {
                frame: axial,
                ..base
            },
        ),
        (
            "axial reversed",
            DrawingProjection {
                frame: reversed,
                ..base
            },
        ),
        (
            "axial hidden",
            DrawingProjection {
                frame: axial,
                hidden_lines: true,
                ..base
            },
        ),
        (
            "side",
            DrawingProjection {
                frame: side,
                ..base
            },
        ),
        (
            "side hidden",
            DrawingProjection {
                frame: side,
                hidden_lines: true,
                ..base
            },
        ),
        (
            "oblique",
            DrawingProjection {
                frame: oblique,
                ..base
            },
        ),
        (
            "oblique hidden",
            DrawingProjection {
                frame: oblique,
                hidden_lines: true,
                ..base
            },
        ),
        ("axial section", section(axial)),
        (
            "axial section-only",
            DrawingProjection {
                section_only: true,
                ..section(axial)
            },
        ),
        (
            "axial finite depth",
            DrawingProjection {
                section_depth_mm: Some(1.),
                ..section(axial)
            },
        ),
        ("side section", section(side)),
        (
            "side section-only",
            DrawingProjection {
                section_only: true,
                ..section(side)
            },
        ),
        (
            "side finite depth",
            DrawingProjection {
                section_depth_mm: Some(1.),
                ..section(side)
            },
        ),
        (
            "side lateral window",
            DrawingProjection {
                section_u_mm: Some([-2., 2.]),
                section_v_mm: Some([-2., 2.]),
                ..section(side)
            },
        ),
        ("oblique section", section(oblique)),
        (
            "oblique section-only",
            DrawingProjection {
                section_only: true,
                ..section(oblique)
            },
        ),
        (
            "axial detail",
            DrawingProjection {
                frame: axial,
                detail: Some(([0., 0.], 3.)),
                ..base
            },
        ),
        (
            "side detail",
            DrawingProjection {
                frame: side,
                detail: Some(([0., 0.], 3.)),
                ..base
            },
        ),
    ]
}

fn hatch() -> Option<DrawingHatch> {
    Some(DrawingHatch {
        angle_rad: 0.6,
        spacing_mm: 0.5,
        line_width_mm: 0.18,
        pattern: multicad_core::materials::HatchPattern::Lines,
    })
}

pub fn center(solid: &Solid) -> [f64; 3] {
    let b = solid.bounds().unwrap();
    std::array::from_fn(|i| (b.min[i] + b.max[i]) / 2.)
}

pub fn outcome(result: &Result<Vec<ProjectedVector>, GeometryError>) -> String {
    match result {
        Ok(vectors) => {
            let count = |role| vectors.iter().filter(|v| v.role == role).count();
            format!(
                "ok v{} h{} s{} x{}",
                count(VectorRole::Visible),
                count(VectorRole::Hidden),
                count(VectorRole::Section),
                count(VectorRole::Hatch)
            )
        }
        Err(e) => format!("ERR {:?}: {}", e.kind, e.message),
    }
}

fn well_formed(vectors: &[ProjectedVector], sources: usize) {
    for vector in vectors {
        assert!(vector.source_index < sources);
        assert!(vector.preview_mm.len() >= 2);
        assert!(vector.preview_mm.iter().flatten().all(|v| v.is_finite()));
        if let VectorCurve::Line { a, b } = vector.curve {
            assert!(a.iter().chain(&b).all(|v| v.is_finite()));
        }
    }
}

#[test]
#[ignore = "inventory report; run explicitly with --ignored --nocapture"]
fn print_rust_drawing_inventory() {
    let families = families();
    for (name, built) in &families {
        let solid = match built {
            Ok(solid) => solid,
            Err(e) => {
                println!("INV|{name}|construction|ERR {:?}: {}", e.kind, e.message);
                continue;
            }
        };
        for (view, options) in views(center(solid)) {
            let result = project_drawing(
                &[DrawingSolid {
                    solid,
                    hatch: hatch(),
                }],
                options,
            );
            if let Ok(vectors) = &result {
                well_formed(vectors, 1);
            }
            println!("INV|{name}|{view}|{}", outcome(&result));
        }
    }
    let policy = GeometryPolicy::default();
    let place = |solid: &Solid, x: f64| {
        solid
            .transformed(RigidTransform {
                rotation: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
                translation_mm: [x, 0., 20.],
            })
            .unwrap()
    };
    let by_name = |n: &str| {
        families
            .iter()
            .find(|(name, _)| *name == n)
            .unwrap()
            .1
            .as_ref()
            .ok()
    };
    for (a, b) in SCENES {
        let (Some(first), Some(second)) = (by_name(a), by_name(b)) else {
            println!("INV|{a} + {b}|construction|ERR missing member");
            continue;
        };
        let _ = policy;
        for (layout, offset) in [("overlapping", 2.), ("separated", 80.)] {
            let near = place(second, offset);
            let scene = [
                DrawingSolid {
                    solid: first,
                    hatch: hatch(),
                },
                DrawingSolid {
                    solid: &near,
                    hatch: hatch(),
                },
            ];
            for (view, options) in views([2., 0., 10.]) {
                let result = project_drawing(&scene, options);
                if let Ok(vectors) = &result {
                    well_formed(vectors, 2);
                }
                println!("INV|{a} + {b} {layout}|{view}|{}", outcome(&result));
            }
        }
    }
}

/// Every inventoried family, view and two-body scene draws well-formed
/// vectors.
#[test]
fn every_family_view_and_scene_is_a_named_regression_or_finding() {
    let mut seen = 0;
    let mut check = |name: &str,
                     view: &str,
                     result: Result<Vec<ProjectedVector>, GeometryError>,
                     sources: usize| {
        let (_, _, expect) = EXPECTED
            .iter()
            .find(|(family, case, _)| *family == name && *case == view)
            .unwrap_or_else(|| panic!("{name} | {view} is not inventoried"));
        seen += 1;
        match (expect, &result) {
            (Drawn, Ok(vectors)) => well_formed(vectors, sources),
            _ => panic!(
                "{name} | {view}: expected {expect:?}, got {}",
                outcome(&result)
            ),
        }
    };
    let families = families();
    for (name, built) in &families {
        let solid = built.as_ref().unwrap_or_else(|e| panic!("{name}: {e:?}"));
        for (view, options) in views(center(solid)) {
            let result = project_drawing(
                &[DrawingSolid {
                    solid,
                    hatch: hatch(),
                }],
                options,
            );
            check(name, view, result, 1);
        }
    }
    let place = |solid: &Solid, x: f64| {
        solid
            .transformed(RigidTransform {
                rotation: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
                translation_mm: [x, 0., 20.],
            })
            .unwrap()
    };
    let by_name = |n: &str| {
        families
            .iter()
            .find(|(name, _)| *name == n)
            .unwrap()
            .1
            .as_ref()
            .unwrap()
    };
    for (a, b) in SCENES {
        for (layout, offset) in [("overlapping", 2.), ("separated", 80.)] {
            let near = place(by_name(b), offset);
            let scene = [
                DrawingSolid {
                    solid: by_name(a),
                    hatch: hatch(),
                },
                DrawingSolid {
                    solid: &near,
                    hatch: hatch(),
                },
            ];
            for (view, options) in views([2., 0., 10.]) {
                let name = format!("{a} + {b} {layout}");
                check(&name, view, project_drawing(&scene, options), 2);
            }
        }
    }
    assert_eq!(seen, EXPECTED.len(), "every inventoried combination ran");
}

/// Negative control: every inventoried family and two-body scene refuses an
/// invalid request with a typed input error and an exhausted budget with a
/// typed limit error, never with partial or empty output.
#[test]
fn invalid_requests_and_exhausted_budgets_refuse_atomically_for_every_family_and_scene() {
    let families = families();
    let by_name = |n: &str| {
        families
            .iter()
            .find(|(name, _)| *name == n)
            .unwrap()
            .1
            .as_ref()
            .unwrap()
    };
    let place = |solid: &Solid, x: f64| {
        solid
            .transformed(RigidTransform {
                rotation: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
                translation_mm: [x, 0., 20.],
            })
            .unwrap()
    };
    let side = frame([1., 0., 0.], [0., 0., 1.]);
    let base = DrawingProjection {
        frame: side,
        hidden_lines: true,
        ..DrawingProjection::default()
    };
    let invalid: [(&str, DrawingProjection); 9] = [
        (
            "parallel frame axes",
            DrawingProjection {
                frame: frame([1., 0., 0.], [2., 0., 0.]),
                ..base
            },
        ),
        (
            "zero frame axis",
            DrawingProjection {
                frame: frame([0., 0., 0.], [0., 0., 1.]),
                ..base
            },
        ),
        (
            "non-finite frame origin",
            DrawingProjection {
                frame: PlaneFrame {
                    origin: [f64::NAN, 0., 0.],
                    ..side
                },
                ..base
            },
        ),
        (
            "zero chord tolerance",
            DrawingProjection {
                chord_tolerance_mm: 0.,
                ..base
            },
        ),
        (
            "empty curve budget",
            DrawingProjection {
                maximum_curves: 0,
                ..base
            },
        ),
        (
            "empty point budget",
            DrawingProjection {
                maximum_points: 0,
                ..base
            },
        ),
        (
            "section-only without a section plane",
            DrawingProjection {
                section_only: true,
                ..base
            },
        ),
        (
            "zero detail radius",
            DrawingProjection {
                detail: Some(([0., 0.], 0.)),
                ..base
            },
        ),
        (
            "inverted lateral window",
            DrawingProjection {
                section_offset_mm: Some(0.),
                section_u_mm: Some([2., -2.]),
                ..base
            },
        ),
    ];
    let mut checked = 0;
    let mut refuse = |name: &str, scene: &[DrawingSolid<'_>]| {
        for (label, options) in invalid {
            match project_drawing(scene, options) {
                Err(error) if error.kind == multicad_kernel::GeometryErrorKind::Input => {}
                other => panic!(
                    "{name} | {label}: expected an input error, got {}",
                    outcome(&other)
                ),
            }
        }
        // One curve or one point fewer than the full drawing needs refuses
        // the whole drawing.
        let full =
            project_drawing(scene, base).unwrap_or_else(|e| panic!("{name} | full drawing: {e:?}"));
        let points: usize = full.iter().map(|v| v.preview_mm.len()).sum();
        let mut budgets = vec![(
            "one point short",
            DrawingProjection {
                maximum_points: u32::try_from(points - 1).unwrap(),
                ..base
            },
        )];
        if full.len() > 1 {
            budgets.push((
                "one curve short",
                DrawingProjection {
                    maximum_curves: u32::try_from(full.len() - 1).unwrap(),
                    ..base
                },
            ));
        }
        for (label, options) in budgets {
            match project_drawing(scene, options) {
                Err(error) if error.kind == multicad_kernel::GeometryErrorKind::Limit => {}
                other => panic!(
                    "{name} | {label}: expected a limit error, got {}",
                    outcome(&other)
                ),
            }
        }
        checked += 1;
    };
    for (name, built) in &families {
        let solid = built.as_ref().unwrap_or_else(|e| panic!("{name}: {e:?}"));
        refuse(
            name,
            &[DrawingSolid {
                solid,
                hatch: hatch(),
            }],
        );
    }
    for (a, b) in SCENES {
        let near = place(by_name(b), 2.);
        refuse(
            &format!("{a} + {b} overlapping"),
            &[
                DrawingSolid {
                    solid: by_name(a),
                    hatch: hatch(),
                },
                DrawingSolid {
                    solid: &near,
                    hatch: hatch(),
                },
            ],
        );
    }
    assert_eq!(checked, families.len() + SCENES.len());
    assert!(
        project_drawing(&[], base).is_err(),
        "an empty scene refuses"
    );
}
