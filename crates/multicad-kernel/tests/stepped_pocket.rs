//! A box pocket entering one cap of a bored normal sweep keeps the bore
//! below its floor, and is refused when it reaches a profile edge or cuts
//! part of a hole.
use multicad_core::{document::EntityId, references::Origin};
use multicad_kernel::{
    BooleanOperation, CurvedEdge2d, CurvedProfile, GeometryPolicy, PlaneFrame, RigidTransform,
    Solid,
};
use std::f64::consts::PI;

fn policy() -> GeometryPolicy {
    GeometryPolicy::default()
}

/// A disk of radius 10 about the origin, swept by `height` (negative sweeps
/// against the frame normal), with a through bore of radius 3.3.
fn bored_disk(height: f64) -> Solid {
    let owner = EntityId::new();
    let profile = CurvedProfile {
        frame: PlaneFrame::default(),
        outer: vec![CurvedEdge2d::Circle {
            center: [0., 0.],
            radius_mm: 10.,
        }],
        holes: vec![],
    };
    let origins = [Origin {
        source: owner,
        role: "curve/0".into(),
    }];
    let disk =
        Solid::rust_extrude_normal_profile_named(&profile, height, owner, &origins, policy())
            .unwrap();
    let [low, high] = [height.min(0.), height.max(0.)];
    let bore = Solid::cylinder(
        PlaneFrame {
            origin: [0., 0., low - 1.],
            ..Default::default()
        },
        3.3,
        high - low + 2.,
        policy(),
    )
    .unwrap();
    disk.boolean(&bore, BooleanOperation::Subtract, policy())
        .unwrap()
        .unwrap()
}

/// An axis-aligned box `size` with its minimum corner at `min`.
fn block(min: [f64; 3], size: [f64; 3]) -> Solid {
    Solid::box_solid(size, policy())
        .unwrap()
        .transformed(RigidTransform {
            translation_mm: min,
            ..Default::default()
        })
        .unwrap()
}

fn near(a: f64, b: f64) {
    assert!((a - b).abs() <= 1e-9 * b.abs().max(1.), "{a} != {b}");
}

/// The disk with its bore, less a square pocket of side 8.4 and `depth`.
fn expected(height: f64, depth: f64) -> f64 {
    let bore = PI * 3.3 * 3.3;
    (PI * 100. - bore) * height - (8.4 * 8.4 - bore) * depth
}

#[test]
fn a_square_pocket_around_a_bore_keeps_the_bore_below_its_floor() {
    for height in [9., -9.] {
        let stock = bored_disk(height);
        // The box enters the upper cap and ends 3.2 mm below it: z = 9 for
        // the positive sweep, z = 0 for the one swept against the normal.
        let top = height.max(0.);
        let tool = block([-4.2, -4.2, top - 3.2], [8.4, 8.4, 4.]);
        let cut = stock
            .boolean(&tool, BooleanOperation::Subtract, policy())
            .unwrap()
            .unwrap();
        assert!(cut.uses_rust_geometry());
        let bounds = cut.bounds().unwrap();
        near(bounds.max[0], 10.);
        near(bounds.min[2], height.min(0.));
        near(bounds.max[2], top);
        let volume = cut.volume_mm3().unwrap();
        near(volume, expected(9., 3.2));
        // The exchange file reads back as the same solid.
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("pocket.step");
        cut.write_step(&path).unwrap();
        let read = Solid::read_step(&path, policy()).unwrap();
        near(read.volume_mm3().unwrap(), volume);
        assert!(
            !cut.manufacturing_mesh_with_receipt(Default::default())
                .unwrap()
                .mesh
                .indices
                .is_empty()
        );
    }
}

#[test]
fn a_pocket_from_the_bottom_cap_and_one_clear_of_a_second_hole() {
    // From below: the box spans z = -1 .. 2.5 of a 9 mm sweep.
    let stock = bored_disk(9.);
    let cut = stock
        .boolean(
            &block([-4.2, -4.2, -1.], [8.4, 8.4, 3.5]),
            BooleanOperation::Subtract,
            policy(),
        )
        .unwrap()
        .unwrap();
    near(cut.volume_mm3().unwrap(), expected(9., 2.5));
    // A second bore away from the pocket stays a through hole.
    let second = Solid::cylinder(
        PlaneFrame {
            origin: [7., 0., -1.],
            ..Default::default()
        },
        1.,
        11.,
        policy(),
    )
    .unwrap();
    let drilled = stock
        .boolean(&second, BooleanOperation::Subtract, policy())
        .unwrap()
        .unwrap();
    let cut = drilled
        .boolean(
            &block([-4.2, -4.2, 5.8], [8.4, 8.4, 4.]),
            BooleanOperation::Subtract,
            policy(),
        )
        .unwrap()
        .unwrap();
    near(cut.volume_mm3().unwrap(), expected(9., 3.2) - PI * 9.);
}

#[test]
fn a_pocket_reaching_the_profile_or_cutting_part_of_a_hole_is_not_this_cut() {
    let stock = bored_disk(9.);
    let expected_whole = expected(9., 3.2);
    for tool in [
        // Its corners reach past the disk's edge.
        block([-7.5, -7.5, 5.8], [15., 15., 4.]),
        // Its side crosses the bore.
        block([-1., -4.2, 5.8], [8.4, 8.4, 4.]),
    ] {
        if let Ok(Some(cut)) = stock.boolean(&tool, BooleanOperation::Subtract, policy()) {
            // Another exact provider may handle it, never with this volume.
            assert!((cut.volume_mm3().unwrap() - expected_whole).abs() > 1e-3);
        }
    }
}
