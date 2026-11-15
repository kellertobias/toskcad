//! The retained boundary provider against analytic oracles and tampered
//! exchange text.
use super::*;
use crate::manufacturing::certify;
use crate::step_read::brep::{read_boundaries, read_solids};
use crate::step_read::{ReadError, read_analytic_compound};
use std::f64::consts::PI;

const TOLERANCE: f64 = 1e-7;

macro_rules! fixture {
    ($name:literal) => {
        include_str!(concat!(
            "../../../multicad-app/tests/fixtures/exchange/",
            $name
        ))
    };
}

const FRUSTUM: &str = fixture!("corpus/occt-cone-frustum.step");
const SPHERE: &str = fixture!("corpus/occt-sphere.step");
const TORUS: &str = fixture!("corpus/occt-torus.step");
const CYLINDER: &str = fixture!("corpus/occt-cylinder.step");
const VOID: &str = fixture!("corpus/occt-box-void.step");
const VASE: &str = fixture!("corpus/occt-revolved-vase.step");
const SPLINE: &str = fixture!("corpus/occt-spline-prism.step");
const DRILLED: &str = fixture!("freecad-drilled-box.step");
const STANDOFF: &str =
    fixture!("corpus/kicad-library/Mounting_Wuerth_WA-SMSE-ExternalM3_H5mm_9771050360.step");
const BUMPS: &str = fixture!("corpus/kicad-library/WLP-4_0.83x0.83mm_P0.4mm.step");

fn solids(text: &str) -> Result<Vec<Solid>> {
    read_solids(text, TOLERANCE, None)
}

fn solid(text: &str) -> Solid {
    let mut all = solids(text).unwrap();
    assert_eq!(all.len(), 1);
    all.remove(0)
}

fn options(deflection_mm: f64) -> mesh::Options {
    mesh::Options {
        deflection_mm,
        angle_rad: 0.3,
        triangle_limit: 60_000,
    }
}

fn close(actual: f64, expected: f64, relative: f64) -> bool {
    (actual - expected).abs() <= relative * expected.abs().max(1.)
}

/// Tessellate, certify the mesh closed and oriented, and bound its volume
/// defect by the deflection swept over the boundary area.
fn certified(solid: &Solid, deflection: f64) -> mesh::Mesh {
    let mesh = mesh::tessellate(solid, options(deflection), None).unwrap();
    let proof = certify(&mesh.display()).unwrap();
    let area: f64 = solid.faces().iter().map(|face| face.area_mm2).sum();
    assert!(mesh.maximum_deviation_mm <= deflection);
    assert!((proof.signed_volume_mm3 - solid.volume_mm3()).abs() <= deflection * area);
    mesh
}

fn sorted_areas(solid: &Solid) -> Vec<f64> {
    let mut areas: Vec<f64> = solid.faces().iter().map(|face| face.area_mm2).collect();
    areas.sort_by(f64::total_cmp);
    areas
}

#[test]
fn frustum_matches_its_closed_forms_and_meshes_within_deflection() {
    let solid = solid(FRUSTUM);
    let exact = PI * 10. / 3. * (64. + 32. + 16.);
    assert!(close(solid.volume_mm3(), exact, 1e-12));
    assert!(solid.volume_error_mm3() < 1e-4);
    let areas = sorted_areas(&solid);
    let expected = [16. * PI, 64. * PI, PI * 12. * 116f64.sqrt()];
    assert!(
        areas.iter().zip(expected).all(|(&a, e)| close(a, e, 1e-12)),
        "{areas:?}"
    );
    let [low, high] = solid.bounds_mm();
    assert!(
        low.iter()
            .zip([-8., -8., 0.])
            .all(|(a, b)| (a - b).abs() < 1e-9),
        "{low:?}"
    );
    assert!(
        high.iter()
            .zip([8., 8., 10.])
            .all(|(a, b)| (a - b).abs() < 1e-9),
        "{high:?}"
    );
    for deflection in [0.2, 0.05, 0.01] {
        certified(&solid, deflection);
    }
}

#[test]
fn complete_sphere_and_torus_use_their_whole_parameter_rectangle() {
    let sphere = solid(SPHERE);
    assert!(close(sphere.volume_mm3(), 4. / 3. * PI * 125., 1e-12));
    assert!(close(sphere.faces()[0].area_mm2, 100. * PI, 1e-12));
    assert!(
        sphere.faces()[0]
            .centroid_mm
            .iter()
            .all(|x| x.abs() < 1e-12)
    );
    let [low, high] = sphere.bounds_mm();
    assert!(low.iter().chain(&high).all(|x| (x.abs() - 5.).abs() < 1e-9));
    let mesh = certified(&sphere, 0.05);
    assert!(mesh.positions.iter().all(|p| (norm(*p) - 5.).abs() < 1e-9));

    let torus = solid(TORUS);
    assert!(close(torus.volume_mm3(), 2. * PI * PI * 12. * 9., 1e-12));
    assert!(close(torus.faces()[0].area_mm2, 4. * PI * PI * 36., 1e-12));
    let [low, high] = torus.bounds_mm();
    assert!(
        (low[0] + 15.).abs() < 1e-9 && (high[2] - 3.).abs() < 1e-9,
        "{low:?} {high:?}"
    );
    certified(&torus, 0.05);
}

#[test]
fn cylinder_void_and_drilled_solids_match_their_closed_forms() {
    let cylinder = solid(CYLINDER);
    assert!(close(cylinder.volume_mm3(), PI * 36. * 15., 1e-12));
    certified(&cylinder, 0.05);
    let hollow = solid(VOID);
    assert!(close(hollow.volume_mm3(), 5760., 1e-12));
    assert_eq!(hollow.boundary().shells.len(), 2);
    assert_eq!(certified(&hollow, 0.05).triangles.len(), 24);
    let drilled = solid(DRILLED);
    assert!(close(drilled.volume_mm3(), 4699.469035085, 1e-9));
    certified(&drilled, 0.05);
}

#[test]
fn revolved_swept_and_blended_solids_mesh_watertight() {
    for (text, faces) in [(VASE, 5), (SPLINE, 4), (STANDOFF, 11), (BUMPS, 12)] {
        let solid = solid(text);
        assert_eq!(solid.faces().len(), faces);
        assert!(solid.volume_error_mm3() <= 1e-6 * solid.volume_mm3().max(1.));
        certified(&solid, 0.05);
        certified(&solid, 0.01);
    }
}

#[test]
fn a_rigid_motion_keeps_volume_and_moves_bounds_and_centroids() {
    let solid = solid(FRUSTUM);
    // Quarter turn about x, then a shift.
    let linear = |p: Point| [p[0], -p[2], p[1]];
    let map = |p: Point| add(linear(p), [100., 20., -3.]);
    let moved = solid.moved(&map, &linear).unwrap();
    assert_eq!(moved.volume_mm3(), solid.volume_mm3());
    let [low, high] = moved.bounds_mm();
    assert!(distance(low, [92., 10., -11.]) < 1e-9 && distance(high, [108., 20., 5.]) < 1e-9);
    assert!(
        distance(
            moved.faces()[0].centroid_mm,
            map(solid.faces()[0].centroid_mm)
        ) < 1e-12
    );
    assert_eq!(moved.tolerance(), TOLERANCE);
    certified(&moved, 0.05);
    let rebuilt = Solid::new(moved.boundary().clone(), TOLERANCE, None).unwrap();
    assert!(close(rebuilt.volume_mm3(), solid.volume_mm3(), 1e-12));
}

#[test]
fn interrupts_and_invalid_options_are_refused() {
    let boundary = read_boundaries(FRUSTUM, TOLERANCE).unwrap().remove(0);
    let stop = || true;
    assert_eq!(
        Solid::new(boundary.clone(), TOLERANCE, Some(&stop)).unwrap_err(),
        Error::Cancelled
    );
    assert_eq!(
        Solid::new(boundary.clone(), 0., None).unwrap_err(),
        Error::Invalid
    );
    let solid = Solid::new(boundary, TOLERANCE, None).unwrap();
    assert_eq!(
        mesh::tessellate(&solid, options(0.05), Some(&stop)).unwrap_err(),
        Error::Cancelled
    );
    assert_eq!(
        mesh::tessellate(&solid, options(0.), None).unwrap_err(),
        Error::Invalid
    );
    let tight = mesh::Options {
        triangle_limit: 10,
        ..options(0.05)
    };
    assert_eq!(
        mesh::tessellate(&solid, tight, None).unwrap_err(),
        Error::Limit
    );
}

/// Replace the first occurrence, which must exist.
fn tampered(text: &str, from: &str, to: &str) -> String {
    assert!(text.contains(from), "{from}");
    text.replacen(from, to, 1)
}

#[test]
fn flipped_faces_missing_faces_and_moved_vertices_are_invalid() {
    // One face sense flipped: its edges now run the same way twice.
    let flipped = tampered(
        FRUSTUM,
        "#54 = ADVANCED_FACE('',(#55),#43,.T.)",
        "#54 = ADVANCED_FACE('',(#55),#43,.F.)",
    );
    assert_eq!(solids(&flipped).unwrap_err(), Error::Invalid);
    // A closed shell without one of its faces is open.
    let open = tampered(
        FRUSTUM,
        "CLOSED_SHELL('',(#17,#54,#109))",
        "CLOSED_SHELL('',(#17,#54))",
    );
    assert_eq!(solids(&open).unwrap_err(), Error::Unsupported);
    // A vertex moved by a micrometre leaves its edges; a hundredth of the
    // tolerance stays within it.
    let from = "#23 = CARTESIAN_POINT('',(8.,0.,0.))";
    let moved = tampered(FRUSTUM, from, "#23 = CARTESIAN_POINT('',(8.000001,0.,0.))");
    assert_eq!(solids(&moved).unwrap_err(), Error::Invalid);
    let within = tampered(
        FRUSTUM,
        from,
        "#23 = CARTESIAN_POINT('',(8.000000001,0.,0.))",
    );
    assert!(close(
        solids(&within).unwrap()[0].volume_mm3(),
        PI * 10. / 3. * 112.,
        1e-9
    ));
}

#[test]
fn exported_boundaries_read_back_with_the_same_geometry() {
    for text in [
        FRUSTUM, SPHERE, TORUS, CYLINDER, VOID, VASE, SPLINE, DRILLED, STANDOFF, BUMPS,
    ] {
        let source = solid(text);
        let mut writer = crate::step::Writer::new();
        writer.append_brep(&source).unwrap();
        let exported = writer.finish().unwrap();
        let back = solid(&exported);
        assert!(close(back.volume_mm3(), source.volume_mm3(), 1e-12));
        assert_eq!(back.faces().len(), source.faces().len());
        assert_eq!(back.boundary().edges.len(), source.boundary().edges.len());
        let sorted = |solid: &Solid| {
            let mut vertices = solid.boundary().vertices.clone();
            vertices.sort_by(|a, b| a.partial_cmp(b).unwrap());
            vertices
        };
        assert_eq!(sorted(&back), sorted(&source));
        let (a, b) = (back.bounds_mm(), source.bounds_mm());
        assert!(distance(a[0], b[0]) < 1e-9 && distance(a[1], b[1]) < 1e-9);
        for (one, other) in back.faces().iter().zip(source.faces()) {
            assert!(close(one.area_mm2, other.area_mm2, 1e-12));
        }
        // The independent readback entry point accepts both texts.
        for text in [text, exported.as_str()] {
            let members = read_analytic_compound(text, TOLERANCE).unwrap();
            let volume: f64 = members.iter().map(|member| member.volume_mm3).sum();
            assert!(close(volume, source.volume_mm3(), 1e-9), "{volume}");
        }
    }
}

#[test]
fn a_refusal_of_both_readers_keeps_the_first_reason() {
    let open = tampered(
        FRUSTUM,
        "CLOSED_SHELL('',(#17,#54,#109))",
        "CLOSED_SHELL('',(#17,#54))",
    );
    let first = crate::step_read::read_analytic_solid(&open, TOLERANCE).unwrap_err();
    assert_eq!(read_analytic_compound(&open, TOLERANCE).unwrap_err(), first);
    assert_eq!(ReadError::from(Error::Cancelled), ReadError::Invalid);
    assert_eq!(ReadError::from(Error::Limit), ReadError::Limit);
}

#[test]
fn sections_lie_on_the_exact_boundary_and_close() {
    // The frustum's radius falls from 8 to 4 over 10 mm: 6 mm at mid-height.
    let frustum = solid(FRUSTUM);
    let loops =
        section::section(&frustum, [0., 0., 5.], [0., 0., 1.], options(0.01), None).unwrap();
    assert_eq!(loops.len(), 1);
    let radial = loops[0]
        .iter()
        .map(|p| (p[0].hypot(p[1]) - 6.).abs())
        .fold(0., f64::max);
    let axial = loops[0]
        .iter()
        .map(|p| (p[2] - 5.).abs())
        .fold(0., f64::max);
    assert!(
        radial < 1e-9 && axial < 1e-12,
        "{radial:e} {axial:e} {}",
        loops[0].len()
    );
    let length = section::length(&loops[0]);
    assert!(length < 12. * PI && length > 12. * PI * 0.999, "{length}");
    // Material is left of the loop seen against the normal: counter-clockwise.
    let turn: f64 = (0..loops[0].len())
        .map(|i| {
            let (a, b) = (loops[0][i], loops[0][(i + 1) % loops[0].len()]);
            a[0] * b[1] - b[0] * a[1]
        })
        .sum();
    assert!(turn > 0.);

    // A box with a void: the outer rectangle and the void's, as four
    // straight sides each.
    let hollow = solid(VOID);
    let mut loops =
        section::section(&hollow, [0., 0., 5.], [0., 0., 1.], options(0.01), None).unwrap();
    loops.sort_by(|a, b| section::length(a).total_cmp(&section::length(b)));
    assert_eq!(loops.iter().map(Vec::len).collect::<Vec<_>>(), [4, 4]);
    assert!(close(section::length(&loops[1]), 100., 1e-12));

    // A plane that misses the solid cuts nothing, and an interrupt stops it.
    assert!(
        section::section(&frustum, [0., 0., 50.], [0., 0., 1.], options(0.01), None)
            .unwrap()
            .is_empty()
    );
    let stop = || true;
    let stopped = section::section(
        &frustum,
        [0., 0., 5.],
        [0., 0., 1.],
        options(0.01),
        Some(&stop),
    );
    assert_eq!(stopped.unwrap_err(), Error::Cancelled);

    // A sphere cut off-centre: a circle of radius 4 at height 3.
    let sphere = solid(SPHERE);
    let loops = section::section(&sphere, [0., 0., 3.], [0., 0., 1.], options(0.01), None).unwrap();
    assert_eq!(loops.len(), 1);
    assert!(
        loops[0]
            .iter()
            .all(|p| (p[0].hypot(p[1]) - 4.).abs() < 1e-9)
    );
}
