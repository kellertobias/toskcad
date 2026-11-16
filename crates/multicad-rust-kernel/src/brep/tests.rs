//! The retained boundary provider against analytic oracles and tampered
//! exchange text.
use super::domain::{Domain, Key, LoopPoint, UvLoop};
use super::*;
use crate::manufacturing::certify;
use crate::step_read::ReadError;
use crate::step_read::brep::{read_boundaries, read_solids, verify_retained};
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

/// KiCad library model, copied byte for byte with its CC-BY-SA 4.0 licence
/// header (KiCad library exception); its two B-spline wedges meet planes
/// at a thin angle along curved edges.
const WEDGES: &str = include_str!("fixtures/D_SOD-923.step");

/// KiCad library model, copied with its licence header like `WEDGES`: a
/// disc whose revolved B-spline faces close near their axis.
const DISC: &str = include_str!("fixtures/RV_Disc_D9mm_W6.1mm_P5mm.step");

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
fn spline_extrusion_matches_its_closed_forms() {
    // The profile is the cubic x = 20 - 60 t^2 + 40 t^3, y = 36 t (1 - t),
    // swept 6 mm: the prism holds 6 mm times the 144 mm^2 cap. The side is
    // 6 mm times the profile's arc length, 30.1537396807876 mm by
    // independent quadrature.
    let prism = solid(SPLINE);
    assert!(close(prism.volume_mm3(), 864., 1e-12));
    let mut areas = sorted_areas(&prism);
    let side = areas.remove(3);
    assert!(close(side, 6. * 30.153_739_680_787_6, 1e-10), "{side}");
    assert!(
        areas
            .iter()
            .zip([120., 144., 144.])
            .all(|(&a, b)| close(a, b, 1e-12)),
        "{areas:?}"
    );
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
    }
}

#[test]
fn read_errors_carry_the_provider_refusal_kind() {
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

#[test]
fn sections_publish_exact_lines_and_ellipses() {
    use crate::warped_union::Curve as SectionCurve;
    use std::f64::consts::TAU;
    let edges = |solid: &Solid, origin: Point, normal: Point| {
        section::section_edges(solid, origin, normal, options(0.01), None).unwrap()
    };
    // The frustum cut across its axis at mid-height: one whole circle of
    // radius 6.
    let frustum = solid(FRUSTUM);
    let loops = edges(&frustum, [0., 0., 5.], [0., 0., 1.]);
    assert_eq!(loops.len(), 1);
    let [edge] = loops[0].as_slice() else {
        panic!("{loops:?}")
    };
    let SectionCurve::Ellipse(arc) = edge.curve else {
        panic!("{edge:?}")
    };
    assert_eq!(arc.sweep_rad(), TAU);
    let ellipse = arc.ellipse();
    assert!(
        ellipse.semiaxes_mm().iter().all(|r| (r - 6.).abs() < 1e-9),
        "{ellipse:?}"
    );
    assert!(distance(ellipse.center_mm(), [0., 0., 5.]) < 1e-9);

    // Tilted through the side and both caps: two chords of the caps and two
    // arcs of the cone's ellipse. The cone x^2 + y^2 = (8 - 0.4 z)^2 meets
    // the plane z = 5 + 0.5 x in the ellipse whose vertices on y = 0 solve
    // x = +-(8 - 0.4 (5 + 0.5 x)).
    let normal = [-0.5, 0., 1.];
    let loops = edges(&frustum, [0., 0., 5.], normal);
    assert_eq!(loops.len(), 1, "{loops:?}");
    let kinds: Vec<_> = loops[0]
        .iter()
        .map(|e| matches!(e.curve, SectionCurve::Line))
        .collect();
    assert_eq!(kinds.iter().filter(|&&line| line).count(), 0, "{kinds:?}");
    let [edge] = loops[0].as_slice() else {
        panic!("{loops:?}")
    };
    let SectionCurve::Ellipse(arc) = edge.curve else {
        panic!("{edge:?}")
    };
    let [right, left] = [6. / 1.2, -6. / 0.8];
    let ellipse = arc.ellipse();
    assert!((ellipse.semiaxes_mm()[0] - 0.5 * (right - left) * 1.25f64.sqrt()).abs() < 1e-9);
    assert!((ellipse.center_mm()[0] - 0.5 * (right + left)).abs() < 1e-9);
    // Every crossing of the mesh lies on that exact arc.
    for ring in section::section(&frustum, [0., 0., 5.], normal, options(0.01), None).unwrap() {
        for p in ring {
            assert!(ellipse.parameter_at(p, 1e-9).is_some(), "{p:?}");
        }
    }

    // Steep enough to reach both caps: two exact chords and two arcs.
    let normal = [-2., 0., 1.];
    let loops = edges(&frustum, [0., 0., 5.], normal);
    assert_eq!(loops.len(), 1, "{loops:?}");
    let lines = loops[0]
        .iter()
        .filter(|e| matches!(e.curve, SectionCurve::Line))
        .count();
    let arcs = loops[0]
        .iter()
        .filter(|e| matches!(e.curve, SectionCurve::Ellipse(_)))
        .count();
    assert_eq!((lines, arcs), (2, 2), "{loops:?}");
    for pair in loops[0].iter().zip(loops[0].iter().cycle().skip(1)) {
        assert_eq!(pair.0.ends_mm[1], pair.1.ends_mm[0]);
    }

    // A sphere's section is a circle; a torus keeps its chords.
    let sphere = solid(SPHERE);
    let loops = edges(&sphere, [0., 0., 3.], [0., 0., 1.]);
    assert!(
        matches!(loops[0].as_slice(), [e] if matches!(e.curve, SectionCurve::Ellipse(a) if (a.ellipse().semiaxes_mm()[0] - 4.).abs() < 1e-9))
    );
    let torus = solid(TORUS);
    let bounds = torus.bounds_mm();
    let middle = std::array::from_fn(|i| 0.5 * (bounds[0][i] + bounds[1][i]));
    let loops = edges(&torus, middle, [1., 0., 0.]);
    assert!(!loops.is_empty());
    assert!(
        loops
            .iter()
            .flatten()
            .all(|e| matches!(e.curve, SectionCurve::Line))
    );
}

const STATED: &str = "#23 = CARTESIAN_POINT('',(8.,0.,0.))";
const DECLARED: &str = "UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-07)";

fn declared(text: &str, uncertainty: &str) -> String {
    tampered(
        text,
        DECLARED,
        &format!("UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE({uncertainty})"),
    )
}

#[test]
fn a_declared_uncertainty_admits_vertices_and_re_derives_them_on_their_supports() {
    let exact = PI * 10. / 3. * 112.;
    let moved = tampered(
        FRUSTUM,
        STATED,
        "#23 = CARTESIAN_POINT('',(8.000005,0.,0.))",
    );
    // The file's own 1e-7 mm declaration does not admit a 5 um gap.
    assert_eq!(solids(&moved).unwrap_err(), Error::Invalid);
    // 2e-5 mm is capped at 1e-5 mm, which does: the vertex returns to the
    // point where its circle, seam, plane and cone actually meet.
    let wider = declared(&moved, "2.E-05");
    let solid = solids(&wider).unwrap().remove(0);
    assert_eq!(
        solid.tolerance(),
        crate::step_read::brep::MAX_DECLARED_TOLERANCE_MM
    );
    let vertex = solid
        .boundary()
        .vertices
        .iter()
        .copied()
        .min_by(|a, b| distance(*a, [8., 0., 0.]).total_cmp(&distance(*b, [8., 0., 0.])))
        .unwrap();
    assert!(distance(vertex, [8., 0., 0.]) < 1e-12, "{vertex:?}");
    assert!(close(solid.volume_mm3(), exact, 1e-9));
    certified(&solid, 0.05);
    // A gap beyond the capped tolerance is never read. Within the file's
    // own declared uncertainty it is unsupported; beyond that, invalid.
    let far = tampered(FRUSTUM, STATED, "#23 = CARTESIAN_POINT('',(8.00002,0.,0.))");
    assert_eq!(
        solids(&declared(&far, "1.E-03")).unwrap_err(),
        Error::Unsupported
    );
    assert_eq!(
        solids(&declared(&far, "1.5E-05")).unwrap_err(),
        Error::Invalid
    );
    // So do a flipped face and an open shell.
    let flipped = tampered(
        &wider,
        "#54 = ADVANCED_FACE('',(#55),#43,.T.)",
        "#54 = ADVANCED_FACE('',(#55),#43,.F.)",
    );
    assert_eq!(solids(&flipped).unwrap_err(), Error::Invalid);
}

#[test]
fn only_positive_length_uncertainties_widen_incidence() {
    let moved = tampered(
        FRUSTUM,
        STATED,
        "#23 = CARTESIAN_POINT('',(8.000005,0.,0.))",
    );
    let source = |text: &str| {
        crate::step_read::brep::Source::parse(text, TOLERANCE).map(|s| s.declared_uncertainty())
    };
    assert_eq!(source(FRUSTUM), Ok(Some(1e-7)));
    let angle = tampered(
        &moved,
        "LENGTH_MEASURE(1.E-07)",
        "PLANE_ANGLE_MEASURE(1.E-02)",
    );
    assert_eq!(source(&angle), Ok(None));
    assert_eq!(solids(&angle).unwrap_err(), Error::Invalid);
    assert_eq!(source(&declared(FRUSTUM, "0.")), Err(ReadError::Invalid));
    assert_eq!(
        source(&declared(FRUSTUM, "1.E-07,2.")),
        Err(ReadError::Invalid)
    );
    let untyped = tampered(
        FRUSTUM,
        "UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-07),",
        "UNCERTAINTY_MEASURE_WITH_UNIT(1.E-07,",
    );
    assert_eq!(source(&untyped), Err(ReadError::Invalid));
    let dangling = tampered(
        FRUSTUM,
        "GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#117))",
        "GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#99999))",
    );
    assert_eq!(source(&dangling), Err(ReadError::Invalid));
    let bare = tampered(
        FRUSTUM,
        "GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#117))",
        "GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#117),#117)",
    );
    assert_eq!(source(&bare), Err(ReadError::Invalid));
}

#[test]
fn re_derivation_refuses_supports_that_do_not_meet_near_the_vertex() {
    let boundary = read_boundaries(FRUSTUM, TOLERANCE).unwrap().remove(0);
    let (same, residual) = vertex::rederive(&boundary, TOLERANCE).unwrap();
    assert!(residual < 1e-9, "{residual:e}");
    assert_eq!(same.faces, boundary.faces);
    let mut shifted = boundary.clone();
    shifted.vertices[0][2] += 1e-3;
    assert_eq!(
        vertex::rederive(&shifted, 1e-5).unwrap_err(),
        Error::Invalid
    );
}

#[test]
fn a_stated_vertex_nearer_to_its_supports_than_their_balance_stands() {
    // A vertex 0.95 tolerances above three faces in z = 0 and an edge in it,
    // and 0.95 below an edge 1.9 above them: the supports have no common
    // point, and their balance at 0.38 lies 1.52 from the high edge. The
    // stated vertex is within the tolerance of all five.
    let tolerance = 1e-5;
    let line = |origin: Point, direction: Point| Curve::Line { origin, direction };
    let stated = [0., 0., 0.95 * tolerance];
    let vertices = vec![stated, [10., 0., 1.9 * tolerance], [0., 10., 0.]];
    let high = line([0., 0., 1.9 * tolerance], [1., 0., 0.]);
    let low = line([0.; 3], [0., 1., 0.]);
    let edges = vec![
        Edge::new(high, [0, 1], [vertices[0], vertices[1]], tolerance).unwrap(),
        Edge::new(low, [0, 2], [vertices[0], vertices[2]], tolerance).unwrap(),
    ];
    let plane = Surface::Plane(Frame {
        origin: [0.; 3],
        x: [1., 0., 0.],
        y: [0., 1., 0.],
        z: [0., 0., 1.],
    });
    let face = Face {
        surface: plane,
        same_sense: true,
        loops: vec![Loop {
            uses: vec![(1, true)],
            vertex: None,
        }],
    };
    let boundary = Boundary {
        vertices,
        edges,
        faces: vec![face.clone(), face.clone(), face],
        shells: vec![],
    };
    let (kept, residual) = vertex::rederive(&boundary, tolerance).unwrap();
    assert_eq!(kept.vertices[0], stated);
    assert!((residual - 0.95 * tolerance).abs() < 1e-12, "{residual:e}");
    // A tolerance below the stated vertex's own distance refuses it.
    assert_eq!(
        vertex::rederive(&boundary, 0.9 * tolerance).unwrap_err(),
        Error::Invalid
    );
}

#[test]
fn thin_wedges_along_curved_edges_mesh_without_folds() {
    // OpenCascade 7.9.3 measures 0.17402361314366 mm3 for this model.
    let solid = solid(WEDGES);
    assert!(close(solid.volume_mm3(), 0.174_023_613_143_661_5, 1e-11));
    assert_eq!(solid.faces().len(), 27);
    let area: f64 = solid.faces().iter().map(|face| face.area_mm2).sum();
    assert!(close(area, 2.053_581_607_713_32, 1e-9), "{area}");
    certified(&solid, 0.05);
    certified(&solid, 0.005);
}

#[test]
fn revolved_faces_closing_near_their_axis_refine_to_the_deflection() {
    // OpenCascade 7.9.3 measures 298.900116969481 mm3 and 234.676869 mm2.
    let solid = solid(DISC);
    assert!(close(solid.volume_mm3(), 298.900_116_969_481, 1e-9));
    let area: f64 = solid.faces().iter().map(|face| face.area_mm2).sum();
    assert!(close(area, 234.676_869_496_623, 1e-7), "{area}");
    certified(&solid, 0.05);
}

#[test]
fn retained_exports_verify_only_against_their_own_boundary() {
    let frustum = solid(FRUSTUM);
    let mut writer = crate::step::Writer::new();
    writer.append_brep(&frustum).unwrap();
    let text = writer.finish().unwrap();
    let volume = verify_retained(&text, &[&frustum], TOLERANCE).unwrap();
    assert!(close(volume, frustum.volume_mm3(), 1e-12));
    let hollow = solid(VOID);
    assert_eq!(
        verify_retained(&text, &[&hollow], TOLERANCE),
        Err(ReadError::Invalid)
    );
    assert_eq!(
        verify_retained(&text, &[&frustum, &frustum], TOLERANCE),
        Err(ReadError::Invalid)
    );
    let shift = |p: Point| add(p, [1e-3, 0., 0.]);
    let moved = frustum.moved(&shift, &|v| v).unwrap();
    assert_eq!(
        verify_retained(&text, &[&moved], TOLERANCE),
        Err(ReadError::Invalid)
    );
    assert_eq!(
        verify_retained("ISO-10303-21;", &[&frustum], TOLERANCE),
        Err(ReadError::Invalid)
    );
}

/// A millimetre/radian Part 21 file around `body`, whose representation
/// holds the solids `items`.
fn part21(items: &str, body: &str) -> String {
    format!(
        "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION((''),'2;1');\nFILE_NAME('t','',(''),(''),'','','');\n\
         FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\nENDSEC;\nDATA;\n\
         #900=ADVANCED_BREP_SHAPE_REPRESENTATION('',({items}),#901);\n\
         #901=(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#905))\
         GLOBAL_UNIT_ASSIGNED_CONTEXT((#902,#903,#904))REPRESENTATION_CONTEXT('',''));\n\
         #902=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.));\n\
         #903=(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT($,.RADIAN.));\n\
         #904=(NAMED_UNIT(*)SI_UNIT($,.STERADIAN.)SOLID_ANGLE_UNIT());\n\
         #905=UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-07),#902,'','');\n\
         {body}ENDSEC;\nEND-ISO-10303-21;\n"
    )
}

/// An elliptic cylinder, semiaxes 3 and 2 mm and 5 mm tall: two planar
/// caps bounded by ellipses and a linear extrusion of the lower ellipse.
const ELLIPTIC: &str = "\
#1=CARTESIAN_POINT('',(0.,0.,0.));#2=CARTESIAN_POINT('',(0.,0.,5.));\
#3=DIRECTION('',(0.,0.,1.));#4=DIRECTION('',(1.,0.,0.));\
#5=AXIS2_PLACEMENT_3D('',#1,#3,#4);#6=AXIS2_PLACEMENT_3D('',#2,#3,#4);\
#7=ELLIPSE('',#5,3.,2.);#8=ELLIPSE('',#6,3.,2.);\
#9=VERTEX_POINT('',#10);#10=CARTESIAN_POINT('',(3.,0.,0.));\
#11=VERTEX_POINT('',#12);#12=CARTESIAN_POINT('',(3.,0.,5.));\
#13=EDGE_CURVE('',#9,#9,#7,.T.);#14=EDGE_CURVE('',#11,#11,#8,.T.);\
#15=VECTOR('',#3,1.);#16=LINE('',#10,#15);#17=EDGE_CURVE('',#9,#11,#16,.T.);\
#20=PLANE('',#5);#21=PLANE('',#6);#22=SURFACE_OF_LINEAR_EXTRUSION('',#7,#15);\
#30=ORIENTED_EDGE('',*,*,#13,.F.);#31=EDGE_LOOP('',(#30));#32=FACE_OUTER_BOUND('',#31,.T.);\
#33=ADVANCED_FACE('',(#32),#20,.F.);\
#40=ORIENTED_EDGE('',*,*,#14,.T.);#41=EDGE_LOOP('',(#40));#42=FACE_OUTER_BOUND('',#41,.T.);\
#43=ADVANCED_FACE('',(#42),#21,.T.);\
#50=ORIENTED_EDGE('',*,*,#13,.T.);#51=ORIENTED_EDGE('',*,*,#17,.T.);\
#52=ORIENTED_EDGE('',*,*,#14,.F.);#53=ORIENTED_EDGE('',*,*,#17,.F.);\
#54=EDGE_LOOP('',(#50,#51,#52,#53));#55=FACE_BOUND('',#54,.T.);\
#56=ADVANCED_FACE('',(#55),#22,.T.);\
#60=CLOSED_SHELL('',(#33,#43,#56));#61=MANIFOLD_SOLID_BREP('',#60);\n";

/// Perimeter of an ellipse by a fine midpoint rule, independent of the
/// provider's quadrature.
fn perimeter(a: f64, b: f64) -> f64 {
    let n = 200_000;
    (0..n)
        .map(|k| {
            let t = (k as f64 + 0.5) / n as f64 * 2. * PI;
            (a * t.sin()).hypot(b * t.cos()) * 2. * PI / n as f64
        })
        .sum()
}

#[test]
fn elliptic_cylinder_matches_its_closed_forms_and_round_trips() {
    let text = part21("#61", ELLIPTIC);
    let elliptic = solid(&text);
    assert!(
        close(elliptic.volume_mm3(), 30. * PI, 1e-12),
        "{}",
        elliptic.volume_mm3()
    );
    let areas = sorted_areas(&elliptic);
    assert!(close(areas[0], 6. * PI, 1e-12) && close(areas[1], 6. * PI, 1e-12));
    assert!(close(areas[2], 5. * perimeter(3., 2.), 1e-9), "{areas:?}");
    let [low, high] = elliptic.bounds_mm();
    assert!(distance(low, [-3., -2., 0.]) < 1e-9 && distance(high, [3., 2., 5.]) < 1e-9);
    certified(&elliptic, 0.01);
    let loops =
        section::section(&elliptic, [0., 0., 2.5], [0., 0., 1.], options(0.001), None).unwrap();
    assert_eq!(loops.len(), 1);
    assert!(
        loops[0]
            .iter()
            .all(|p| ((p[0] / 3.).hypot(p[1] / 2.) - 1.).abs() < 1e-9)
    );
    let mut writer = crate::step::Writer::new();
    writer.append_brep(&elliptic).unwrap();
    let back = solid(&writer.finish().unwrap());
    assert!(close(back.volume_mm3(), 30. * PI, 1e-12));
    let turned = elliptic
        .moved(&|p| [p[1], -p[0], p[2]], &|p| [p[1], -p[0], p[2]])
        .unwrap();
    assert!(distance(turned.bounds_mm()[0], [-2., -3., 0.]) < 1e-9);
}

/// The authored elliptic cylinder with `from` replaced by `to`.
fn elliptic(from: &str, to: &str) -> String {
    tampered(&part21("#61", ELLIPTIC), from, to)
}

#[test]
fn malformed_and_foreign_structure_is_refused_with_its_kind() {
    let refused = |text: &str| solids(text).map(|solids| solids.len()).unwrap_err();
    let extra = |entity: &str| elliptic("#60=", &format!("{entity}#60="));
    // Shapes that are not solids keep their meaning.
    assert_eq!(
        refused(&extra("#990=OPEN_SHELL('',(#33));")),
        Error::Unsupported
    );
    // Product structure belongs to the assembly layer.
    assert_eq!(
        refused(&extra("#990=MAPPED_ITEM('',#5,#5);")),
        Error::Unsupported
    );
    // A face outside every solid would be dropped.
    assert_eq!(
        refused(&extra("#990=ADVANCED_FACE('',(#32),#20,.F.);")),
        Error::Unsupported
    );
    // No solid, or no units.
    assert_eq!(
        refused(&elliptic("#61=MANIFOLD_SOLID_BREP('',#60);", "")),
        Error::Unsupported
    );
    assert_eq!(
        refused(&elliptic(
            "GLOBAL_UNIT_ASSIGNED_CONTEXT((#902,#903,#904))",
            ""
        )),
        Error::Unsupported
    );
    assert_eq!(
        refused(&elliptic(
            "GLOBAL_UNIT_ASSIGNED_CONTEXT((#902,#903,#904))",
            "GLOBAL_UNIT_ASSIGNED_CONTEXT((#902),#903)"
        )),
        Error::Invalid
    );
    assert_eq!(
        refused(&elliptic(
            "GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#905))",
            "GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#905),#905)"
        )),
        Error::Invalid
    );
    // Inches are not millimetres; a second length unit is ambiguous.
    assert_eq!(
        refused(&elliptic("SI_UNIT(.MILLI.,.METRE.)", "SI_UNIT($,.METRE.)")),
        Error::Unsupported
    );
    assert_eq!(
        refused(&elliptic("((#902,#903,#904))", "((#902,#902,#903,#904))")),
        Error::Unsupported
    );
    assert_eq!(
        refused(&elliptic("((#902,#903,#904))", "((#902,#903,#903,#904))")),
        Error::Unsupported
    );
    assert_eq!(
        refused(&elliptic("((#902,#903,#904))", "((#902,#904))")),
        Error::Unsupported
    );
    assert_eq!(
        refused(&elliptic("((#902,#903,#904))", "((#902,#903,#999))")),
        Error::Invalid
    );
    // Wrong arities and kinds.
    assert_eq!(
        refused(&elliptic(
            "MANIFOLD_SOLID_BREP('',#60)",
            "MANIFOLD_SOLID_BREP('',#60,#60)"
        )),
        Error::Invalid
    );
    assert_eq!(
        refused(&elliptic(
            "MANIFOLD_SOLID_BREP('',#60)",
            "MANIFOLD_SOLID_BREP('',#33)"
        )),
        Error::Unsupported
    );
    assert_eq!(
        refused(&elliptic(
            "CLOSED_SHELL('',(#33,#43,#56))",
            "CLOSED_SHELL('',(#33,#43,#56),#1)"
        )),
        Error::Invalid
    );
    assert_eq!(
        refused(&elliptic(
            "CLOSED_SHELL('',(#33,#43,#56))",
            "CLOSED_SHELL('',(#33,#33,#43,#56))"
        )),
        Error::Invalid
    );
    assert_eq!(
        refused(&elliptic("#15=VECTOR('',#3,1.)", "#15=VECTOR('',#3,0.)")),
        Error::Invalid
    );
    assert_eq!(
        refused(&elliptic("#15=VECTOR('',#3,1.)", "#15=VECTOR('',#3)")),
        Error::Invalid
    );
    assert_eq!(
        refused(&elliptic("#15=VECTOR('',#3,1.)", "#15=VECTOR('',#3,1.E13)")),
        Error::Limit
    );
    assert_eq!(
        refused(&elliptic(
            "#3=DIRECTION('',(0.,0.,1.))",
            "#3=DIRECTION('',(0.,0.,1.),1)"
        )),
        Error::Invalid
    );
    assert_eq!(
        refused(&elliptic(
            "#3=DIRECTION('',(0.,0.,1.))",
            "#3=DIRECTION('',(0.,0.,0.))"
        )),
        Error::Invalid
    );
    assert_eq!(
        refused(&elliptic(
            "#1=CARTESIAN_POINT('',(0.,0.,0.))",
            "#1=CARTESIAN_POINT('',(0.,0.,0.),1)"
        )),
        Error::Invalid
    );
    assert_eq!(
        refused(&elliptic(
            "#5=AXIS2_PLACEMENT_3D('',#1,#3,#4)",
            "#5=AXIS2_PLACEMENT_3D('',#1,#3)"
        )),
        Error::Invalid
    );
    assert_eq!(
        refused(&elliptic(
            "#5=AXIS2_PLACEMENT_3D('',#1,#3,#4)",
            "#5=AXIS2_PLACEMENT_3D('',#1,#3,#3)"
        )),
        Error::Invalid
    );
    assert_eq!(
        refused(&elliptic(
            "#7=ELLIPSE('',#5,3.,2.)",
            "#7=ELLIPSE('',#5,3.,0.)"
        )),
        Error::Invalid
    );
    assert_eq!(
        refused(&elliptic(
            "#7=ELLIPSE('',#5,3.,2.)",
            "#7=OFFSET_CURVE_3D('',#5,3.,2.)"
        )),
        Error::Unsupported
    );
    assert_eq!(
        refused(&elliptic(
            "#16=LINE('',#10,#15)",
            "#16=(LINE('',#10,#15)CURVE())"
        )),
        Error::Unsupported
    );
    assert_eq!(
        refused(&elliptic(
            "#20=PLANE('',#5)",
            "#20=CYLINDRICAL_SURFACE('',#5,0.)"
        )),
        Error::Invalid
    );
    assert_eq!(
        refused(&elliptic(
            "#20=PLANE('',#5)",
            "#20=OFFSET_SURFACE('',#21,1.,.F.)"
        )),
        Error::Unsupported
    );
    assert_eq!(
        refused(&elliptic("#20=PLANE('',#5)", "#20=(PLANE('',#5)SURFACE())")),
        Error::Unsupported
    );
    assert_eq!(
        refused(&elliptic(
            "#9=VERTEX_POINT('',#10)",
            "#9=VERTEX_POINT('',#10,#10)"
        )),
        Error::Invalid
    );
    assert_eq!(
        refused(&elliptic(
            "#13=EDGE_CURVE('',#9,#9,#7,.T.)",
            "#13=EDGE_CURVE('',#9,#9,#7)"
        )),
        Error::Invalid
    );
    assert_eq!(
        refused(&elliptic(
            "#13=EDGE_CURVE('',#9,#9,#7,.T.)",
            "#13=EDGE_CURVE('',#9,#9,#7,.U.)"
        )),
        Error::Invalid
    );
    assert_eq!(
        refused(&elliptic(
            "#32=FACE_OUTER_BOUND('',#31,.T.)",
            "#32=FACE_OUTER_BOUND('',#31)"
        )),
        Error::Invalid
    );
    assert_eq!(
        refused(&elliptic(
            "#32=FACE_OUTER_BOUND('',#31,.T.)",
            "#32=POLY_LOOP('',#31,.T.)"
        )),
        Error::Invalid
    );
    assert_eq!(
        refused(&elliptic("#31=EDGE_LOOP('',(#30))", "#31=EDGE_LOOP('',())")),
        Error::Invalid
    );
    assert_eq!(
        refused(&elliptic(
            "#31=EDGE_LOOP('',(#30))",
            "#31=POLY_LOOP('',(#10))"
        )),
        Error::Unsupported
    );
    assert_eq!(
        refused(&elliptic(
            "#30=ORIENTED_EDGE('',*,*,#13,.F.)",
            "#30=ORIENTED_EDGE('',*,*,#13)"
        )),
        Error::Invalid
    );
    assert_eq!(
        refused(&elliptic(
            "#33=ADVANCED_FACE('',(#32),#20,.F.)",
            "#33=ADVANCED_FACE('',(#32),#20)"
        )),
        Error::Invalid
    );
    // Outside the model's coordinate range.
    assert_eq!(
        refused(&elliptic(
            "#7=ELLIPSE('',#5,3.,2.)",
            "#7=ELLIPSE('',#5,3.,'2')"
        )),
        Error::Invalid
    );
    // A loop that passes a vertex of another edge is discontinuous.
    assert_eq!(
        refused(&elliptic(
            "#51=ORIENTED_EDGE('',*,*,#17,.T.)",
            "#51=ORIENTED_EDGE('',*,*,#17,.F.)"
        )),
        Error::Invalid
    );
}

#[test]
fn optional_placement_axes_take_their_standard_defaults() {
    // Omitted axes default to z and x: the same cylinder.
    let text = elliptic(
        "#5=AXIS2_PLACEMENT_3D('',#1,#3,#4)",
        "#5=AXIS2_PLACEMENT_3D('',#1,$,$)",
    );
    assert!(close(solid(&text).volume_mm3(), 30. * PI, 1e-12));
    // A void must face into its cavity: a plain closed shell, or an oriented
    // one kept in its stated sense, encloses material instead.
    for (from, to) in [
        (
            "BREP_WITH_VOIDS('',#16,(#345))",
            "BREP_WITH_VOIDS('',#16,(#346))",
        ),
        (
            "ORIENTED_CLOSED_SHELL('',*,#346,.F.)",
            "ORIENTED_CLOSED_SHELL('',*,#346,.T.)",
        ),
    ] {
        assert_eq!(
            solids(&tampered(VOID, from, to)).unwrap_err(),
            Error::Invalid
        );
    }
}

#[test]
fn point_triangle_distances_cover_every_closest_feature() {
    let triangle = [[0., 0., 0.], [2., 0., 0.], [0., 2., 0.]];
    for (point, expected) in [
        ([-1., -1., 0.], 2f64.sqrt()),
        ([3., -1., 0.], 2f64.sqrt()),
        ([-1., 3., 0.], 2f64.sqrt()),
        ([1., -1., 0.], 1.),
        ([-1., 1., 0.], 1.),
        ([2., 2., 0.], 2f64.sqrt()),
        ([0.5, 0.5, 3.], 3.),
    ] {
        assert!(
            (point_triangle(point, triangle) - expected).abs() < 1e-12,
            "{point:?}"
        );
    }
    // A ray in the triangle's own plane never crosses it.
    let ray = [
        0.574_352_934_563_6,
        0.577_919_283_817_1,
        0.579_868_302_102_6,
    ];
    assert!(!crosses_ray([0.; 3], &[[0.; 3], ray, [1., 0., 0.]]));
    assert!(crosses_ray(
        [-0.5, -0.5, -1.],
        &[[0., 0., 0.], [2., 0., 0.], [0., 2., 0.]]
    ));
    assert!(!crosses_ray(
        [-0.5, -0.5, 1.],
        &[[0., 0., 0.], [2., 0., 0.], [0., 2., 0.]]
    ));
}

/// The box with its void moved by `map` (`replace`), or with a moved copy
/// of the void as a second void.
fn void_copy(map: &dyn Fn(Point) -> Point, replace: bool) -> Boundary {
    let source = read_boundaries(VOID, TOLERANCE).unwrap().remove(0);
    let moved = source.mapped(map, &|v| v);
    let inner = &source.shells[1];
    let mut boundary = source.clone();
    if replace {
        for &face in &inner.faces {
            boundary.faces[face] = moved.faces[face].clone();
            for trim in &moved.faces[face].loops {
                for &(edge, _) in &trim.uses {
                    boundary.edges[edge] = moved.edges[edge].clone();
                    for vertex in moved.edges[edge].vertices {
                        boundary.vertices[vertex] = moved.vertices[vertex];
                    }
                }
            }
        }
        return boundary;
    }
    let (vertices, edges, faces) = (
        boundary.vertices.len(),
        boundary.edges.len(),
        boundary.faces.len(),
    );
    let mut copy = Shell {
        faces: vec![],
        reversed: inner.reversed,
    };
    let mut edge_ids = std::collections::BTreeMap::new();
    let mut vertex_ids = std::collections::BTreeMap::new();
    for &face in &inner.faces {
        let mut face_copy = moved.faces[face].clone();
        for trim in &mut face_copy.loops {
            for used in &mut trim.uses {
                let next = edges + edge_ids.len();
                used.0 = *edge_ids.entry(used.0).or_insert(next);
            }
        }
        copy.faces.push(faces + copy.faces.len());
        boundary.faces.push(face_copy);
    }
    let mut ordered: Vec<_> = edge_ids.into_iter().collect();
    ordered.sort_by_key(|&(_, new)| new);
    for (old, _) in ordered {
        let mut edge = moved.edges[old].clone();
        for vertex in &mut edge.vertices {
            let next = vertices + vertex_ids.len();
            let id = *vertex_ids.entry(*vertex).or_insert(next);
            if id == vertices + boundary.vertices.len() - vertices {
                boundary.vertices.push(moved.vertices[*vertex]);
            }
            *vertex = id;
        }
        boundary.edges.push(edge);
    }
    boundary.shells.push(copy);
    boundary
}

#[test]
fn voids_must_lie_inside_the_exterior_and_apart() {
    let shift = |by: Point| move |p: Point| add(p, by);
    // The 240 mm3 void moved within the box, and a copy of it 0.5 mm
    // from both the void and the wall.
    let same = Solid::new(void_copy(&shift([-5., 1., 0.]), true), TOLERANCE, None).unwrap();
    assert!(close(same.volume_mm3(), 5760., 1e-12));
    let two = Solid::new(void_copy(&shift([0., -6.5, 0.]), false), TOLERANCE, None).unwrap();
    assert!(
        close(two.volume_mm3(), 5520., 1e-12),
        "{}",
        two.volume_mm3()
    );
    // A void outside the box, through its wall, or inside the first void
    // is invalid; one touching the wall cannot be resolved and is refused.
    for (by, error) in [
        ([25., 0., 0.], Error::Invalid),
        ([0., 0., 5.], Error::Invalid),
        ([0., 0., 3.], Error::Unsupported),
    ] {
        let moved = Solid::new(void_copy(&shift(by), true), TOLERANCE, None);
        assert_eq!(moved.map(|s| s.volume_mm3()).unwrap_err(), error, "{by:?}");
    }
    // A copy overlapping the first void has vertices inside it.
    let inside = Solid::new(void_copy(&shift([2., 1., 0.5]), false), TOLERANCE, None);
    assert_eq!(inside.map(|s| s.volume_mm3()).unwrap_err(), Error::Invalid);
}

#[test]
fn every_solid_is_one_item_of_one_representation_with_its_own_records() {
    let refused = |items: &str, extra: &str| {
        solids(&part21(items, &format!("{ELLIPTIC}{extra}")))
            .map(|s| s.len())
            .unwrap_err()
    };
    assert_eq!(refused("#5", ""), Error::Unsupported);
    assert_eq!(refused("#61,#61", ""), Error::Unsupported);
    assert_eq!(
        refused("#61,#62", "#62=MANIFOLD_SOLID_BREP('',#60);\n"),
        Error::Invalid
    );
    assert_eq!(solids(&part21("#61,#5", ELLIPTIC)).unwrap().len(), 1);
}

/// A cylinder of radius 2 and height 3 whose rims are closed rational
/// quadratic B-spline circles and whose side is bounded by the two rims
/// alone, without a seam.
fn seamless_cylinder() -> String {
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let square = [
        (2., 0.),
        (2., 2.),
        (0., 2.),
        (-2., 2.),
        (-2., 0.),
        (-2., -2.),
        (0., -2.),
        (2., -2.),
        (2., 0.),
    ];
    let mut body = String::from(
        "#1=CARTESIAN_POINT('',(0.,0.,0.));#2=CARTESIAN_POINT('',(0.,0.,3.));\
         #3=DIRECTION('',(0.,0.,1.));#4=DIRECTION('',(1.,0.,0.));\
         #5=AXIS2_PLACEMENT_3D('',#1,#3,#4);#6=AXIS2_PLACEMENT_3D('',#2,#3,#4);\n",
    );
    for (row, z) in [(100, 0.), (110, 3.)] {
        for (k, (x, y)) in square.iter().enumerate() {
            body += &format!("#{}=CARTESIAN_POINT('',({x:?},{y:?},{z:?}));", row + k + 1);
        }
        let poles: Vec<String> = (1..=9).map(|k| format!("#{}", row + k)).collect();
        body += &format!(
            "#{}=(BOUNDED_CURVE()B_SPLINE_CURVE(2,({}),.CIRCULAR_ARC.,.T.,.F.)\
             B_SPLINE_CURVE_WITH_KNOTS((3,2,2,2,3),(0.,0.25,0.5,0.75,1.),.UNSPECIFIED.)CURVE()\
             GEOMETRIC_REPRESENTATION_ITEM()RATIONAL_B_SPLINE_CURVE((1.,{w:?},1.,{w:?},1.,{w:?},1.,{w:?},1.))\
             REPRESENTATION_ITEM(''));\n",
            row + 20,
            poles.join(",")
        );
    }
    body += "#135=VERTEX_POINT('',#101);#136=VERTEX_POINT('',#111);\
             #140=EDGE_CURVE('',#135,#135,#120,.T.);#141=EDGE_CURVE('',#136,#136,#130,.T.);\
             #150=PLANE('',#5);#151=PLANE('',#6);#152=CYLINDRICAL_SURFACE('',#5,2.);\
             #160=ORIENTED_EDGE('',*,*,#140,.F.);#161=EDGE_LOOP('',(#160));#162=FACE_OUTER_BOUND('',#161,.T.);\
             #163=ADVANCED_FACE('',(#162),#150,.F.);\
             #170=ORIENTED_EDGE('',*,*,#141,.T.);#171=EDGE_LOOP('',(#170));#172=FACE_OUTER_BOUND('',#171,.T.);\
             #173=ADVANCED_FACE('',(#172),#151,.T.);\
             #180=ORIENTED_EDGE('',*,*,#140,.T.);#181=EDGE_LOOP('',(#180));#182=FACE_BOUND('',#181,.T.);\
             #183=ORIENTED_EDGE('',*,*,#141,.F.);#184=EDGE_LOOP('',(#183));#185=FACE_BOUND('',#184,.T.);\
             #186=ADVANCED_FACE('',(#182,#185),#152,.T.);\
             #190=CLOSED_SHELL('',(#163,#173,#186));#191=MANIFOLD_SOLID_BREP('',#190);\n";
    part21("#191", &body)
}

/// A quarter torus elbow, major radius 10 and minor radius 2, closed by
/// two discs on its meridian circles: the toroidal face is a band between
/// two loops that wind along v.
const ELBOW: &str = "\
#1=CARTESIAN_POINT('',(0.,0.,0.));#2=DIRECTION('',(0.,0.,1.));#3=DIRECTION('',(1.,0.,0.));\
#4=AXIS2_PLACEMENT_3D('',#1,#2,#3);#5=TOROIDAL_SURFACE('',#4,10.,2.);\
#10=CARTESIAN_POINT('',(10.,0.,0.));#11=DIRECTION('',(0.,-1.,0.));#12=AXIS2_PLACEMENT_3D('',#10,#11,#3);\
#13=CIRCLE('',#12,2.);#14=PLANE('',#12);#15=CARTESIAN_POINT('',(12.,0.,0.));#16=VERTEX_POINT('',#15);\
#17=EDGE_CURVE('',#16,#16,#13,.T.);\
#20=CARTESIAN_POINT('',(0.,10.,0.));#21=DIRECTION('',(1.,0.,0.));#22=DIRECTION('',(0.,1.,0.));\
#23=AXIS2_PLACEMENT_3D('',#20,#21,#22);#24=CIRCLE('',#23,2.);#25=CARTESIAN_POINT('',(0.,12.,0.));\
#26=VERTEX_POINT('',#25);#27=EDGE_CURVE('',#26,#26,#24,.T.);\
#28=DIRECTION('',(-1.,0.,0.));#29=AXIS2_PLACEMENT_3D('',#20,#28,#22);#30=PLANE('',#29);\
#40=ORIENTED_EDGE('',*,*,#17,.F.);#41=EDGE_LOOP('',(#40));#42=FACE_BOUND('',#41,.T.);\
#43=ORIENTED_EDGE('',*,*,#27,.T.);#44=EDGE_LOOP('',(#43));#45=FACE_BOUND('',#44,.T.);\
#46=ADVANCED_FACE('',(#42,#45),#5,.T.);\
#50=ORIENTED_EDGE('',*,*,#17,.T.);#51=EDGE_LOOP('',(#50));#52=FACE_OUTER_BOUND('',#51,.T.);\
#53=ADVANCED_FACE('',(#52),#14,.T.);\
#60=ORIENTED_EDGE('',*,*,#27,.F.);#61=EDGE_LOOP('',(#60));#62=FACE_OUTER_BOUND('',#61,.T.);\
#63=ADVANCED_FACE('',(#62),#30,.T.);\
#70=CLOSED_SHELL('',(#46,#53,#63));#71=MANIFOLD_SOLID_BREP('',#70);\n";

/// A sphere of radius 3 bounded by one meridian seam from pole to pole.
const MERIDIAN: &str = "\
#1=CARTESIAN_POINT('',(0.,0.,0.));#2=DIRECTION('',(0.,0.,1.));#3=DIRECTION('',(1.,0.,0.));\
#4=AXIS2_PLACEMENT_3D('',#1,#2,#3);#5=SPHERICAL_SURFACE('',#4,3.);\
#6=DIRECTION('',(0.,-1.,0.));#7=AXIS2_PLACEMENT_3D('',#1,#6,#3);#8=CIRCLE('',#7,3.);\
#10=CARTESIAN_POINT('',(0.,0.,-3.));#11=VERTEX_POINT('',#10);\
#12=CARTESIAN_POINT('',(0.,0.,3.));#13=VERTEX_POINT('',#12);\
#14=EDGE_CURVE('',#11,#13,#8,.T.);\
#20=ORIENTED_EDGE('',*,*,#14,.T.);#21=ORIENTED_EDGE('',*,*,#14,.F.);\
#22=EDGE_LOOP('',(#20,#21));#23=FACE_OUTER_BOUND('',#22,.T.);#24=ADVANCED_FACE('',(#23),#5,.T.);\
#30=CLOSED_SHELL('',(#24));#31=MANIFOLD_SOLID_BREP('',#30);\n";

#[test]
fn seamless_bands_and_pole_to_pole_seams_match_their_closed_forms() {
    let cylinder = solid(&seamless_cylinder());
    assert!(
        close(cylinder.volume_mm3(), 12. * PI, 1e-12),
        "{}",
        cylinder.volume_mm3()
    );
    let areas = sorted_areas(&cylinder);
    assert!(
        close(areas[0], 4. * PI, 1e-12) && close(areas[2], 12. * PI, 1e-12),
        "{areas:?}"
    );
    certified(&cylinder, 0.01);
    let elbow = solid(&part21("#71", ELBOW));
    assert!(
        close(elbow.volume_mm3(), 20. * PI * PI, 1e-12),
        "{}",
        elbow.volume_mm3()
    );
    let areas = sorted_areas(&elbow);
    assert!(close(areas[2], 20. * PI * PI, 1e-12), "{areas:?}");
    let [low, high] = elbow.bounds_mm();
    assert!(
        distance(low, [0., 0., -2.]) < 1e-9 && distance(high, [12., 12., 2.]) < 1e-9,
        "{low:?} {high:?}"
    );
    certified(&elbow, 0.01);
    let sphere = solid(&part21("#31", MERIDIAN));
    assert!(
        close(sphere.volume_mm3(), 36. * PI, 1e-12),
        "{}",
        sphere.volume_mm3()
    );
    assert!(close(sphere.faces()[0].area_mm2, 36. * PI, 1e-12));
    certified(&sphere, 0.01);
    for text in [
        seamless_cylinder(),
        part21("#71", ELBOW),
        part21("#31", MERIDIAN),
    ] {
        let source = solid(&text);
        let mut writer = crate::step::Writer::new();
        writer.append_brep(&source).unwrap();
        assert!(close(
            solid(&writer.finish().unwrap()).volume_mm3(),
            source.volume_mm3(),
            1e-12
        ));
    }
}

#[test]
fn bezier_and_rational_supports_and_degree_units_read_as_written() {
    let text = elliptic(
        "#20=PLANE('',#5)",
        "#20=(BOUNDED_SURFACE()B_SPLINE_SURFACE(1,1,((#70,#71),(#72,#73)),.UNSPECIFIED.,.F.,.F.,.F.)B_SPLINE_SURFACE_WITH_KNOTS((2,2),(2,2),(0.,1.),(0.,1.),.UNSPECIFIED.)GEOMETRIC_REPRESENTATION_ITEM()RATIONAL_B_SPLINE_SURFACE(((1.,1.),(1.,1.)))REPRESENTATION_ITEM('')SURFACE());#70=CARTESIAN_POINT('',(-4.,-3.,0.));#71=CARTESIAN_POINT('',(-4.,3.,0.));#72=CARTESIAN_POINT('',(4.,-3.,0.));#73=CARTESIAN_POINT('',(4.,3.,0.))",
    );
    let text = tampered(
        &text,
        "#21=PLANE('',#6)",
        "#21=BEZIER_SURFACE('',1,1,((#74,#75),(#76,#77)),.UNSPECIFIED.,.F.,.F.,.F.);#74=CARTESIAN_POINT('',(-4.,-3.,5.));#75=CARTESIAN_POINT('',(-4.,3.,5.));#76=CARTESIAN_POINT('',(4.,-3.,5.));#77=CARTESIAN_POINT('',(4.,3.,5.))",
    );
    let text = tampered(
        &text,
        "#16=LINE('',#10,#15)",
        "#16=BEZIER_CURVE('',1,(#10,#12),.UNSPECIFIED.,.F.,.F.)",
    );
    let read = solid(&text);
    assert!(
        close(read.volume_mm3(), 30. * PI, 1e-12),
        "{}",
        read.volume_mm3()
    );
    certified(&read, 0.01);
    let mut writer = crate::step::Writer::new();
    writer.append_brep(&read).unwrap();
    assert!(close(
        solid(&writer.finish().unwrap()).volume_mm3(),
        30. * PI,
        1e-12
    ));
    // The frustum in degrees: 0.380506377112 rad is 21.8014094863518 deg.
    let degrees = tampered(
        FRUSTUM,
        "#115 = ( NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.) );",
        "#115 = ( CONVERSION_BASED_UNIT('DEGREE',#119) NAMED_UNIT(*) PLANE_ANGLE_UNIT() );\n#119 = PLANE_ANGLE_MEASURE_WITH_UNIT(PLANE_ANGLE_MEASURE(1.745329251994E-02),#120);\n#120 = ( NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.) );",
    );
    let degrees = tampered(&degrees, "8.,0.380506377112)", "8.,21.8014094863518)");
    assert!(close(
        solid(&degrees).volume_mm3(),
        PI * 10. / 3. * 112.,
        1e-9
    ));
    let gradians = tampered(
        &degrees,
        "PLANE_ANGLE_MEASURE(1.745329251994E-02)",
        "PLANE_ANGLE_MEASURE(1.570796326795E-02)",
    );
    let gradians = tampered(&gradians, "8.,21.8014094863518)", "8.,24.2237883181687)");
    assert!(close(
        solid(&gradians).volume_mm3(),
        PI * 10. / 3. * 112.,
        1e-9
    ));
}

/// Half of a sphere of radius 3, the lune 0 <= u <= pi, closed by the
/// disc in the plane y = 0: its spherical face runs from pole to pole
/// along two meridians.
const LUNE: &str = "\
#1=CARTESIAN_POINT('',(0.,0.,0.));#2=DIRECTION('',(0.,0.,1.));#3=DIRECTION('',(1.,0.,0.));\
#4=AXIS2_PLACEMENT_3D('',#1,#2,#3);#5=SPHERICAL_SURFACE('',#4,3.);\
#6=DIRECTION('',(0.,-1.,0.));#7=AXIS2_PLACEMENT_3D('',#1,#6,#3);#8=CIRCLE('',#7,3.);\
#9=DIRECTION('',(0.,1.,0.));#15=DIRECTION('',(-1.,0.,0.));#16=AXIS2_PLACEMENT_3D('',#1,#9,#15);\
#17=CIRCLE('',#16,3.);#18=PLANE('',#7);\
#10=CARTESIAN_POINT('',(0.,0.,-3.));#11=VERTEX_POINT('',#10);\
#12=CARTESIAN_POINT('',(0.,0.,3.));#13=VERTEX_POINT('',#12);\
#14=EDGE_CURVE('',#11,#13,#8,.T.);#19=EDGE_CURVE('',#11,#13,#17,.T.);\
#20=ORIENTED_EDGE('',*,*,#19,.T.);#21=ORIENTED_EDGE('',*,*,#14,.F.);\
#22=EDGE_LOOP('',(#20,#21));#23=FACE_OUTER_BOUND('',#22,.T.);#24=ADVANCED_FACE('',(#23),#5,.T.);\
#25=ORIENTED_EDGE('',*,*,#14,.T.);#26=ORIENTED_EDGE('',*,*,#19,.F.);\
#27=EDGE_LOOP('',(#25,#26));#28=FACE_OUTER_BOUND('',#27,.T.);#29=ADVANCED_FACE('',(#28),#18,.T.);\
#30=CLOSED_SHELL('',(#24,#29));#31=MANIFOLD_SOLID_BREP('',#30);\n";

#[test]
fn a_lune_between_two_meridians_closes_at_both_poles() {
    let lune = solid(&part21("#31", LUNE));
    assert!(
        close(lune.volume_mm3(), 18. * PI, 1e-12),
        "{}",
        lune.volume_mm3()
    );
    let areas = sorted_areas(&lune);
    assert!(
        close(areas[0], 9. * PI, 1e-12) && close(areas[1], 18. * PI, 1e-12),
        "{areas:?}"
    );
    let [low, high] = lune.bounds_mm();
    assert!(
        distance(low, [-3., 0., -3.]) < 1e-9 && distance(high, [3., 3., 3.]) < 1e-9,
        "{low:?} {high:?}"
    );
    certified(&lune, 0.01);
}

#[test]
fn spline_constructors_refuse_inconsistent_nets_with_their_kind() {
    let line = vec![[0., 0., 0.], [1., 0., 0.]];
    let curve =
        |degree, points: Vec<Point>, weights, knots: Vec<f64>, multiplicities: Vec<usize>| {
            NurbsCurve::new(degree, points, weights, knots, multiplicities).map(|_| ())
        };
    assert_eq!(
        curve(1, line.clone(), None, vec![0., 1.], vec![2, 2]),
        Ok(())
    );
    assert_eq!(
        curve(0, line.clone(), None, vec![0., 1.], vec![1, 1]),
        Err(Error::Invalid)
    );
    assert_eq!(
        curve(1, line.clone(), None, vec![0., 1.], vec![2]),
        Err(Error::Invalid)
    );
    assert_eq!(
        curve(2, line.clone(), None, vec![0., 1.], vec![3, 3]),
        Err(Error::Limit)
    );
    assert_eq!(
        curve(1, line.clone(), None, vec![1., 0.], vec![2, 2]),
        Err(Error::Invalid)
    );
    assert_eq!(
        curve(1, line.clone(), None, vec![0., 1.], vec![3, 2]),
        Err(Error::Invalid)
    );
    assert_eq!(
        curve(1, line.clone(), None, vec![0., 1., 2.], vec![2, 1, 2]),
        Err(Error::Invalid)
    );
    assert_eq!(
        curve(1, line.clone(), Some(vec![1.]), vec![0., 1.], vec![2, 2]),
        Err(Error::Invalid)
    );
    assert_eq!(
        curve(
            1,
            line.clone(),
            Some(vec![1., 0.]),
            vec![0., 1.],
            vec![2, 2]
        ),
        Err(Error::Invalid)
    );
    let surface = |counts, points: Vec<Point>| {
        NurbsSurface::new(
            [1, 1],
            counts,
            points,
            None,
            [vec![0., 1.], vec![0., 1.]],
            [vec![2, 2], vec![2, 2]],
        )
        .map(|_| ())
    };
    let square = vec![[0., 0., 0.], [0., 1., 0.], [1., 0., 0.], [1., 1., 0.]];
    assert_eq!(surface([2, 2], square.clone()), Ok(()));
    assert_eq!(surface([2, 3], square.clone()), Err(Error::Invalid));
    assert_eq!(
        surface([300, 300], vec![[0.; 3]; 90_000]),
        Err(Error::Limit)
    );
    // A degree above the admitted maximum.
    assert_eq!(
        curve(
            nurbs::MAX_DEGREE + 1,
            vec![[0.; 3]; 40],
            None,
            vec![0., 1.],
            vec![34, 34]
        ),
        Err(Error::Invalid)
    );
}

#[test]
fn every_support_family_moves_rigidly() {
    // A quarter turn about z and a shift: volume and areas are kept, and
    // the moved boundary validates and integrates to the same volume.
    let linear = |p: Point| [-p[1], p[0], p[2]];
    let map = move |p: Point| add(linear(p), [7., -3., 2.]);
    let texts = [
        FRUSTUM.to_owned(),
        SPHERE.to_owned(),
        TORUS.to_owned(),
        VASE.to_owned(),
        SPLINE.to_owned(),
        WEDGES.to_owned(),
        DISC.to_owned(),
        seamless_cylinder(),
        part21("#71", ELBOW),
        part21("#61", ELLIPTIC),
    ];
    for text in texts {
        for source in solids(&text).unwrap() {
            let moved = source.moved(&map, &linear).unwrap();
            let again = Solid::new(moved.boundary().clone(), source.tolerance(), None).unwrap();
            assert!(
                close(again.volume_mm3(), source.volume_mm3(), 1e-10),
                "{}",
                source.volume_mm3()
            );
            let [low, high] = source.bounds_mm();
            let [a, b] = [map(low), map(high)];
            let expected = [
                [a[0].min(b[0]), a[1].min(b[1]), a[2]],
                [a[0].max(b[0]), a[1].max(b[1]), b[2]],
            ];
            let [ml, mh] = again.bounds_mm();
            assert!(
                distance(ml, expected[0]) < 1e-7 && distance(mh, expected[1]) < 1e-7,
                "{ml:?} {expected:?}"
            );
        }
    }
}

#[test]
fn surface_inversion_handles_poles_apexes_and_both_nappes() {
    let frame = Frame {
        origin: [0.; 3],
        x: [1., 0., 0.],
        y: [0., 1., 0.],
        z: [0., 0., 1.],
    };
    let sphere = Surface::Sphere { frame, radius: 2. };
    // At a pole u is free: the seed decides it.
    assert_eq!(
        sphere.invert([0., 0., 2.], Some([1.25, 1.5])),
        [1.25, std::f64::consts::FRAC_PI_2]
    );
    let cone = Surface::Cone {
        frame,
        radius: 1.,
        half_angle: std::f64::consts::FRAC_PI_4,
    };
    // Below the apex the cone's radius is negative and u turns by pi.
    let below = cone.evaluate([0.5, -3.]).0;
    let uv = cone.invert(below, None);
    assert!(
        distance(cone.evaluate(uv).0, below) < 1e-12 && (uv[0] - 0.5).abs() < 1e-12,
        "{uv:?}"
    );
    assert_eq!(cone.pole(true, 0.), None);
    assert!(cone.pole(false, 0.).is_some_and(|v| (v + 1.).abs() < 1e-12));
    let plane = Surface::Plane(frame);
    assert_eq!(plane.bounded(), [None, None]);
    assert!(plane.candidates([1., 2., 3.]).is_empty());
    assert_eq!(plane.complete(), None);
}

/// The upper half of a sphere of radius 3 closed by the disc z = 0: its
/// spherical face is bounded by the equator alone and closes at the pole.
const DOME: &str = "\
#1=CARTESIAN_POINT('',(0.,0.,0.));#2=DIRECTION('',(0.,0.,1.));#3=DIRECTION('',(1.,0.,0.));\
#4=AXIS2_PLACEMENT_3D('',#1,#2,#3);#5=SPHERICAL_SURFACE('',#4,3.);#6=CIRCLE('',#4,3.);\
#7=PLANE('',#4);#8=CARTESIAN_POINT('',(3.,0.,0.));#9=VERTEX_POINT('',#8);\
#10=EDGE_CURVE('',#9,#9,#6,.T.);\
#11=ORIENTED_EDGE('',*,*,#10,.T.);#12=EDGE_LOOP('',(#11));#13=FACE_BOUND('',#12,.T.);\
#14=ADVANCED_FACE('',(#13),#5,.T.);\
#15=ORIENTED_EDGE('',*,*,#10,.F.);#16=EDGE_LOOP('',(#15));#17=FACE_OUTER_BOUND('',#16,.T.);\
#18=ADVANCED_FACE('',(#17),#7,.F.);\
#19=CLOSED_SHELL('',(#14,#18));#20=MANIFOLD_SOLID_BREP('',#19);\n";

#[test]
fn caps_off_seam_meridians_and_off_surface_rims() {
    let dome = solid(&part21("#20", DOME));
    assert!(
        close(dome.volume_mm3(), 18. * PI, 1e-12),
        "{}",
        dome.volume_mm3()
    );
    let areas = sorted_areas(&dome);
    assert!(
        close(areas[0], 9. * PI, 1e-12) && close(areas[1], 18. * PI, 1e-12),
        "{areas:?}"
    );
    assert!(distance(dome.bounds_mm()[1], [3., 3., 3.]) < 1e-9);
    // A cap closed at its pole by a loop without a seam edge meshes
    // closed, every vertex on its sphere or its plane.
    for deflection in [0.05, 0.01] {
        let mesh = certified(&dome, deflection);
        assert!(
            mesh.positions
                .iter()
                .all(|p| (norm(*p) - 3.).abs() < 1e-9 || p[2].abs() < 1e-12)
        );
    }
    // The meridian away from the sphere's seam: up one side and back.
    let turned = tampered(
        &part21("#31", MERIDIAN),
        "#4=AXIS2_PLACEMENT_3D('',#1,#2,#3)",
        "#4=AXIS2_PLACEMENT_3D('',#1,#2,#40);#40=DIRECTION('',(0.,1.,0.))",
    );
    let sphere = solid(&turned);
    assert!(
        close(sphere.volume_mm3(), 36. * PI, 1e-12),
        "{}",
        sphere.volume_mm3()
    );
    certified(&sphere, 0.01);
    // Rims 1 um off their cylinder do not bound it.
    let off = tampered(
        &seamless_cylinder(),
        "CYLINDRICAL_SURFACE('',#5,2.)",
        "CYLINDRICAL_SURFACE('',#5,2.000001)",
    );
    assert_eq!(solids(&off).unwrap_err(), Error::Invalid);
}

#[test]
fn open_conics_evaluate_invert_and_move_exactly() {
    let frame = Frame {
        origin: [1., 2., 3.],
        x: [0., 1., 0.],
        y: [-1., 0., 0.],
        z: [0., 0., 1.],
    };
    let hyperbola = Curve::Hyperbola {
        frame,
        semiaxes: [2., 3.],
    };
    let parabola = Curve::Parabola { frame, focal: 0.5 };
    for curve in [&hyperbola, &parabola] {
        assert_eq!(curve.period(), None);
        for t in [-1.5, -0.2, 0., 0.7, 2.] {
            let (p, d) = curve.evaluate(t);
            assert!((curve.invert(p, None) - t).abs() < 1e-12, "{t}");
            let h = 1e-6;
            let difference = scale(
                sub(curve.evaluate(t + h).0, curve.evaluate(t - h).0),
                0.5 / h,
            );
            assert!(distance(difference, d) < 1e-6, "{t}");
        }
        let breaks = curve.breaks(-1., 1.);
        assert_eq!(breaks.len(), 9);
        let moved = curve.mapped(&|p| add(p, [1., 0., 0.]), &|v| v);
        assert!(
            distance(
                moved.evaluate(0.3).0,
                add(curve.evaluate(0.3).0, [1., 0., 0.])
            ) < 1e-12
        );
        let samples = curve.sample(-1., 1., 1e-3, 0.3).unwrap();
        assert!(samples.len() > 9 && samples[0] == -1. && *samples.last().unwrap() == 1.);
    }
    // a cosh t along x, b sinh t along y of the frame.
    assert!(distance(hyperbola.evaluate(0.).0, [1., 4., 3.]) < 1e-12);
    assert!(distance(parabola.evaluate(1.).0, [0., 2.5, 3.]) < 1e-12);
}

/// `a` and `b` in one boundary: `b`'s shells beside `a`'s, or its faces
/// joined into `a`'s exterior shell.
fn beside(a: &Boundary, b: &Boundary, joined: bool) -> Boundary {
    let (vertices, edges, faces) = (a.vertices.len(), a.edges.len(), a.faces.len());
    let mut out = a.clone();
    out.vertices.extend(&b.vertices);
    out.edges.extend(b.edges.iter().map(|edge| Edge {
        vertices: edge.vertices.map(|v| v + vertices),
        ..edge.clone()
    }));
    out.faces.extend(b.faces.iter().map(|face| {
        Face {
            loops: face
                .loops
                .iter()
                .map(|trim| Loop {
                    uses: trim.uses.iter().map(|&(e, f)| (e + edges, f)).collect(),
                    vertex: trim.vertex.map(|v| v + vertices),
                })
                .collect(),
            ..face.clone()
        }
    }));
    let shells = b.shells.iter().map(|shell| Shell {
        faces: shell.faces.iter().map(|f| f + faces).collect(),
        reversed: shell.reversed,
    });
    if joined {
        let joined: Vec<usize> = shells.flat_map(|shell| shell.faces).collect();
        out.shells[0].faces.extend(joined);
    } else {
        out.shells.extend(shells);
    }
    out
}

#[test]
fn constructed_boundaries_with_shared_or_stray_records_are_invalid() {
    let elliptic = read_boundaries(&part21("#61", ELLIPTIC), TOLERANCE)
        .unwrap()
        .remove(0);
    let shifted = elliptic.mapped(&|p| add(p, [10., 0., 0.]), &|v| v);
    let refused = |boundary: Boundary| Solid::new(boundary, TOLERANCE, None).unwrap_err();
    // Two cylinders side by side are two solids, not one shell.
    assert_eq!(refused(beside(&elliptic, &shifted, true)), Error::Invalid);
    // One face claimed by two shells.
    let mut twice = elliptic.clone();
    twice.shells.push(Shell {
        faces: vec![0],
        reversed: false,
    });
    assert_eq!(refused(twice), Error::Invalid);
    // A vertex loop naming no vertex.
    let mut stray = elliptic.clone();
    stray.faces[0].loops.push(Loop {
        uses: vec![],
        vertex: Some(99),
    });
    assert_eq!(refused(stray), Error::Invalid);
    // Both bottom discs as one planar face with two outer loops: a face is
    // one connected region.
    let mut discs = beside(&elliptic, &shifted, true);
    let moved = discs.faces.remove(3);
    discs.faces[0].loops.extend(moved.loops);
    discs.shells[0].faces.retain(|&face| face != 3);
    for face in &mut discs.shells[0].faces {
        *face -= usize::from(*face > 3);
    }
    assert_eq!(refused(discs), Error::Invalid);
}

/// A prism 5 mm tall over the region between an open conic from `a` to
/// `b` and the chord back to `a`. `conic` names the curve with `{}` for
/// its placement; the curve opens along +x, so its own normal points into
/// the region.
fn conic_prism(conic: &str, a: Point, b: Point) -> String {
    let at = |p: Point, z: f64| format!("({:?},{:?},{z:?})", p[0], p[1]);
    let bottom = conic.replace("{}", "#5");
    let top = conic.replace("{}", "#6");
    let body = format!(
        "#1=CARTESIAN_POINT('',(0.,0.,0.));#2=CARTESIAN_POINT('',(0.,0.,5.));\
         #3=DIRECTION('',(0.,0.,1.));#4=DIRECTION('',(1.,0.,0.));\
         #5=AXIS2_PLACEMENT_3D('',#1,#3,#4);#6=AXIS2_PLACEMENT_3D('',#2,#3,#4);\
         #7={bottom};#8={top};\
         #10=CARTESIAN_POINT('',{});#11=CARTESIAN_POINT('',{});\
         #12=CARTESIAN_POINT('',{});#13=CARTESIAN_POINT('',{});\
         #14=VERTEX_POINT('',#10);#15=VERTEX_POINT('',#11);\
         #16=VERTEX_POINT('',#12);#17=VERTEX_POINT('',#13);\
         #18=DIRECTION('',(0.,-1.,0.));#19=VECTOR('',#18,1.);#20=VECTOR('',#3,1.);\
         #21=EDGE_CURVE('',#14,#15,#7,.T.);#22=LINE('',#11,#19);#23=EDGE_CURVE('',#15,#14,#22,.T.);\
         #24=EDGE_CURVE('',#16,#17,#8,.T.);#25=LINE('',#13,#19);#26=EDGE_CURVE('',#17,#16,#25,.T.);\
         #27=LINE('',#10,#20);#28=EDGE_CURVE('',#14,#16,#27,.T.);\
         #29=LINE('',#11,#20);#30=EDGE_CURVE('',#15,#17,#29,.T.);\
         #31=PLANE('',#5);#32=PLANE('',#6);#33=SURFACE_OF_LINEAR_EXTRUSION('',#7,#20);\
         #34=AXIS2_PLACEMENT_3D('',#11,#4,#18);#35=PLANE('',#34);\
         #40=ORIENTED_EDGE('',*,*,#21,.T.);#41=ORIENTED_EDGE('',*,*,#23,.T.);\
         #42=EDGE_LOOP('',(#40,#41));#43=FACE_OUTER_BOUND('',#42,.T.);#44=ADVANCED_FACE('',(#43),#31,.F.);\
         #45=ORIENTED_EDGE('',*,*,#26,.F.);#46=ORIENTED_EDGE('',*,*,#24,.F.);\
         #47=EDGE_LOOP('',(#45,#46));#48=FACE_OUTER_BOUND('',#47,.T.);#49=ADVANCED_FACE('',(#48),#32,.T.);\
         #50=ORIENTED_EDGE('',*,*,#21,.F.);#51=ORIENTED_EDGE('',*,*,#28,.T.);\
         #52=ORIENTED_EDGE('',*,*,#24,.T.);#53=ORIENTED_EDGE('',*,*,#30,.F.);\
         #54=EDGE_LOOP('',(#50,#51,#52,#53));#55=FACE_BOUND('',#54,.T.);#56=ADVANCED_FACE('',(#55),#33,.F.);\
         #57=ORIENTED_EDGE('',*,*,#23,.F.);#58=ORIENTED_EDGE('',*,*,#30,.T.);\
         #59=ORIENTED_EDGE('',*,*,#26,.T.);#60=ORIENTED_EDGE('',*,*,#28,.F.);\
         #61=EDGE_LOOP('',(#57,#58,#59,#60));#62=FACE_BOUND('',#61,.T.);#63=ADVANCED_FACE('',(#62),#35,.T.);\
         #70=CLOSED_SHELL('',(#44,#49,#56,#63));#71=MANIFOLD_SOLID_BREP('',#70);\n",
        at(a, 0.),
        at(b, 0.),
        at(a, 5.),
        at(b, 5.),
    );
    part21("#71", &body)
}

#[test]
fn parabolic_and_hyperbolic_prisms_match_their_closed_forms_and_round_trip() {
    // x = t^2, y = 2 t for -1 <= t <= 1 against the chord x = 1: the
    // region holds 8/3 mm^2.
    let parabola = conic_prism("PARABOLA('',{},1.)", [1., -2., 0.], [1., 2., 0.]);
    // x = cosh t, y = sinh t: the region holds sinh(2) / 2 - 1.
    let (c, h) = (1f64.cosh(), 1f64.sinh());
    let hyperbola = conic_prism("HYPERBOLA('',{},1.,1.)", [c, -h, 0.], [c, h, 0.]);
    for (text, cap) in [(parabola, 8. / 3.), (hyperbola, 2f64.sinh() / 2. - 1.)] {
        let prism = solid(&text);
        assert!(
            close(prism.volume_mm3(), 5. * cap, 1e-12),
            "{} {}",
            prism.volume_mm3(),
            5. * cap
        );
        let areas = sorted_areas(&prism);
        assert!(close(areas[0], cap, 1e-12) && close(areas[1], cap, 1e-12));
        certified(&prism, 0.01);
        let mut writer = crate::step::Writer::new();
        writer.append_brep(&prism).unwrap();
        let back = solid(&writer.finish().unwrap());
        assert!(close(back.volume_mm3(), prism.volume_mm3(), 1e-12));
        assert_eq!(back.boundary().edges, prism.boundary().edges);
    }
}

/// The cone of half-angle 30 degrees over the equator disc of `DOME`: its
/// side, 18 pi mm^2 like the hemisphere's, closes at the apex without a
/// vertex.
fn comet() -> String {
    tampered(
        &part21("#20", DOME),
        "#5=SPHERICAL_SURFACE('',#4,3.)",
        "#5=CONICAL_SURFACE('',#41,3.,0.5235987755982988);\
         #41=AXIS2_PLACEMENT_3D('',#1,#42,#3);#42=DIRECTION('',(0.,0.,-1.))",
    )
}

#[test]
fn seam_free_loops_close_at_apexes_and_split_seams_cover_the_torus() {
    // The cone is 27^0.5 mm tall.
    let cone = solid(&comet());
    let height = 27f64.sqrt();
    assert!(
        close(cone.volume_mm3(), 3. * PI * height, 1e-12),
        "{}",
        cone.volume_mm3()
    );
    let areas = sorted_areas(&cone);
    assert!(close(areas[0], 9. * PI, 1e-12) && close(areas[1], 18. * PI, 1e-12));
    assert!((cone.bounds_mm()[1][2] - height).abs() < 1e-9);
    for deflection in [0.05, 0.01] {
        certified(&cone, deflection);
    }
    // The same cone with its apex stated as a vertex bounding the side.
    let apex = tampered(
        &comet(),
        "#14=ADVANCED_FACE('',(#13),#5,.T.)",
        &format!(
            "#14=ADVANCED_FACE('',(#13,#45),#5,.T.);#43=CARTESIAN_POINT('',(0.,0.,{height:?}));\
             #44=VERTEX_POINT('',#43);#46=VERTEX_LOOP('',#44);#45=FACE_BOUND('',#46,.T.)"
        ),
    );
    let pointed = solid(&apex);
    assert!(close(pointed.volume_mm3(), cone.volume_mm3(), 1e-12));
    certified(&pointed, 0.01);
    // The torus with its two seams as two loops, each running along both
    // seams once: neither loop is anchored, and together they bound the
    // whole surface.
    let split = tampered(
        TORUS,
        "#17 = ADVANCED_FACE('',(#18),#31,.T.);",
        "#17 = ADVANCED_FACE('',(#18,#90),#31,.T.);\n\
         #90 = FACE_BOUND('',#91,.T.);\n#91 = EDGE_LOOP('',(#71,#72));",
    );
    let split = tampered(
        &split,
        "#19 = EDGE_LOOP('',(#20,#49,#71,#72));",
        "#19 = EDGE_LOOP('',(#20,#49));",
    );
    let torus = solid(&split);
    assert!(close(torus.volume_mm3(), 2. * PI * PI * 12. * 9., 1e-12));
    certified(&torus, 0.05);
}

/// A sphere of radius 3 bounded by its pole-to-pole seam, less the cap
/// beyond the plane y = -2.6, which a disc closes.
const NICKED: &str = "\
#1=CARTESIAN_POINT('',(0.,0.,0.));#2=DIRECTION('',(0.,0.,1.));#3=DIRECTION('',(1.,0.,0.));\
#4=AXIS2_PLACEMENT_3D('',#1,#2,#3);#5=SPHERICAL_SURFACE('',#4,3.);\
#6=DIRECTION('',(0.,-1.,0.));#7=AXIS2_PLACEMENT_3D('',#1,#6,#3);#8=CIRCLE('',#7,3.);\
#10=CARTESIAN_POINT('',(0.,0.,-3.));#11=VERTEX_POINT('',#10);\
#12=CARTESIAN_POINT('',(0.,0.,3.));#13=VERTEX_POINT('',#12);\
#14=EDGE_CURVE('',#11,#13,#8,.T.);\
#15=CARTESIAN_POINT('',(0.,-2.6,0.));#16=AXIS2_PLACEMENT_3D('',#15,#6,#3);\
#17=CIRCLE('',#16,1.4966629547095765);#18=PLANE('',#16);\
#19=CARTESIAN_POINT('',(1.4966629547095765,-2.6,0.));#9=VERTEX_POINT('',#19);\
#25=EDGE_CURVE('',#9,#9,#17,.T.);\
#20=ORIENTED_EDGE('',*,*,#14,.T.);#21=ORIENTED_EDGE('',*,*,#14,.F.);\
#22=EDGE_LOOP('',(#20,#21));#23=FACE_OUTER_BOUND('',#22,.T.);\
#26=ORIENTED_EDGE('',*,*,#25,.T.);#27=EDGE_LOOP('',(#26));#28=FACE_BOUND('',#27,.T.);\
#24=ADVANCED_FACE('',(#23,#28),#5,.T.);\
#29=ORIENTED_EDGE('',*,*,#25,.F.);#32=EDGE_LOOP('',(#29));#33=FACE_OUTER_BOUND('',#32,.T.);\
#34=ADVANCED_FACE('',(#33),#18,.F.);\
#30=CLOSED_SHELL('',(#24,#34));#31=MANIFOLD_SOLID_BREP('',#30);\n";

#[test]
fn faces_whose_loops_do_not_place_their_region_are_refused() {
    // A seam loop from pole to pole and back does not say where the hole
    // beside it cuts the sphere's parameter plane.
    assert_eq!(
        solids(&part21("#31", NICKED)).unwrap_err(),
        Error::Unsupported
    );
    // A reversed sphere bounded only by its pole could be either side.
    assert_eq!(
        solids(&tampered(
            SPHERE,
            "#17 = ADVANCED_FACE('',(#18),#22,.T.)",
            "#17 = ADVANCED_FACE('',(#18),#22,.F.)"
        ))
        .unwrap_err(),
        Error::Unsupported
    );
    // Both rims of the seamless side wind the same way once the top disc
    // faces into the cylinder: the band between them is not a region.
    let crossed = tampered(
        &seamless_cylinder(),
        "#185=FACE_BOUND('',#184,.T.)",
        "#185=FACE_BOUND('',#184,.F.)",
    );
    let crossed = tampered(
        &crossed,
        "#173=ADVANCED_FACE('',(#172),#151,.T.)",
        "#173=ADVANCED_FACE('',(#172),#151,.F.)",
    );
    assert_eq!(solids(&crossed).unwrap_err(), Error::Invalid);
}

/// A 10 mm cube whose front top edge bends up 1e-9 mm at its middle, at
/// (5, 0, 10 + 1e-9): the triangle from the bend down to the straight top
/// edge is a fin on the plane y = 0 that encloses 5e-9 mm².
const SLIVER: &str = "\
#1=CARTESIAN_POINT('',(0.0,0.0,0.0));#2=VERTEX_POINT('',#1);\
#3=CARTESIAN_POINT('',(10.0,0.0,0.0));#4=VERTEX_POINT('',#3);\
#5=CARTESIAN_POINT('',(10.0,10.0,0.0));#6=VERTEX_POINT('',#5);\
#7=CARTESIAN_POINT('',(0.0,10.0,0.0));#8=VERTEX_POINT('',#7);\
#9=CARTESIAN_POINT('',(0.0,0.0,10.0));#10=VERTEX_POINT('',#9);\
#11=CARTESIAN_POINT('',(10.0,0.0,10.0));#12=VERTEX_POINT('',#11);\
#13=CARTESIAN_POINT('',(10.0,10.0,10.0));#14=VERTEX_POINT('',#13);\
#15=CARTESIAN_POINT('',(0.0,10.0,10.0));#16=VERTEX_POINT('',#15);\
#17=CARTESIAN_POINT('',(5.0,0.0,10.000000001));#18=VERTEX_POINT('',#17);\
#19=CARTESIAN_POINT('',(0.0,0.0,0.0));#20=DIRECTION('',(1.0,0.0,0.0));\
#21=VECTOR('',#20,10.0);#22=LINE('',#19,#21);#23=EDGE_CURVE('',#2,#4,#22,.T.);\
#24=CARTESIAN_POINT('',(10.0,0.0,0.0));#25=DIRECTION('',(0.0,1.0,0.0));\
#26=VECTOR('',#25,10.0);#27=LINE('',#24,#26);#28=EDGE_CURVE('',#4,#6,#27,.T.);\
#29=CARTESIAN_POINT('',(10.0,10.0,0.0));#30=DIRECTION('',(-1.0,0.0,0.0));\
#31=VECTOR('',#30,10.0);#32=LINE('',#29,#31);#33=EDGE_CURVE('',#6,#8,#32,.T.);\
#34=CARTESIAN_POINT('',(0.0,10.0,0.0));#35=DIRECTION('',(0.0,-1.0,0.0));\
#36=VECTOR('',#35,10.0);#37=LINE('',#34,#36);#38=EDGE_CURVE('',#8,#2,#37,.T.);\
#39=CARTESIAN_POINT('',(0.0,0.0,10.0));#40=DIRECTION('',(1.0,0.0,0.0));\
#41=VECTOR('',#40,10.0);#42=LINE('',#39,#41);#43=EDGE_CURVE('',#10,#12,#42,.T.);\
#44=CARTESIAN_POINT('',(10.0,0.0,10.0));#45=DIRECTION('',(0.0,1.0,0.0));\
#46=VECTOR('',#45,10.0);#47=LINE('',#44,#46);#48=EDGE_CURVE('',#12,#14,#47,.T.);\
#49=CARTESIAN_POINT('',(10.0,10.0,10.0));#50=DIRECTION('',(-1.0,0.0,0.0));\
#51=VECTOR('',#50,10.0);#52=LINE('',#49,#51);#53=EDGE_CURVE('',#14,#16,#52,.T.);\
#54=CARTESIAN_POINT('',(0.0,10.0,10.0));#55=DIRECTION('',(0.0,-1.0,0.0));\
#56=VECTOR('',#55,10.0);#57=LINE('',#54,#56);#58=EDGE_CURVE('',#16,#10,#57,.T.);\
#59=CARTESIAN_POINT('',(0.0,0.0,0.0));#60=DIRECTION('',(0.0,0.0,1.0));\
#61=VECTOR('',#60,10.0);#62=LINE('',#59,#61);#63=EDGE_CURVE('',#2,#10,#62,.T.);\
#64=CARTESIAN_POINT('',(10.0,0.0,0.0));#65=DIRECTION('',(0.0,0.0,1.0));\
#66=VECTOR('',#65,10.0);#67=LINE('',#64,#66);#68=EDGE_CURVE('',#4,#12,#67,.T.);\
#69=CARTESIAN_POINT('',(10.0,10.0,0.0));#70=DIRECTION('',(0.0,0.0,1.0));\
#71=VECTOR('',#70,10.0);#72=LINE('',#69,#71);#73=EDGE_CURVE('',#6,#14,#72,.T.);\
#74=CARTESIAN_POINT('',(0.0,10.0,0.0));#75=DIRECTION('',(0.0,0.0,1.0));\
#76=VECTOR('',#75,10.0);#77=LINE('',#74,#76);#78=EDGE_CURVE('',#8,#16,#77,.T.);\
#79=CARTESIAN_POINT('',(10.0,0.0,10.0));#80=DIRECTION('',(-1.0,0.0,2.000000165480742E-10));\
#81=VECTOR('',#80,5.0);#82=LINE('',#79,#81);#83=EDGE_CURVE('',#12,#18,#82,.T.);\
#84=CARTESIAN_POINT('',(5.0,0.0,10.000000001));\
#85=DIRECTION('',(-1.0,0.0,-2.000000165480742E-10));#86=VECTOR('',#85,5.0);\
#87=LINE('',#84,#86);#88=EDGE_CURVE('',#18,#10,#87,.T.);#89=ORIENTED_EDGE('',*,*,#38,.F.);\
#90=ORIENTED_EDGE('',*,*,#33,.F.);#91=ORIENTED_EDGE('',*,*,#28,.F.);\
#92=ORIENTED_EDGE('',*,*,#23,.F.);#93=EDGE_LOOP('',(#89,#90,#91,#92));\
#94=FACE_OUTER_BOUND('',#93,.T.);#95=CARTESIAN_POINT('',(0.0,0.0,0.0));\
#96=DIRECTION('',(0.0,0.0,-1.0));#97=DIRECTION('',(1.0,0.0,0.0));\
#98=AXIS2_PLACEMENT_3D('',#95,#96,#97);#99=PLANE('',#98);\
#100=ADVANCED_FACE('',(#94),#99,.T.);#101=ORIENTED_EDGE('',*,*,#43,.T.);\
#102=ORIENTED_EDGE('',*,*,#48,.T.);#103=ORIENTED_EDGE('',*,*,#53,.T.);\
#104=ORIENTED_EDGE('',*,*,#58,.T.);#105=EDGE_LOOP('',(#101,#102,#103,#104));\
#106=FACE_OUTER_BOUND('',#105,.T.);#107=CARTESIAN_POINT('',(0.0,0.0,10.0));\
#108=DIRECTION('',(0.0,0.0,1.0));#109=DIRECTION('',(1.0,0.0,0.0));\
#110=AXIS2_PLACEMENT_3D('',#107,#108,#109);#111=PLANE('',#110);\
#112=ADVANCED_FACE('',(#106),#111,.T.);#113=ORIENTED_EDGE('',*,*,#23,.T.);\
#114=ORIENTED_EDGE('',*,*,#68,.T.);#115=ORIENTED_EDGE('',*,*,#83,.T.);\
#116=ORIENTED_EDGE('',*,*,#88,.T.);#117=ORIENTED_EDGE('',*,*,#63,.F.);\
#118=EDGE_LOOP('',(#113,#114,#115,#116,#117));#119=FACE_OUTER_BOUND('',#118,.T.);\
#120=CARTESIAN_POINT('',(0.0,0.0,0.0));#121=DIRECTION('',(0.0,-1.0,0.0));\
#122=DIRECTION('',(1.0,0.0,0.0));#123=AXIS2_PLACEMENT_3D('',#120,#121,#122);\
#124=PLANE('',#123);#125=ADVANCED_FACE('',(#119),#124,.T.);\
#126=ORIENTED_EDGE('',*,*,#28,.T.);#127=ORIENTED_EDGE('',*,*,#73,.T.);\
#128=ORIENTED_EDGE('',*,*,#48,.F.);#129=ORIENTED_EDGE('',*,*,#68,.F.);\
#130=EDGE_LOOP('',(#126,#127,#128,#129));#131=FACE_OUTER_BOUND('',#130,.T.);\
#132=CARTESIAN_POINT('',(10.0,0.0,0.0));#133=DIRECTION('',(1.0,0.0,0.0));\
#134=DIRECTION('',(0.0,1.0,0.0));#135=AXIS2_PLACEMENT_3D('',#132,#133,#134);\
#136=PLANE('',#135);#137=ADVANCED_FACE('',(#131),#136,.T.);\
#138=ORIENTED_EDGE('',*,*,#33,.T.);#139=ORIENTED_EDGE('',*,*,#78,.T.);\
#140=ORIENTED_EDGE('',*,*,#53,.F.);#141=ORIENTED_EDGE('',*,*,#73,.F.);\
#142=EDGE_LOOP('',(#138,#139,#140,#141));#143=FACE_OUTER_BOUND('',#142,.T.);\
#144=CARTESIAN_POINT('',(0.0,10.0,0.0));#145=DIRECTION('',(0.0,1.0,0.0));\
#146=DIRECTION('',(1.0,0.0,0.0));#147=AXIS2_PLACEMENT_3D('',#144,#145,#146);\
#148=PLANE('',#147);#149=ADVANCED_FACE('',(#143),#148,.T.);\
#150=ORIENTED_EDGE('',*,*,#38,.T.);#151=ORIENTED_EDGE('',*,*,#63,.T.);\
#152=ORIENTED_EDGE('',*,*,#58,.F.);#153=ORIENTED_EDGE('',*,*,#78,.F.);\
#154=EDGE_LOOP('',(#150,#151,#152,#153));#155=FACE_OUTER_BOUND('',#154,.T.);\
#156=CARTESIAN_POINT('',(0.0,0.0,0.0));#157=DIRECTION('',(-1.0,0.0,0.0));\
#158=DIRECTION('',(0.0,1.0,0.0));#159=AXIS2_PLACEMENT_3D('',#156,#157,#158);\
#160=PLANE('',#159);#161=ADVANCED_FACE('',(#155),#160,.T.);\
#162=ORIENTED_EDGE('',*,*,#43,.F.);#163=ORIENTED_EDGE('',*,*,#88,.F.);\
#164=ORIENTED_EDGE('',*,*,#83,.F.);#165=EDGE_LOOP('',(#162,#163,#164));\
#166=FACE_OUTER_BOUND('',#165,.T.);#167=CARTESIAN_POINT('',(0.0,0.0,10.0));\
#168=DIRECTION('',(0.0,1.0,0.0));#169=DIRECTION('',(1.0,0.0,0.0));\
#170=AXIS2_PLACEMENT_3D('',#167,#168,#169);#171=PLANE('',#170);\
#172=ADVANCED_FACE('',(#166),#171,.T.);\
#173=CLOSED_SHELL('',(#100,#112,#125,#137,#149,#161,#172));\
#174=MANIFOLD_SOLID_BREP('',#173);\n";

#[test]
fn faces_narrower_than_the_tolerance_never_turn_on_their_rounding() {
    // Stated in agreement with its loop, the fin reads as one more face.
    let cube = solid(&part21("#174", SLIVER));
    assert!(close(cube.volume_mm3(), 1000., 1e-9));
    // Stated against its loop by an area within the tolerance, the fin's
    // direction is rounding: it is not turned to fit, nor refused as
    // invalid, but cannot be placed.
    let against = tampered(
        SLIVER,
        "#172=ADVANCED_FACE('',(#166),#171,.T.)",
        "#172=ADVANCED_FACE('',(#166),#171,.F.)",
    );
    assert_eq!(
        solids(&part21("#174", &against)).unwrap_err(),
        Error::Unsupported
    );
    // With its bend flattened, the fin bounds nothing at all.
    let flat = SLIVER
        .replace("(5.0,0.0,10.000000001)", "(5.0,0.0,10.0)")
        .replace("2.000000165480742E-10", "0.0");
    assert_eq!(
        solids(&part21("#174", &flat)).unwrap_err(),
        Error::Unsupported
    );
}

/// A crescent prism 1 mm tall between two arcs over the chord from
/// (-1, 0) to (1, 0): a radius 6 circle, sampled at its ends alone, and a
/// flatter radius 6.66 circle, of 0.9 times the first one's
/// height, written as a four-span rational spline,
/// sampled at its knots. The chords of the crescent enclose a negative
/// area, its arcs 0.0112 mm2.
const CRESCENT: &str = "\
#1=CARTESIAN_POINT('',(-1.0,0.0,0.0));#2=VERTEX_POINT('',#1);\
#3=CARTESIAN_POINT('',(1.0,0.0,0.0));#4=VERTEX_POINT('',#3);\
#5=CARTESIAN_POINT('',(-1.0,0.0,1.0));#6=VERTEX_POINT('',#5);\
#7=CARTESIAN_POINT('',(1.0,0.0,1.0));#8=VERTEX_POINT('',#7);\
#9=CARTESIAN_POINT('',(0.0,-5.916079783099616,0.0));#10=DIRECTION('',(0.0,0.0,-1.0));\
#11=DIRECTION('',(1.0,0.0,0.0));#12=AXIS2_PLACEMENT_3D('',#9,#10,#11);\
#13=CIRCLE('',#12,6.0);#14=EDGE_CURVE('',#2,#4,#13,.T.);\
#15=CARTESIAN_POINT('',(-1.0,0.0,0.0));\
#16=CARTESIAN_POINT('',(-0.7517796174689421,0.037710394270029646,0.0));\
#17=CARTESIAN_POINT('',(-0.5014240990099426,0.05661929474018734,0.0));\
#18=CARTESIAN_POINT('',(-0.251068580550943,0.07552819521034593,0.0));\
#19=CARTESIAN_POINT('',(4.0767318863933883E-16,0.07552819521034593,0.0));\
#20=CARTESIAN_POINT('',(0.25106858055094383,0.07552819521034593,0.0));\
#21=CARTESIAN_POINT('',(0.5014240990099433,0.05661929474018734,0.0));\
#22=CARTESIAN_POINT('',(0.751779617468943,0.037710394270029646,0.0));\
#23=CARTESIAN_POINT('',(1.0,0.0,0.0));\
#24=(BOUNDED_CURVE()B_SPLINE_CURVE(2,(#15,#16,#17,#18,#19,#20,#21,#22,#23),.CIRCULAR_ARC.,.F.,.F.)B_SPLINE_CURVE_WITH_KNOTS((3,2,2,2,3),(0.,1.,2.,3.,4.),.UNSPECIFIED.)CURVE()GEOMETRIC_REPRESENTATION_ITEM()RATIONAL_B_SPLINE_CURVE((1.0,0.9992897205446609,1.0,0.9992897205446609,1.0,0.9992897205446609,1.0,0.9992897205446609,1.0))REPRESENTATION_ITEM(''));\
#25=EDGE_CURVE('',#2,#4,#24,.T.);#26=CARTESIAN_POINT('',(0.0,-5.916079783099616,1.0));\
#27=DIRECTION('',(0.0,0.0,-1.0));#28=DIRECTION('',(1.0,0.0,0.0));\
#29=AXIS2_PLACEMENT_3D('',#26,#27,#28);#30=CIRCLE('',#29,6.0);\
#31=EDGE_CURVE('',#6,#8,#30,.T.);#32=CARTESIAN_POINT('',(-1.0,0.0,1.0));\
#33=CARTESIAN_POINT('',(-0.7517796174689421,0.037710394270029646,1.0));\
#34=CARTESIAN_POINT('',(-0.5014240990099426,0.05661929474018734,1.0));\
#35=CARTESIAN_POINT('',(-0.251068580550943,0.07552819521034593,1.0));\
#36=CARTESIAN_POINT('',(4.0767318863933883E-16,0.07552819521034593,1.0));\
#37=CARTESIAN_POINT('',(0.25106858055094383,0.07552819521034593,1.0));\
#38=CARTESIAN_POINT('',(0.5014240990099433,0.05661929474018734,1.0));\
#39=CARTESIAN_POINT('',(0.751779617468943,0.037710394270029646,1.0));\
#40=CARTESIAN_POINT('',(1.0,0.0,1.0));\
#41=(BOUNDED_CURVE()B_SPLINE_CURVE(2,(#32,#33,#34,#35,#36,#37,#38,#39,#40),.CIRCULAR_ARC.,.F.,.F.)B_SPLINE_CURVE_WITH_KNOTS((3,2,2,2,3),(0.,1.,2.,3.,4.),.UNSPECIFIED.)CURVE()GEOMETRIC_REPRESENTATION_ITEM()RATIONAL_B_SPLINE_CURVE((1.0,0.9992897205446609,1.0,0.9992897205446609,1.0,0.9992897205446609,1.0,0.9992897205446609,1.0))REPRESENTATION_ITEM(''));\
#42=EDGE_CURVE('',#6,#8,#41,.T.);#43=CARTESIAN_POINT('',(-1.0,0.0,0.0));\
#44=DIRECTION('',(0.0,0.0,1.0));#45=VECTOR('',#44,1.);#46=LINE('',#43,#45);\
#47=EDGE_CURVE('',#2,#6,#46,.T.);#48=CARTESIAN_POINT('',(1.0,0.0,0.0));\
#49=DIRECTION('',(0.0,0.0,1.0));#50=VECTOR('',#49,1.);#51=LINE('',#48,#50);\
#52=EDGE_CURVE('',#4,#8,#51,.T.);#53=CARTESIAN_POINT('',(0.0,0.0,1.0));\
#54=DIRECTION('',(0.0,0.0,1.0));#55=DIRECTION('',(1.0,0.0,0.0));\
#56=AXIS2_PLACEMENT_3D('',#53,#54,#55);#57=PLANE('',#56);#58=ORIENTED_EDGE('',*,*,#42,.T.);\
#59=ORIENTED_EDGE('',*,*,#31,.F.);#60=EDGE_LOOP('',(#58,#59));\
#61=FACE_OUTER_BOUND('',#60,.T.);#62=ADVANCED_FACE('',(#61),#57,.T.);\
#63=CARTESIAN_POINT('',(0.0,0.0,0.0));#64=DIRECTION('',(0.0,0.0,-1.0));\
#65=DIRECTION('',(1.0,0.0,0.0));#66=AXIS2_PLACEMENT_3D('',#63,#64,#65);#67=PLANE('',#66);\
#68=ORIENTED_EDGE('',*,*,#14,.T.);#69=ORIENTED_EDGE('',*,*,#25,.F.);\
#70=EDGE_LOOP('',(#68,#69));#71=FACE_OUTER_BOUND('',#70,.T.);\
#72=ADVANCED_FACE('',(#71),#67,.T.);#73=CARTESIAN_POINT('',(0.0,-5.916079783099616,0.0));\
#74=DIRECTION('',(0.0,0.0,1.0));#75=DIRECTION('',(1.0,0.0,0.0));\
#76=AXIS2_PLACEMENT_3D('',#73,#74,#75);#77=CYLINDRICAL_SURFACE('',#76,6.0);\
#78=ORIENTED_EDGE('',*,*,#31,.T.);#79=ORIENTED_EDGE('',*,*,#52,.F.);\
#80=ORIENTED_EDGE('',*,*,#14,.F.);#81=ORIENTED_EDGE('',*,*,#47,.T.);\
#82=EDGE_LOOP('',(#78,#79,#80,#81));#83=FACE_OUTER_BOUND('',#82,.T.);\
#84=ADVANCED_FACE('',(#83),#77,.T.);#85=CARTESIAN_POINT('',(0.0,-6.582280226339064,0.0));\
#86=DIRECTION('',(0.0,0.0,1.0));#87=DIRECTION('',(1.0,0.0,0.0));\
#88=AXIS2_PLACEMENT_3D('',#85,#86,#87);#89=CYLINDRICAL_SURFACE('',#88,6.65780842154941);\
#90=ORIENTED_EDGE('',*,*,#25,.T.);#91=ORIENTED_EDGE('',*,*,#52,.T.);\
#92=ORIENTED_EDGE('',*,*,#42,.F.);#93=ORIENTED_EDGE('',*,*,#47,.F.);\
#94=EDGE_LOOP('',(#90,#91,#92,#93));#95=FACE_OUTER_BOUND('',#94,.T.);\
#96=ADVANCED_FACE('',(#95),#89,.F.);#97=CLOSED_SHELL('',(#62,#72,#84,#96));\
#98=MANIFOLD_SOLID_BREP('',#97);\n";

#[test]
fn thin_faces_between_near_arcs_take_their_sense_from_the_arcs() {
    // Each arc's segment is its sector less the triangle on the chord.
    let segment = |r: f64| r * r * (1. / r).asin() - (r * r - 1.).sqrt();
    let flat = 0.9 * (6. - 35f64.sqrt());
    let area = segment(6.) - segment((1. + flat * flat) / (2. * flat));
    let text = part21("#98", CRESCENT);
    let boundary = read_boundaries(&text, TOLERANCE).unwrap().remove(0);
    let samples = boundary.sample(f64::INFINITY, f64::INFINITY).unwrap();
    let top = &boundary.faces[0];
    let points = domain::raw(&boundary, &top.loops[0].uses, &samples, false);
    let unrolled = domain::unroll(&top.surface, points, TOLERANCE)
        .unwrap()
        .unwrap();
    let n = unrolled.points.len();
    let chords: f64 = (0..n)
        .map(|i| {
            let (a, b) = (unrolled.points[i].uv, unrolled.points[(i + 1) % n].uv);
            0.5 * (a[0] * b[1] - b[0] * a[1])
        })
        .sum();
    assert!(chords < -0.05, "{chords}");
    assert!(unrolled.exact_area(&boundary, &top.surface) > 0.);
    let crescent = solid(&text);
    assert!(
        close(crescent.volume_mm3(), area, 1e-6),
        "{} {area}",
        crescent.volume_mm3()
    );
}

/// A square picture frame: a closed bilinear spline net, 20 mm across with
/// a 0.4 by 4 mm square section scaled about its axis, whose every cell
/// is a flat quadrilateral.
fn frame_net() -> (NurbsSurface, impl Fn(usize, usize) -> Point) {
    let at = |i: usize, j: usize| {
        let corner = [[10., 10.], [-10., 10.], [-10., -10.], [10., -10.]][i % 4];
        let (k, z) = [(0.8, -2.), (1.2, -2.), (1.2, 2.), (0.8, 2.)][j % 4];
        [corner[0] * k, corner[1] * k, z]
    };
    let points = (0..5).flat_map(|i| (0..5).map(move |j| at(i, j))).collect();
    let knots = vec![0., 1., 2., 3., 4.];
    let multiplicities = vec![2, 1, 1, 1, 2];
    let net = NurbsSurface::new(
        [1, 1],
        [5, 5],
        points,
        None,
        [knots.clone(), knots],
        [multiplicities.clone(), multiplicities],
    )
    .unwrap();
    (net, at)
}

#[test]
fn loops_winding_both_ways_round_a_closed_net_are_refused() {
    let (net, at) = frame_net();
    let surface = Surface::nurbs(net, TOLERANCE);
    assert_eq!(surface.periods(), [Some(4.), Some(4.)]);
    // A closed net repeats past its knot range.
    let (inside, _) = surface.evaluate([0.5, 0.3]);
    assert!(distance(surface.evaluate([4.5, 4.3]).0, inside) < 1e-12);
    // Two diagonals across the flat cells each wind once along u and once
    // along v; together they cut the frame into two strips.
    let diagonal = |offset: usize| {
        let points = (0..5).map(|k| at(k, k + offset)).collect();
        let curve = NurbsCurve::new(
            1,
            points,
            None,
            vec![0., 1., 2., 3., 4.],
            vec![2, 1, 1, 1, 2],
        )
        .unwrap();
        Curve::nurbs(curve, TOLERANCE)
    };
    let vertices = vec![at(0, 0), at(0, 2)];
    let edges = vec![
        Edge::new(diagonal(0), [0, 0], [vertices[0]; 2], TOLERANCE).unwrap(),
        Edge::new(diagonal(2), [1, 1], [vertices[1]; 2], TOLERANCE).unwrap(),
    ];
    let strip = |first: usize, second: usize| Face {
        surface: surface.clone(),
        same_sense: true,
        loops: vec![
            Loop {
                uses: vec![(first, true)],
                vertex: None,
            },
            Loop {
                uses: vec![(second, false)],
                vertex: None,
            },
        ],
    };
    let boundary = Boundary {
        vertices,
        edges,
        faces: vec![strip(0, 1), strip(1, 0)],
        shells: vec![Shell {
            faces: vec![0, 1],
            reversed: false,
        }],
    };
    assert_eq!(
        Solid::new(boundary, TOLERANCE, None).unwrap_err(),
        Error::Unsupported
    );
}

fn flat(points: &[[f64; 2]], first: usize) -> UvLoop {
    UvLoop {
        points: points
            .iter()
            .enumerate()
            .map(|(k, &uv)| LoopPoint {
                uv,
                pole: false,
                key: Key::Extra(0, first + k),
            })
            .collect(),
        segments: vec![None; points.len()],
        winding: [0, 0],
    }
}

fn circle(radius: f64, clockwise: bool) -> Vec<[f64; 2]> {
    (0..12)
        .map(|k| {
            let t = f64::from(k) * PI / 6. * if clockwise { -1. } else { 1. };
            [radius * t.cos(), radius * t.sin()]
        })
        .collect()
}

fn plane() -> Surface {
    Surface::Plane(Frame {
        origin: [0.; 3],
        x: [1., 0., 0.],
        y: [0., 1., 0.],
        z: [0., 0., 1.],
    })
}

fn covered(triangles: &[polygon::Triangle]) -> f64 {
    triangles
        .iter()
        .map(|[a, b, c]| {
            0.5 * ((b.uv[0] - a.uv[0]) * (c.uv[1] - a.uv[1])
                - (b.uv[1] - a.uv[1]) * (c.uv[0] - a.uv[0]))
        })
        .sum()
}

#[test]
fn parameter_regions_triangulate_or_refuse_what_they_cannot_cover() {
    let square = [[-4., -4.], [4., -4.], [4., 4.], [-4., 4.]];
    let region = |loops: Vec<UvLoop>| Domain::Trimmed { loops, pole: None };
    // A hole inside another hole cannot be bridged to the outer loop.
    let nested = region(vec![
        flat(&square, 0),
        flat(&circle(2., true), 10),
        flat(&circle(1., true), 30),
    ]);
    assert_eq!(
        polygon::triangulate(&nested, &plane(), 0, &[]).unwrap_err(),
        Error::Unsupported
    );
    // A triangle whose corners repeat a boundary point has no extent.
    let mut pinched = flat(&[[0., 0.], [1., 0.], [0., 1.]], 0);
    pinched.points[2].key = pinched.points[0].key;
    assert_eq!(
        polygon::triangulate(&region(vec![pinched]), &plane(), 0, &[]).unwrap_err(),
        Error::Unsupported
    );
    // A strip whose ends repeat the same two boundary points, as a band
    // cut open along a seam: the second diagonal between them is refused
    // once the first is laid, and the strip is still covered once.
    let mut band = flat(
        &[[0., 0.], [1., 0.], [2., 0.], [2., 1.], [1., 1.], [0., 1.]],
        0,
    );
    band.points[2].key = band.points[0].key;
    band.points[5].key = band.points[3].key;
    let triangles = polygon::triangulate(&region(vec![band]), &plane(), 0, &[]).unwrap();
    assert!(close(covered(&triangles), 2., 1e-12), "{triangles:?}");
    // A boundary point on the diagonal of an ear keeps that ear closed:
    // the triangle beside it takes the point as a corner.
    let notch = flat(&[[0., 0.], [2., 0.], [2., 2.], [1., 1.]], 0);
    let triangles = polygon::triangulate(&region(vec![notch]), &plane(), 0, &[]).unwrap();
    assert!(close(covered(&triangles), 2., 1e-12));
    assert!(
        triangles
            .iter()
            .flatten()
            .any(|c| c.key == Key::Extra(0, 3))
    );
    // Two squares touching at one point, as one ring through it twice:
    // once both are covered the pinch is a sliver without area.
    let mut pinch = flat(
        &[
            [0., 0.],
            [1., 0.],
            [1., 1.],
            [2., 1.],
            [2., 2.],
            [1., 2.],
            [1., 1.],
            [0., 1.],
        ],
        0,
    );
    pinch.points[6].key = pinch.points[2].key;
    let triangles = polygon::triangulate(&region(vec![pinch]), &plane(), 0, &[]).unwrap();
    assert!(close(covered(&triangles), 2., 1e-12));
    // A ring that meets one boundary point twice at different parameters:
    // the first diagonal laid to that point forbids a second, which an
    // ear judged before must notice when it is taken.
    let points = [
        [0.75, 0.25],
        [0.25, 0.75],
        [-0.25, 1.5],
        [-1.25, 0.25],
        [-0.75, -0.5],
        [-0.25, -0.75],
        [1., -1.],
    ];
    let mut twice = flat(&points, 0);
    twice.points[3].key = twice.points[5].key;
    let triangles = polygon::triangulate(&region(vec![twice]), &plane(), 0, &[]).unwrap();
    assert!(close(covered(&triangles), 3.0625, 1e-12));
    // Points on one straight side all become corners.
    let run = flat(
        &[[0., 0.], [1., 0.], [2., 0.], [3., 0.], [3., 1.], [0., 1.]],
        0,
    );
    let triangles = polygon::triangulate(&region(vec![run]), &plane(), 0, &[]).unwrap();
    assert!(close(covered(&triangles), 3., 1e-12));
    for k in 0..6 {
        assert!(
            triangles
                .iter()
                .flatten()
                .any(|c| c.key == Key::Extra(0, k))
        );
    }
}

#[test]
fn degenerate_parameters_seeds_and_budgets_behave() {
    // A cone made by revolving a line through its axis: at the apex the
    // u derivative vanishes and Newton steps along v alone.
    let cone = Surface::Revolution {
        curve: Curve::Line {
            origin: [0.; 3],
            direction: [1., 0., 1.],
        },
        origin: [0.; 3],
        axis: [0., 0., 1.],
    };
    let target = cone.evaluate([0.7, 0.5]).0;
    for seed in [Some([0.7, 0.]), None] {
        let uv = cone.invert(target, seed);
        assert!(distance(cone.evaluate(uv).0, target) < 1e-12, "{seed:?}");
    }
    assert_eq!(plane().pole(true, 0.), None);
    // A seeded spline inversion settles on the foot point.
    let spline = Curve::nurbs(
        NurbsCurve::new(
            2,
            vec![[0., 0., 0.], [1., 2., 0.], [3., 2., 0.], [4., 0., 0.]],
            None,
            vec![0., 0.5, 1.],
            vec![3, 1, 3],
        )
        .unwrap(),
        TOLERANCE,
    );
    let (point, _) = spline.evaluate(0.37);
    assert!((spline.invert(point, Some(0.2)) - 0.37).abs() < 1e-9);
    // An edge's polyline runs along its curve from end to end.
    let elliptic = solid(&part21("#61", ELLIPTIC));
    for edge in &elliptic.boundary().edges {
        let line = edge.polyline(0.01, 0.3).unwrap();
        assert!(line.len() >= 2);
        assert_eq!(line[0], edge.curve.evaluate(edge.range[0]).0);
        assert_eq!(*line.last().unwrap(), edge.curve.evaluate(edge.range[1]).0);
        for p in &line {
            let t = edge.curve.invert(*p, None);
            assert!(distance(edge.curve.evaluate(t).0, *p) < 1e-9);
        }
    }
    // Refinement stops at the triangle budget.
    let sphere = solid(SPHERE);
    let tight = mesh::Options {
        triangle_limit: 300,
        ..options(0.001)
    };
    assert_eq!(
        mesh::tessellate(&sphere, tight, None).unwrap_err(),
        Error::Limit
    );
    // Coordinates beyond the exchange range are not written.
    let far = elliptic.moved(&|p| add(p, [2e12, 0., 0.]), &|v| v).unwrap();
    assert_eq!(
        crate::step::Writer::new().append_brep(&far).unwrap_err(),
        crate::step::StepError::InvalidPlacement
    );
}

#[test]
fn spline_records_units_and_uncertainties_are_read_to_their_arity() {
    let refused = |from: &str, to: &str| solids(&elliptic(from, to)).map(|s| s.len()).unwrap_err();
    let line = "#16=LINE('',#10,#15)";
    let rational = |base: &str, knots: &str, weights: &str| {
        format!(
            "#16=(BOUNDED_CURVE(){base}{knots}CURVE()GEOMETRIC_REPRESENTATION_ITEM()\
             {weights}REPRESENTATION_ITEM(''))"
        )
    };
    let base = "B_SPLINE_CURVE(1,(#10,#12),.UNSPECIFIED.,.F.,.F.)";
    let knots = "B_SPLINE_CURVE_WITH_KNOTS((2,2),(0.,1.),.UNSPECIFIED.)";
    let weights = "RATIONAL_B_SPLINE_CURVE((1.,1.))";
    // The rational, plain and Bezier lines written as complex records read
    // as the line.
    for to in [
        rational(base, knots, weights),
        rational(base, knots, ""),
        rational(base, "BEZIER_CURVE()", ""),
    ] {
        let read = solids(&elliptic(line, &to)).unwrap();
        assert!(close(read[0].volume_mm3(), 30. * PI, 1e-12), "{to}");
    }
    for (to, error) in [
        (
            "#16=B_SPLINE_CURVE_WITH_KNOTS('',1,(#10,#12))".to_owned(),
            Error::Unsupported,
        ),
        (
            rational(base, "B_SPLINE_CURVE_WITH_KNOTS((2,2),(0.,1.))", weights),
            Error::Invalid,
        ),
        (
            rational(base, knots, "RATIONAL_B_SPLINE_CURVE((1.,1.),1.)"),
            Error::Invalid,
        ),
        (
            rational("B_SPLINE_CURVE(1)", knots, weights),
            Error::Invalid,
        ),
        // Without knots a spline is a Bezier curve, which must say so.
        (rational(base, "", weights), Error::Unsupported),
    ] {
        assert_eq!(refused(line, &to), error, "{to}");
    }
    let plane = "#20=PLANE('',#5)";
    let points = "#70=CARTESIAN_POINT('',(-4.,-3.,0.));#71=CARTESIAN_POINT('',(-4.,3.,0.));\
                  #72=CARTESIAN_POINT('',(4.,-3.,0.));#73=CARTESIAN_POINT('',(4.,3.,0.))";
    let net = |base: &str, knots: &str, weights: &str| {
        format!(
            "#20=(BOUNDED_SURFACE(){base}{knots}GEOMETRIC_REPRESENTATION_ITEM()\
             {weights}REPRESENTATION_ITEM('')SURFACE());{points}"
        )
    };
    let base = "B_SPLINE_SURFACE(1,1,((#70,#71),(#72,#73)),.UNSPECIFIED.,.F.,.F.,.F.)";
    let knots = "B_SPLINE_SURFACE_WITH_KNOTS((2,2),(2,2),(0.,1.),(0.,1.),.UNSPECIFIED.)";
    let weights = "RATIONAL_B_SPLINE_SURFACE(((1.,1.),(1.,1.)))";
    for to in [net(base, knots, ""), net(base, "BEZIER_SURFACE()", "")] {
        let read = solids(&elliptic(plane, &to)).unwrap();
        assert!(close(read[0].volume_mm3(), 30. * PI, 1e-12), "{to}");
    }
    for (to, error) in [
        (net(base, "", weights), Error::Unsupported),
        (
            format!("#20=B_SPLINE_SURFACE_WITH_KNOTS('',1,1);{points}"),
            Error::Unsupported,
        ),
        (
            net(
                base,
                "B_SPLINE_SURFACE_WITH_KNOTS((2,2),(2,2),(0.,1.),(0.,1.))",
                weights,
            ),
            Error::Invalid,
        ),
        (
            net(
                base,
                knots,
                "RATIONAL_B_SPLINE_SURFACE(((1.,1.),(1.,1.)),1.)",
            ),
            Error::Invalid,
        ),
        (net("B_SPLINE_SURFACE(1,1)", knots, weights), Error::Invalid),
        // Knot multiplicities must be a list.
        (
            net(
                base,
                "B_SPLINE_SURFACE_WITH_KNOTS(2,(2,2),(0.,1.),(0.,1.),.UNSPECIFIED.)",
                weights,
            ),
            Error::Invalid,
        ),
    ] {
        assert_eq!(refused(plane, &to), error, "{to}");
    }
    // A vertex written as a complex record is outside the boundary kinds.
    assert_eq!(
        refused(
            "#9=VERTEX_POINT('',#10)",
            "#9=(VERTEX_POINT('',#10)REPRESENTATION_ITEM(''))"
        ),
        Error::Unsupported
    );
    // An edge stated against its curve reads as the same edge.
    let reversed = elliptic(
        "#17=EDGE_CURVE('',#9,#11,#16,.T.)",
        "#17=EDGE_CURVE('',#11,#9,#16,.F.)",
    );
    let reversed = tampered(
        &reversed,
        "#51=ORIENTED_EDGE('',*,*,#17,.T.)",
        "#51=ORIENTED_EDGE('',*,*,#17,.F.)",
    );
    let reversed = tampered(
        &reversed,
        "#53=ORIENTED_EDGE('',*,*,#17,.F.)",
        "#53=ORIENTED_EDGE('',*,*,#17,.T.)",
    );
    assert!(close(solid(&reversed).volume_mm3(), 30. * PI, 1e-12));
    // Uncertainties: one length per measure, and only of a simple record.
    assert_eq!(
        refused("LENGTH_MEASURE(1.E-07)", "LENGTH_MEASURE(1.E-07,1.)"),
        Error::Invalid
    );
    assert_eq!(
        refused(
            "GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#905))",
            "GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#902))"
        ),
        Error::Unsupported
    );
    // A second context in degrees disagrees with the first in radians.
    let degrees = "#990=(GEOMETRIC_REPRESENTATION_CONTEXT(3)\
                   GLOBAL_UNIT_ASSIGNED_CONTEXT((#902,#991,#904))REPRESENTATION_CONTEXT('',''));\
                   #991=(CONVERSION_BASED_UNIT('DEGREE',#992)NAMED_UNIT(*)PLANE_ANGLE_UNIT());\
                   #992=PLANE_ANGLE_MEASURE_WITH_UNIT(PLANE_ANGLE_MEASURE(1.745329251994E-02),#903);";
    assert_eq!(
        refused("#60=", &format!("{degrees}#60=")),
        Error::Unsupported
    );
    // A conversion names its measure, a simple record holding one value.
    let converted = tampered(
        FRUSTUM,
        "#115 = ( NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.) );",
        "#115 = ( CONVERSION_BASED_UNIT('DEGREE',#119) NAMED_UNIT(*) PLANE_ANGLE_UNIT() );\n\
         #119 = PLANE_ANGLE_MEASURE_WITH_UNIT(PLANE_ANGLE_MEASURE(1.745329251994E-02),#120);\n\
         #120 = ( NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.) );",
    );
    for (from, to, error) in [
        ("'DEGREE',#119)", "'DEGREE')", Error::Invalid),
        ("'DEGREE',#119)", "'DEGREE',#120)", Error::Unsupported),
        (
            "WITH_UNIT(PLANE_ANGLE_MEASURE(1.745329251994E-02),#120)",
            "WITH_UNIT(1.745329251994E-02,#120)",
            Error::Invalid,
        ),
        (
            "PLANE_ANGLE_MEASURE(1.745329251994E-02)",
            "PLANE_ANGLE_MEASURE(1.745329251994E-02,1.)",
            Error::Invalid,
        ),
    ] {
        let text = tampered(&converted, from, to);
        assert_eq!(solids(&text).map(|s| s.len()).unwrap_err(), error, "{to}");
    }
}

#[test]
fn sources_list_their_solids_and_verification_compares_every_measure() {
    let text = part21("#61", ELLIPTIC);
    let source = crate::step_read::brep::Source::parse(&text, TOLERANCE).unwrap();
    assert_eq!(source.solids(), &[61]);
    assert!(!source.is_assembly());
    // The cone over the equator disc has the hemisphere's vertex, counts
    // and face areas but not its volume.
    let dome = part21("#20", DOME);
    let cone = solid(&comet());
    assert!(verify_retained(&comet(), &[&cone], TOLERANCE).is_ok());
    assert_eq!(
        verify_retained(&dome, &[&cone], TOLERANCE),
        Err(ReadError::Invalid)
    );
    // A steeper cone differs in area as well.
    let steeper = solid(&tampered(&comet(), "3.,0.5235987755982988)", "3.,0.6)"));
    assert_eq!(
        verify_retained(&dome, &[&steeper], TOLERANCE),
        Err(ReadError::Invalid)
    );
    // A source the provider does not read keeps its refusal kind.
    let open = elliptic("#60=", "#990=OPEN_SHELL('',(#33));#60=");
    let elliptic = solid(&text);
    assert_eq!(
        verify_retained(&open, &[&elliptic], TOLERANCE),
        Err(ReadError::Unsupported)
    );
    // A revolved boundary writes its axis and profile back.
    assert_eq!(
        solids(&tampered(
            DISC,
            "AXIS1_PLACEMENT('',#150,#151)",
            "AXIS1_PLACEMENT('',#150)"
        ))
        .unwrap_err(),
        Error::Invalid
    );
    let disc = solid(DISC);
    let mut writer = crate::step::Writer::new();
    writer.append_brep(&disc).unwrap();
    let written = writer.finish().unwrap();
    assert!(written.contains("SURFACE_OF_REVOLUTION"));
    assert!(verify_retained(&written, &[&disc], TOLERANCE).is_ok());
    // A sphere whose pole vertex is stated 2 um off its surface reads to
    // its declared 5 um uncertainty with the vertex re-derived onto it.
    let loose = tampered(SPHERE, "(3.061616996623E-16,0.,-5.)", "(0.,0.,-5.000002)");
    assert_eq!(solids(&loose).unwrap_err(), Error::Invalid);
    let declared = tampered(&loose, "LENGTH_MEASURE(1.E-07)", "LENGTH_MEASURE(5.E-06)");
    let sphere = solid(&declared);
    assert!(close(sphere.volume_mm3(), 4. / 3. * PI * 125., 1e-12));
    assert!(distance(sphere.boundary().vertices[0], [0., 0., -5.]) < 1e-12);
    // Of several declared length uncertainties the smallest governs.
    let both = |second: &str| {
        let text = tampered(
            &declared,
            "GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#31))",
            "GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#31,#90))",
        );
        tampered(
            &text,
            "#32 = ",
            &format!(
                "#90 = UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE({second}),#28,'','');\n#32 = "
            ),
        )
    };
    assert!(solids(&both("3.E-06")).is_ok());
    assert_eq!(solids(&both("1.E-06")).unwrap_err(), Error::Invalid);
}

#[test]
fn caps_closing_below_their_loop_and_parallel_parameters_integrate_or_refuse() {
    // The lower hemisphere: the equator runs the other way round the
    // sphere, so the cap closes at the south pole below its loop.
    let bowl = tampered(
        &part21("#20", DOME),
        "#6=CIRCLE('',#4,3.)",
        "#6=CIRCLE('',#40,3.);#40=AXIS2_PLACEMENT_3D('',#1,#41,#3);#41=DIRECTION('',(0.,0.,-1.))",
    );
    let bowl = tampered(
        &bowl,
        "#18=ADVANCED_FACE('',(#17),#7,.F.)",
        "#18=ADVANCED_FACE('',(#17),#7,.T.)",
    );
    let bowl = solid(&bowl);
    assert!(
        close(bowl.volume_mm3(), 18. * PI, 1e-12),
        "{}",
        bowl.volume_mm3()
    );
    assert!(distance(bowl.bounds_mm()[0], [-3., -3., -3.]) < 1e-9);
    certified(&bowl, 0.01);
    // A face swept along its own profile has no area element: the rate of
    // its parameters along the boundary is undefined and is refused.
    let line = Curve::Line {
        origin: [0.; 3],
        direction: [1., 0., 0.],
    };
    let boundary = Boundary {
        vertices: vec![[0.; 3], [1., 0., 0.]],
        edges: vec![Edge::new(line.clone(), [0, 1], [[0.; 3], [1., 0., 0.]], TOLERANCE).unwrap()],
        faces: vec![],
        shells: vec![],
    };
    let face = Face {
        surface: Surface::Extrusion {
            curve: line,
            direction: [1., 0., 0.],
        },
        same_sense: true,
        loops: vec![],
    };
    let point = |vertex: usize, t: [f64; 2]| {
        let uv = face.surface.invert(boundary.vertices[vertex], None);
        (
            LoopPoint {
                uv,
                pole: false,
                key: Key::Vertex(vertex),
            },
            Some(domain::Segment { edge: 0, t }),
        )
    };
    let (points, segments) = [point(0, [0., 1.]), point(1, [1., 0.])].into_iter().unzip();
    let domain = Domain::Trimmed {
        loops: vec![UvLoop {
            points,
            segments,
            winding: [0, 0],
        }],
        pole: None,
    };
    assert_eq!(
        mass::integrate(&boundary, &face, &domain, [0.; 3], 1, TOLERANCE).unwrap_err(),
        Error::Unsupported
    );
}

#[test]
fn interior_extrema_are_sought_only_where_newton_settles() {
    // The graph z = u^4 / 4 - u^2 + 2 u + y^2 over 0 <= u <= 1 and
    // -1 <= y <= 1 as one Bezier patch. Its height has no stationary point
    // over the face; from the corners u = 0 and u = 1 Newton's method on
    // the height's gradient cycles between them and is abandoned.
    // The patch spans -1 <= u <= 2 and -2 <= y <= 2 so that the face's
    // corners lie inside it, where the height's derivatives are exact.
    let height = [-2.75, -0.5, 2.5, -0.5, 4.];
    let square = [4., -4., 4.];
    let points = (0..5)
        .flat_map(|i| {
            (0..3).map(move |j| {
                [
                    -1. + 3. * i as f64 / 4.,
                    2. * j as f64 - 2.,
                    height[i] + square[j],
                ]
            })
        })
        .collect();
    let patch = NurbsSurface::new(
        [4, 2],
        [5, 3],
        points,
        None,
        [vec![0., 1.], vec![0., 1.]],
        [vec![5, 5], vec![3, 3]],
    )
    .unwrap();
    let surface = Surface::nurbs(patch, TOLERANCE);
    let (top, _) = surface.evaluate([2. / 3., 0.5]);
    assert!(distance(top, [1., 0., 1.25]) < 1e-12);
    let boundary = Boundary {
        faces: vec![Face {
            surface,
            same_sense: true,
            loops: vec![],
        }],
        ..Boundary::default()
    };
    let (a, b) = (1. / 3., 2. / 3.);
    let domain = Domain::Trimmed {
        loops: vec![flat(&[[a, 0.25], [b, 0.25], [b, 0.75], [a, 0.75]], 0)],
        pole: None,
    };
    let [low, high] = extent::bounds(&boundary, &[domain]).unwrap();
    assert_eq!((low, high), ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]));
}

/// A horn torus (major radius equal to the minor radius, 2 mm) without the
/// inner upper quarter that runs from its top circle down to the point
/// where it touches its axis; the plane through that circle closes it.
const HORN_RING: &str = "\
#1=CARTESIAN_POINT('',(0.,0.,0.));#2=DIRECTION('',(0.,0.,1.));#3=DIRECTION('',(1.,0.,0.));\
#4=AXIS2_PLACEMENT_3D('',#1,#2,#3);#5=TOROIDAL_SURFACE('',#4,2.,2.);\
#21=CARTESIAN_POINT('',(0.,0.,2.));#22=AXIS2_PLACEMENT_3D('',#21,#2,#3);\
#6=CIRCLE('',#22,2.);#7=PLANE('',#22);#8=CARTESIAN_POINT('',(2.,0.,2.));#9=VERTEX_POINT('',#8);\
#10=EDGE_CURVE('',#9,#9,#6,.T.);\
#11=ORIENTED_EDGE('',*,*,#10,.F.);#12=EDGE_LOOP('',(#11));#13=FACE_BOUND('',#12,.T.);\
#14=ADVANCED_FACE('',(#13),#5,.T.);\
#15=ORIENTED_EDGE('',*,*,#10,.T.);#16=EDGE_LOOP('',(#15));#17=FACE_OUTER_BOUND('',#16,.T.);\
#18=ADVANCED_FACE('',(#17),#7,.T.);\
#19=CLOSED_SHELL('',(#14,#18));#20=MANIFOLD_SOLID_BREP('',#19);\n";

#[test]
fn horn_torus_faces_close_at_the_point_on_their_axis() {
    let text = part21("#20", HORN_RING);
    let ring = solid(&text);
    // The horn torus encloses 2 pi^2 R^3; the cup between the plane z = R
    // and the removed funnel rho = R - sqrt(R^2 - z^2) adds
    // pi R^3 (5/3 - pi/2). The funnel's area 2 pi R^2 (pi/2 - 1) leaves
    // the torus's 4 pi^2 R^2.
    let volume = 16. * PI * PI + 8. * PI * (5. / 3. - PI / 2.);
    assert!(
        close(ring.volume_mm3(), volume, 1e-10),
        "{}",
        ring.volume_mm3()
    );
    let areas = sorted_areas(&ring);
    let torus = 16. * PI * PI - 8. * PI * (PI / 2. - 1.);
    assert!(
        close(areas[0], 4. * PI, 1e-12) && close(areas[1], torus, 1e-10),
        "{areas:?}"
    );
    let bounds = ring.bounds_mm();
    assert!(distance(bounds[0], [-4., -4., -2.]) < 1e-9, "{bounds:?}");
    assert!(distance(bounds[1], [4., 4., 2.]) < 1e-9, "{bounds:?}");
    for deflection in [0.05, 0.01] {
        let mesh = certified(&ring, deflection);
        // Every vertex lies on the torus or the plane, and one is the point
        // where the torus touches its axis.
        assert!(mesh.positions.iter().all(|p| {
            let rho = p[0].hypot(p[1]);
            ((rho - 2.).hypot(p[2]) - 2.).abs() < 1e-9 || (p[2] - 2.).abs() < 1e-12
        }));
        assert!(mesh.positions.iter().any(|p| norm(*p) < 1e-12));
    }
    // Through the axis the section is the two tube circles, each without
    // its inner upper quarter, joined by the plane's chord and meeting at
    // the axis point: length 4 + 6 pi, area 2 (3 pi + 4).
    let loops = section::section(&ring, [0.; 3], [0., 1., 0.], options(0.01), None).unwrap();
    assert_eq!(loops.len(), 1, "{loops:?}");
    let shoelace: f64 = (0..loops[0].len())
        .map(|i| {
            let (a, b) = (loops[0][i], loops[0][(i + 1) % loops[0].len()]);
            0.5 * (a[2] * b[0] - b[2] * a[0])
        })
        .sum();
    let length = section::length(&loops[0]);
    // Chords run inside the curved arcs, every point on the exact section.
    assert!(
        length <= 4. + 6. * PI && length > 4. + 6. * PI - 0.05,
        "{length}"
    );
    assert!(loops[0].iter().all(|p| {
        p[1].abs() < 1e-12
            && ((p[0].abs() - 2.).hypot(p[2]) - 2.)
                .abs()
                .min((p[2] - 2.).abs())
                < 1e-9
    }));
    assert!(
        (shoelace.abs() - 2. * (3. * PI + 4.)).abs() < 0.01 * length,
        "{shoelace}"
    );
    // The loop passes the axis point within the merge distance: beside it
    // the plane's height grows only with the square of the distance.
    let nearest = loops[0]
        .iter()
        .map(|p| norm(*p))
        .fold(f64::INFINITY, f64::min);
    assert!(nearest < 1e-5, "{nearest}");
    // Across the removed funnel the section is one disk of radius 2 + sqrt 3.
    let loops = section::section(&ring, [0., 0., 1.], [0., 0., 1.], options(0.01), None).unwrap();
    assert_eq!(loops.len(), 1, "{loops:?}");
    assert!(
        loops[0]
            .iter()
            .all(|p| (p[0].hypot(p[1]) - (2. + 3f64.sqrt())).abs() < 1e-9)
    );
    // The removed quarter with the plane is the cup alone. Its torus face is
    // reversed and its lone loop bounds a region on either side, both
    // closing at the axis point, so its side is not decided.
    let cup = tampered(
        &text,
        "#14=ADVANCED_FACE('',(#13),#5,.T.)",
        "#14=ADVANCED_FACE('',(#13),#5,.F.)",
    );
    assert_eq!(solids(&cup).unwrap_err(), Error::Unsupported);
    // A spindle torus, its minor radius past the major one, stays refused.
    let spindle = tampered(
        &text,
        "TOROIDAL_SURFACE('',#4,2.,2.)",
        "TOROIDAL_SURFACE('',#4,2.,2.5)",
    );
    assert_eq!(solids(&spindle).unwrap_err(), Error::Unsupported);
}
