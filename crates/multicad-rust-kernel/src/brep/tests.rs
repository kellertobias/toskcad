//! The retained boundary provider against analytic oracles and tampered
//! exchange text.
use super::*;
use crate::manufacturing::certify;
use crate::step_read::ReadError;
use crate::step_read::brep::{read_boundaries, read_solids};
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
    // A gap beyond the capped tolerance stays refused, whatever is declared.
    let far = tampered(FRUSTUM, STATED, "#23 = CARTESIAN_POINT('',(8.00002,0.,0.))");
    assert_eq!(
        solids(&declared(&far, "1.E-03")).unwrap_err(),
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
    use crate::step_read::brep::verify_retained;
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
            "#7=HYPERBOLA('',#5,3.,2.)"
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
