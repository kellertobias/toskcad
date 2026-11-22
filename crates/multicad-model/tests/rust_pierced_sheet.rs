use multicad_core::{
    commands::{Actor, Channel, Session},
    document::{Document, EntityId},
    exchange::LengthUnit,
    model::{self, Cutout, CutoutFace, CutoutShape, length},
    project_file,
};
use multicad_kernel::{GeometryPolicy, Solid};
use multicad_model::{
    Compiler, GeometryBackend,
    export::sheet_dxf::{SheetDxfLayers, encode_developed},
};
use std::{
    f64::consts::{FRAC_PI_2, PI},
    io::Cursor,
};

fn front_panel() -> (Document, EntityId, EntityId) {
    let document = project_file::read(
        &mut include_bytes!("../../../examples/front-panel.tcd").as_slice(),
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap()
    .project
    .document;
    let part = document
        .entities()
        .find(|entity| entity.name == "Front panel sheet")
        .unwrap()
        .id;
    let feature = model::sheet_feature(&document, part).unwrap();
    (document, part, feature)
}

fn compiler() -> Compiler {
    Compiler {
        omit_display: true,
        ..Default::default()
    }
}

fn near(a: f64, b: f64) {
    assert!((a - b).abs() <= 1e-7 * b.abs().max(1.), "{a} != {b}");
}

#[test]
fn saved_front_panel_keeps_all_folded_openings_and_develops_them_into_dxf() {
    let (document, part, feature) = front_panel();
    let compiler = compiler();
    let evaluate = |backend| {
        compiler
            .compile_sources_with_backend(&document, 0, &[part], backend, None)
            .unwrap()
    };
    let rust = evaluate(GeometryBackend::Rust);
    let solid = rust.solid(feature).unwrap();
    assert!(solid.uses_rust_geometry());
    assert_eq!(solid.solid_count().unwrap(), 1);
    // Independent sheet expectation from the saved sheet definition, its
    // parameters and its sheet stock: a rectangular base with two 90° down
    // flanges on edges 0 and 2, circular and rectangular cutouts and
    // fastener clearance holes, all in flat material.
    let sheet = SavedSheet::read(&document, feature);
    let SavedSheet {
        t, r, k, w, h, leg, ..
    } = sheet;
    let allowance = FRAC_PI_2 * (r + k * t);
    let blank = w * (h + 2. * (leg + allowance));
    let flat_area = blank - sheet.holes.iter().map(|hole| hole.area).sum::<f64>();
    let perimeter = 2. * (w + h + 2. * (leg + allowance))
        + sheet.holes.iter().map(|hole| hole.perimeter).sum::<f64>();
    // Folding moves the neutral layer from r + k t to the mid-surface r + t/2:
    // each 90° bend of width w gains θ t² (1/2 - k) w of volume, and its faces
    // and two side sectors change area by θ t (w (1 - 2k) + 2 t (1/2 - k)).
    let volume = t * flat_area + 2. * FRAC_PI_2 * t * t * (0.5 - k) * w;
    let surface = 2. * flat_area
        + t * perimeter
        + 2. * FRAC_PI_2 * t * (w * (1. - 2. * k) + 2. * t * (0.5 - k));
    near(solid.volume_mm3().unwrap(), volume);
    let surface_area = |s: &multicad_kernel::Solid| {
        s.faces()
            .unwrap()
            .into_iter()
            .map(|face| face.area_mm2)
            .sum()
    };
    near(surface_area(solid), surface);
    let bounds = solid.bounds().unwrap();
    for (axis, [min, max]) in [[0., w], [-t, r + leg], [-(r + t), h + r + t]]
        .into_iter()
        .enumerate()
    {
        near(bounds.min[axis], min);
        near(bounds.max[axis], max);
    }
    let receipt = solid
        .manufacturing_mesh_with_receipt(Default::default())
        .unwrap();
    assert!(!receipt.mesh.indices.is_empty());
    let temporary = tempfile::tempdir().unwrap();
    let path = temporary.path().join("pierced-front-panel.step");
    solid.write_step(&path).unwrap();
    let policy = GeometryPolicy::default();
    for imported in [
        Solid::read_step(&path, policy).unwrap(),
        Solid::import_rust_step(&std::fs::read_to_string(&path).unwrap(), feature, policy).unwrap(),
    ] {
        assert_eq!(
            imported.faces().unwrap().len(),
            solid.faces().unwrap().len()
        );
        near(imported.volume_mm3().unwrap(), solid.volume_mm3().unwrap());
    }

    let developed = rust.unfold(&document, part).unwrap();
    near(developed.thickness_mm, t);
    assert_eq!(developed.bends.len(), 2);
    for (index, bend) in developed.bends.iter().enumerate() {
        assert_eq!(bend.id, format!("flange.{index}.bend.0"));
        assert_eq!(bend.order, index);
        let record = &bend.record;
        assert_eq!(
            (record.flange, record.bend, record.edge, record.direction),
            (index, 0, 2 * index, "down")
        );
        near(record.angle_deg, 90.);
        near(record.inside_radius_mm, r);
        near(record.allowance_mm, allowance);
        near(record.zone_mm[0], 0.);
        near(record.zone_mm[1], allowance);
        near(record.length_mm, w);
        let y = if index == 0 {
            -allowance / 2.
        } else {
            h + allowance / 2.
        };
        let ends = if index == 0 { [0., w] } else { [w, 0.] };
        for (point, x) in record.line.iter().zip(ends) {
            near(point[0], x);
            near(point[1], y);
        }
    }
    // Circular cutouts keep their base coordinates; flange holes are
    // developed beyond the bend zone along each counterclockwise edge.
    let drills: Vec<(String, [f64; 2], f64)> = sheet
        .holes
        .iter()
        .filter_map(|hole| {
            let center = match hole.face {
                None => hole.center,
                Some(0) => [hole.center[0], -(allowance + hole.center[1])],
                Some(_) => [w - hole.center[0], h + allowance + hole.center[1]],
            };
            Some((hole.drill.clone()?, center, hole.diameter?))
        })
        .collect();
    assert_eq!(developed.drills.len(), drills.len());
    for (drill, (id, center, diameter)) in developed.drills.iter().zip(&drills) {
        assert_eq!(&drill.id, id);
        near(drill.center[0], center[0]);
        near(drill.center[1], center[1]);
        near(drill.diameter_mm, *diameter);
    }
    assert_eq!(developed.cut.loops.len(), 16);
    near(developed.area_mm2, flat_area);
    let mut areas: Vec<_> = developed
        .cut
        .loops
        .iter()
        .map(|wire| wire.signed_area())
        .collect();
    areas.sort_by(f64::total_cmp);
    let mut expected_areas = vec![blank];
    expected_areas.extend(sheet.holes.iter().map(|hole| -hole.area));
    expected_areas.sort_by(f64::total_cmp);
    for (a, b) in areas.into_iter().zip(expected_areas) {
        near(a, b);
    }
    let bytes = encode_developed(
        &developed,
        LengthUnit::Millimeter,
        &SheetDxfLayers::default(),
        compiler.policy.linear_tolerance_mm,
    )
    .unwrap();
    let dxf = ::dxf::Drawing::load(&mut Cursor::new(bytes)).unwrap();
    assert_eq!(
        dxf.header.default_drawing_units,
        ::dxf::enums::Units::Millimeters
    );
    for (layer, count) in [
        ("CUT", 1),
        ("CUT_INNER", 15),
        ("BEND_DOWN", 2),
        ("DRILL", 14),
    ] {
        assert_eq!(
            dxf.entities().filter(|e| e.common.layer == layer).count(),
            count,
            "{layer}"
        );
    }
}

#[test]
fn formed_sheet_refuses_overlapping_openings_and_a_cut_crossing_a_bend() {
    let (document, part, feature) = front_panel();
    let session = Session::new(document).unwrap();
    let source = model::read_sheet(session.document(), feature).unwrap();
    let duplicate = source.cutouts[0].clone();
    let across_bend = Cutout {
        face: CutoutFace::Flange { index: 0 },
        shape: CutoutShape::Circle {
            center: [length(50.), length(0.5)],
            diameter: length(4.),
        },
    };
    for opening in [duplicate, across_bend] {
        let mut edited = source.clone();
        edited.cutouts.push(opening);
        let prepared =
            model::prepare_sheet_edit(&session, 0, Actor::local(Channel::Ui), feature, &edited)
                .unwrap();
        let error = compiler()
            .compile_sources_with_backend(
                prepared.prepared.document(),
                1,
                &[part],
                GeometryBackend::Rust,
                None,
            )
            .err()
            .expect("the opening must stay in disjoint material inside one face");
        assert_eq!(error.entity, Some(feature), "{error}");
    }
}

struct Hole {
    /// None on the base, else the flange index.
    face: Option<usize>,
    center: [f64; 2],
    diameter: Option<f64>,
    drill: Option<String>,
    area: f64,
    perimeter: f64,
}

/// The saved front-panel sheet definition with every dimension evaluated.
struct SavedSheet {
    t: f64,
    r: f64,
    k: f64,
    w: f64,
    h: f64,
    leg: f64,
    holes: Vec<Hole>,
}

impl SavedSheet {
    fn read(document: &Document, feature: EntityId) -> Self {
        use multicad_core::{
            constraint_sketch::Dimension,
            model::{BendDirection, FlangeKind, HardwarePreset, PanelPlane, SheetBase, Thread},
            parameters::{self, Unit},
        };
        let values = parameters::evaluate(document).unwrap();
        let mm = |d: &Dimension| d.value_in(Unit::Millimeter, Some(&values)).unwrap();
        let sheet = model::read_sheet(document, feature).unwrap();
        let (t, _, _) =
            model::sheet_stock_values(&model::read_sheet_stock(document, sheet.stock).unwrap())
                .unwrap();
        let SheetBase::Rectangle {
            plane: PanelPlane::Xz,
            origin,
            length,
            width,
        } = &sheet.base
        else {
            panic!("Expected a rectangular XZ base");
        };
        assert!(origin.iter().all(|d| mm(d) == 0.));
        assert!(!sheet.reverse);
        assert_eq!(sheet.flanges.len(), 2);
        let mut legs = vec![];
        for (flange, edge) in sheet.flanges.iter().zip([0, 2]) {
            assert_eq!(flange.edge, edge);
            assert_eq!(flange.direction, BendDirection::Down);
            let FlangeKind::Bend { angle } = &flange.kind else {
                panic!("Expected a bend flange");
            };
            assert_eq!(angle.value_in(Unit::Degree, Some(&values)).unwrap(), 90.);
            assert!(flange.inside_radius.is_none() && flange.relief.is_none());
            assert!(flange.start_inset.is_none() && flange.end_inset.is_none());
            legs.push(mm(&flange.length));
        }
        assert_eq!(legs[0], legs[1]);
        let mut holes = vec![];
        for (i, cutout) in sheet.cutouts.iter().enumerate() {
            assert_eq!(cutout.face, model::CutoutFace::Base);
            holes.push(match &cutout.shape {
                model::CutoutShape::Circle { center, diameter } => {
                    let d = mm(diameter);
                    Hole {
                        face: None,
                        center: center.each_ref().map(mm),
                        diameter: Some(d),
                        drill: Some(format!("cutout.{i}")),
                        area: PI * d * d / 4.,
                        perimeter: PI * d,
                    }
                }
                model::CutoutShape::Rectangle { center, size } => {
                    let [a, b] = size.each_ref().map(mm);
                    Hole {
                        face: None,
                        center: center.each_ref().map(mm),
                        diameter: None,
                        drill: None,
                        area: a * b,
                        perimeter: 2. * (a + b),
                    }
                }
            });
        }
        for (p, pattern) in sheet.fasteners.iter().enumerate() {
            let thread = match &pattern.preset {
                HardwarePreset::Screw { thread, .. } | HardwarePreset::Standoff { thread, .. } => {
                    *thread
                }
                other => panic!("Unexpected hardware {other:?}"),
            };
            // ISO 273 medium clearance for M3.
            assert_eq!(thread, Thread::M3);
            let d = 3.4;
            let face = match pattern.face {
                model::CutoutFace::Base => None,
                model::CutoutFace::Flange { index } => Some(index as usize),
            };
            for (j, point) in pattern.points.iter().enumerate() {
                holes.push(Hole {
                    face,
                    center: point.each_ref().map(mm),
                    diameter: Some(d),
                    drill: Some(format!("fastener.{p}.{j}")),
                    area: PI * d * d / 4.,
                    perimeter: PI * d,
                });
            }
        }
        let k = sheet.k_factor.as_ref().unwrap();
        Self {
            t,
            r: mm(sheet.inside_radius.as_ref().unwrap()),
            k: k.value_in(Unit::Ratio, Some(&values)).unwrap(),
            w: mm(length),
            h: mm(width),
            leg: legs[0],
            holes,
        }
    }
}
