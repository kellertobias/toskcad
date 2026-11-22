//! CAD-346 acceptance through the normal desktop drawing worker for the body
//! families that draw through the exact-curve scene engine: twisted, internal,
//! bevel and worm gears, an overlapping gear and box scene, an oblique bore
//! and the thread, screw and pipe treatments. Every family is authored with
//! the ordinary model commands, drawn in side, hidden-line and through-axis
//! section views, edited, cancelled, superseded, saved and reopened, and a
//! page of these views is exported to PDF, DXF and SVG and read back.
use super::tests::preview;
use super::*;
use multicad_core::{
    commands::{Actor, Channel, Session},
    constraint_sketch::Dimension,
    document::Document,
    drawings::{self, VectorRole, ViewEvaluation},
    materials::*,
    model::{
        self, BevelDimensions, BevelToothSystem, BodyInput, GearFamily, GearHand, GearModuleSystem,
        GearParameters, GearRequest, MechanicalPartParameters, MechanicalPartRequest,
        PrimitiveFrame, PrimitiveParameters, PrimitiveRequest, ThreadEndTreatment,
        ThreadEndTreatments, ThreadFamily, ThreadFit, ThreadMode, ThreadParameters, ThreadRequest,
        ThreadStartTreatment, WormDimensions, WormMember, WormToothSystem, angle, length, ratio,
    },
    references::ReferenceFrame,
};
use multicad_model::{Compiler, GeometryBackend};
use multicad_precision::PortableMath as _;
use std::time::{Duration, Instant};

fn actor() -> Actor {
    Actor::local(Channel::Ui)
}

fn count(n: f64) -> Dimension {
    Dimension::Literal {
        value: n,
        unit: "count".into(),
    }
}

fn at(x: f64) -> PrimitiveFrame {
    PrimitiveFrame {
        origin: [length(x), length(0.), length(0.)],
        ..Default::default()
    }
}

/// Horizontal spacing between families so that only the deliberate scene
/// overlaps.
const PITCH: f64 = 100.;

/// What a family proves in its views beyond well-formed vectors.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Check {
    /// Rotary gear about +Z: side outline spans ±tip, through-axis section
    /// cuts the bore walls at ±bore along the whole face width.
    Gear { tip: f64, bore: f64, face: f64 },
    /// The side outline spans exactly `u` across the view frame.
    Outline { u: [f64; 2] },
    /// A fitted thread whose fitted crest lies within its nominal radius.
    Within { radius: f64 },
    /// No closed-form extents; well-formed and source-attributed only.
    Generic,
}

struct Family {
    name: &'static str,
    /// The drawn bodies; the first one defines the view frame.
    sources: Vec<EntityId>,
    check: Check,
    /// Whether the side view has hidden lines: rear sharp edges such as
    /// bore rims, flat thread crests and treatment rims. Fully rounded thread
    /// forms have only tangent-continuous rear faces and so none.
    hidden: bool,
    hatched: bool,
}

struct Fixture {
    session: Session,
    families: Vec<Family>,
    helical: GearRequest,
    relief: ThreadRequest,
    book: EntityId,
    page: drawings::Page,
    export_page: drawings::Page,
}

fn commit_gear(session: &mut Session, name: &str, parameters: GearParameters) -> GearRequest {
    let request = GearRequest::new(name, parameters);
    let prepared = model::prepare_gear(session, session.revision(), actor(), &request)
        .unwrap_or_else(|e| panic!("{name}: {e:?}"));
    session.commit(prepared.prepared).unwrap();
    request
}

fn commit_thread(session: &mut Session, request: &ThreadRequest) {
    let prepared = model::prepare_thread(session, session.revision(), actor(), request)
        .unwrap_or_else(|e| panic!("{}: {e:?}", request.name));
    session.commit(prepared.prepared).unwrap();
}

fn commit_primitive(session: &mut Session, request: &PrimitiveRequest) -> BodyInput {
    let prepared = model::prepare_primitive(session, session.revision(), actor(), request)
        .unwrap_or_else(|e| panic!("{}: {e:?}", request.name));
    session.commit(prepared.prepared).unwrap();
    BodyInput {
        body: request.body,
        producer: request.feature,
    }
}

fn gear(family: GearFamily, x: f64) -> GearParameters {
    let twisted = family != GearFamily::Spur;
    GearParameters {
        family,
        face_width: length(if twisted { 8. } else { 4. }),
        helix_angle: angle(if twisted { 15. } else { 0. }),
        hand: if twisted {
            GearHand::Right
        } else {
            GearHand::None
        },
        center_gap: length(if family == GearFamily::DoubleHelical {
            2.
        } else {
            0.
        }),
        bore_diameter: length(6.),
        frame: at(x),
        ..Default::default()
    }
}

fn worm(wheel: bool, x: f64) -> GearParameters {
    GearParameters {
        family: GearFamily::Worm,
        worm: Some(WormDimensions {
            role: if wheel {
                WormMember::Wheel
            } else {
                WormMember::Worm
            },
            tooth_system: WormToothSystem::ZaAxialStraight,
            mate_teeth: count(if wheel { 1. } else { 32. }),
            diameter_quotient: ratio(10.),
        }),
        module_system: GearModuleSystem::Axial,
        module: length(1.),
        teeth: count(if wheel { 32. } else { 1. }),
        hand: GearHand::Right,
        face_width: length(if wheel { 2. } else { 20. }),
        bore_diameter: length(if wheel { 6. } else { 2. }),
        root_fillet_ratio: ratio(0.),
        frame: at(x),
        ..Default::default()
    }
}

fn metric(x: f64, ends: ThreadEndTreatments) -> ThreadParameters {
    ThreadParameters {
        major_diameter: length(6.),
        pitch: length(1.),
        length: length(6.6),
        fit: ThreadFit::Metric6gExternal,
        end_treatments: Some(ends),
        frame: at(x),
        ..Default::default()
    }
}

fn relief(width: f64) -> ThreadEndTreatments {
    ThreadEndTreatments {
        start: ThreadStartTreatment::FlatCut,
        end: ThreadEndTreatment::CylindricalRelief {
            length: length(width),
            radial_depth: length(0.1),
        },
    }
}

/// A section style so the view shows that material hatching follows the
/// canonical assignment.
fn assign_hatched_material(session: &mut Session, bodies: &[EntityId]) {
    let material = multicad_core::document::Entity::new(
        "Steel",
        multicad_core::document::EntityData::Material {
            density_kg_m3: None,
            color_rgba: [120, 120, 130, 255],
            section_hatch: None,
        },
    );
    let hatch = LibraryDefinition::Hatch(HatchStyle {
        color_rgb: None,
        pattern: HatchPattern::Lines,
        angle_deg: 45.,
        spacing_mm: 1.5,
        line_width_mm: 0.18,
    })
    .entity(EntityId::new(), "Section")
    .unwrap();
    let prepared = prepare_library(
        session,
        session.revision(),
        actor(),
        &[material.clone(), hatch.clone()],
    )
    .unwrap();
    session.commit(prepared).unwrap();
    let style = LibraryDefinition::MaterialStyle(MaterialStyle {
        material: material.id,
        appearance: None,
        hatch: Some(LibraryRef::capture(session.document(), hatch.id).unwrap()),
    })
    .entity(EntityId::new(), "Defaults")
    .unwrap();
    let prepared = prepare_library(session, session.revision(), actor(), &[style]).unwrap();
    session.commit(prepared).unwrap();
    let prepared = prepare_assignment(
        session,
        session.revision(),
        actor(),
        bodies,
        &Assignment {
            material: Some(material.id),
            stock: None,
            appearance: None,
            hatch: None,
        },
    )
    .unwrap();
    session.commit(prepared).unwrap();
}

/// The three views drawn for every family: a side projection, the same with
/// hidden lines and a retained section through the family's axis.
const VIEWS: [&str; 3] = ["side", "side hidden", "axis section"];

fn view(
    name: String,
    sources: Vec<EntityId>,
    center: [f64; 3],
    kind: usize,
) -> drawings::DrawingView {
    drawings::DrawingView {
        id: EntityId::new(),
        name,
        source: drawings::ViewSource::Model { entities: sources },
        frame: ReferenceFrame {
            origin_mm: center,
            u: [1., 0., 0.],
            v: [0., 0., 1.],
        },
        kind: if kind == 2 {
            drawings::ViewKind::Section {
                offset_mm: 0.,
                range: drawings::SectionRange::Unlimited,
                window: Default::default(),
            }
        } else {
            drawings::ViewKind::Projection
        },
        position_mm: [0.; 2],
        scale: drawings::Scale {
            numerator: 1,
            denominator: 1,
        },
        hidden_lines: kind == 1,
        tangent_edges: false,
        paper_chord_tolerance_mm: 0.01,
        guide: None,
        equipment_categories: None,
    }
}

fn fixture() -> Fixture {
    let mut session = Session::new(Document::new("Scene families").unwrap()).unwrap();
    let mut families = Vec::new();
    let mut x = 0.;
    let mut next = || {
        x += PITCH;
        x
    };

    // Spur gear overlapped by a box: one scene with two sources.
    let spur_x = next();
    let spur = commit_gear(&mut session, "Spur", gear(GearFamily::Spur, spur_x));
    let block = commit_primitive(
        &mut session,
        &PrimitiveRequest::new(
            "Overlapping block",
            PrimitiveParameters::Box {
                size: [length(20.), length(10.), length(6.)],
                frame: PrimitiveFrame {
                    origin: [length(spur_x + 18.), length(-5.), length(-1.)],
                    ..Default::default()
                },
            },
        ),
    );
    families.push(Family {
        name: "spur gear + box overlapping",
        sources: vec![spur.body, block.body],
        // Module 2, 24 teeth: tip radius 26; the block ends 38 mm right of
        // the gear axis.
        check: Check::Outline { u: [-26., 38.] },
        hidden: true,
        hatched: false,
    });

    let ring = commit_gear(
        &mut session,
        "Ring",
        GearParameters {
            family: GearFamily::InternalRing,
            module: length(1.),
            teeth: count(48.),
            face_width: length(3.),
            bore_diameter: length(0.),
            outer_diameter: Some(length(56.)),
            frame: at(next()),
            ..Default::default()
        },
    );
    families.push(Family {
        name: "internal ring",
        sources: vec![ring.body],
        check: Check::Outline { u: [-28., 28.] },
        hidden: true,
        hatched: false,
    });

    let helical = commit_gear(&mut session, "Helical", gear(GearFamily::Helical, next()));
    let herringbone = commit_gear(
        &mut session,
        "Herringbone",
        gear(GearFamily::Herringbone, next()),
    );
    let double = commit_gear(
        &mut session,
        "Double helical",
        gear(GearFamily::DoubleHelical, next()),
    );
    // Normal module 2, 24 teeth, 15 degree helix: transverse pitch radius
    // 24 / cos 15°, tip one module beyond it.
    let tip = 24. / 15_f64.to_radians().portable_cos() + 2.;
    for (name, request) in [
        ("helical gear", &helical),
        ("herringbone gear", &herringbone),
        ("double-helical gear", &double),
    ] {
        families.push(Family {
            name,
            sources: vec![request.body],
            check: Check::Gear {
                tip,
                bore: 3.,
                face: 8.,
            },
            hidden: true,
            hatched: name == "helical gear",
        });
    }
    let bevel = commit_gear(
        &mut session,
        "Bevel",
        GearParameters {
            family: GearFamily::Bevel,
            teeth: count(24.),
            face_width: length(6.),
            bevel: Some(BevelDimensions {
                tooth_system: BevelToothSystem::SphericalInvoluteStraight,
                mate_teeth: count(32.),
                shaft_angle: angle(90.),
            }),
            frame: at(next()),
            ..Default::default()
        },
    );
    families.push(Family {
        name: "bevel gear",
        sources: vec![bevel.body],
        // Outer tip circle: pitch radius 24 plus one module times cos δ,
        // with tan δ = 24 / 32.
        check: Check::Outline { u: [-25.6, 25.6] },
        hidden: true,
        hatched: false,
    });
    let worm_gear = commit_gear(&mut session, "Worm", worm(false, next()));
    let wheel = commit_gear(&mut session, "Worm wheel", worm(true, next()));
    families.push(Family {
        name: "za worm",
        sources: vec![worm_gear.body],
        // Axial module 1, quotient 10: tip radius 5 + 1.
        check: Check::Outline { u: [-6., 6.] },
        hidden: true,
        hatched: false,
    });
    families.push(Family {
        name: "worm wheel",
        sources: vec![wheel.body],
        check: Check::Generic,
        hidden: true,
        hatched: false,
    });

    // A 10 mm cube with a 0.5 mm bore along its main diagonal.
    let cube_x = next();
    let cube = commit_primitive(
        &mut session,
        &PrimitiveRequest::new(
            "Cube",
            PrimitiveParameters::Box {
                size: [length(10.), length(10.), length(10.)],
                frame: at(cube_x),
            },
        ),
    );
    let (sx, ty) = (0.5_f64.sqrt(), (1_f64 / 6.).sqrt());
    let start = 1. / 3. - 3_f64.sqrt().recip();
    let drill = PrimitiveRequest::new(
        "Diagonal drill",
        PrimitiveParameters::Cylinder {
            radius: length(0.5),
            height: length(20.),
            frame: PrimitiveFrame {
                origin: [length(cube_x + start), length(start), length(start)],
                u: [sx, -sx, 0.],
                v: [ty, ty, -2. * ty],
            },
        },
    );
    let tool = commit_primitive(&mut session, &drill);
    let cut = model::BooleanRequest::new(model::Operation::Cut, cube.clone(), vec![tool]);
    let prepared = model::prepare_boolean(&session, session.revision(), actor(), &cut).unwrap();
    session.commit(prepared.prepared).unwrap();
    families.push(Family {
        name: "oblique bore",
        sources: vec![cube.body],
        check: Check::Outline { u: [-5., 5.] },
        hidden: true,
        hatched: false,
    });

    let relief_thread = ThreadRequest::new("M6 relieved shaft", metric(next(), relief(0.6)));
    commit_thread(&mut session, &relief_thread);
    families.push(Family {
        name: "thread end relief",
        sources: vec![relief_thread.body],
        // The 6g major diameter is drawn at its mid-tolerance 5.884 mm.
        check: Check::Outline { u: [-2.942, 2.942] },
        hidden: true,
        hatched: false,
    });

    let bsw = model::whitworth_thread_presets()
        .into_iter()
        .find(|p| p.family == ThreadFamily::Bsw)
        .unwrap();
    let rounded = ThreadRequest::new(
        "BSW rounded shaft",
        ThreadParameters {
            family: ThreadFamily::Bsw,
            major_diameter: length(bsw.diameter_mm),
            pitch: length(bsw.pitch_mm),
            length: length(bsw.pitch_mm * 2.37),
            fit: ThreadFit::WhitworthMediumExternal,
            frame: at(next()),
            ..Default::default()
        },
    );
    commit_thread(&mut session, &rounded);
    families.push(Family {
        name: "rounded thread",
        sources: vec![rounded.body],
        check: Check::Within {
            radius: bsw.diameter_mm / 2.,
        },
        hidden: false,
        hatched: false,
    });

    let chamfered = ThreadRequest::new(
        "M6 conical entry",
        metric(
            next(),
            ThreadEndTreatments {
                start: ThreadStartTreatment::Chamfer45 {
                    radial_depth: length(0.3),
                },
                end: ThreadEndTreatment::FlatCut,
            },
        ),
    );
    commit_thread(&mut session, &chamfered);
    families.push(Family {
        name: "thread entry chamfer",
        sources: vec![chamfered.body],
        check: Check::Outline { u: [-2.942, 2.942] },
        hidden: true,
        hatched: false,
    });

    let mut screw_parameters = model::default_machine_screw();
    let MechanicalPartParameters::MachineScrew {
        underhead_length,
        thread,
        ..
    } = &mut screw_parameters
    else {
        unreachable!()
    };
    *underhead_length = Some(length(6.));
    thread.length = length(6.);
    thread.frame = at(next());
    let screw = MechanicalPartRequest::new("Machine screw", screw_parameters);
    let prepared =
        model::prepare_mechanical_part(&session, session.revision(), actor(), &screw).unwrap();
    session.commit(prepared.prepared).unwrap();
    families.push(Family {
        name: "machine screw",
        sources: vec![screw.body],
        // The 5.68 mm socket-cap head is the widest part.
        check: Check::Outline { u: [-2.84, 2.84] },
        hidden: true,
        hatched: false,
    });

    let stock_x = next();
    let stock = commit_primitive(
        &mut session,
        &PrimitiveRequest::new(
            "Cylinder stock",
            PrimitiveParameters::Cylinder {
                radius: length(2.5),
                height: length(4.),
                frame: at(stock_x),
            },
        ),
    );
    let internal = ThreadRequest::receiving(
        "Treated M3 internal",
        ThreadMode::CutInternal,
        stock.clone(),
        ThreadParameters {
            major_diameter: length(3.),
            pitch: length(0.5),
            length: length(4.),
            fit: ThreadFit::Metric6HInternal,
            minimum_wall: length(0.5),
            end_treatments: Some(ThreadEndTreatments {
                start: ThreadStartTreatment::Chamfer45 {
                    radial_depth: length(0.5),
                },
                end: ThreadEndTreatment::CylindricalRelief {
                    length: length(0.5),
                    radial_depth: length(0.125),
                },
            }),
            frame: at(stock_x),
            ..Default::default()
        },
    )
    .unwrap();
    commit_thread(&mut session, &internal);
    families.push(Family {
        name: "internal thread treatment",
        sources: vec![stock.body],
        check: Check::Outline { u: [-2.5, 2.5] },
        hidden: true,
        hatched: true,
    });

    // An ISO 7 R 1/8 taper pipe thread. Posed pipe constructions are not
    // certified for stock clipping, so it keeps its authored frame at the
    // origin, clear of the other families.
    let pipe_thread = model::PipeThreadIntent {
        form: model::PipeThreadForm::R,
        size: model::PipeThreadSize::OneEighth,
        gauge: model::PipeGaugeConvention::Iso7Nominal1994,
    }
    .parameters();
    let pipe = ThreadRequest::new("R 1/8 taper", pipe_thread);
    commit_thread(&mut session, &pipe);
    families.push(Family {
        name: "axial pipe thread",
        sources: vec![pipe.body],
        check: Check::Generic,
        hidden: false,
        hatched: false,
    });

    let hatched: Vec<_> = families
        .iter()
        .filter(|f| f.hatched)
        .map(|f| f.sources[0])
        .collect();
    assign_hatched_material(&mut session, &hatched);

    // View frames sit on each family's axis at mid-height, from the exact
    // model bounds of its first source.
    let model = Compiler {
        omit_display: true,
        ..Default::default()
    }
    .compile_with_backend(
        session.document(),
        session.revision(),
        GeometryBackend::Rust,
        None,
    )
    .unwrap();
    // An A0 sheet at 1:1 with every family view on a 120 mm grid.
    let mut page = drawings::Page::new("1");
    page.settings.paper = drawings::PaperSize::A0;
    page.settings.orientation = drawings::Orientation::Landscape;
    for family in &families {
        let bounds = model.solid(family.sources[0]).unwrap().bounds().unwrap();
        let center = [0, 1, 2].map(|i| (bounds.min[i] + bounds.max[i]) / 2.);
        for (kind, label) in VIEWS.iter().enumerate() {
            let mut view = view(
                format!("{} {label}", family.name),
                family.sources.clone(),
                center,
                kind,
            );
            let cell = page.views.len();
            view.position_mm = [
                80. + 120. * (cell % 9) as f64,
                80. + 120. * (cell / 9) as f64,
            ];
            page.views.push(view);
        }
    }
    let export_page = export_layout(&page);
    let mut books = Vec::new();
    for (name, page) in [("Scene families", &page), ("Scene export", &export_page)] {
        let (book, prepared) = drawings::prepare_create(
            &session,
            session.revision(),
            actor(),
            name,
            Some(page.clone()),
        )
        .unwrap();
        session.commit(prepared).unwrap();
        books.push(book);
    }
    Fixture {
        session,
        families,
        helical,
        relief: relief_thread,
        book: books[1],
        page,
        export_page,
    }
}

/// An A3 landscape sheet at 1:1 holding the helical gear (side, hidden-line
/// side, hatched section and axial view), the worm wheel side view, the
/// oblique bore section and the hatched treated internal thread section.
fn export_layout(page: &drawings::Page) -> drawings::Page {
    let mut sheet = drawings::Page::new("E-1");
    sheet.settings.paper = drawings::PaperSize::A3;
    sheet.settings.orientation = drawings::Orientation::Landscape;
    for (name, position) in [
        ("helical gear side", [60., 60.]),
        ("helical gear side hidden", [140., 60.]),
        ("helical gear axis section", [220., 60.]),
        ("worm wheel side", [60., 130.]),
        ("oblique bore axis section", [140., 130.]),
        ("internal thread treatment axis section", [220., 130.]),
    ] {
        let mut view = page
            .views
            .iter()
            .find(|v| v.name == name)
            .unwrap_or_else(|| panic!("{name}"))
            .clone();
        view.id = EntityId::new();
        view.position_mm = position;
        sheet.views.push(view);
    }
    // An axial projection of the helical gear shows its bore circle.
    let mut axial = sheet.views[0].clone();
    axial.id = EntityId::new();
    axial.name = "helical gear axial".into();
    axial.frame.v = [0., 1., 0.];
    axial.position_mm = [320., 90.];
    sheet.views.push(axial);
    sheet
}

fn ready(result: &DrawingEvaluation, view: &drawings::DrawingView) -> Vec<drawings::DrawingVector> {
    match &result.views[&view.id] {
        ViewEvaluation::Ready { vectors, .. } => vectors.clone(),
        other => panic!("{}: {other:?}", view.name),
    }
}

/// Paper extents of the drawn preview points of `roles`.
fn extents(vectors: &[drawings::DrawingVector], roles: &[VectorRole]) -> [f64; 4] {
    let mut b = [
        f64::INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NEG_INFINITY,
    ];
    for p in vectors
        .iter()
        .filter(|v| roles.contains(&v.role))
        .flat_map(|v| &v.preview_mm)
    {
        b = [
            b[0].min(p[0]),
            b[1].min(p[1]),
            b[2].max(p[0]),
            b[3].max(p[1]),
        ];
    }
    b
}

/// Through-axis section points on a bore wall at `u = ±radius`, by side.
fn bore_walls(vectors: &[drawings::DrawingVector], radius: f64) -> [[f64; 2]; 2] {
    [-radius, radius].map(|u| {
        let mut span = [f64::INFINITY, f64::NEG_INFINITY];
        for p in vectors
            .iter()
            .filter(|v| v.role == VectorRole::Section)
            .flat_map(|v| &v.preview_mm)
            .filter(|p| (p[0] - u).abs() < 1e-6)
        {
            span = [span[0].min(p[1]), span[1].max(p[1])];
        }
        span
    })
}

fn check(result: &DrawingEvaluation, f: &Fixture) {
    for (family, views) in f.families.iter().zip(f.page.views.chunks(VIEWS.len())) {
        let [side, hidden, section] = [0, 1, 2].map(|i| ready(result, &views[i]));
        for (vectors, view) in [&side, &hidden, &section].into_iter().zip(views) {
            assert!(!vectors.is_empty(), "{}", view.name);
            for vector in vectors.iter() {
                assert!(family.sources.contains(&vector.source), "{}", view.name);
                assert!(
                    vector.preview_mm.len() >= 2
                        && vector.preview_mm.iter().flatten().all(|c| c.is_finite()),
                    "{}: {vector:?}",
                    view.name
                );
            }
            // Every source of a scene is drawn.
            for source in &family.sources {
                assert!(vectors.iter().any(|v| v.source == *source), "{}", view.name);
            }
        }
        let roles =
            |vectors: &[drawings::DrawingVector], role| vectors.iter().any(|v| v.role == role);
        for role in [VectorRole::Hidden, VectorRole::Section, VectorRole::Hatch] {
            assert!(!roles(&side, role), "{} side {role:?}", family.name);
        }
        assert!(roles(&side, VectorRole::Visible), "{}", family.name);
        assert_eq!(
            roles(&hidden, VectorRole::Hidden),
            family.hidden,
            "{} hidden lines",
            family.name
        );
        assert!(
            roles(&section, VectorRole::Section),
            "{} section",
            family.name
        );
        assert_eq!(
            roles(&section, VectorRole::Hatch),
            family.hatched,
            "{} hatch",
            family.name
        );
        // The hidden-line view adds hidden lines to the same visible outline.
        assert_eq!(
            extents(&side, &[VectorRole::Visible]),
            extents(&hidden, &[VectorRole::Visible]),
            "{}",
            family.name
        );
        let tolerance = views[0].paper_chord_tolerance_mm + 1e-6;
        let outline = extents(&side, &[VectorRole::Visible]);
        match family.check {
            Check::Gear { tip, bore, face } => {
                assert!(
                    (outline[0] + tip).abs() <= tolerance && (outline[2] - tip).abs() <= tolerance,
                    "{} tip {outline:?} vs {tip}",
                    family.name
                );
                assert!(
                    (outline[3] - outline[1] - face).abs() < 1e-6,
                    "{} face {outline:?}",
                    family.name
                );
                for wall in bore_walls(&section, bore) {
                    assert!(
                        (wall[1] - wall[0] - face).abs() < 1e-6,
                        "{} bore wall {wall:?}",
                        family.name
                    );
                }
            }
            Check::Outline { u } => {
                assert!(
                    (outline[0] - u[0]).abs() <= tolerance
                        && (outline[2] - u[1]).abs() <= tolerance,
                    "{} outline {outline:?} vs {u:?}",
                    family.name
                );
            }
            Check::Within { radius } => {
                assert!(
                    outline[0] >= -radius
                        && outline[2] <= radius
                        && outline[2] - outline[0] >= 2. * radius * 0.95,
                    "{} outline {outline:?} vs ±{radius}",
                    family.name
                );
            }
            Check::Generic => {}
        }
    }
}

fn stamp(session: &Session) -> DrawingStamp {
    DrawingStamp {
        document: session.document().id(),
        revision: session.revision(),
        source_hash: multicad_core::project_file::ContentHash::document(session.document())
            .unwrap(),
    }
}

#[test]
fn scene_families_draw_through_the_worker_edits_cancellation_reopen_and_exports() {
    use multicad_core::project_file::ProjectFile;
    use multicad_ui::{application::DocumentStorage, storage::SaveMode};

    let directory = tempfile::tempdir().unwrap();
    let cache = cache::Cache::isolated(directory.path().join("drawing-cache"));
    let mut f = fixture();
    let mut worker = Backend::with_cache_and_backend(Some(cache.clone()), GeometryBackend::Rust);
    let started = Instant::now();
    let original = preview(&mut worker, &f.session, f.page.id);
    let cold = started.elapsed();
    check(&original, &f);
    let started = Instant::now();
    assert_eq!(
        preview(&mut worker, &f.session, f.page.id).views,
        original.views
    );
    let warm = started.elapsed();

    // A cancelled request publishes nothing, and the worker accepts the next one.
    worker
        .request_page(f.session.document(), stamp(&f.session), f.page.id)
        .unwrap();
    worker.cancel();
    let deadline = Instant::now() + Duration::from_millis(500);
    while Instant::now() < deadline {
        assert!(worker.poll().is_none(), "a cancelled page must stay silent");
        std::thread::sleep(Duration::from_millis(5));
    }

    // A request superseded by source edits is never delivered: widen the
    // helical bore and lengthen the thread end relief.
    worker
        .request_page(f.session.document(), stamp(&f.session), f.page.id)
        .unwrap();
    let stale = stamp(&f.session);
    let mut helical = f.helical.parameters.clone();
    helical.bore_diameter = length(8.);
    let prepared = model::prepare_gear_edit(
        &f.session,
        f.session.revision(),
        actor(),
        f.helical.feature,
        &helical,
    )
    .unwrap();
    f.session.commit(prepared.prepared).unwrap();
    let mut relief_parameters = f.relief.parameters.clone();
    relief_parameters.end_treatments = Some(relief(1.1));
    let prepared = model::prepare_thread_edit(
        &f.session,
        f.session.revision(),
        actor(),
        f.relief.feature,
        &relief_parameters,
    )
    .unwrap();
    f.session.commit(prepared.prepared).unwrap();
    let edited = preview(&mut worker, &f.session, f.page.id);
    assert_ne!(edited.stamp, stale);
    assert_eq!(edited.stamp.revision, f.session.revision());
    assert_ne!(edited.stamp.source_hash, original.stamp.source_hash);
    let helical_index = f
        .families
        .iter()
        .position(|family| family.name == "helical gear")
        .unwrap();
    f.families[helical_index].check = match f.families[helical_index].check {
        Check::Gear { tip, face, .. } => Check::Gear {
            tip,
            bore: 4.,
            face,
        },
        other => other,
    };
    check(&edited, &f);
    // The bore is invisible in the plain side view of the helical gear; its
    // hidden lines and section change. The longer relief changes all three.
    let edited_family = |view: &drawings::DrawingView| {
        (view.name.starts_with("helical gear") && view.name != "helical gear side")
            || view.name.starts_with("thread end relief")
    };
    for view in &f.page.views {
        if edited_family(view) {
            assert_ne!(
                edited.views[&view.id], original.views[&view.id],
                "{}",
                view.name
            );
        } else {
            assert_eq!(
                edited.views[&view.id], original.views[&view.id],
                "{}",
                view.name
            );
        }
    }
    // The widened bore no longer cuts at the original radius.
    let section = ready(&edited, &f.page.views[helical_index * VIEWS.len() + 2]);
    assert!(bore_walls(&section, 3.).iter().all(|w| w[0] > w[1]));
    // The longer relief leaves less threaded material on the side outline.
    let relief_index = f
        .families
        .iter()
        .position(|family| family.name == "thread end relief")
        .unwrap();
    let length_of = |result: &DrawingEvaluation| {
        ready(result, &f.page.views[relief_index * VIEWS.len()])
            .iter()
            .map(|v| {
                v.preview_mm
                    .windows(2)
                    .map(|w| (w[1][0] - w[0][0]).portable_hypot(w[1][1] - w[0][1]))
                    .sum::<f64>()
            })
            .sum::<f64>()
    };
    assert!(length_of(&edited) < length_of(&original));
    drop(worker);

    // Save, reopen from disk and restore the edited page from the cache.
    let path = directory.path().join("scene-families.tcd");
    let mut store = crate::files::ProjectStore::default();
    store
        .save_project(
            &path,
            &ProjectFile::new(f.session.document().clone()).unwrap(),
            SaveMode::Create,
        )
        .unwrap();
    let reopened = Session::new(store.open_project(&path).unwrap().file.document).unwrap();
    assert_eq!(reopened.document(), f.session.document());
    let mut cached = Backend::with_cache_and_backend(Some(cache), GeometryBackend::Rust);
    let started = Instant::now();
    let restored = preview(&mut cached, &reopened, f.page.id);
    let cached_reopen = started.elapsed();
    assert_eq!(restored.views, edited.views);
    assert!(
        cached
            .runtime
            .as_ref()
            .unwrap()
            .last_compile_stats()
            .is_none()
    );
    // A fresh worker without a cache rebuilds the identical page.
    let mut fresh = Backend::with_cache_and_backend(None, GeometryBackend::Rust);
    let started = Instant::now();
    assert_eq!(
        preview(&mut fresh, &reopened, f.page.id).views,
        edited.views
    );
    let uncached_reopen = started.elapsed();
    exports(reopened, &f);
    let views = f.page.views.len() as f64;
    println!(
        "{{\"case\":\"CAD-346 scene families page\",\"views\":{views},\"cold_ms\":{},\"warm_ms\":{},\"cached_reopen_ms\":{},\"uncached_reopen_ms\":{},\"cold_ms_per_view\":{}}}",
        cold.as_secs_f64() * 1e3,
        warm.as_secs_f64() * 1e3,
        cached_reopen.as_secs_f64() * 1e3,
        uncached_reopen.as_secs_f64() * 1e3,
        cold.as_secs_f64() * 1e3 / views,
    );
}

/// One drawn paper segment read back from an exported file, in paper
/// millimeters with y downward, and its layer where the format keeps one.
#[derive(Clone, Debug)]
struct Segment {
    layer: String,
    a: [f64; 2],
    b: [f64; 2],
}

/// Paper rectangle `[x0, y0, x1, y1]` of a view's drawn vectors of `roles`.
fn paper_box(
    view: &drawings::DrawingView,
    vectors: &[drawings::DrawingVector],
    roles: &[VectorRole],
) -> [f64; 4] {
    let scale = view.scale.numerator as f64 / view.scale.denominator as f64;
    let [u0, v0, u1, v1] = extents(vectors, roles);
    let [x, y] = view.position_mm;
    [
        x + u0 * scale,
        y - v1 * scale,
        x + u1 * scale,
        y - v0 * scale,
    ]
}

fn inside(r: [f64; 4], p: [f64; 2], slack: f64) -> bool {
    p[0] >= r[0] - slack && p[0] <= r[2] + slack && p[1] >= r[1] - slack && p[1] <= r[3] + slack
}

/// Extents of the file segments on `layers` that lie within `region`.
fn file_box(segments: &[Segment], layers: &[&str], region: [f64; 4]) -> ([f64; 4], usize) {
    let mut b = [
        f64::INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NEG_INFINITY,
    ];
    let mut n = 0;
    for s in segments
        .iter()
        .filter(|s| layers.contains(&s.layer.as_str()))
    {
        if inside(region, s.a, 0.01) && inside(region, s.b, 0.01) {
            n += 1;
            for p in [s.a, s.b] {
                b = [
                    b[0].min(p[0]),
                    b[1].min(p[1]),
                    b[2].max(p[0]),
                    b[3].max(p[1]),
                ];
            }
        }
    }
    (b, n)
}

fn segments_of(points: &[[f64; 2]], closed: bool, layer: &str, out: &mut Vec<Segment>) {
    for w in points.windows(2) {
        out.push(Segment {
            layer: layer.into(),
            a: w[0],
            b: w[1],
        });
    }
    if closed && points.len() > 2 {
        out.push(Segment {
            layer: layer.into(),
            a: points[points.len() - 1],
            b: points[0],
        });
    }
}

fn svg_segments(svg: &str) -> Vec<Segment> {
    let mut out = Vec::new();
    for element in svg
        .lines()
        .filter(|l| l.starts_with("<polyline") || l.starts_with("<polygon"))
    {
        let attribute = |name: &str| {
            element
                .split(&format!(" {name}=\""))
                .nth(1)
                .and_then(|rest| rest.split('"').next())
                .unwrap()
        };
        let points: Vec<[f64; 2]> = attribute("points")
            .split_whitespace()
            .map(|pair| {
                let (x, y) = pair.split_once(',').unwrap();
                [x.parse().unwrap(), y.parse().unwrap()]
            })
            .collect();
        segments_of(
            &points,
            element.starts_with("<polygon"),
            attribute("class"),
            &mut out,
        );
    }
    out
}

/// Stroked PDF paths `w [dash] 0 d x y m x y l … S|s`, converted back to
/// paper millimeters. PDF keeps no layers: dashed strokes are reported as
/// "dashed", others by their stroke width in millimeters.
fn pdf_segments(text: &str, height: f64) -> Vec<Segment> {
    let to_mm = |pt: f64| pt * 25.4 / 72.;
    let mut out = Vec::new();
    for line in text.lines() {
        let tokens: Vec<_> = line.split_whitespace().collect();
        let Some(end) = tokens.iter().position(|t| *t == "S" || *t == "s") else {
            continue;
        };
        let Some(width) = tokens.iter().position(|t| *t == "w") else {
            continue;
        };
        let dashed = tokens
            .get(width + 1)
            .is_some_and(|t| t.starts_with('[') && *t != "[]");
        let mut points = Vec::new();
        for i in 2..end {
            if (tokens[i] == "m" || tokens[i] == "l")
                && let (Ok(x), Ok(y)) = (tokens[i - 2].parse::<f64>(), tokens[i - 1].parse::<f64>())
            {
                points.push([to_mm(x), height - to_mm(y)]);
            }
        }
        let layer = if dashed {
            "dashed".to_owned()
        } else {
            format!("{:.3}", to_mm(tokens[width - 1].parse().unwrap()))
        };
        segments_of(&points, tokens[end] == "s", &layer, &mut out);
    }
    out
}

fn dxf_segments(drawing: &::dxf::Drawing, height: f64) -> Vec<Segment> {
    drawing
        .entities()
        .filter_map(|e| match &e.specific {
            ::dxf::entities::EntityType::Line(l) => Some(Segment {
                layer: e.common.layer.clone(),
                a: [l.p1.x, height - l.p1.y],
                b: [l.p2.x, height - l.p2.y],
            }),
            _ => None,
        })
        .collect()
}

/// Render the export sheet through the desktop export backend (PDF and
/// DXF, written by the export workspace's native worker) and the SVG page
/// writer, then read each file back.
fn exports(mut session: Session, f: &Fixture) {
    use multicad_core::exchange::*;
    use multicad_ui::exchange::ExportBackend;
    let sheet = &f.export_page;
    let set = EntityId::new();
    let definition = SetDefinition {
        pages: Pages::Selected {
            pages: vec![PageRef {
                book: f.book,
                page: sheet.id,
            }],
        },
        recipes: vec![
            Recipe {
                id: EntityId::new(),
                name: "PDF".into(),
                output: RecipeOutput::Pdf,
            },
            Recipe {
                id: EntityId::new(),
                name: "DXF".into(),
                output: RecipeOutput::DrawingDxf {
                    units: LengthUnit::Millimeter,
                    layers: Layers::default(),
                },
            },
        ],
        ..Default::default()
    };
    let prepared = prepare_set(
        &session,
        session.revision(),
        actor(),
        set,
        "Scene sheet",
        &definition,
    )
    .unwrap();
    session.commit(prepared).unwrap();

    // The worker evaluation of the sheet the SVG is written from.
    let mut worker = Backend::with_cache_and_backend(None, GeometryBackend::Rust);
    let evaluation = preview(&mut worker, &session, sheet.id);
    let document = session.document();
    let revision = session.revision();
    let svg = drawings::svg::page_svg(document, revision, f.book, sheet.id, &evaluation).unwrap();
    let render =
        drawings::render::render_page(document, revision, f.book, sheet.id, &evaluation).unwrap();

    let mut backend = crate::exchange::Backend::with_geometry_backend(GeometryBackend::Rust);
    let wait = |backend: &mut crate::exchange::Backend| {
        let deadline = Instant::now() + Duration::from_secs(600);
        loop {
            if let Some(reply) = backend.poll() {
                return reply.result;
            }
            assert!(Instant::now() < deadline, "Drawing export did not reply");
            std::thread::sleep(Duration::from_millis(2));
        }
    };
    backend.preview(document, revision, set, actor()).unwrap();
    let ExportResult::Preview(Ok(export_preview)) = wait(&mut backend) else {
        panic!("Rust drawing export preview failed");
    };
    assert!(export_preview.can_run(), "{:?}", export_preview.recipes);
    let root = tempfile::tempdir().unwrap();
    backend
        .run(
            document,
            revision,
            &export_preview.key,
            root.path(),
            actor(),
        )
        .unwrap();
    let ExportResult::Run(Ok(receipt)) = wait(&mut backend) else {
        panic!("Rust drawing export failed");
    };
    let directory = std::path::Path::new(&receipt.directory);
    let file = |extension: &str| {
        let written = receipt
            .files
            .iter()
            .find(|file| file.artifact.output.filename.ends_with(extension))
            .unwrap_or_else(|| panic!("{extension}: {receipt:?}"));
        directory.join(&written.artifact.output.filename)
    };
    let pdf_path = file(".pdf");
    let dxf_path = file(".dxf");
    let svg_path = root.path().join("scene-sheet.svg");
    std::fs::write(&svg_path, &svg).unwrap();
    if let Some(keep) = std::env::var_os("CAD346_ARTIFACTS") {
        std::fs::create_dir_all(&keep).unwrap();
        for (path, name) in [
            (&pdf_path, "scene-sheet.pdf"),
            (&dxf_path, "scene-sheet.dxf"),
            (&svg_path, "scene-sheet.svg"),
        ] {
            std::fs::copy(path, std::path::Path::new(&keep).join(name)).unwrap();
        }
    }

    // Sheet size: A3 landscape at 1:1 in every file.
    let [width, height] = render.size_mm;
    assert_eq!([width, height], [420., 297.]);
    let svg = std::fs::read_to_string(&svg_path).unwrap();
    assert!(svg.contains(&format!(r#"width="{width}mm" height="{height}mm""#)));
    let pdf = std::fs::read(&pdf_path).unwrap();
    assert!(pdf.starts_with(b"%PDF-"));
    let pdf = String::from_utf8_lossy(&pdf).into_owned();
    let media: Vec<f64> = pdf
        .split("/MediaBox [")
        .nth(1)
        .and_then(|rest| rest.split(']').next())
        .unwrap()
        .split_whitespace()
        .map(|v| v.parse().unwrap())
        .collect();
    let to_mm = |pt: f64| pt * 25.4 / 72.;
    assert!(
        (to_mm(media[2]) - width).abs() < 1e-3 && (to_mm(media[3]) - height).abs() < 1e-3,
        "{media:?}"
    );
    assert_eq!(
        pdf.matches("/Type /Page ").count() + pdf.matches("/Type /Page\n").count(),
        1
    );
    // The export worker compiles the sheet independently of the drawing
    // worker; both publish the identical page.
    assert_eq!(
        std::fs::read_to_string(&dxf_path).unwrap(),
        drawings::dxf::render_dxf(&render, 1., 4)
    );
    assert_eq!(
        std::fs::read(&pdf_path).unwrap(),
        drawings::pdf::render_pdf(std::slice::from_ref(&render), document.name())
    );
    let dxf = ::dxf::Drawing::load_file(&dxf_path).unwrap();
    assert_eq!(dxf.header.maximum_drawing_extents.x, width);
    assert_eq!(dxf.header.maximum_drawing_extents.y, height);

    let files = [
        ("SVG", svg_segments(&svg), 1e-3),
        ("PDF", pdf_segments(&pdf, height), 1e-3),
        ("DXF", dxf_segments(&dxf, height), 1e-6),
    ];
    // The render's own segments, per layer, are what every file must carry.
    let mut expected = Vec::new();
    for line in &render.polylines {
        segments_of(&line.points, line.closed, line.layer.name(), &mut expected);
    }
    let pdf_layer = |layer: &str| -> String {
        let line = render
            .polylines
            .iter()
            .find(|l| l.layer.name() == layer)
            .unwrap_or_else(|| panic!("{layer}"));
        if line.dashed {
            "dashed".into()
        } else {
            format!("{:.3}", line.width_mm)
        }
    };
    let dxf_layer = |layer: &str| -> &'static str {
        match layer {
            "visible" => "VIEW_VISIBLE",
            "hidden" => "VIEW_HIDDEN",
            "section" => "VIEW_SECTION",
            "hatch" => "VIEW_HATCH",
            other => panic!("{other}"),
        }
    };
    let file_layer = |format: &str, layer: &str| -> String {
        match format {
            "SVG" => layer.into(),
            "PDF" => pdf_layer(layer),
            _ => dxf_layer(layer).into(),
        }
    };

    let hatched = [
        "helical gear axis section",
        "internal thread treatment axis section",
    ];
    for view in &sheet.views {
        let ViewEvaluation::Ready { vectors, .. } = &evaluation.views[&view.id] else {
            panic!("{}", view.name);
        };
        let all = [
            VectorRole::Visible,
            VectorRole::Hidden,
            VectorRole::Section,
            VectorRole::Hatch,
        ];
        let region = paper_box(view, vectors, &all);
        for (role, layer) in [
            (VectorRole::Visible, "visible"),
            (VectorRole::Hidden, "hidden"),
            (VectorRole::Section, "section"),
            (VectorRole::Hatch, "hatch"),
        ] {
            let drawn = vectors.iter().any(|v| v.role == role);
            let (want, want_count) = file_box(&expected, &[layer], region);
            assert_eq!(drawn, want_count > 0, "{} {layer}: render", view.name);
            if !drawn {
                continue;
            }
            // At 1:1 the paper extents of each role equal the model extents.
            let model = paper_box(view, vectors, &[role]);
            for i in 0..4 {
                assert!(
                    (want[i] - model[i]).abs() < 1e-9,
                    "{} {layer} {want:?} {model:?}",
                    view.name
                );
            }
            for (format, segments, tolerance) in &files {
                let layer = file_layer(format, layer);
                let (got, count) = file_box(segments, &[layer.as_str()], region);
                for i in 0..4 {
                    assert!(
                        (got[i] - want[i]).abs() <= *tolerance,
                        "{format} {} {layer}: {got:?} vs {want:?}",
                        view.name
                    );
                }
                // Every render segment is written once: neither missing
                // nor duplicated. PDF strokes of equal width share a key,
                // so its count is checked on the role's own key only.
                if *format != "PDF" {
                    assert_eq!(count, want_count, "{format} {} {layer}", view.name);
                }
            }
        }
        assert_eq!(
            vectors.iter().any(|v| v.role == VectorRole::Hatch),
            hatched.contains(&view.name.as_str()),
            "{} hatch",
            view.name
        );
    }
    // PDF: the whole sheet's strokes of each width or dash equal the render's.
    let pdf_segments = &files[1].1;
    let mut keys: Vec<_> = expected.iter().map(|s| pdf_layer(&s.layer)).collect();
    keys.sort();
    keys.dedup();
    for key in keys {
        let want = render
            .polylines
            .iter()
            .filter(|l| {
                (if l.dashed {
                    "dashed".into()
                } else {
                    format!("{:.3}", l.width_mm)
                }) == key
            })
            .map(|l| l.points.len() - 1 + usize::from(l.closed && l.points.len() > 2))
            .sum::<usize>();
        // Arcs (annotations and sheet frame) are stroked too; views have none.
        let got = pdf_segments.iter().filter(|s| s.layer == key).count();
        assert!(got >= want, "PDF {key}: {got} < {want}");
    }

    // Closed-form content at 1:1: the helical gear's tip extents, its 8 mm
    // bore as hidden and section walls and as a circle in the axial view,
    // and the hatched sections.
    let tip = 24. / 15_f64.to_radians().portable_cos() + 2.;
    let tolerance = 0.01 + 2e-3;
    let by_name = |name: &str| sheet.views.iter().find(|v| v.name == name).unwrap();
    for (format, segments, _) in &files {
        let side = by_name("helical gear side");
        let region = [
            side.position_mm[0] - 30.,
            side.position_mm[1] - 6.,
            side.position_mm[0] + 30.,
            side.position_mm[1] + 6.,
        ];
        let (outline, _) = file_box(segments, &[file_layer(format, "visible").as_str()], region);
        assert!(
            (outline[0] - (side.position_mm[0] - tip)).abs() <= tolerance
                && (outline[2] - (side.position_mm[0] + tip)).abs() <= tolerance
                && (outline[3] - outline[1] - 8.).abs() <= 2e-3,
            "{format} helical tip {outline:?}"
        );
        for (name, layer) in [
            ("helical gear side hidden", "hidden"),
            ("helical gear axis section", "section"),
        ] {
            let view = by_name(name);
            let layer = file_layer(format, layer);
            for x in [view.position_mm[0] - 4., view.position_mm[0] + 4.] {
                let wall: f64 = segments
                    .iter()
                    .filter(|s| s.layer == layer)
                    .filter(|s| (s.a[0] - x).abs() < 2e-3 && (s.b[0] - x).abs() < 2e-3)
                    .map(|s| (s.b[1] - s.a[1]).abs())
                    .sum();
                assert!(
                    (wall - 8.).abs() < 4e-3,
                    "{format} {name} bore wall at {x}: {wall}"
                );
            }
        }
        let axial = by_name("helical gear axial");
        let center = axial.position_mm;
        let ring = segments
            .iter()
            .filter(|s| s.layer == file_layer(format, "visible"))
            .filter(|s| {
                [s.a, s.b].iter().all(|p| {
                    ((p[0] - center[0]).portable_hypot(p[1] - center[1]) - 4.).abs() <= tolerance
                })
            })
            .count();
        assert!(ring > 16, "{format} bore circle: {ring} segments");
        for name in hatched {
            let view = by_name(name);
            let region = [
                view.position_mm[0] - 30.,
                view.position_mm[1] - 6.,
                view.position_mm[0] + 30.,
                view.position_mm[1] + 6.,
            ];
            let (_, count) = file_box(segments, &[file_layer(format, "hatch").as_str()], region);
            assert!(count > 0, "{format} {name} hatch");
        }
    }

    // Every view is labelled exactly once in each file.
    let captions: Vec<_> = render
        .texts
        .iter()
        .filter(|t| t.layer == drawings::render::Layer::ViewLabel)
        .map(|t| t.text.clone())
        .collect();
    assert_eq!(captions.len(), sheet.views.len(), "{captions:?}");
    let dxf_texts: Vec<_> = dxf
        .entities()
        .filter_map(|e| match &e.specific {
            ::dxf::entities::EntityType::Text(t) if e.common.layer == "VIEW_LABELS" => {
                Some(t.value.clone())
            }
            _ => None,
        })
        .collect();
    // Views whose caption starts with `name` followed by a space.
    let prefixed = |name: &str| {
        sheet
            .views
            .iter()
            .filter(|v| v.name == name || v.name.starts_with(&format!("{name} ")))
            .count()
    };
    for (caption, view) in captions.iter().zip(&sheet.views) {
        assert!(caption.starts_with(&view.name), "{caption}");
        // DXF R12 text is ASCII: the caption's dash becomes a hyphen.
        let ascii = caption.replace('—', "-");
        assert_eq!(
            dxf_texts.iter().filter(|t| **t == ascii).count(),
            1,
            "DXF {caption}"
        );
        let shared = |haystack: &str| haystack.matches(caption.as_str()).count();
        assert_eq!(
            shared(&svg),
            captions
                .iter()
                .filter(|c| c.contains(caption.as_str()))
                .count(),
            "SVG {caption}"
        );
        // PDF literals are WinAnsi; the ASCII view name opens the literal once.
        assert_eq!(
            pdf.matches(&format!("({} ", view.name)).count(),
            prefixed(&view.name),
            "PDF {caption}"
        );
    }
    if let Some(text) = std::process::Command::new("pdftotext")
        .args(["-layout", pdf_path.to_str().unwrap(), "-"])
        .output()
        .ok()
        .filter(|o| o.status.success())
    {
        let text = String::from_utf8_lossy(&text.stdout);
        for view in &sheet.views {
            assert_eq!(
                text.matches(&format!("{} ", view.name)).count(),
                prefixed(&view.name),
                "pdftotext {}",
                view.name
            );
        }
    }
    println!(
        "{{\"case\":\"CAD-346 scene export\",\"views\":{},\"svg_segments\":{},\"pdf_segments\":{},\"dxf_lines\":{},\"render_segments\":{}}}",
        sheet.views.len(),
        files[0].1.len(),
        files[1].1.len(),
        files[2].1.len(),
        expected.len()
    );
}

/// Per-view cost through the worker: each view on its own page, requested
/// from one uncached worker so that each request compiles only the view's
/// own sources. Run in release with `--ignored --nocapture`.
#[test]
#[ignore = "timing report"]
fn scene_family_view_timings() {
    let mut f = fixture();
    let mut pages = Vec::new();
    for view in &f.page.views {
        let mut page = drawings::Page::new(view.name.clone());
        page.settings = f.page.settings.clone();
        let mut single = view.clone();
        single.id = EntityId::new();
        page.views.push(single);
        let (_, prepared) = drawings::prepare_create(
            &f.session,
            f.session.revision(),
            actor(),
            &view.name,
            Some(page.clone()),
        )
        .unwrap();
        f.session.commit(prepared).unwrap();
        pages.push(page);
    }
    let mut worker = Backend::with_cache_and_backend(None, GeometryBackend::Rust);
    let mut slowest = (0., String::new());
    for page in &pages {
        let started = Instant::now();
        let result = preview(&mut worker, &f.session, page.id);
        let ms = started.elapsed().as_secs_f64() * 1e3;
        ready(&result, &page.views[0]);
        println!(
            "{{\"case\":\"CAD-346 scene view\",\"view\":\"{}\",\"ms\":{ms}}}",
            page.views[0].name
        );
        if ms > slowest.0 {
            slowest = (ms, page.views[0].name.clone());
        }
    }
    println!(
        "{{\"case\":\"CAD-346 slowest scene view\",\"view\":\"{}\",\"ms\":{}}}",
        slowest.1, slowest.0
    );
}
