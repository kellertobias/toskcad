use multicad_core::{
    commands::{Actor, Channel, Session},
    document::{EntityData, EntityId},
    model::{self, OpeningKind, SlabHole, VerticalCutout, WallFace, WallOpeningChange, length},
    project_file,
    references::TopologyKind,
};
use multicad_model::{Compiler, GeometryBackend, Model};
use multicad_precision::PortableMath as _;

fn apartment() -> Session {
    let mut bytes = include_bytes!("../../../examples/apartment.tcd").as_slice();
    let opened = project_file::read(
        &mut bytes,
        Default::default(),
        &project_file::supported_capabilities(),
    )
    .unwrap();
    Session::new(opened.project.document).unwrap()
}

fn architecture(session: &Session, type_id: &str) -> Vec<EntityId> {
    session
        .document()
        .entities()
        .filter_map(|entity| {
            let EntityData::Feature { definition, .. } = &entity.data else {
                return None;
            };
            (definition.type_id == type_id).then_some(entity.id)
        })
        .collect()
}

fn compile(session: &Session, sources: &[EntityId], backend: GeometryBackend) -> Model {
    Compiler::default()
        .compile_sources_with_backend(
            session.document(),
            session.revision(),
            sources,
            backend,
            None,
        )
        .unwrap()
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

fn face_has_role(model: &Model, feature: EntityId, role: &str) -> bool {
    model.snapshot().shapes[&feature]
        .topology
        .elements
        .iter()
        .filter(|element| element.key.kind == TopologyKind::Face)
        .flat_map(|element| &element.origins)
        .any(|origin| origin.source == feature && origin.role == role)
}

/// Independent prism expectation for an axis-aligned, center-justified wall
/// without junctions: footprint times height, less every full-depth opening.
fn expected_wall(session: &Session, feature: EntityId) -> (f64, [[f64; 3]; 2]) {
    let source = model::read_wall(session.document(), feature).unwrap();
    let wall = model::resolve_wall(session.document(), session.parameters(), &source).unwrap();
    assert!(wall.junctions.is_empty() && wall.cutouts.is_empty());
    assert_eq!(wall.justification, model::Justification::Center);
    let (t, h) = (wall.thickness_mm, wall.height_mm);
    let n = wall.points.len();
    let segments = if wall.closed { n } else { n - 1 };
    let mut min = [f64::INFINITY; 2];
    let mut max = [f64::NEG_INFINITY; 2];
    let mut centerline = 0.;
    for i in 0..segments {
        let [a, b] = [wall.points[i], wall.points[(i + 1) % n]];
        assert!(a[0] == b[0] || a[1] == b[1], "axis-aligned segment");
        centerline += (b[0] - a[0]).abs() + (b[1] - a[1]).abs();
        for p in [a, b] {
            for axis in 0..2 {
                min[axis] = min[axis].min(p[axis] - t / 2.);
                max[axis] = max[axis].max(p[axis] + t / 2.);
            }
        }
    }
    if !wall.closed {
        assert_eq!(n, 2, "one straight open segment");
        // Square ends: the open path is not widened along its own axis.
        let [a, b] = [wall.points[0], wall.points[n - 1]];
        for axis in 0..2 {
            if a[1 - axis] == b[1 - axis] {
                min[axis] = a[axis].min(b[axis]);
                max[axis] = a[axis].max(b[axis]);
            }
        }
    }
    // Mitered closed corners and square open ends both enclose exactly the
    // centerline length times the thickness.
    let openings: f64 = wall
        .openings
        .iter()
        .map(|opening| opening.width_mm * opening.height_mm)
        .sum();
    (
        (centerline * h - openings) * t,
        [
            [min[0], min[1], wall.base_z_mm],
            [max[0], max[1], wall.base_z_mm + h],
        ],
    )
}

fn expected_slab(session: &Session, feature: EntityId) -> (f64, [[f64; 3]; 2]) {
    let source = model::read_slab(session.document(), feature).unwrap();
    let slab = model::resolve_slab(session.document(), session.parameters(), &source).unwrap();
    let n = slab.outline.len();
    let area = (0..n)
        .map(|i| {
            let [a, b] = [slab.outline[i], slab.outline[(i + 1) % n]];
            a[0] * b[1] - b[0] * a[1]
        })
        .sum::<f64>()
        .abs()
        / 2.;
    let holes: f64 = slab
        .holes
        .iter()
        .map(|hole| {
            let c = hole.corners;
            (0..4)
                .map(|i| c[i][0] * c[(i + 1) % 4][1] - c[(i + 1) % 4][0] * c[i][1])
                .sum::<f64>()
                .abs()
                / 2.
        })
        .sum();
    let fold = |pick: fn(f64, f64) -> f64, axis: usize| {
        slab.outline.iter().map(|p| p[axis]).reduce(pick).unwrap()
    };
    (
        (area - holes) * (slab.top_mm - slab.bottom_mm),
        [
            [fold(f64::min, 0), fold(f64::min, 1), slab.bottom_mm],
            [fold(f64::max, 0), fold(f64::max, 1), slab.top_mm],
        ],
    )
}

#[test]
fn shipped_apartment_walls_and_slab_rebuild_with_openings() {
    let session = apartment();
    let walls = architecture(&session, model::WALL_TYPE);
    let slabs = architecture(&session, model::SLAB_TYPE);
    assert_eq!(walls.len(), 4);
    assert_eq!(slabs.len(), 1);
    let sources: Vec<_> = walls.iter().chain(&slabs).copied().collect();
    let rust = compile(&session, &sources, GeometryBackend::Rust);
    for feature in sources {
        assert!(rust.solid(feature).unwrap().uses_rust_geometry());
        let shape = &rust.snapshot().shapes[&feature];
        let (volume, [min, max]) = if slabs.contains(&feature) {
            expected_slab(&session, feature)
        } else {
            expected_wall(&session, feature)
        };
        assert_eq!(shape.solid_count, 1);
        near(shape.volume_mm3, volume, 1e-3);
        for axis in 0..3 {
            near(shape.bounds.min[axis], min[axis], 1e-6);
            near(shape.bounds.max[axis], max[axis], 1e-6);
        }
        assert!(!shape.mesh.indices.is_empty());
    }
    let mut openings = 0;
    for feature in walls {
        let source = model::read_wall(session.document(), feature).unwrap();
        for opening in source.openings {
            for side in ["head", "jamb.start", "jamb.end"] {
                let role = format!("opening.{}.{side}", opening.id);
                assert!(face_has_role(&rust, feature, &role), "missing {role}");
            }
            let opening_type =
                model::read_opening_type(session.document(), opening.opening_type).unwrap();
            if opening_type.kind == OpeningKind::Window {
                let role = format!("opening.{}.sill", opening.id);
                assert!(face_has_role(&rust, feature, &role), "missing {role}");
            }
            openings += 1;
        }
    }
    assert_eq!(openings, 10);
}

#[test]
fn slab_shaft_edit_removes_material_and_keeps_the_four_sides_named() {
    let mut session = apartment();
    let [feature] = architecture(&session, model::SLAB_TYPE).try_into().unwrap();
    let baseline = compile(&session, &[feature], GeometryBackend::Rust);
    let mut slab = model::read_slab(session.document(), feature).unwrap();
    let shaft = EntityId::new();
    slab.holes.push(SlabHole {
        id: shaft,
        center: [length(500.), length(500.)],
        width: length(500.),
        depth: length(800.),
        rotation_deg: 30.,
    });
    let edit = model::prepare_slab_edit(
        &session,
        session.revision(),
        Actor::local(Channel::Ui),
        feature,
        &slab,
    )
    .unwrap();
    session.commit(edit.prepared).unwrap();
    let updated = compile(&session, &[feature], GeometryBackend::Rust);
    let old = &baseline.snapshot().shapes[&feature];
    let new = &updated.snapshot().shapes[&feature];
    near(new.volume_mm3, expected_slab(&session, feature).0, 1e-3);
    assert_eq!(old.bounds, new.bounds);
    let height = old.bounds.max[2] - old.bounds.min[2];
    near(old.volume_mm3 - new.volume_mm3, 500. * 800. * height, 1e-3);
    for side in 0..4 {
        assert!(face_has_role(
            &updated,
            feature,
            &format!("cutout.{shaft}.side.{side}")
        ));
    }
}

#[test]
fn edited_wall_can_cut_a_chase_and_shorten_its_end() {
    let mut session = apartment();
    let (feature, wall) = uncut_partition(&session);
    let baseline = compile(&session, &[feature], GeometryBackend::Rust);
    let chase = EntityId::new();
    let end = EntityId::new();
    for cutout in [
        VerticalCutout {
            id: chase,
            segment: 0,
            offset: length(600.),
            width: length(200.),
            depth: Some(length(40.)),
            face: WallFace::Right,
        },
        VerticalCutout {
            id: end,
            segment: 0,
            offset: length(3_090.),
            width: length(300.),
            depth: None,
            face: WallFace::Left,
        },
    ] {
        let edit = model::prepare_wall_opening(
            &session,
            session.revision(),
            Actor::local(Channel::Ui),
            feature,
            &[],
            &WallOpeningChange::AddCutout(cutout),
        )
        .unwrap();
        session.commit(edit.prepared).unwrap();
    }
    let updated = compile(&session, &[feature], GeometryBackend::Rust);
    let old = &baseline.snapshot().shapes[&feature];
    let new = &updated.snapshot().shapes[&feature];
    assert_eq!(new.solid_count, 1);
    near(new.bounds.max[1], old.bounds.max[1] - 300., 1e-6);
    let wall = model::resolve_wall(session.document(), session.parameters(), &wall).unwrap();
    near(
        old.volume_mm3 - new.volume_mm3,
        (200. * 40. + 300. * wall.thickness_mm) * wall.height_mm,
        1e-3,
    );
    assert!(face_has_role(
        &updated,
        feature,
        &format!("cutout.{chase}.back")
    ));
    assert!(face_has_role(
        &updated,
        feature,
        &format!("cutout.{end}.jamb.start")
    ));
}

fn uncut_partition(session: &Session) -> (EntityId, model::WallDefinition) {
    architecture(session, model::WALL_TYPE)
        .into_iter()
        .find_map(|feature| {
            let wall = model::read_wall(session.document(), feature).unwrap();
            (wall.path.len() == 2 && wall.openings.is_empty()).then_some((feature, wall))
        })
        .unwrap()
}

#[test]
fn edited_diagonal_wall_matches_the_analytic_volume_and_boundary() {
    let mut session = apartment();
    let (feature, mut wall) = uncut_partition(&session);
    wall.path[1][0] = length(4700.);
    let edit = model::prepare_wall_edit(
        &session,
        session.revision(),
        Actor::local(Channel::Ui),
        feature,
        &wall,
    )
    .unwrap();
    session.commit(edit.prepared).unwrap();
    let rust = compile(&session, &[feature], GeometryBackend::Rust);
    let actual = &rust.snapshot().shapes[&feature];
    // A square-ended straight wall is a rectangle swept vertically: its
    // corners are the path ends offset by half the thickness along the normal.
    let resolved = model::resolve_wall(session.document(), session.parameters(), &wall).unwrap();
    assert!(resolved.junctions.is_empty());
    assert_eq!(resolved.justification, model::Justification::Center);
    let [a, b] = [resolved.points[0], resolved.points[1]];
    let length = (b[0] - a[0]).portable_hypot(b[1] - a[1]);
    let normal = [(a[1] - b[1]) / length, (b[0] - a[0]) / length];
    let half = resolved.thickness_mm / 2.;
    let corners = [a, b].into_iter().flat_map(|p| {
        [-half, half].map(|offset| [p[0] + offset * normal[0], p[1] + offset * normal[1]])
    });
    let mut min = [f64::INFINITY; 3];
    let mut max = [f64::NEG_INFINITY; 3];
    for corner in corners {
        for axis in 0..2 {
            min[axis] = min[axis].min(corner[axis]);
            max[axis] = max[axis].max(corner[axis]);
        }
    }
    [min[2], max[2]] = [resolved.base_z_mm, resolved.base_z_mm + resolved.height_mm];
    assert_eq!(actual.solid_count, 1);
    assert!(rust.solid(feature).unwrap().uses_rust_geometry());
    near(
        actual.volume_mm3,
        length * resolved.thickness_mm * resolved.height_mm,
        1e-3,
    );
    for axis in 0..3 {
        near(actual.bounds.min[axis], min[axis], 1e-6);
        near(actual.bounds.max[axis], max[axis], 1e-6);
    }
    assert!(face_has_role(&rust, feature, "wall.seg.0.left"));
}
