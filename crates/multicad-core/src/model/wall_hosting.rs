//! The host relationship between a wall and the doors, windows and vertical
//! cutouts it carries. An opening is hosted by exactly the wall whose source
//! lists it, so this index is read from canonical wall definitions, never from
//! geometry. Moving a door or window to another wall is one transaction on
//! both walls, so the UI tree, MCP and file reopen see the same host.
use super::{
    CreatedBody, ModelPrepared,
    common::{error, revision},
    wall::{Opening, OpeningKind, read_opening_type, read_wall, resolve_wall, walls},
};
use crate::{
    commands::{Actor, CommandError, Request, Session},
    document::{Document, DocumentError, EntityData, EntityId, Mutation},
};

/// What a wall-hosted opening is.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HostedKind {
    Door,
    Window,
    VerticalCutout,
}
impl HostedKind {
    pub fn label(self) -> &'static str {
        match self {
            Self::Door => "Door",
            Self::Window => "Window",
            Self::VerticalCutout => "Vertical cutout",
        }
    }
}

/// One opening and the wall that hosts it.
#[derive(Clone, Debug, PartialEq)]
pub struct HostedOpening {
    pub id: EntityId,
    /// The hosting wall feature.
    pub wall: EntityId,
    /// The hosting wall's output Part.
    pub part: EntityId,
    pub kind: HostedKind,
    /// The opening type's name, or the cutout kind.
    pub name: String,
    /// Zero-based wall segment.
    pub segment: usize,
}

fn kind_of(document: &Document, opening: &Opening) -> HostedKind {
    match read_opening_type(document, opening.opening_type).map(|t| t.kind) {
        Ok(OpeningKind::Window) => HostedKind::Window,
        // A door, or a type that no longer reads: the wall itself reports
        // the broken type when it is evaluated.
        _ => HostedKind::Door,
    }
}

/// The openings hosted by one wall feature, in source order: doors and
/// windows first, then vertical cutouts.
pub fn wall_openings(
    document: &Document,
    wall: EntityId,
) -> Result<Vec<HostedOpening>, DocumentError> {
    let source = read_wall(document, wall)?;
    let part = source.output.body;
    let mut out = Vec::with_capacity(source.openings.len() + source.cutouts.len());
    for opening in &source.openings {
        let name = document
            .get(opening.opening_type)
            .map(|t| t.name.clone())
            .unwrap_or_else(|| "Missing opening type".into());
        out.push(HostedOpening {
            id: opening.id,
            wall,
            part,
            kind: kind_of(document, opening),
            name,
            segment: opening.segment,
        });
    }
    for cutout in &source.cutouts {
        out.push(HostedOpening {
            id: cutout.id,
            wall,
            part,
            kind: HostedKind::VerticalCutout,
            name: HostedKind::VerticalCutout.label().into(),
            segment: cutout.segment,
        });
    }
    Ok(out)
}

/// Every wall-hosted opening in the document. Walls that do not read are
/// skipped; their own evaluation reports why.
pub fn hosted_openings(document: &Document) -> Vec<HostedOpening> {
    walls(document)
        .filter_map(|wall| wall_openings(document, wall.id).ok())
        .flatten()
        .collect()
}

/// The wall hosting the opening or cutout with this identity, if any.
pub fn opening_host(document: &Document, id: EntityId) -> Option<HostedOpening> {
    if id.is_nil() {
        return None;
    }
    hosted_openings(document).into_iter().find(|o| o.id == id)
}

/// Stage moving a door or window from wall `from` to wall `to` as one
/// transaction on both walls. `opening` keeps the identity, type and size
/// overrides and carries its placement (segment and offset) on the new wall.
/// Both walls are planned and their openings checked before anything is
/// staged; a native evaluator still builds both solids.
pub fn prepare_rehost_opening(
    session: &Session,
    expected_revision: u64,
    actor: Actor,
    from: EntityId,
    to: EntityId,
    opening: &Opening,
) -> Result<ModelPrepared, CommandError> {
    revision(session, expected_revision)?;
    let document = session.document();
    if from == to {
        return Err(error(
            "Choose a different wall; edit the opening to move it along this wall",
        ));
    }
    let mut source = read_wall(document, from)?;
    let mut target = read_wall(document, to)?;
    let before = source.openings.len();
    source.openings.retain(|o| o.id != opening.id);
    if source.openings.len() == before {
        return Err(if source.cutouts.iter().any(|c| c.id == opening.id) {
            error("A vertical cutout cannot move to another wall; remove it and place a new one")
        } else {
            error("This wall hosts no door or window with that identity")
        });
    }
    if opening.id.is_nil() || target.has_opening(opening.id) {
        return Err(error("The target wall already uses that opening identity"));
    }
    if source.level != target.level {
        return Err(error(
            "The target wall is on another level; move openings between walls on one level",
        ));
    }
    target.openings.push(opening.clone());
    let mut edits = vec![];
    for (feature, wall) in [(from, &source), (to, &target)] {
        wall.validate()?;
        let mut entity = document
            .get(feature)
            .ok_or_else(|| error("Wall is missing"))?
            .clone();
        let EntityData::Feature { suppressed, .. } = entity.data else {
            return Err(error("Wall is not a feature"));
        };
        entity.data = EntityData::Feature {
            definition: wall.definition()?,
            inputs: wall.inputs(),
            suppressed,
        };
        edits.push(Mutation::Replace(entity));
    }
    let candidate = document.edited(&edits)?;
    let parameters = crate::parameters::evaluate(&candidate)?;
    for wall in [&source, &target] {
        let resolved = resolve_wall(&candidate, &parameters, wall)?;
        resolved.check_openings(&resolved.plan()?)?;
    }
    let noun = match kind_of(document, opening) {
        HostedKind::Window => "window",
        _ => "door",
    };
    let target_name = document
        .get(target.output.body)
        .map(|part| part.name.as_str())
        .unwrap_or("wall");
    let prepared = session.prepare(Request {
        expected_revision,
        actor,
        label: format!("Move {noun} to {target_name}"),
        edits,
    })?;
    Ok(ModelPrepared {
        prepared,
        outputs: vec![
            CreatedBody {
                feature: from,
                body: source.output.body,
            },
            CreatedBody {
                feature: to,
                body: target.output.body,
            },
        ],
    })
}

/// The default placement of a moved opening on its new wall: the first
/// segment, centered along it. The caller may adjust it before staging.
pub fn centered_on_first_segment(
    document: &Document,
    to: EntityId,
    width_mm: f64,
) -> Result<f64, CommandError> {
    let wall = read_wall(document, to)?;
    let parameters = crate::parameters::evaluate(document)?;
    let resolved = resolve_wall(document, &parameters, &wall)?;
    let plan = resolved.plan()?;
    let segment = plan
        .segments
        .first()
        .ok_or_else(|| error("The target wall has no segment"))?;
    Ok(((segment.length_mm - width_mm) / 2.).max(0.))
}

/// The plan footprint of a hosted opening or cutout in world X/Y
/// millimetres: the rectangle it cuts through the wall band, counter-
/// clockwise from the start jamb on the right face. `None` when no current
/// wall hosts it or the wall does not resolve.
pub fn opening_footprint(
    document: &Document,
    parameters: &crate::parameters::EvaluatedParameters,
    id: EntityId,
) -> Option<[[f64; 2]; 4]> {
    let host = opening_host(document, id)?;
    let wall = resolve_wall(document, parameters, &read_wall(document, host.wall).ok()?).ok()?;
    let plan = wall.plan().ok()?;
    let ([a, b], [low, high]) = if let Some(o) = wall.openings.iter().find(|o| o.id == id) {
        (
            [o.offset_mm, o.offset_mm + o.width_mm],
            [plan.right_mm, plan.left_mm],
        )
    } else {
        let c = wall.cutouts.iter().find(|c| c.id == id)?;
        wall.cutout_span(&plan, c)
    };
    let seg = plan.segments.get(host.segment)?;
    let at = |along: f64, across: f64| {
        [
            seg.start[0] + seg.direction[0] * along + seg.normal[0] * across,
            seg.start[1] + seg.direction[1] * along + seg.normal[1] * across,
        ]
    };
    Some([at(a, low), at(b, low), at(b, high), at(a, high)])
}

#[cfg(test)]
mod tests {
    use super::super::{
        BodyOutput, Hinge, Justification, NewType, NewTypeKind, OpeningTypeDefinition, Swing,
        VerticalCutout, WallChange, WallDefinition, WallFace, WallOpeningChange, WallRequest,
        WallTypeDefinition, length, prepare_wall_change, prepare_wall_opening,
    };
    use super::*;
    use crate::commands::Channel;

    fn actor() -> Actor {
        Actor::local(Channel::Ui)
    }
    struct Fixture {
        session: Session,
        walls: [EntityId; 2],
        door: EntityId,
        window: EntityId,
    }
    /// Two parallel 6 m walls, 3 m apart, and a door and window type.
    fn fixture() -> Fixture {
        let mut session = Session::new(Document::new("Hosting").unwrap()).unwrap();
        let ty = NewType::new(
            "Block 200",
            NewTypeKind::Wall(WallTypeDefinition {
                thickness: length(200.),
                height: length(2500.),
                justification: Justification::Center,
                material: None,
            }),
        );
        let opening = |name: &str, kind, w, h, s| {
            NewType::new(
                name,
                NewTypeKind::Opening(OpeningTypeDefinition {
                    kind,
                    width: length(w),
                    height: length(h),
                    sill: length(s),
                }),
            )
        };
        let door = opening(
            "Door 900",
            OpeningKind::Door {
                hinge: Hinge::Start,
                swing: Swing::Left,
            },
            900.,
            2100.,
            0.,
        );
        let window = opening("Window 1200", OpeningKind::Window, 1200., 1200., 900.);
        let mut types = vec![ty.clone(), door.clone(), window.clone()];
        let mut ids = vec![];
        for (name, y) in [("North wall", 3000.), ("South wall", 0.)] {
            let request = WallRequest::new(
                name,
                WallDefinition {
                    layer_placement: None,
                    output: BodyOutput {
                        body: EntityId::new(),
                    },
                    wall_type: ty.id,
                    path: vec![[length(0.), length(y)], [length(6000.), length(y)]],
                    closed: false,
                    level: None,
                    base_z: length(0.),
                    thickness: None,
                    height: None,
                    justification: None,
                    openings: vec![],
                    cutouts: vec![],
                },
            );
            let prepared = prepare_wall_change(
                &session,
                session.revision(),
                actor(),
                &types,
                &WallChange::Create(&request),
            )
            .unwrap();
            session.commit(prepared.prepared).unwrap();
            types.clear();
            ids.push(request.feature);
        }
        Fixture {
            session,
            walls: [ids[0], ids[1]],
            door: door.id,
            window: window.id,
        }
    }
    fn opening(ty: EntityId, offset: f64) -> Opening {
        Opening {
            id: EntityId::new(),
            opening_type: ty,
            segment: 0,
            offset: length(offset),
            width: None,
            height: None,
            sill: None,
        }
    }
    fn add(f: &mut Fixture, wall: EntityId, change: WallOpeningChange) {
        let prepared = prepare_wall_opening(
            &f.session,
            f.session.revision(),
            actor(),
            wall,
            &[],
            &change,
        )
        .unwrap();
        f.session.commit(prepared.prepared).unwrap();
    }

    #[test]
    fn hosted_openings_follow_the_wall_source_through_add_remove_and_rehost() {
        let mut f = fixture();
        let [north, south] = f.walls;
        assert!(hosted_openings(f.session.document()).is_empty());
        let door = opening(f.door, 1000.);
        let window = opening(f.window, 3000.);
        let cutout = VerticalCutout {
            id: EntityId::new(),
            segment: 0,
            offset: length(5000.),
            width: length(300.),
            depth: Some(length(50.)),
            face: WallFace::Left,
        };
        add(&mut f, north, WallOpeningChange::AddOpening(door.clone()));
        add(&mut f, north, WallOpeningChange::AddOpening(window.clone()));
        add(&mut f, south, WallOpeningChange::AddCutout(cutout.clone()));
        let document = f.session.document();
        let north_part = read_wall(document, north).unwrap().output.body;
        let listed = wall_openings(document, north).unwrap();
        assert_eq!(
            listed
                .iter()
                .map(|o| (o.id, o.kind, o.name.as_str(), o.part))
                .collect::<Vec<_>>(),
            vec![
                (door.id, HostedKind::Door, "Door 900", north_part),
                (window.id, HostedKind::Window, "Window 1200", north_part),
            ]
        );
        let host = opening_host(document, cutout.id).unwrap();
        assert_eq!(
            (host.wall, host.kind, host.name.as_str()),
            (south, HostedKind::VerticalCutout, "Vertical cutout")
        );
        assert_eq!(opening_host(document, EntityId::new()), None);
        // The north wall runs along +X at y = 3000 and is 200 mm thick; its
        // left face is +Y. The door spans x 1000–1900 across the band.
        let parameters = crate::parameters::evaluate(document).unwrap();
        assert_eq!(
            opening_footprint(document, &parameters, door.id),
            Some([
                [1000., 2900.],
                [1900., 2900.],
                [1900., 3100.],
                [1000., 3100.]
            ])
        );
        // A 50 mm chase on the left face of the south wall (y = 0).
        assert_eq!(
            opening_footprint(document, &parameters, cutout.id),
            Some([[5000., 50.], [5300., 50.], [5300., 100.], [5000., 100.]])
        );
        assert_eq!(
            opening_footprint(document, &parameters, EntityId::new()),
            None
        );

        // Rehost the door: one transaction on both walls, identity kept.
        let mut moved = door.clone();
        moved.offset = length(centered_on_first_segment(document, south, 900.).unwrap());
        assert_eq!(moved.offset, length(2550.));
        let revision = f.session.revision();
        let prepared =
            prepare_rehost_opening(&f.session, revision, actor(), north, south, &moved).unwrap();
        assert_eq!(
            prepared
                .outputs
                .iter()
                .map(|o| o.feature)
                .collect::<Vec<_>>(),
            vec![north, south]
        );
        f.session.commit(prepared.prepared).unwrap();
        assert_eq!(f.session.revision(), revision + 1);
        assert_eq!(
            f.session.activity().last().map(|a| a.label.as_str()),
            Some("Move door to South wall")
        );
        let document = f.session.document();
        assert_eq!(opening_host(document, door.id).unwrap().wall, south);
        assert_eq!(
            wall_openings(document, north)
                .unwrap()
                .iter()
                .map(|o| o.id)
                .collect::<Vec<_>>(),
            vec![window.id]
        );
        assert_eq!(
            read_wall(document, south).unwrap().openings,
            vec![moved.clone()]
        );

        // Undo restores the original host in one step.
        f.session.undo(f.session.revision(), actor()).unwrap();
        assert_eq!(
            opening_host(f.session.document(), door.id).unwrap().wall,
            north
        );

        // Removing the window leaves only the door listed on the north wall.
        add(&mut f, north, WallOpeningChange::Remove(window.id));
        assert_eq!(opening_host(f.session.document(), window.id), None);

        // The host relationship is file source: save and reopen keep it.
        let hosted = hosted_openings(f.session.document());
        assert_eq!(hosted.len(), 2);
        let file = crate::project_file::ProjectFile::new(f.session.document().clone()).unwrap();
        let mut bytes = vec![];
        crate::project_file::write(&file, &mut bytes, Default::default()).unwrap();
        let reopened = crate::project_file::read(
            &mut bytes.as_slice(),
            Default::default(),
            &crate::project_file::supported_capabilities(),
        )
        .unwrap()
        .project
        .document;
        assert_eq!(hosted_openings(&reopened), hosted);
    }

    #[test]
    fn refused_rehosts_stage_nothing() {
        let mut f = fixture();
        let [north, south] = f.walls;
        let door = opening(f.door, 1000.);
        add(&mut f, north, WallOpeningChange::AddOpening(door.clone()));
        let cutout = VerticalCutout {
            id: EntityId::new(),
            segment: 0,
            offset: length(4000.),
            width: length(300.),
            depth: Some(length(50.)),
            face: WallFace::Left,
        };
        add(&mut f, north, WallOpeningChange::AddCutout(cutout.clone()));
        let before = f.session.document().clone();
        let revision = f.session.revision();
        let refused = |from, to, o: &Opening| {
            prepare_rehost_opening(&f.session, revision, actor(), from, to, o)
                .unwrap_err()
                .to_string()
        };
        assert!(refused(north, north, &door).contains("different wall"));
        assert!(refused(south, north, &door).contains("no door or window"));
        let as_opening = Opening {
            id: cutout.id,
            ..door.clone()
        };
        assert!(refused(north, south, &as_opening).contains("vertical cutout"));
        let mut misfit = door.clone();
        misfit.offset = length(5500.);
        assert!(!refused(north, south, &misfit).is_empty());
        assert_eq!(f.session.document(), &before);
    }
}
