//! Evaluated attachment and revision guards precede typed drawing and solving.
use crate::{EditError, InputError, SketchEdit};
use multicad_core::{
    commands::{Actor, CommandError, ErrorCode, Prepared, Session},
    constraint_sketch::source::Source,
    document::{EntityId, Mutation},
    references::ReferenceReport,
    sketch_draw::{self, DrawRequest},
};

pub struct DrawEdit {
    pub edit: SketchEdit,
    pub created_curves: Vec<EntityId>,
}
/// Numerical conflicts retain their diagnostics but produce no commit-capable
/// preview. Legacy sketches require an explicit source upgrade before drawing.
pub fn prepare_draw(
    session: &Session,
    expected_revision: u64,
    actor: Actor,
    sketch: EntityId,
    report: Option<&ReferenceReport>,
    request: &DrawRequest,
) -> Result<DrawEdit, EditError> {
    prepare_draw_after(
        session,
        expected_revision,
        actor,
        sketch,
        report,
        request,
        Vec::new(),
        session.parameters(),
        |_, _| Ok(()),
    )
}
/// Draws onto an uncommitted candidate, such as a captured external
/// reference whose output point the drawing reuses, and returns ONE edit that
/// contains both: one commit, one Undo step (CAD-383). `report` must describe
/// the candidate document. A candidate from another session or revision is
/// refused.
pub fn prepare_draw_on(
    session: &Session,
    expected_revision: u64,
    actor: Actor,
    sketch: EntityId,
    base: &Prepared,
    report: Option<&ReferenceReport>,
    request: &DrawRequest,
) -> Result<DrawEdit, EditError> {
    prepare_draw_on_with(
        session,
        expected_revision,
        actor,
        sketch,
        base,
        report,
        request,
        |_, _| Ok(()),
    )
}
/// [`prepare_draw_on`] whose `finish` may add constraints, such as entered
/// sizes, to the drawn source before it is solved, so they join the same
/// transaction (CAD-383). `finish` sees the created curves in drawing order;
/// an error from it prepares nothing.
#[allow(clippy::too_many_arguments)]
pub fn prepare_draw_on_with(
    session: &Session,
    expected_revision: u64,
    actor: Actor,
    sketch: EntityId,
    base: &Prepared,
    report: Option<&ReferenceReport>,
    request: &DrawRequest,
    finish: impl FnOnce(&mut Source, &[EntityId]) -> Result<(), String>,
) -> Result<DrawEdit, EditError> {
    session.validate_preview(base).map_err(EditError::Command)?;
    if base.base_revision() != expected_revision {
        return Err(stale(sketch));
    }
    let prefix = session.document().mutations_to(base.document());
    prepare_draw_after(
        session,
        expected_revision,
        actor,
        sketch,
        report,
        request,
        prefix,
        base.parameters(),
        finish,
    )
}
fn stale(sketch: EntityId) -> EditError {
    EditError::Command(CommandError {
        code: ErrorCode::StaleRevision,
        message: "Drawing is based on a stale document revision".into(),
        involved: vec![sketch],
        parameter: None,
    })
}
#[allow(clippy::too_many_arguments)]
fn prepare_draw_after(
    session: &Session,
    expected_revision: u64,
    actor: Actor,
    sketch: EntityId,
    report: Option<&ReferenceReport>,
    request: &DrawRequest,
    prefix: Vec<Mutation>,
    parameters: &multicad_core::parameters::EvaluatedParameters,
    finish: impl FnOnce(&mut Source, &[EntityId]) -> Result<(), String>,
) -> Result<DrawEdit, EditError> {
    if session.revision() != expected_revision {
        return Err(stale(sketch));
    }
    let input = |e: multicad_core::document::DocumentError| {
        EditError::Input(InputError {
            source: Some(sketch),
            message: e.to_string(),
        })
    };
    let finished = |message: String| {
        EditError::Input(InputError {
            source: Some(sketch),
            message,
        })
    };
    let base = if prefix.is_empty() {
        session.document().clone()
    } else {
        session.document().edited(&prefix).map_err(input)?
    };
    let (source, _, _) =
        Source::read_for_evaluation(&base, sketch, expected_revision, report).map_err(input)?;
    let semantic = source.minimum_schema_version() >= 2
        || matches!(&base.get(sketch).expect("read sketch").data,multicad_core::document::EntityData::Sketch {definition,..} if (2..=4).contains(&definition.schema_version));
    if semantic {
        let evaluated = crate::evaluate_document_sketch_with_parameters(
            &base,
            expected_revision,
            sketch,
            report,
            parameters,
        )
        .map_err(EditError::Input)?;
        let solved = evaluated.solved().ok_or_else(|| {
            EditError::Input(InputError {
                source: Some(sketch),
                message: "Resolve the current sketch before drawing".into(),
            })
        })?;
        let drawing_source = Source {
            u: source.u,
            sketch: solved.clone(),
        };
        let draft = sketch_draw::draft(&drawing_source, request).map_err(input)?;
        let mut canonical = source.clone();
        canonical.sketch.points.extend(
            draft
                .source
                .sketch
                .points
                .iter()
                .filter(|p| !source.sketch.points.iter().any(|old| old.id == p.id))
                .cloned(),
        );
        canonical.sketch.curves.extend(
            draft
                .source
                .sketch
                .curves
                .iter()
                .filter(|c| !source.sketch.curves.iter().any(|old| old.id == c.id))
                .cloned(),
        );
        canonical.sketch.constraints.extend(
            draft
                .source
                .sketch
                .constraints
                .iter()
                .filter(|c| !source.sketch.constraints.iter().any(|old| old.id == c.id))
                .cloned(),
        );
        canonical.sketch.shapes.extend(
            draft
                .source
                .sketch
                .shapes
                .iter()
                .filter(|s| !source.sketch.shapes.iter().any(|old| old.id == s.id))
                .cloned(),
        );
        finish(&mut canonical, &draft.created_curves).map_err(finished)?;
        let edit = crate::transaction::prepare_edit_after(
            session,
            expected_revision,
            actor,
            sketch,
            canonical,
            report,
            Some(&draft.source.sketch),
            None,
            prefix,
            Vec::new(),
        )?;
        return Ok(DrawEdit {
            edit,
            created_curves: draft.created_curves,
        });
    }
    let mut draft = sketch_draw::draft(&source, request).map_err(input)?;
    finish(&mut draft.source, &draft.created_curves).map_err(finished)?;
    let edit = crate::transaction::prepare_edit_after(
        session,
        expected_revision,
        actor,
        sketch,
        draft.source,
        None,
        None,
        None,
        prefix,
        Vec::new(),
    )?;
    Ok(DrawEdit {
        edit,
        created_curves: draft.created_curves,
    })
}
