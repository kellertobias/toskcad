//! Size and position edits of a recorded drawn shape (CAD-398). Sizes are
//! canonical dimension constraints on the shape's own curves; position is a
//! Fix of the shape's anchor point. Nothing here solves: callers publish the
//! returned source through the solver transaction, which refuses conflicts.
use super::{
    Constraint, ConstraintKind, CurveKind, Dimension, Shape, ShapeKind, Sketch,
    binding::{DimensionInput, draft_dimension},
    source::Source,
};
use crate::document::{Document, DocumentError, EntityId};
use std::collections::BTreeSet;

/// One editable size of a shape.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ShapeSize {
    /// Rectangle first edge.
    Width,
    /// Rectangle second edge.
    Height,
    /// Line length.
    Length,
    Diameter,
    Radius,
    /// Chamfer setback from the original corner (CAD-414).
    Distance,
}
impl ShapeSize {
    pub fn label(self) -> &'static str {
        match self {
            Self::Width => "Width",
            Self::Height => "Height",
            Self::Length => "Length",
            Self::Diameter => "Diameter",
            Self::Radius => "Radius",
            Self::Distance => "Distance",
        }
    }
}
fn invalid(message: impl Into<String>) -> DocumentError {
    DocumentError::Invalid(message.into())
}
/// Sizes a shape offers, in inspector order.
pub fn sizes(shape: &Shape) -> &'static [ShapeSize] {
    match shape.kind {
        ShapeKind::Rectangle if shape.curves.len() == 4 => &[ShapeSize::Width, ShapeSize::Height],
        ShapeKind::Line => &[ShapeSize::Length],
        ShapeKind::Circle => &[ShapeSize::Diameter],
        ShapeKind::Arc | ShapeKind::Fillet => &[ShapeSize::Radius],
        ShapeKind::Chamfer => &[ShapeSize::Distance],
        _ => &[],
    }
}
/// The curve a size measures.
pub fn size_curve(shape: &Shape, size: ShapeSize) -> Option<EntityId> {
    if !sizes(shape).contains(&size) {
        return None;
    }
    match size {
        ShapeSize::Height => shape.curves.get(1).copied(),
        _ => shape.curves.first().copied(),
    }
}
fn matches_size(kind: &ConstraintKind, curve: EntityId, size: ShapeSize) -> bool {
    match (kind, size) {
        (
            ConstraintKind::Length { line, .. },
            ShapeSize::Width | ShapeSize::Height | ShapeSize::Length,
        ) => *line == curve,
        (ConstraintKind::Diameter { curve: c, .. }, ShapeSize::Diameter)
        | (ConstraintKind::Radius { curve: c, .. }, ShapeSize::Radius) => *c == curve,
        _ => false,
    }
}
/// The driving dimension constraints that size a fillet or chamfer
/// (CAD-414): a fillet's Radius on its arc; a chamfer's setback on each side,
/// a Length on the construction line from a chamfer end to the virtual sharp
/// (line sides) or a point Distance from a chamfer end (arc sides).
pub fn corner_constraints(sketch: &Sketch, shape: &Shape) -> Vec<EntityId> {
    let Some(curve) = shape.curves.first().copied() else {
        return vec![];
    };
    match shape.kind {
        ShapeKind::Fillet => sketch
            .constraints
            .iter()
            .filter(|c| {
                !c.reference
                    && matches!(c.kind, ConstraintKind::Radius { curve: arc, .. } if arc == curve)
            })
            .map(|c| c.id)
            .take(1)
            .collect(),
        ShapeKind::Chamfer => {
            let Some(CurveKind::Line { start, end }) = sketch
                .curves
                .iter()
                .find(|c| c.id == curve)
                .map(|c| c.kind.clone())
            else {
                return vec![];
            };
            let ends = [start, end];
            let helper_from_end = |line: EntityId| {
                sketch.curves.iter().any(|c| {
                    c.id == line
                        && c.construction
                        && matches!(c.kind, CurveKind::Line { start, .. } if ends.contains(&start))
                })
            };
            sketch
                .constraints
                .iter()
                .filter(|c| !c.reference)
                .filter(|c| match c.kind {
                    ConstraintKind::Length { line, .. } => helper_from_end(line),
                    ConstraintKind::Distance { a, .. } => ends.contains(&a),
                    _ => false,
                })
                .map(|c| c.id)
                .collect()
        }
        _ => vec![],
    }
}
/// The canonical radius or setback of a fillet or chamfer (CAD-414).
pub fn corner_size(sketch: &Sketch, shape: &Shape) -> Option<Dimension> {
    let id = corner_constraints(sketch, shape).into_iter().next()?;
    sketch
        .constraints
        .iter()
        .find(|c| c.id == id)
        .and_then(|c| c.kind.dimension().cloned())
}
/// List label of a fillet or chamfer with its size, e.g. "Fillet R10" or
/// "Chamfer 2.5" (millimetres, trailing zeros trimmed); an expression size is
/// shown as written. Other shapes keep their name.
pub fn corner_label(sketch: &Sketch, shape: &Shape) -> String {
    let prefix = match shape.kind {
        ShapeKind::Fillet => "Fillet R",
        ShapeKind::Chamfer => "Chamfer ",
        _ => return shape.name.clone(),
    };
    match corner_size(sketch, shape) {
        Some(Dimension::Literal { value, unit }) => {
            let text = format!("{value:.3}");
            let text = text.trim_end_matches('0').trim_end_matches('.');
            if unit == "mm" {
                format!("{prefix}{text}")
            } else {
                format!("{prefix}{text} {unit}")
            }
        }
        Some(Dimension::Expression { expression }) => format!("{prefix}{}", expression.text),
        None => shape.kind.label().into(),
    }
}
/// The driving (non-reference) dimension constraint for `size`, if any.
pub fn size_constraint<'a>(
    sketch: &'a Sketch,
    shape: &Shape,
    size: ShapeSize,
) -> Option<&'a Constraint> {
    if shape.kind.corner() {
        let id = corner_constraints(sketch, shape).into_iter().next()?;
        return sketch.constraints.iter().find(|c| c.id == id);
    }
    let curve = size_curve(shape, size)?;
    sketch
        .constraints
        .iter()
        .find(|c| !c.reference && matches_size(&c.kind, curve, size))
}
/// The point that positions the shape: a rectangle's or line's first corner,
/// a circle's or arc's center.
pub fn anchor(sketch: &Sketch, shape: &Shape) -> Option<EntityId> {
    let curve = sketch
        .curves
        .iter()
        .find(|c| Some(c.id) == shape.curves.first().copied())?;
    match (&curve.kind, shape.kind) {
        (
            CurveKind::Line { start, .. },
            ShapeKind::Rectangle | ShapeKind::Line | ShapeKind::Polyline,
        ) => Some(*start),
        (CurveKind::Circle { center, .. }, ShapeKind::Circle)
        | (CurveKind::Arc { center, .. }, ShapeKind::Arc) => Some(*center),
        _ => None,
    }
}
/// Every point the shape's curves use, in first-use order.
pub fn points(sketch: &Sketch, shape: &Shape) -> Vec<EntityId> {
    let mut seen = BTreeSet::new();
    let mut out = vec![];
    for id in &shape.curves {
        let Some(curve) = sketch.curves.iter().find(|c| c.id == *id) else {
            continue;
        };
        let refs = match &curve.kind {
            CurveKind::Line { start, end } => vec![*start, *end],
            CurveKind::Arc {
                center, start, end, ..
            } => vec![*center, *start, *end],
            CurveKind::Circle { center, .. } => vec![*center],
            CurveKind::BSpline { poles, .. } => poles.clone(),
        };
        out.extend(refs.into_iter().filter(|p| seen.insert(*p)));
    }
    out
}
fn find_shape(source: &Source, shape: EntityId) -> Result<Shape, DocumentError> {
    source
        .sketch
        .shapes
        .iter()
        .find(|s| s.id == shape)
        .cloned()
        .ok_or_else(|| invalid("The sketch shape no longer exists"))
}
/// Set one size as a canonical dimension: the existing driving dimension is
/// edited, otherwise a new one is added on the shape's own curve. Expressions
/// (`=…`) bind parameters exactly as other sketch dimensions do.
pub fn set_size(
    document: &Document,
    source: &Source,
    shape: EntityId,
    size: ShapeSize,
    input: &DimensionInput,
) -> Result<Source, DocumentError> {
    let found = find_shape(source, shape)?;
    let curve = size_curve(&found, size).ok_or_else(|| {
        invalid(format!(
            "{} has no {}",
            found.name,
            size.label().to_lowercase()
        ))
    })?;
    if found.kind.corner() {
        // Every side of a chamfer keeps the same setback; one transaction.
        let constraints = corner_constraints(&source.sketch, &found);
        if constraints.is_empty() {
            return Err(invalid(format!(
                "{} has lost its {} dimension",
                found.name,
                size.label().to_lowercase()
            )));
        }
        let mut draft = source.clone();
        for constraint in constraints {
            draft = draft_dimension(document, &draft, constraint, input, false)?;
        }
        return Ok(draft);
    }
    let mut draft = source.clone();
    let constraint = match size_constraint(&draft.sketch, &found, size) {
        Some(existing) => existing.id,
        None => {
            let dimension = Dimension::Literal {
                value: 1.,
                unit: "mm".into(),
            };
            draft.add_constraint(match size {
                ShapeSize::Diameter => ConstraintKind::Diameter { curve, dimension },
                ShapeSize::Radius => ConstraintKind::Radius { curve, dimension },
                _ => ConstraintKind::Length {
                    line: curve,
                    dimension,
                },
            })
        }
    };
    draft_dimension(document, &draft, constraint, input, false)
}
/// Move the shape so its anchor lies at `target` (sketch millimeters) and fix
/// the anchor there. `solved` supplies current positions; every point of the
/// shape, including points shared with other shapes, moves by the same offset
/// so the shape keeps its form and attached shapes stay attached.
pub fn set_position(
    source: &Source,
    solved: &Sketch,
    shape: EntityId,
    target: [f64; 2],
) -> Result<Source, DocumentError> {
    if target.iter().any(|v| !v.is_finite() || v.abs() > 1e9) {
        return Err(invalid("Position must be finite millimeters within ±1e9"));
    }
    let found = find_shape(source, shape)?;
    let anchor = anchor(&source.sketch, &found)
        .ok_or_else(|| invalid(format!("{} has no position point", found.name)))?;
    let at = |id: EntityId| {
        solved
            .points
            .iter()
            .find(|p| p.id == id)
            .map(|p| p.position_mm)
            .ok_or_else(|| invalid("Wait for the current sketch solution before moving a shape"))
    };
    let origin = at(anchor)?;
    let delta = [target[0] - origin[0], target[1] - origin[1]];
    let mut draft = source.clone();
    for id in points(&source.sketch, &found) {
        let p = at(id)?;
        if let Some(point) = draft.sketch.points.iter_mut().find(|q| q.id == id) {
            point.position_mm = [p[0] + delta[0], p[1] + delta[1]];
        }
    }
    if !draft
        .sketch
        .constraints
        .iter()
        .any(|c| matches!(c.kind, ConstraintKind::Fix { point } if point == anchor))
    {
        draft.add_constraint(ConstraintKind::Fix { point: anchor });
    }
    draft.validate()?;
    Ok(draft)
}
