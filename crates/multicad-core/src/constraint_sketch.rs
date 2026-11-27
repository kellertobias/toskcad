//! Solver-independent sketch source. Distances are explicitly unit-bearing;
//! coordinates and radius guesses use millimeters in the sketch plane.
use crate::document::EntityId;
use serde::{Deserialize, Serialize};

pub mod binding;
pub mod shape_edit;
pub mod source;
mod spline;
pub mod text;
pub use spline::{MAX_SPLINE_POLES, SplineBasis, SplineSample};

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Sketch {
    pub points: Vec<Point>,
    pub curves: Vec<Curve>,
    pub constraints: Vec<Constraint>,
    /// Parametric text elements (CAD-142). Their exact glyph contours are
    /// derived during evaluation and are not solver points or curves.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub texts: Vec<text::SketchText>,
    /// Drawn shapes as recorded by the draw request (CAD-398): which curves
    /// one Rectangle, Circle, … gesture created. Identity is canonical, never
    /// guessed from geometry. Sources saved before shapes were recorded have
    /// none; their curves are listed individually. Curves removed by later
    /// edits are pruned from their shape when the source is written.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub shapes: Vec<Shape>,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Shape {
    pub id: EntityId,
    pub kind: ShapeKind,
    /// Readable, stable name such as "Rectangle 1".
    pub name: String,
    /// Curves in authored drawing order; a rectangle's first curve is its
    /// width edge and its second the height edge. A fillet or chamfer lists
    /// only its new arc or line.
    pub curves: Vec<EntityId>,
    /// For a fillet or chamfer, the shape whose corner it rounds or cuts
    /// (CAD-414); `None` for every other kind, and when that corner's curves
    /// belong to no recorded shape.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub modifies: Option<EntityId>,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ShapeKind {
    Line,
    Rectangle,
    Circle,
    Arc,
    Polyline,
    Spline,
    /// A tangent arc replacing a corner (CAD-414).
    Fillet,
    /// A straight cut replacing a corner (CAD-414).
    Chamfer,
}
impl ShapeKind {
    pub fn label(self) -> &'static str {
        match self {
            Self::Line => "Line",
            Self::Rectangle => "Rectangle",
            Self::Circle => "Circle",
            Self::Arc => "Arc",
            Self::Polyline => "Polyline",
            Self::Spline => "Spline",
            Self::Fillet => "Fillet",
            Self::Chamfer => "Chamfer",
        }
    }
    /// Fillets and chamfers modify another shape's corner.
    pub fn corner(self) -> bool {
        matches!(self, Self::Fillet | Self::Chamfer)
    }
}
/// Longest accepted shape name.
pub const MAX_SHAPE_NAME: usize = 128;
impl Sketch {
    /// Record `curves` as one drawn shape named "<Kind> N", where N is one
    /// past the highest number already used by that kind. Returns its id.
    pub fn record_shape(&mut self, kind: ShapeKind, curves: Vec<EntityId>) -> EntityId {
        self.record(kind, curves, None)
    }
    /// Record a fillet or chamfer `curve` that treats a corner of the shape
    /// `modifies` (CAD-414), named "Fillet N" or "Chamfer N".
    pub fn record_corner(
        &mut self,
        kind: ShapeKind,
        curve: EntityId,
        modifies: Option<EntityId>,
    ) -> EntityId {
        debug_assert!(kind.corner());
        self.record(kind, vec![curve], modifies)
    }
    fn record(
        &mut self,
        kind: ShapeKind,
        curves: Vec<EntityId>,
        modifies: Option<EntityId>,
    ) -> EntityId {
        let prefix = format!("{} ", kind.label());
        let next = self
            .shapes
            .iter()
            .filter_map(|s| s.name.strip_prefix(&prefix)?.parse::<u32>().ok())
            .max()
            .unwrap_or(0)
            + 1;
        let id = EntityId::new();
        self.shapes.push(Shape {
            id,
            kind,
            name: format!("{prefix}{next}"),
            curves,
            modifies,
        });
        id
    }
    /// The recorded shape that owns `curve`, if any.
    pub fn shape_of(&self, curve: EntityId) -> Option<&Shape> {
        self.shapes.iter().find(|s| s.curves.contains(&curve))
    }
    /// Drop curve references that no longer exist and shapes left empty.
    pub fn prune_shapes(&mut self) {
        let curves: std::collections::BTreeSet<_> = self.curves.iter().map(|c| c.id).collect();
        for shape in &mut self.shapes {
            shape.curves.retain(|id| curves.contains(id));
        }
        self.shapes.retain(|s| !s.curves.is_empty());
        // A treatment whose modified shape is gone stands on its own.
        let shapes: std::collections::BTreeSet<_> = self.shapes.iter().map(|s| s.id).collect();
        for shape in &mut self.shapes {
            if shape.modifies.is_some_and(|id| !shapes.contains(&id)) {
                shape.modifies = None;
            }
        }
    }
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Point {
    pub id: EntityId,
    pub position_mm: [f64; 2],
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Curve {
    pub id: EntityId,
    pub construction: bool,
    pub kind: CurveKind,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum CurveKind {
    BSpline {
        poles: Vec<EntityId>,
        basis: SplineBasis,
    },
    Line {
        start: EntityId,
        end: EntityId,
    },
    Arc {
        center: EntityId,
        start: EntityId,
        end: EntityId,
        clockwise: bool,
    },
    Circle {
        center: EntityId,
        radius_mm: f64,
    },
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(untagged, deny_unknown_fields)]
pub enum Dimension {
    Literal {
        value: f64,
        unit: String,
    },
    Expression {
        expression: crate::parameters::BoundExpression,
    },
}
impl Dimension {
    pub fn unit(&self) -> &str {
        match self {
            Self::Literal { unit, .. } => unit,
            Self::Expression { expression } => &expression.unit,
        }
    }
    pub fn expression(&self) -> Option<&crate::parameters::BoundExpression> {
        match self {
            Self::Expression { expression } => Some(expression),
            _ => None,
        }
    }
    pub fn value_in(
        &self,
        target: crate::parameters::Unit,
        parameters: Option<&crate::parameters::EvaluatedParameters>,
    ) -> Result<f64, String> {
        let unit =
            crate::parameters::Unit::parse(self.unit()).ok_or("Unknown sketch dimension unit")?;
        if unit.dimension() != target.dimension() {
            return Err("Incompatible sketch dimension units".into());
        }
        let value = match self {
            Self::Literal { value, .. }
                if unit == target && unit.dimension() == crate::parameters::Dimension::LENGTH =>
            {
                *value
            }
            Self::Literal { value, .. } => {
                let exact = (unit.dimension() == crate::parameters::Dimension::LENGTH)
                    .then(|| crate::parameters::literal_length_in(*value, unit, target))
                    .and_then(Result::ok);
                exact.map_or_else(
                    || *value * unit.scale() / target.scale(),
                    |value| value.as_f64(),
                )
            }
            Self::Expression { expression } => {
                let parameters = parameters
                    .ok_or("Parameter-driven sketch requires an evaluated parameter context")?;
                let exact = (unit.dimension() == crate::parameters::Dimension::LENGTH)
                    .then(|| crate::parameters::exact_length_in(self, target, parameters))
                    .and_then(Result::ok);
                match exact {
                    Some(value) => value.as_f64(),
                    None => expression
                        .evaluate(parameters)
                        .map_err(|e| e.to_string())?
                        .in_unit(target)
                        .ok_or("Incompatible expression dimension")?,
                }
            }
        };
        if !value.is_finite() || value.abs() > 1e9 {
            return Err("Sketch dimension result is nonfinite or exceeds limits".into());
        }
        Ok(value)
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Axis {
    Direct,
    X,
    Y,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Constraint {
    pub id: EntityId,
    /// Reference dimensions measure the solved geometry and add no equations.
    pub reference: bool,
    pub kind: ConstraintKind,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum ConstraintKind {
    Fix {
        point: EntityId,
    },
    Coincident {
        a: EntityId,
        b: EntityId,
    },
    Horizontal {
        line: EntityId,
    },
    Vertical {
        line: EntityId,
    },
    Parallel {
        a: EntityId,
        b: EntityId,
    },
    Perpendicular {
        a: EntityId,
        b: EntityId,
    },
    Tangent {
        a: EntityId,
        b: EntityId,
    },
    Equal {
        a: EntityId,
        b: EntityId,
    },
    Symmetric {
        a: EntityId,
        b: EntityId,
        about: EntityId,
    },
    Distance {
        a: EntityId,
        b: EntityId,
        axis: Axis,
        dimension: Dimension,
    },
    Length {
        line: EntityId,
        dimension: Dimension,
    },
    Radius {
        curve: EntityId,
        dimension: Dimension,
    },
    Diameter {
        curve: EntityId,
        dimension: Dimension,
    },
    Angle {
        from: EntityId,
        to: EntityId,
        dimension: Dimension,
    },
}
impl ConstraintKind {
    pub fn dimension(&self) -> Option<&Dimension> {
        match self {
            Self::Distance { dimension, .. }
            | Self::Length { dimension, .. }
            | Self::Radius { dimension, .. }
            | Self::Diameter { dimension, .. }
            | Self::Angle { dimension, .. } => Some(dimension),
            _ => None,
        }
    }
    pub fn references(&self) -> Vec<EntityId> {
        match self {
            Self::Fix { point } => vec![*point],
            Self::Horizontal { line } | Self::Vertical { line } | Self::Length { line, .. } => {
                vec![*line]
            }
            Self::Radius { curve, .. } | Self::Diameter { curve, .. } => vec![*curve],
            Self::Coincident { a, b }
            | Self::Parallel { a, b }
            | Self::Perpendicular { a, b }
            | Self::Tangent { a, b }
            | Self::Equal { a, b }
            | Self::Distance { a, b, .. } => vec![*a, *b],
            Self::Symmetric { a, b, about } => vec![*a, *b, *about],
            Self::Angle { from, to, .. } => vec![*from, *to],
        }
    }
    pub fn dimension_mut(&mut self) -> Option<&mut Dimension> {
        match self {
            Self::Distance { dimension, .. }
            | Self::Length { dimension, .. }
            | Self::Radius { dimension, .. }
            | Self::Diameter { dimension, .. }
            | Self::Angle { dimension, .. } => Some(dimension),
            _ => None,
        }
    }
}
