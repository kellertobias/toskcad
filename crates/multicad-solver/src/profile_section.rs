//! Solve and validate the owned section sketch of profile stock (CAD-385).
//! The UI section modal and MCP both publish only a solved candidate whose
//! section resolves; the stored positions are then canonical geometry.
use crate::evaluate;
use multicad_core::{
    constraint_sketch::{Sketch, source::Source},
    model::{SketchSection, resolve_sketch_section},
    sketch_draw::{DrawRequest, draft},
};

/// A solved section sketch and its resolved section.
#[derive(Clone, Debug, PartialEq)]
pub struct SolvedSection {
    pub sketch: Sketch,
    pub section: SketchSection,
}

/// Solve `sketch` with its constraints and dimensions. Returns the solved
/// candidate even when it does not yet form a valid section, so an editor
/// can keep drawing; `resolve` reports why it is not a section.
pub fn solve_section_sketch(sketch: &Sketch) -> Result<Sketch, String> {
    let evaluation = evaluate(sketch).map_err(|e| e.to_string())?;
    evaluation.candidate.ok_or_else(|| {
        format!(
            "The section sketch does not solve ({:?}); remove or change the conflicting constraint",
            evaluation.diagnostics.state
        )
    })
}

/// Solve, then require a valid section: the publishing boundary.
pub fn solve_profile_section(sketch: &Sketch) -> Result<SolvedSection, String> {
    let sketch = solve_section_sketch(sketch)?;
    let section = resolve_sketch_section(&sketch).map_err(|e| e.to_string())?;
    Ok(SolvedSection { sketch, section })
}

/// Add drawn geometry with the same typed draw command the sketch editor
/// uses, then solve.
pub fn draw_section(sketch: &Sketch, request: &DrawRequest) -> Result<Sketch, String> {
    let source = Source {
        u: [1., 0., 0.],
        sketch: sketch.clone(),
    };
    let drafted = draft(&source, request).map_err(|e| e.to_string())?;
    solve_section_sketch(&drafted.source.sketch)
}

#[cfg(test)]
mod tests {
    use super::*;
    use multicad_core::{
        constraint_sketch::{Constraint, ConstraintKind, CurveKind, Dimension},
        document::EntityId,
        sketch_draw::{DrawInput, PointInput},
    };

    fn draw(sketch: &Sketch, input: DrawInput) -> Sketch {
        draw_section(
            sketch,
            &DrawRequest {
                input,
                construction: false,
            },
        )
        .unwrap()
    }

    #[test]
    fn drawn_hollow_section_solves_and_a_dimension_resizes_it() {
        let rect = draw(
            &Sketch::default(),
            DrawInput::Rectangle {
                a: PointInput::Position([0., 0.]),
                b: PointInput::Position([40., 30.]),
            },
        );
        let tube = draw(
            &rect,
            DrawInput::Circle {
                center: PointInput::Position([20., 15.]),
                radius_mm: 8.,
            },
        );
        let solved = solve_profile_section(&tube).unwrap();
        assert_eq!(solved.section.size_mm, [40., 30.]);
        assert_eq!(solved.section.holes.len(), 1);
        // Dimensioning the bottom edge to 60 mm moves the solved geometry,
        // and the section follows the solved positions.
        let mut sized = solved.sketch.clone();
        let bottom = sized
            .curves
            .iter()
            .find_map(|c| match c.kind {
                CurveKind::Line { start, end } => {
                    let y = |id| {
                        sized
                            .points
                            .iter()
                            .find(|p| p.id == id)
                            .unwrap()
                            .position_mm[1]
                    };
                    (y(start) == 0. && y(end) == 0.).then_some(c.id)
                }
                _ => None,
            })
            .unwrap();
        sized.constraints.push(Constraint {
            id: EntityId::new(),
            reference: false,
            kind: ConstraintKind::Length {
                line: bottom,
                dimension: Dimension::Literal {
                    value: 60.,
                    unit: "mm".into(),
                },
            },
        });
        let resized = solve_profile_section(&sized).unwrap();
        assert!((resized.section.size_mm[0] - 60.).abs() < 1e-6);
        assert!((resized.section.size_mm[1] - 30.).abs() < 1e-6);
    }

    #[test]
    fn open_drawings_are_refused_at_the_publishing_boundary() {
        let open = draw(
            &Sketch::default(),
            DrawInput::Polyline {
                vertices: vec![
                    PointInput::Position([0., 0.]),
                    PointInput::Position([10., 0.]),
                    PointInput::Position([10., 10.]),
                ],
                closed: false,
            },
        );
        let error = solve_profile_section(&open).unwrap_err();
        assert!(error.contains("open"), "{error}");
    }
}
