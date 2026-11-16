//! Sections of solids that retain an imported boundary.
use super::*;

/// The retained members of a solid, or `None` when any part of it is not
/// a retained imported boundary.
fn members(solid: &Solid) -> Option<Vec<&crate::geometry::imported_brep_body::ImportedBrepBody>> {
    if let Some(body) = solid.imported_brep_body() {
        return Some(vec![body]);
    }
    solid
        .rust_compound()?
        .members()
        .map(Solid::imported_brep_body)
        .collect()
}

/// Cut every solid with the section plane on its exact boundary and publish
/// the loops with their hatch. `None` leaves scenes with other bodies to the
/// general projector.
pub(super) fn section(
    solids: &[DrawingSolid<'_>],
    options: DrawingProjection,
    offset: f64,
) -> Option<Result<Vec<ProjectedVector>, GeometryError>> {
    let bodies = solids
        .iter()
        .map(|input| members(input.solid))
        .collect::<Option<Vec<_>>>()?;
    let [_, _, normal] = planar::view_axes(options.frame);
    let origin: [f64; 3] = std::array::from_fn(|i| options.frame.origin[i] + normal[i] * offset);
    let mut vectors = Vec::new();
    // One hatch budget for the whole section, as in the other section
    // projectors; every published curve and point counts against the
    // drawing's output limits.
    let mut hatch_budget = multicad_rust_kernel::section_hatch::Budget::default();
    for (source, (input, bodies)) in solids.iter().zip(bodies).enumerate() {
        let policy = input.solid.geometry_policy();
        // Each solid of a compound is hatched as its own region: solids may
        // touch, and their loops must not be read as one even-odd area.
        for body in bodies {
            let deflection = options.chord_tolerance_mm.max(1e-3);
            let published = body
                .section_edges(origin, normal, deflection)
                .and_then(|loops| {
                    section::warped::publish_world(
                        loops,
                        options,
                        policy,
                        source,
                        input.hatch.is_some(),
                    )
                })
                .and_then(|(mut output, boundary)| {
                    if let Some(hatch) = input.hatch {
                        let mut curves = options.maximum_curves as usize;
                        let mut points = options.maximum_points as usize;
                        output.extend(section::publish_hatch(
                            source,
                            hatch,
                            &boundary,
                            policy.linear_tolerance_mm,
                            &mut curves,
                            &mut points,
                            &mut hatch_budget,
                        )?);
                    }
                    Ok(output)
                });
            match published {
                Ok(published) => vectors.extend(published),
                Err(error) => return Some(Err(error)),
            }
        }
    }
    let points: usize = vectors.iter().map(|vector| vector.preview_mm.len()).sum();
    if vectors.len() > options.maximum_curves as usize || points > options.maximum_points as usize {
        return Some(Err(err(
            GeometryErrorKind::Limit,
            "Rust drawing section exceeded its output budget",
        )));
    }
    Some(Ok(vectors))
}

/// Section solids through their exported boundary: each solid is written
/// to STEP with its exact supports and read back as a retained boundary.
/// `None` when a solid has no exact export or the export does not read
/// back, so the family refusal stands.
pub(super) fn section_through_export(
    solids: &[DrawingSolid<'_>],
    options: DrawingProjection,
    offset: f64,
) -> Option<Result<Vec<ProjectedVector>, GeometryError>> {
    let retained = solids
        .iter()
        .map(|input| input.solid.retained_boundary())
        .collect::<Option<Vec<_>>>()?;
    let inputs: Vec<_> = solids
        .iter()
        .zip(&retained)
        .map(|(input, solid)| DrawingSolid {
            solid,
            hatch: input.hatch,
        })
        .collect();
    section(&inputs, options, offset)
}
