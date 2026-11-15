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
    let [u, v, normal] = planar::view_axes(options.frame);
    let origin: [f64; 3] = std::array::from_fn(|i| options.frame.origin[i] + normal[i] * offset);
    let flat = |p: [f64; 3]| {
        let d: [f64; 3] = std::array::from_fn(|i| p[i] - options.frame.origin[i]);
        let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        [dot(d, u), dot(d, v)]
    };
    let mut vectors = Vec::new();
    for (source, (input, bodies)) in solids.iter().zip(bodies).enumerate() {
        // Each solid of a compound is hatched as its own region: solids may
        // touch, and their loops must not be read as one even-odd area.
        for body in bodies {
            let deflection = options.chord_tolerance_mm.max(1e-3);
            let rings = match body.section(origin, normal, deflection) {
                Ok(loops) => loops
                    .into_iter()
                    .map(|ring| ring.into_iter().map(flat).collect())
                    .collect(),
                Err(error) => return Some(Err(error)),
            };
            match section::sampled::publish(input, rings, options) {
                Ok(published) => {
                    vectors.extend(published.into_iter().map(|vector| ProjectedVector {
                        source_index: source,
                        ..vector
                    }))
                }
                Err(error) => return Some(Err(error)),
            }
        }
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
