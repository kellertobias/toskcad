//! Comparison of Rust vectors with an independent native drawing.
use multicad_kernel::{DrawingProjection, DrawingSolid, Solid, project_drawing};

/// Compare Rust vectors with the drawing of independently built native
/// solids, with the same rules as `assert_matches_native_with`.
pub fn assert_matches_reference(
    rust: &[multicad_kernel::ProjectedVector],
    native_solids: &[&Solid],
    options: DrawingProjection,
    references: &[DrawingProjection],
    native_to_rust: bool,
    label: &str,
) {
    let project = |solids: &[&Solid], options| {
        let inputs: Vec<_> = solids
            .iter()
            .map(|&solid| DrawingSolid { solid, hatch: None })
            .collect();
        project_drawing(&inputs, options).unwrap()
    };
    assert!(
        native_solids
            .iter()
            .all(|solid| !solid.uses_rust_geometry())
    );
    let extra: Vec<_> = references
        .iter()
        .flat_map(|&reference| project(native_solids, reference))
        .collect();
    let native = project(native_solids, options);
    let rust = rust.to_vec();
    let mut worst: f64 = 0.;
    let directions = [("Rust", &rust, &native), ("OpenCascade", &native, &rust)];
    for (from_label, from, to) in directions.into_iter().take(1 + usize::from(native_to_rust)) {
        for vector in from {
            // As in the established worm comparison, a native curve may lie
            // on a Rust curve of any role: rear edges that coincide with a
            // visible edge are published once by Rust and twice natively.
            let spans: Vec<_> = to
                .iter()
                .filter(|other| {
                    other.source_index == vector.source_index
                        && (other.role == vector.role || from_label == "OpenCascade")
                })
                .chain(extra.iter().filter(|other| {
                    from_label == "Rust" && other.source_index == vector.source_index
                }))
                .flat_map(|other| other.preview_mm.windows(2).map(|w| [w[0], w[1]]))
                .collect();
            for pair in vector.preview_mm.windows(2) {
                for [x, y] in [
                    pair[0],
                    [
                        pair[0][0].midpoint(pair[1][0]),
                        pair[0][1].midpoint(pair[1][1]),
                    ],
                    pair[1],
                ] {
                    let distance = spans
                        .iter()
                        .map(|&[[a, b], [c, d]]| {
                            let [u, v] = [c - a, d - b];
                            let along = (((x - a) * u + (y - b) * v) / (u * u + v * v).max(1e-300))
                                .clamp(0., 1.);
                            (x - a - along * u).hypot(y - b - along * v)
                        })
                        .fold(f64::INFINITY, f64::min);
                    worst = worst.max(distance);
                    assert!(
                        distance <= 3. * options.chord_tolerance_mm,
                        "{label}: {from_label} {:?} source {} at [{x}, {y}]: {distance} mm",
                        vector.role,
                        vector.source_index
                    );
                }
            }
        }
    }
    println!("{{\"case\":\"{label}\",\"worst_preview_distance_mm\":{worst}}}");
}
