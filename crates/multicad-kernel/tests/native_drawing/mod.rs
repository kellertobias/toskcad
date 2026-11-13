//! Independent OpenCascade references for Rust drawing regressions.
use multicad_kernel::{DrawingProjection, DrawingSolid, Solid, project_drawing};

/// Every Rust preview point lies on a native curve of the same role and
/// source, and vice versa, within three chord tolerances. OpenCascade reads
/// the same finite solids from the Rust STEP export, so the comparison
/// isolates drawing parity from different construction histories.
pub fn assert_matches_native(solids: &[&Solid], options: DrawingProjection, label: &str) {
    assert_matches_native_with(solids, options, &[], true, label);
}

/// As `assert_matches_native`, but a Rust curve may also lie on any curve
/// of the native `references` projections when OpenCascade's own result for
/// `options` is known to be incomplete. With `native_to_rust`, native curves
/// must also lie on Rust curves of the same role; callers disable it only
/// where the native result is measured to leave the true geometry.
pub fn assert_matches_native_with(
    solids: &[&Solid],
    options: DrawingProjection,
    references: &[DrawingProjection],
    native_to_rust: bool,
    label: &str,
) {
    let project = |solids: &[&Solid]| {
        let inputs: Vec<_> = solids
            .iter()
            .map(|&solid| DrawingSolid { solid, hatch: None })
            .collect();
        project_drawing(&inputs, options).unwrap()
    };
    let rust = project(solids);
    let directory = tempfile::tempdir().unwrap();
    let native: Vec<_> = solids
        .iter()
        .enumerate()
        .map(|(index, solid)| {
            let path = directory.path().join(format!("solid-{index}.step"));
            solid.write_step(&path).unwrap();
            Solid::read_step(&path, solid.geometry_policy()).unwrap()
        })
        .collect();
    assert!(native.iter().all(|solid| !solid.uses_rust_geometry()));
    let native_solids: Vec<_> = native.iter().collect();
    let extra: Vec<_> = references
        .iter()
        .flat_map(|&reference| {
            let inputs: Vec<_> = native_solids
                .iter()
                .map(|&solid| DrawingSolid { solid, hatch: None })
                .collect();
            project_drawing(&inputs, reference).unwrap()
        })
        .collect();
    let native = project(&native_solids);
    let mut worst: f64 = 0.;
    let directions = [("Rust", &rust, &native), ("OpenCascade", &native, &rust)];
    for (from_label, from, to) in directions.into_iter().take(1 + usize::from(native_to_rust)) {
        for vector in from {
            let spans: Vec<_> = to
                .iter()
                .filter(|other| {
                    other.source_index == vector.source_index && other.role == vector.role
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
