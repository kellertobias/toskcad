//! Independent OpenCascade references for Rust drawing regressions.
use multicad_kernel::{DrawingProjection, DrawingSolid, Solid, project_drawing};

mod reference;
pub use reference::assert_matches_reference;

/// Every Rust preview point lies on a native curve of the same role and
/// source, and every native point on a Rust curve of the same source,
/// within three chord tolerances. OpenCascade reads
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
    let native: Vec<_> = native.iter().collect();
    assert_matches_reference(&rust, &native, options, references, native_to_rust, label);
}
