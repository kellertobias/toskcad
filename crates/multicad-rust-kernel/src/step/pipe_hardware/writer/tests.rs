// The codec's fresh valid Writer normally dominates these transport errors.
// Exercise the adapter's own contract using real public exporter refusals;
// this does not claim those failures are reachable through valid propose().
#[test]
fn real_public_exporter_refusals_keep_typed_transport_and_atomic_members() {
    use super::AppendError;
    use crate::step::{MAX_COMPOUND_SOLIDS, StepError, Writer};
    use crate::{Rational, prism::Prism};
    let origin = [0.; 3];
    let x = [1., 0., 0.];
    let y = [0., 1., 0.];
    let q = |n| Rational::new(n, 1).unwrap();
    let size = [q(2), q(3), q(4)];
    let check = |actual: StepError, expected_step, expected| {
        assert_eq!(actual, expected_step);
        let AppendError(typed) = actual.into();
        assert_eq!(typed, expected);
    };
    check(
        Writer::new().finish().unwrap_err(),
        StepError::InvalidGeometry,
        crate::axial_pipe_thread::Error::Invalid,
    );
    let mut writer = Writer::new();
    writer.append_box(size, origin, x, y).unwrap();
    check(
        writer.append_box(size, origin, x, x).unwrap_err(),
        StepError::InvalidPlacement,
        crate::axial_pipe_thread::Error::Invalid,
    );
    let rounded = Prism::extrude(1, q(3), q(2), q(4)).unwrap();
    check(
        writer
            .append_rectangular_prism(&rounded, origin, x, y)
            .unwrap_err(),
        StepError::Unsupported,
        crate::axial_pipe_thread::Error::Unsupported,
    );
    // Both rejected appends leave the real accepted body byte-for-byte intact.
    let actual = writer.finish().unwrap();
    let mut retry = Writer::new();
    retry.append_box(size, origin, x, y).unwrap();
    assert_eq!(actual, retry.finish().unwrap());
    // Fill the public member frontier without forging the private ledger.
    let mut full = Writer::new();
    for _ in 0..MAX_COMPOUND_SOLIDS {
        full.append_box(size, origin, x, y).unwrap();
    }
    check(
        full.append_box(size, origin, x, y).unwrap_err(),
        StepError::Limit,
        crate::axial_pipe_thread::Error::Limit,
    );
    let accepted = full.finish().unwrap();
    assert_eq!(
        accepted.matches("MANIFOLD_SOLID_BREP(").count(),
        MAX_COMPOUND_SOLIDS
    );
}
