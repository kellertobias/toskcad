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

#[test]
fn public_preflight_spline_and_tensor_metadata_cancellation_keeps_real_source_reusable() {
    use super::{Options, propose};
    use crate::{
        axial_pipe_thread::Error,
        radial_tensor::{Curve, Surface},
    };
    let source = crate::step::pipe_hardware::test_fixture::engineering();
    let input = source.shell().input();
    use crate::axial_pipe_thread::authored::hardware::shell::Curve as ShellCurve;
    // Independent count of actual preflight interrupt checkpoints on the
    // real retained operator tables: scalar charge+size,3 scalars+point size.
    // This identifies the metadata boundary after a COMPLETE original net,
    // rather than guessing a fraction of the whole encoder run.
    let mut prefix = input.vertices.len() * 8;
    let mut spline_metadata = None;
    for curve in &input.curves {
        let checkpoints = match curve {
            ShellCurve::Line(_) => 7,
            ShellCurve::Original(Curve::Spline { poles, .. }) => {
                let n = poles.len() * 7 + 1;
                spline_metadata.get_or_insert(prefix + n);
                n
            }
            ShellCurve::Original(Curve::Circle { .. }) => 23,
        };
        prefix += checkpoints;
    }
    let mut tensor_metadata = None;
    for surface in &input.surfaces {
        let checkpoints = match surface {
            Surface::Tensor(t) => {
                let n = t.poles().len() * 7 + 1;
                tensor_metadata.get_or_insert(prefix + n);
                n
            }
            Surface::Plane(_) => 21,
            Surface::Cone { .. } => 25,
        };
        prefix += checkpoints;
    }
    let stops = [spline_metadata.unwrap(), tensor_metadata.unwrap()];
    assert!(stops[0] < stops[1] && stops[1] <= prefix);
    let observed = std::cell::Cell::new(0usize);
    let count = || {
        observed.set(observed.get() + 1);
        false
    };
    let baseline = propose(&source, Options::default(), Some(&count))
        .unwrap()
        .into_string();
    assert!(observed.get() > prefix);
    for stop in stops {
        assert!(stop < observed.get());
        let visited = std::cell::Cell::new(0usize);
        let cancel = || {
            visited.set(visited.get() + 1);
            visited.get() == stop
        };
        let options = Options {
            additional_resident_bytes: baseline.capacity(),
            ..Default::default()
        };
        assert!(matches!(
            propose(&source, options, Some(&cancel)),
            Err(Error::Cancelled)
        ));
        assert_eq!(visited.get(), stop);
        println!(
            "ACTUAL_PREFLIGHT_CANCEL observed={} stop={stop} visited={}",
            observed.get(),
            visited.get()
        );
        let retry = propose(&source, options, None).unwrap();
        assert_eq!(retry.as_str(), baseline);
        assert_eq!(
            source.capped_thread().classify([0., 0., 3.], None).unwrap(),
            crate::radial_tensor::Location::Inside
        );
    }
}
