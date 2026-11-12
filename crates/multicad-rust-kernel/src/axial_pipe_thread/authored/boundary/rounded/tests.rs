use super::*;
#[test]
fn actual_periodic_rails_cannot_be_cancelled_or_perturbed_into_an_owned_graph() {
    let thread =
        AuthoredThread::new(super::super::super::tests::bspt_spec(1. / 32.), 1e-7, 1e-6).unwrap();
    let raw = thread.raw_boundary_ownership(1e-7, None).unwrap();
    let original = raw.ownership.input().clone();
    let nets = || {
        original
            .surfaces
            .iter()
            .filter_map(|s| {
                if let Surface::Tensor(t) = s {
                    Some(Net {
                        tensor: t.clone(),
                        error: raw.tensor_error_mm,
                    })
                } else {
                    None
                }
            })
            .collect::<Vec<_>>()
    };
    assert_eq!(
        build(&thread, nets(), raw.tensor_error_mm, 1e-7, Some(&|| true)).err(),
        Some(Error::Cancelled)
    );
    let calls = std::cell::Cell::new(0);
    let cancel = || {
        calls.set(calls.get() + 1);
        calls.get() >= 5
    };
    assert_eq!(
        build(&thread, nets(), raw.tensor_error_mm, 1e-7, Some(&cancel)).err(),
        Some(Error::Cancelled)
    );
    let mut wrong = nets();
    let last = wrong.last_mut().unwrap();
    let spans = last.tensor.spans();
    let [w, h] = spans.map(|s| DEGREE * s + 1);
    let mut points = last.tensor.poles().to_vec();
    points[(w - 1) * h + h / 2][0] = points[(w - 1) * h + h / 2][0].next_up();
    last.tensor = RetainedTensor::new(spans, points).unwrap();
    assert_eq!(
        build(&thread, wrong, raw.tensor_error_mm, 1e-7, None).err(),
        Some(Error::Unresolved)
    );
    let retry = build(&thread, nets(), raw.tensor_error_mm, 1e-7, None).unwrap();
    assert_eq!(retry.ownership.input(), &original);
}
