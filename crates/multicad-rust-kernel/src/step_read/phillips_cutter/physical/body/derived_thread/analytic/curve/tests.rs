use super::*;
#[test]
fn numeric_world_transfer_guard_refuses_overflow_before_publishing_a_bound() {
    // Pure numeric helper frontier, not an admitted Body or a file fixture.
    // Preserve the defensive refusal for future callers with different scales.
    let g = Factors {
        a: vec![[1e-8; 2]],
        au: vec![1.],
        b: [1e-8; 2],
        bv: 1.,
        z: [0.; 2],
        major: 1.,
    };
    let mut work = Work {
        remaining: WORK,
        interrupt: None,
        cancelled: Cell::new(false),
    };
    let frame = [[0.; 3], [f64::MAX, 0., 0.], [0., 1., 0.]];
    assert!(world_transfer(0, &Target::raw(frame), [0.; 3], &g, &mut work).is_err());
    work.remaining = 0;
    assert!(world_transfer(0, &Target::raw([[0.; 3]; 3]), [0.; 3], &g, &mut work).is_err());
}
