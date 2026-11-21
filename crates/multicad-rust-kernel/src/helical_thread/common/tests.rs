use super::*;
use crate::helical_thread::Spec;

fn thread(major: f64, depth: f64, crest: f64, root: Option<f64>) -> Thread {
    Thread::new(Spec {
        major_radius_mm: major,
        pitch_mm: 0.5,
        length_mm: 8.,
        depth_mm: depth,
        crest_width_mm: crest,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 1,
        left_hand: false,
        root_radius_mm: root,
    })
    .unwrap()
}

/// External M3-like thread and a nut void with radial clearance.
fn pair() -> (Thread, Thread) {
    let external = thread(1.48, 0.25, 0.06, None);
    let cavity = thread(1.53, 0.25, 0.1, None);
    (external, cavity)
}

/// Shoelace area of a closed chain of straight meridian pieces.
fn area(lens: &Lens) -> f64 {
    lens.pieces
        .iter()
        .map(|p| (p.phase_mm[0] * p.radius_mm[1] - p.phase_mm[1] * p.radius_mm[0]) / 2.)
        .sum()
}

#[test]
fn a_clear_phase_has_no_lens_and_a_wrong_phase_has_one_closed_lens_per_period() {
    let (external, cavity) = pair();
    // The crests sit in the opposite grooves: every phase is clear.
    assert_eq!(external.common_lenses(&cavity, 0., 1e-7), Ok(vec![]));
    for shift in [0.25, 0.21, 0.29, -0.25, 1.75] {
        let lenses = external.common_lenses(&cavity, shift, 1e-7).unwrap();
        assert_eq!(lenses.len(), 1, "shift {shift}");
        let lens = &lenses[0];
        for (piece, next) in lens.pieces.iter().zip(lens.pieces.iter().cycle().skip(1)) {
            assert_eq!(piece.phase_mm[1], next.phase_mm[0]);
            assert_eq!(piece.radius_mm[1], next.radius_mm[0]);
        }
        // An independent midpoint quadrature of the radial overlap.
        let steps = 400_000;
        let h = 0.5 / f64::from(steps);
        let expected: f64 = (0..steps)
            .map(|i| {
                let phase = (f64::from(i) + 0.5) * h;
                let gap = external.radius_at_phase(phase) - cavity.radius_at_phase(phase - shift);
                gap.max(0.) * h
            })
            .sum();
        assert!(
            (area(lens) - expected).abs() < 1e-8,
            "{} vs {expected}",
            area(lens)
        );
        assert!(area(lens) > 1e-3);
    }
}

#[test]
fn multi_start_threads_repeat_the_lens_once_per_start() {
    let spec = |major, crest| Spec {
        major_radius_mm: major,
        pitch_mm: 0.5,
        length_mm: 8.,
        depth_mm: 0.25,
        crest_width_mm: crest,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 3,
        left_hand: true,
        root_radius_mm: None,
    };
    let external = Thread::new(spec(1.48, 0.06)).unwrap();
    let cavity = Thread::new(spec(1.53, 0.1)).unwrap();
    let lenses = external.common_lenses(&cavity, 0.25, 1e-7).unwrap();
    assert_eq!(lenses.len(), 3);
    for (k, lens) in lenses.iter().enumerate() {
        let shift = lens.pieces[0].phase_mm[0] - lenses[0].pieces[0].phase_mm[0];
        assert!((shift - 0.5 * k as f64).abs() < 1e-12);
    }
}

#[test]
fn contacts_vertex_crossings_rounded_roots_and_mismatched_laws_are_refused() {
    let (external, cavity) = pair();
    let other = Thread::new(Spec {
        pitch_mm: 0.6,
        ..cavity.spec()
    })
    .unwrap();
    assert_eq!(
        external.common_lenses(&other, 0.25, 1e-7),
        Err(CommonError::Invalid)
    );
    assert_eq!(
        external.common_lenses(&cavity, f64::NAN, 1e-7),
        Err(CommonError::Invalid)
    );
    assert_eq!(
        external.common_lenses(&cavity, 0.25, 0.),
        Err(CommonError::Invalid)
    );
    // A crest corner of one law exactly on the other law.
    let corner = external.profile_breaks()[1];
    let g = |shift: f64| external.radius_at_phase(corner) - cavity.radius_at_phase(corner - shift);
    let grid: Vec<f64> = (0..=2000).map(|i| f64::from(i) * 0.5 / 2000.).collect();
    let pair = grid
        .windows(2)
        .find(|pair| g(pair[0]) < 0. && g(pair[1]) > 0.)
        .unwrap();
    let (mut low, mut high) = (pair[0], pair[1]);
    for _ in 0..200 {
        let middle = low.midpoint(high);
        if g(middle) < 0. {
            low = middle;
        } else {
            high = middle;
        }
    }
    assert!(g(low).abs() < 1e-12);
    assert_eq!(
        external.common_lenses(&cavity, low, 1e-7),
        Err(CommonError::Unsupported)
    );
    // A cavity whose minor crest reaches a rounded external root. The
    // root depth follows from the flank intersection and fillet radius.
    let (radius, half) = (0.04, 30_f64.to_radians());
    let depth = (0.5 - 0.06) / (2. * half.tan()) - radius * (1. / half.sin() - 1.);
    let rounded = thread(1.48, depth, 0.06, Some(radius));
    let deep = thread(1.53, 0.43, 0.003, None);
    assert!(rounded.radius_at_phase(0.25) > deep.radius_at_phase(0.25));
    assert_eq!(
        rounded.common_lenses(&deep, 0., 1e-7),
        Err(CommonError::Unsupported)
    );
}
