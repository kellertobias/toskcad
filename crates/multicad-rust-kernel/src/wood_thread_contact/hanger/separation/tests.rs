use super::*;
use crate::{
    Rational,
    helical_thread::{Spec, Thread},
    orthogonal::{self, Box3, Operation},
};
fn identity() -> Frame {
    Frame {
        origin_mm: [0.; 3],
        axes: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
    }
}
fn stock() -> OrthogonalSolid {
    orthogonal::boolean(
        &[Box3 {
            min: [0.; 3],
            max: [1.; 3],
        }],
        &[],
        Operation::Union,
        0.,
    )
    .unwrap()
    .unwrap()
}
fn nut() -> Nut {
    Nut::new(
        Thread::new(Spec {
            major_radius_mm: 1.5,
            pitch_mm: 0.5,
            length_mm: 2.4,
            depth_mm: 0.2,
            crest_width_mm: 0.1,
            lower_flank_deg: 30.,
            upper_flank_deg: 30.,
            starts: 1,
            left_hand: false,
            root_radius_mm: None,
        })
        .unwrap(),
        5.5,
    )
    .unwrap()
}
#[test]
fn every_actual_component_source_earns_far_separation_without_classifying_overlapping_envelopes() {
    let hanger = Hanger::new(
        crate::authored_hanger_bolt::tests::fixture(true, true),
        None,
    )
    .unwrap();
    let nut = nut();
    let ring = Annulus::new(1.7, 3.5, 0.2).unwrap();
    let stock = stock();
    let prism = Prism::rectangular(7, [Rational::new(21, 10).unwrap(); 3]).unwrap();
    let sources = [
        HangerComponent::Hanger(&hanger),
        HangerComponent::Nut(&nut),
        HangerComponent::Annulus(&ring),
        HangerComponent::Stock(&stock),
        HangerComponent::Prism(&prism),
    ];
    for source in sources {
        for posed in [false, true] {
            let frame = if posed {
                Frame {
                    origin_mm: [11., -7., 5.],
                    axes: [
                        [0.8, 0.6, 0.],
                        [-0.36, 0.48, 0.8],
                        [0.48, -0.6400000000000001, 0.6],
                    ],
                }
            } else {
                identity()
            };
            let a = PlacedHangerComponent { source, frame };
            let b = PlacedHangerComponent {
                source: HangerComponent::Stock(&stock),
                frame: Frame {
                    origin_mm: [100.; 3],
                    ..identity()
                },
            };
            assert_eq!(
                hanger_component_separation(a, b, &|| false).unwrap(),
                Some(true)
            );
            assert_eq!(
                hanger_component_separation(a, a, &|| false).unwrap(),
                Some(false)
            );
        }
    }
    // This stock is inside the annulus' empty bore. Overlapping outer
    // envelopes remain unresolved, even for these genuinely disjoint solids.
    let inside = orthogonal::boolean(
        &[Box3 {
            min: [-0.5, -0.5, 0.],
            max: [0.5, 0.5, 0.2],
        }],
        &[],
        Operation::Union,
        0.,
    )
    .unwrap()
    .unwrap();
    assert!(2. * 0.5_f64.powi(2) < ring.inner_radius_mm().powi(2));
    assert_eq!(
        hanger_component_separation(
            PlacedHangerComponent {
                source: HangerComponent::Annulus(&ring),
                frame: identity()
            },
            PlacedHangerComponent {
                source: HangerComponent::Stock(&inside),
                frame: identity()
            },
            &|| false
        )
        .unwrap(),
        Some(false)
    );
}
#[test]
fn directed_raw_maps_contain_actual_vertices_and_original_rational_prism_corners() {
    let nut = nut();
    let hanger = Hanger::new(
        crate::authored_hanger_bolt::tests::fixture(false, false),
        None,
    )
    .unwrap();
    let source = HangerComponent::Nut(&nut);
    for origin in [0.3, 1e6 + 0.3] {
        let frame = Frame {
            origin_mm: [origin, -7., 5.],
            axes: [[0.8, 0.6, 0.], [-0.6, 0.8, 0.], [0., 0., 1.]],
        };
        let boxes = [
            (
                source,
                crate::headed_screw::hex_vertices(nut.across_flats_mm())
                    .into_iter()
                    .flat_map(|p| [0., 2.4].map(|z| [p[0], p[1], z]))
                    .collect::<Vec<_>>(),
            ),
            (HangerComponent::Hanger(&hanger), {
                let b = hanger.bounds_mm();
                (0..8)
                    .map(|c| std::array::from_fn(|i| b[i][usize::from(c & (1 << i) != 0)]))
                    .collect()
            }),
        ];
        for (source, points) in boxes {
            let b = bounds(PlacedHangerComponent { source, frame })
                .unwrap()
                .unwrap();
            for p in points {
                for (i, b) in b.iter().enumerate() {
                    let actual = (0..3).fold(q(frame.origin_mm[i]).unwrap(), |s, j| {
                        s + q(frame.axes[j][i]).unwrap() * q(p[j]).unwrap()
                    });
                    assert!(q(b.lo).unwrap() <= actual && q(b.hi).unwrap() >= actual);
                }
            }
        }
        let prism = Prism::rectangular(7, [Rational::new(21, 10).unwrap(); 3]).unwrap();
        let b = bounds(PlacedHangerComponent {
            source: HangerComponent::Prism(&prism),
            frame,
        })
        .unwrap()
        .unwrap();
        for c in 0..8 {
            for (i, b) in b.iter().enumerate() {
                let actual = (0..3).fold(q(frame.origin_mm[i]).unwrap(), |s, j| {
                    s + q(frame.axes[j][i]).unwrap()
                        * if c & (1 << j) == 0 {
                            Q::zero()
                        } else {
                            Q::new(21.into(), 10.into())
                        }
                });
                assert!(q(b.lo).unwrap() <= actual && q(b.hi).unwrap() >= actual);
            }
        }
    }
}
#[test]
fn touching_ulp_uncertainty_unknown_profiles_and_cancellation_never_become_disjointness() {
    use std::cell::Cell;
    let ring = Annulus::new(1.7, 3.5, 0.2).unwrap();
    let stock = stock();
    let a = PlacedHangerComponent {
        source: HangerComponent::Annulus(&ring),
        frame: identity(),
    };
    for x in [3.5, 3.5_f64.next_up()] {
        let b = PlacedHangerComponent {
            source: HangerComponent::Stock(&stock),
            frame: Frame {
                origin_mm: [x, 0., 0.],
                ..identity()
            },
        };
        assert_eq!(
            hanger_component_separation(a, b, &|| false).unwrap(),
            Some(false)
        );
    }
    let b = PlacedHangerComponent {
        source: HangerComponent::Stock(&stock),
        frame: Frame {
            origin_mm: [3.5 + 1e-8, 0., 0.],
            ..identity()
        },
    };
    assert_eq!(
        hanger_component_separation(a, b, &|| false).unwrap(),
        Some(true)
    );
    let rounded = Prism::extrude(
        7,
        Rational::new(10, 1).unwrap(),
        Rational::new(2, 1).unwrap(),
        Rational::new(2, 1).unwrap(),
    )
    .unwrap();
    let unknown = PlacedHangerComponent {
        source: HangerComponent::Prism(&rounded),
        frame: identity(),
    };
    assert_eq!(
        hanger_component_separation(a, unknown, &|| false).unwrap(),
        None
    );
    for stop in [0, 1] {
        let at = Cell::new(0);
        assert!(matches!(
            hanger_component_separation(a, unknown, &|| {
                let n = at.get();
                at.set(n + 1);
                n == stop
            }),
            Err(Error::Cancelled)
        ));
    }
    let mut bad = a;
    bad.frame.axes[0][0] = 1.01;
    assert!(matches!(
        hanger_component_separation(bad, b, &|| false),
        Err(Error::Unsupported)
    ));
}
