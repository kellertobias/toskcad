use super::*;
use crate::{
    orthogonal::{self, Box3, Operation},
    orthogonal_pockets::{OrthogonalPocket, OrthogonalThroughBore},
};
fn identity() -> Frame {
    Frame {
        origin_mm: [0.; 3],
        axes: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
    }
}
struct Fixture {
    wood_frame: Frame,
    clamp_frame: Frame,
    hanger_frame: Frame,
    wood: OrthogonalSolid,
    pilot: OrthogonalPocket,
    clamp: OrthogonalSolid,
    bore: OrthogonalThroughBore,
    hanger: Hanger,
    nut: Nut,
    nut_frame: Frame,
    ring: Annulus,
    ring_frame: Frame,
}
fn fixture(left: bool) -> Fixture {
    let solid = |min, max| {
        orthogonal::boolean(&[Box3 { min, max }], &[], Operation::Union, 0.)
            .unwrap()
            .unwrap()
    };
    let wood = solid([-10., -10., 0.], [10., 10., 20.]);
    let boundary = orthogonal_boundary::build(wood.quads()).unwrap();
    let pilot = orthogonal_pockets::identify_pocket(
        &wood,
        &boundary,
        [0., 0., -0.001],
        [0., 0., 16.],
        2.,
        1e-7,
        1e-5,
    )
    .unwrap();
    let clamp = solid([-10., -10., -4.1], [10., 10., 0.]);
    let boundary = orthogonal_boundary::build(clamp.quads()).unwrap();
    let bore = orthogonal_pockets::identify_through_bore(
        &clamp,
        &boundary,
        [0., 0., -5.],
        [0., 0., 1.],
        1.7,
        1e-7,
        1e-5,
    )
    .unwrap();
    let hanger = Hanger::new(
        crate::authored_hanger_bolt::tests::fixture(left, false),
        None,
    )
    .unwrap();
    let nut = Nut::new(
        crate::helical_thread::Thread::new(crate::helical_thread::Spec {
            major_radius_mm: 1.5,
            pitch_mm: 0.5,
            length_mm: 2.4,
            depth_mm: (3. - (2.459 + 2.599) / 2.) / 2.,
            crest_width_mm: 0.25 - (3. - (2.675 + 2.775) / 2.) / 3_f64.sqrt(),
            lower_flank_deg: 30.,
            upper_flank_deg: 30.,
            starts: 1,
            left_hand: left,
            root_radius_mm: None,
        })
        .unwrap(),
        5.5,
    )
    .unwrap();
    let station = -6.7;
    let start = -hanger.spec().plain_length_mm - hanger.metric().full_length_mm();
    let angle = std::f64::consts::TAU * (station - start) / 0.5 * if left { -1. } else { 1. };
    let (s, c) = angle.portable_sin_cos();
    let nut_frame = Frame {
        origin_mm: [0., 0., station],
        axes: [[c, s, 0.], [-s, c, 0.], [0., 0., c * c + s * s]],
    };
    let ring = Annulus::new(1.7, 3.5, 0.2).unwrap();
    let ring_frame = Frame {
        origin_mm: [0., 0., -4.3],
        ..identity()
    };
    Fixture {
        wood_frame: identity(),
        clamp_frame: identity(),
        hanger_frame: identity(),
        wood,
        pilot,
        clamp,
        bore,
        hanger,
        nut,
        nut_frame,
        ring,
        ring_frame,
    }
}
fn requirements() -> StackRequirements {
    StackRequirements {
        wood: Requirements {
            minimum_wall_mm: 1.,
            minimum_floor_mm: 1.,
            minimum_engagement_mm: 4.,
            tip_clearance_mm: 1.,
            linear_tolerance_mm: 1e-7,
            minimum_feature_mm: 1e-5,
        },
        minimum_clamp_wall_mm: 1.,
        passage: PassageRequirements {
            minimum_radial_clearance_mm: 1e-5,
            linear_tolerance_mm: 1e-7,
            minimum_feature_mm: 1e-5,
        },
        mating: HangerMatingRequirements {
            minimum_engagement_mm: 2.,
            ..Default::default()
        },
        minimum_bearing_land_mm: 1e-5,
        maximum_bearing_gap_mm: 1e-7,
    }
}
fn nut_station(f: &mut Fixture, station: f64) {
    let start = -f.hanger.spec().plain_length_mm - f.hanger.metric().full_length_mm();
    let angle = std::f64::consts::TAU * (station - start) / f.hanger.metric().spec().pitch_mm
        * if f.hanger.metric().spec().left_hand {
            -1.
        } else {
            1.
        };
    let (s, c) = angle.portable_sin_cos();
    f.nut_frame = Frame {
        origin_mm: [0., 0., station],
        axes: [[c, s, 0.], [-s, c, 0.], [0., 0., c * c + s * s]],
    };
}
fn transport(a: Frame, b: Frame) -> Frame {
    Frame {
        origin_mm: std::array::from_fn(|i| {
            a.origin_mm[i] + (0..3).map(|j| a.axes[j][i] * b.origin_mm[j]).sum::<f64>()
        }),
        axes: std::array::from_fn(|k| {
            std::array::from_fn(|i| (0..3).map(|j| a.axes[j][i] * b.axes[k][j]).sum())
        }),
    }
}
fn query(
    f: &Fixture,
    rings: &[PlacedAnnulus<'_>],
    req: StackRequirements,
    interrupt: &dyn Fn() -> bool,
) -> R<HangerStack> {
    certify_hanger_stack(
        StackSources {
            wood: Stock {
                material: &f.wood,
                pilot: f.pilot,
                frame: f.wood_frame,
            },
            clamp: ThroughStock {
                material: &f.clamp,
                bore: f.bore,
                frame: f.clamp_frame,
            },
            hanger: &f.hanger,
            hanger_frame: f.hanger_frame,
            washers: rings,
            nut: &f.nut,
            nut_frame: f.nut_frame,
        },
        req,
        interrupt,
    )
}

#[test]
fn positive_gaps_are_bounded_but_any_actual_intrusion_is_rejected() {
    let f = fixture(false);
    for offset in [-1e-8, 1e-8, -1e-6, 1e-6] {
        let ring = PlacedAnnulus {
            source: &f.ring,
            frame: Frame {
                origin_mm: [0., 0., -4.3 + offset],
                ..identity()
            },
        };
        // Moving this washer either way intrudes into one of its two actual
        // neighbours; tolerance cannot turn either penetration into contact.
        assert!(matches!(
            query(&f, &[ring], requirements(), &|| false),
            Err(Error::Unsupported)
        ));
    }
    let narrow = Annulus::new(1.7, 3.5, 0.19999998).unwrap();
    let ring = PlacedAnnulus {
        source: &narrow,
        frame: f.ring_frame,
    };
    let result = query(&f, &[ring], requirements(), &|| false).unwrap();
    assert!(result.bearings[1].correspondence_gap_upper_mm > 0.);
    let strict = StackRequirements {
        maximum_bearing_gap_mm: 1e-9,
        ..requirements()
    };
    assert!(matches!(
        query(&f, &[ring], strict, &|| false),
        Err(Error::Unsupported)
    ));
    let wide_land = StackRequirements {
        minimum_bearing_land_mm: 1.1,
        ..requirements()
    };
    assert!(matches!(
        query(&f, &[ring], wide_land, &|| false),
        Err(Error::Unsupported)
    ));
}

#[test]
fn cancellation_at_every_publication_prefix_is_atomic_and_retries_keep_sources() {
    use std::cell::Cell;
    let f = fixture(false);
    let ring = PlacedAnnulus {
        source: &f.ring,
        frame: f.ring_frame,
    };
    let count = Cell::new(0);
    let original = query(&f, &[ring], requirements(), &|| {
        count.set(count.get() + 1);
        false
    })
    .unwrap();
    let total = count.get();
    assert!(total > 20);
    for stop in [0, 1, total / 3, total / 2, total - 2, total - 1] {
        let at = Cell::new(0);
        assert!(matches!(
            query(&f, &[ring], requirements(), &|| {
                let n = at.get();
                at.set(n + 1);
                n == stop
            }),
            Err(Error::Cancelled)
        ));
    }
    let retry = query(&f, &[ring], requirements(), &|| false).unwrap();
    assert_eq!(format!("{original:?}"), format!("{retry:?}"));
    for value in [f64::NAN, -1., 0., 1e8] {
        assert!(matches!(
            query(
                &f,
                &[ring],
                StackRequirements {
                    minimum_bearing_land_mm: value,
                    ..requirements()
                },
                &|| false
            ),
            Err(Error::Invalid)
        ));
    }
    for value in [f64::NAN, -1., 1e8] {
        assert!(matches!(
            query(
                &f,
                &[ring],
                StackRequirements {
                    maximum_bearing_gap_mm: value,
                    ..requirements()
                },
                &|| false
            ),
            Err(Error::Invalid)
        ));
    }
    assert!(matches!(
        query(&f, &vec![ring; 63], requirements(), &|| false),
        Err(Error::Limit)
    ));
}

#[test]
fn real_caps_and_full_stock_corners_control_bearing_not_only_the_hole_axis() {
    let f = fixture(false);
    let wood = Placed::new(Material::Cells(&f.wood), identity()).unwrap();
    let face = cap(&wood, 2, [Q::zero(), Q::zero(), Q::zero()], -1, [2., 3.]);
    let angle = 0.001_f64;
    let (s, c) = angle.portable_sin_cos();
    let tilted = Placed::new(
        Material::Cells(&f.clamp),
        Frame {
            origin_mm: [0., 0., -0.004],
            axes: [[1., 0., 0.], [0., c, s], [0., -s, c]],
        },
    )
    .unwrap();
    // The full radius-three wood envelope is clear of this clamp's cap,
    // while far corners of the real twenty-wide plate cross the wood plane.
    certify_hanger_through_stock(
        ThroughStock {
            material: &f.clamp,
            bore: f.bore,
            frame: tilted.frame,
        },
        &f.hanger,
        identity(),
        1.,
        requirements().passage,
        &|| false,
    )
    .unwrap();
    assert!(matches!(
        outside(&face, &tilted, &|| false),
        Err(Error::Unsupported)
    ));
    let ring = Placed::new(Material::Ring(&f.ring), identity()).unwrap();
    let map = map(identity(), identity(), &|| false).unwrap();
    // A transverse actual cap cannot be assigned as an axial bearing plane.
    let transverse = cap(&wood, 0, [Q::zero(), Q::zero(), Q::zero()], 1, [2., 3.]);
    assert!(matches!(
        member_near_cap(&transverse, &ring, &map, &|| false),
        Err(Error::Unsupported)
    ));
    assert!(matches!(
        member_near_cap(&face, &wood, &map, &|| false),
        Err(Error::Invalid)
    ));
    let same_facing = ring_cap(&ring, &f.ring, false).unwrap();
    assert!(matches!(
        bearing(&face, &same_facing, &ring, &map, requirements(), &|| false),
        Err(Error::Unsupported)
    ));
    let narrow = Nut::new(*f.nut.thread(), 3_f64.next_up()).unwrap();
    assert!(matches!(
        nut_land(&narrow, &|| false),
        Err(Error::Unsupported)
    ));
    let near = member_near_cap(
        &Face {
            outward: 1,
            ..face.clone()
        },
        &Placed::new(Material::Nut(&f.nut), identity()).unwrap(),
        &map,
        &|| false,
    )
    .unwrap();
    assert_eq!(near.center[2], Q::zero());
}

#[test]
fn source_policies_and_passage_errors_are_preserved_instead_of_publishing_partial_stack() {
    let f = fixture(false);
    let bad = Annulus::new(1.4, 3.5, 0.2).unwrap();
    assert!(matches!(
        query(
            &f,
            &[PlacedAnnulus {
                source: &bad,
                frame: f.ring_frame
            }],
            requirements(),
            &|| false
        ),
        Err(Error::Unsupported)
    ));
    assert!(matches!(
        query(
            &f,
            &[],
            StackRequirements {
                minimum_clamp_wall_mm: f64::NAN,
                ..requirements()
            },
            &|| false
        ),
        Err(Error::Invalid)
    ));
    let mut flipped = fixture(true);
    flipped.bore.centers_mm.swap(0, 1);
    flipped.bore.faces.swap(0, 1);
    let ring = PlacedAnnulus {
        source: &flipped.ring,
        frame: flipped.ring_frame,
    };
    assert_eq!(
        query(&flipped, &[ring], requirements(), &|| false)
            .unwrap()
            .bearings
            .len(),
        3
    );
}

#[test]
fn empty_and_multiple_washer_chains_and_proper_nonbinary_sources_retain_actual_phase() {
    for (left, nonbinary) in [(false, false), (true, true)] {
        let mut f = fixture(left);
        if nonbinary {
            f.hanger = Hanger::new(
                crate::authored_hanger_bolt::tests::fixture(left, true),
                None,
            )
            .unwrap();
        }
        nut_station(&mut f, -6.5 - 1e-10);
        assert_eq!(
            query(&f, &[], requirements(), &|| false)
                .unwrap()
                .bearings
                .len(),
            2
        );
        nut_station(&mut f, -6.7 - 3e-10);
        let pose = Frame {
            origin_mm: [11., -7., 5.],
            axes: [[0.36, 0.48, 0.8], [-0.8, 0.6, 0.], [-0.48, -0.64, 0.6]],
        };
        f.wood_frame = pose;
        f.clamp_frame = pose;
        f.hanger_frame = pose;
        f.nut_frame = transport(pose, f.nut_frame);
        let half = Annulus::new(1.7, 3.5, 0.1).unwrap();
        let rings = [-4.2 - 1e-10, -4.3 - 2e-10].map(|z| PlacedAnnulus {
            source: &half,
            frame: transport(
                pose,
                Frame {
                    origin_mm: [0., 0., z],
                    ..identity()
                },
            ),
        });
        let result = query(&f, &rings, requirements(), &|| false).unwrap();
        assert_eq!(result.bearings.len(), 4);
        assert_eq!(result.washer_passages.len(), 2);
    }
    let f = fixture(false);
    let ring = PlacedAnnulus {
        source: &f.ring,
        frame: f.ring_frame,
    };
    assert!(matches!(
        query(
            &f,
            &[ring],
            StackRequirements {
                minimum_clamp_wall_mm: f64::from_bits(1),
                ..requirements()
            },
            &|| false
        ),
        Err(Error::Unsupported)
    ));
}

#[test]
fn sealed_binary_decoding_retains_exact_source_values_and_world_overflow_stays_fallible() {
    for value in [
        0.,
        -0.,
        1.,
        -2.4,
        f64::MIN_POSITIVE,
        f64::from_bits(1),
        -f64::from_bits(1),
        f64::MAX,
        -f64::MAX,
    ] {
        assert_eq!(sealed_binary(value), q(value).unwrap());
    }
    let f = fixture(false);
    let stock = orthogonal::boolean(
        &[Box3 {
            min: [-10., -10., f64::MAX.next_down()],
            max: [10., 10., f64::MAX],
        }],
        &[],
        Operation::Union,
        0.,
    )
    .unwrap()
    .unwrap();
    assert!(stock.volume_mm3().is_finite());
    let s = 1.00000001;
    let strained = Frame {
        axes: [[s, 0., 0.], [0., s, 0.], [0., 0., s * s]],
        ..identity()
    };
    let a = Placed::new(Material::Cells(&f.clamp), strained).unwrap();
    let b = Placed::new(Material::Cells(&stock), strained).unwrap();
    let near = cap(&a, 2, [Q::zero(), Q::zero(), Q::zero()], 1, [2., 3.]);
    let far = cap(
        &b,
        2,
        [Q::zero(), Q::zero(), q(stock.bounds().min[2]).unwrap()],
        -1,
        [2., 3.],
    );
    let map = outside(&near, &b, &|| false).unwrap();
    // This finite source cap really separates, but its raw affine world
    // distance exceeds binary64 range. No Native bounds/MFG claim is made.
    let direct = passage::directed(&far.center[2]).unwrap() * b.metric.maximum_scale;
    assert!(!direct.valid());
    assert!(matches!(
        bearing(&near, &far, &b, &map, requirements(), &|| false),
        Err(Error::Unsupported)
    ));
}
#[test]
fn actual_related_material_stack_has_positive_lands_without_waiving_another_pair() {
    for left in [false, true] {
        let f = fixture(left);
        let rings = [PlacedAnnulus {
            source: &f.ring,
            frame: f.ring_frame,
        }];
        let result = certify_hanger_stack(
            StackSources {
                wood: Stock {
                    material: &f.wood,
                    pilot: f.pilot,
                    frame: identity(),
                },
                clamp: ThroughStock {
                    material: &f.clamp,
                    bore: f.bore,
                    frame: identity(),
                },
                hanger: &f.hanger,
                hanger_frame: identity(),
                washers: &rings,
                nut: &f.nut,
                nut_frame: f.nut_frame,
            },
            requirements(),
            &|| false,
        )
        .unwrap();
        assert_eq!(result.bearings.len(), 3);
        assert_eq!(result.washer_passages.len(), 1);
        assert!(
            result
                .bearings
                .iter()
                .all(|b| b.radial_land_lower_mm >= 1e-5 && b.correspondence_gap_upper_mm <= 1e-7)
        );
    }
}
