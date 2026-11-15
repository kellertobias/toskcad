use super::*;
use std::rc::Rc;
pub(super) fn source(both_ends: bool) -> Seats {
    use crate::{
        orthogonal::{self, Box3, Operation},
        orthogonal_bore_seats::{Options as SourceOptions, Seat, Seats},
        orthogonal_boundary, orthogonal_pockets,
    };
    let stock = Rc::new(
        orthogonal::boolean(
            &[
                Box3 {
                    min: [0.; 3],
                    max: [40., 40., 5.],
                },
                Box3 {
                    min: [0.; 3],
                    max: [40., 5., 40.],
                },
            ],
            &[],
            Operation::Union,
            1e-7,
        )
        .unwrap()
        .unwrap(),
    );
    let boundary = orthogonal_boundary::build(stock.quads()).unwrap();
    let bores = Rc::from(
        [
            [[20., 20., 5.1], [20., 20., -0.1]],
            [[20., 5.1, 20.], [20., -0.1, 20.]],
        ]
        .map(|[a, b]| {
            orthogonal_pockets::identify_through_bore(&stock, &boundary, a, b, 4.9, 1e-7, 1e-5)
                .unwrap()
        }),
    );
    let mut seats = [0, 1]
        .map(|bore| Seat {
            bore,
            end: false,
            radius_mm: 8.02,
            depth_mm: 3.12,
        })
        .to_vec();
    if both_ends {
        seats.push(Seat {
            bore: 0,
            end: true,
            radius_mm: 5.8,
            depth_mm: 0.9,
        });
    }
    Seats::new(
        stock,
        bores,
        &seats,
        SourceOptions {
            clearance_mm: 1e-7,
            minimum_feature_mm: 1e-5,
            maximum_work: 16_000_000,
            maximum_bytes: 64 * 1024 * 1024,
        },
        &|| false,
    )
    .unwrap()
}
#[test]
fn original_perpendicular_seats_emit_candidate_with_shared_closed_rims() {
    let source = source(false);
    let placement =
        StockClipPlacement::new([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]], [0.; 3]).unwrap();
    let candidate = candidate(&source, placement, Options::default(), &|| false).unwrap();
    assert!(std::ptr::eq(candidate.source(), &source));
    assert_eq!(
        candidate.face_entities().len(),
        source.boundary().faces.len() + 4
    );
    assert_eq!(candidate.as_str().matches("CONICAL_SURFACE(").count(), 2);
    let proof = crate::step_read::read_analytic_solid(candidate.as_str(), 1e-7).unwrap();
    let volume = source.volume_bounds_mm3();
    assert!((proof.volume_mm3 - volume[0].midpoint(volume[1])).abs() < 1e-7);
}

#[test]
fn actual_placement_and_both_finite_seats_keep_shared_shell() {
    let source = source(true);
    let placement = StockClipPlacement::new(
        [[0.8, -0.6, 0.], [0.6, 0.8, 0.], [0., 0., 1.]],
        [11., -7., 5.],
    )
    .unwrap();
    let candidate = candidate(&source, placement, Options::default(), &|| false).unwrap();
    assert_eq!(
        candidate.face_entities().len(),
        source.boundary().faces.len() + 5
    );
    assert_eq!(candidate.as_str().matches("CONICAL_SURFACE(").count(), 3);
    let proof = crate::step_read::read_analytic_solid(candidate.as_str(), 1e-7).unwrap();
    assert!((proof.volume_mm3 - source.volume_bounds_mm3()[0]).abs() < 1e-7);
}
#[test]
fn candidate_work_storage_and_actual_cancellation_refuse_atomically() {
    use std::cell::Cell;
    let source = source(false);
    let placement =
        StockClipPlacement::new([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]], [0.; 3]).unwrap();
    let earned = candidate(&source, placement, Options::default(), &|| false).unwrap();
    let exact = Options {
        maximum_work: earned.encoding_work(),
        maximum_bytes: earned.accounted_peak_bytes(),
    };
    assert!(candidate(&source, placement, exact, &|| false).is_ok());
    assert!(matches!(
        candidate(
            &source,
            placement,
            Options {
                maximum_work: exact.maximum_work - 1,
                ..exact
            },
            &|| false
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        candidate(
            &source,
            placement,
            Options {
                maximum_bytes: 1,
                ..exact
            },
            &|| false
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        candidate(
            &source,
            placement,
            Options {
                maximum_work: 0,
                ..exact
            },
            &|| false
        ),
        Err(Error::Invalid)
    ));
    let calls = Cell::new(0usize);
    let stop = || {
        calls.set(calls.get() + 1);
        calls.get() > 8
    };
    assert!(matches!(
        candidate(&source, placement, Options::default(), &stop),
        Err(Error::Cancelled)
    ));
    assert!(calls.get() > 8);
    assert_eq!(
        candidate(&source, placement, Options::default(), &|| false)
            .unwrap()
            .as_str(),
        earned.as_str()
    );
}
