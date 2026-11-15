use super::*;
mod fixture;
fn with_head(f: &mut dyn FnMut(&ClosedRepairedHead<'_, '_, '_, '_, '_>)) {
    let text = fixture::actual_text();
    super::super::super::common_mouth::tests::derived_thread_probe::with_head(&text, f);
}
#[test]
fn actual_m2_degree11_packet_retains_owners_complete_nets_and_measured_transfer() {
    with_head(&mut |head| {
        let candidate = head
            .degree11_thread_carrier(Default::default(), None)
            .unwrap();
        assert!(std::ptr::eq(candidate.body(), head.body()));
        assert_eq!(candidate.chart(), head.chart());
        assert_eq!(candidate.faces().len(), 4);
        assert_eq!(candidate.inherited_head_work_used(), head.work_used());
        assert_eq!(
            candidate.total_work_used(),
            candidate.new_work_used() + head.work_used()
        );
        assert!(candidate.new_work_used() < WORK);
        assert!(candidate.total_work_used() > WORK);
        assert!(candidate.retained_budget_bytes() < BYTES);
        assert!(candidate.cap_stations_mm()[0] < candidate.cap_stations_mm()[1]);
        let mut total = 0;
        for (i, face) in candidate.faces().iter().enumerate() {
            let next = &candidate.faces()[(i + 1) % 4];
            let transfer = face.transfer();
            assert_eq!(transfer.source_face_index, 2 + i);
            assert_eq!(
                transfer.original_rail_entities[1],
                next.transfer().original_rail_entities[0]
            );
            assert_eq!(
                face.original_control_entities().len(),
                face.source().poles().len()
            );
            assert_eq!(transfer.original_cells, transfer.source_cells);
            assert!(transfer.distance_mm < 1e-7);
            assert!(transfer.derivative_mm.into_iter().all(|x| x < 1e-7));
            let [width, height] = face.source().dimensions();
            assert_eq!(
                face.source().poles()[(width - 1) * height..],
                next.source().poles()[..height]
            );
            for row in face.source().poles().chunks_exact(height) {
                assert_eq!([row[0][2], row[height - 1][2]], candidate.cap_stations_mm());
            }
            let cell = face.source().cell([0, 0]).unwrap();
            assert!(cell.derivative_bounds().is_ok());
            let _ = face.forward();
            total += face.source().poles().len();
        }
        assert_eq!(total, 41_232);
        eprintln!(
            "ACTUAL D11 newwork{} inherited{} retained{}",
            candidate.new_work_used(),
            candidate.inherited_head_work_used(),
            candidate.retained_budget_bytes()
        );
    });
}
#[test]
fn actual_packet_limits_precision_and_cancellation_leave_body_unchanged() {
    with_head(&mut |head| {
        let options = Degree11ThreadOptions::default();
        let baseline = head.degree11_thread_carrier(options, None).unwrap();
        let expected: Vec<_> = baseline
            .faces()
            .iter()
            .map(|f| f.source().clone())
            .collect();
        let used = baseline.new_work_used();
        let controls = expected.iter().map(|n| n.poles().len()).sum();
        let check = |o| {
            assert!(head.degree11_thread_carrier(o, None).is_err());
        };
        check(Degree11ThreadOptions {
            maximum_distance_mm: 0.,
            ..options
        });
        check(Degree11ThreadOptions {
            maximum_cell_derivative_mm: f64::NAN,
            ..options
        });
        check(Degree11ThreadOptions {
            maximum_controls: MAX_POINTS + 1,
            ..options
        });
        check(Degree11ThreadOptions {
            maximum_work: WORK + 1,
            ..options
        });
        check(Degree11ThreadOptions {
            maximum_bytes: BYTES + 1,
            ..options
        });
        for o in [
            Degree11ThreadOptions {
                maximum_work: 0,
                ..options
            },
            Degree11ThreadOptions {
                maximum_work: used - 1,
                ..options
            },
            Degree11ThreadOptions {
                maximum_bytes: 0,
                ..options
            },
            Degree11ThreadOptions {
                maximum_bytes: head.retained_budget_bytes() + SCRATCH,
                ..options
            },
            Degree11ThreadOptions {
                maximum_bytes: head.retained_budget_bytes() + SCRATCH + 2 * 1024 * 1024,
                ..options
            },
            Degree11ThreadOptions {
                maximum_controls: controls - 1,
                ..options
            },
        ] {
            assert!(matches!(
                head.degree11_thread_carrier(o, None),
                Err(PhillipsReadError::Limit)
            ));
        }
        let (mut low, mut high) = (0, BYTES);
        while low < high {
            let middle = low + (high - low) / 2;
            match head.degree11_thread_carrier(
                Degree11ThreadOptions {
                    maximum_bytes: middle,
                    ..options
                },
                None,
            ) {
                Ok(_) => high = middle,
                Err(PhillipsReadError::Limit) => low = middle + 1,
                Err(_) => panic!("unexpected byte frontier error"),
            }
        }
        assert!(
            head.degree11_thread_carrier(
                Degree11ThreadOptions {
                    maximum_bytes: high,
                    ..options
                },
                None
            )
            .is_ok()
        );
        assert!(matches!(
            head.degree11_thread_carrier(
                Degree11ThreadOptions {
                    maximum_bytes: high - 1,
                    ..options
                },
                None
            ),
            Err(PhillipsReadError::Limit)
        ));
        for o in [
            Degree11ThreadOptions {
                maximum_distance_mm: 1e-18,
                ..options
            },
            Degree11ThreadOptions {
                maximum_cell_derivative_mm: 1e-18,
                ..options
            },
        ] {
            assert!(matches!(
                head.degree11_thread_carrier(o, None),
                Err(PhillipsReadError::Unsupported)
            ));
        }
        assert!(matches!(
            head.degree11_thread_carrier(options, Some(&|| true)),
            Err(PhillipsReadError::Cancelled)
        ));
        let calls = Cell::new(0usize);
        let interrupt = || {
            calls.set(calls.get() + 1);
            false
        };
        head.degree11_thread_carrier(options, Some(&interrupt))
            .unwrap();
        for checkpoint in [
            1,
            2,
            3,
            100,
            1000,
            calls.get() / 2,
            calls.get() - 1,
            calls.get(),
        ] {
            let n = Cell::new(0usize);
            let cancel = || {
                n.set(n.get() + 1);
                n.get() == checkpoint
            };
            assert!(
                matches!(
                    head.degree11_thread_carrier(options, Some(&cancel)),
                    Err(PhillipsReadError::Cancelled)
                ),
                "checkpoint{checkpoint}"
            );
        }
        let exact = head
            .degree11_thread_carrier(
                Degree11ThreadOptions {
                    maximum_work: used,
                    maximum_controls: controls,
                    ..options
                },
                None,
            )
            .unwrap();
        assert_eq!(exact.new_work_used(), used);
        for (face, expected) in exact.faces().iter().zip(&expected) {
            assert_eq!(face.source(), expected);
        }
    });
}

#[test]
fn independent_scaled_and_proper_pose_packets_retain_full_degree11_cells() {
    for scale in [1., 1.137] {
        for posed in [false, true] {
            let text =
                super::super::super::super::branches::representative::scaled_text(scale, posed);
            super::super::super::common_mouth::tests::derived_thread_probe::with_head(
                &text,
                &mut |head| {
                    let candidate = head
                        .degree11_thread_carrier(Default::default(), None)
                        .unwrap();
                    assert!(std::ptr::eq(candidate.body(), head.body()));
                    assert_eq!(candidate.chart(), head.chart());
                    assert_eq!(candidate.faces().len(), 4);
                    for f in candidate.faces() {
                        let t = f.transfer();
                        assert!(
                            t.distance_mm <= 1e-7 && t.derivative_mm.into_iter().all(|d| d <= 1e-7)
                        );
                        assert_eq!(t.source_cells, f.source().spans());
                    }
                },
            );
        }
    }
}
