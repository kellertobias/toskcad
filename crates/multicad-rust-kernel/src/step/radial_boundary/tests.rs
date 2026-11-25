use super::*;
#[test]
fn actual_retained_writer_append_preserves_source_and_refuses_bad_precision_transactionally() {
    let body = crate::radial_tensor::BoundaryOwnership::new(
        crate::radial_tensor::tests::cam::wedge(&[1, 2]),
        1e-7,
        None,
    )
    .unwrap()
    .bind_surfaces(1e-7, None)
    .unwrap()
    .certify_axial(1e-7, 1e-6, None)
    .unwrap()
    .certify_caps(None)
    .unwrap()
    .certify_ordering(None)
    .unwrap()
    .certify_material(1e-7, None)
    .unwrap();
    let mut writer = Writer::new();
    for tolerance in [0., f64::NAN, -1., f64::INFINITY] {
        assert_eq!(
            write_radial_boundary(&body, tolerance),
            Err(StepError::InvalidGeometry)
        );
        assert_eq!(
            writer.append_radial_boundary(&body, tolerance),
            Err(StepError::InvalidGeometry)
        );
    }
    writer.append_radial_boundary(&body, 1e-7).unwrap();
    let text = writer.finish().unwrap();
    verify_candidate(&text, &body, 1e-7).unwrap();
    assert_eq!(text, write_radial_boundary(&body, 1e-7).unwrap());
}
#[test]
fn retained_export_requires_the_complete_unchanged_finite_source() {
    let body = crate::radial_tensor::tests::material();
    let text = write_radial_boundary(&body, 1e-7).unwrap();
    verify_candidate(&text, &body, 1e-7).unwrap();
    assert_eq!(
        verify_candidate("not STEP", &body, 1e-7),
        Err(StepError::Unsupported)
    );
    let oversized = " ".repeat(crate::step_read::MAX_PLANAR_STEP_BYTES + 1);
    assert_eq!(
        verify_candidate(&oversized, &body, 1e-7),
        Err(StepError::Limit)
    );
    let tensor = body
        .input()
        .surfaces
        .iter()
        .find_map(|s| {
            if let Surface::Tensor(t) = s {
                Some(t)
            } else {
                None
            }
        })
        .unwrap();
    let p = tensor.poles()[tensor.dimensions()[1] + 1];
    let old = Coordinates(p).to_string();
    let mut p = p;
    p[0] += 0.001;
    let modified = text.replacen(&old, &Coordinates(p).to_string(), 1);
    assert_ne!(modified, text);
    // Abstract ownership still holds; the independent finite-source proof
    // must detect the changed interior pole, not merely compare metrics.
    crate::step_read::read_radial_boundary(&modified, 1e-7).unwrap();
    assert_eq!(
        verify_candidate(&modified, &body, 1e-7),
        Err(StepError::Unsupported)
    );
    verify_candidate(&text, &body, 1e-7).unwrap();
}

/// Internal 1/8-27 NPT cavity, mouth at z=0 and negative taper.
pub(crate) fn internal_npt_cavity() -> RadialBoundary {
    use crate::axial_pipe_thread::authored::{AuthoredSection, AuthoredSpec, AuthoredThread};
    let p = 25.4 / 27.;
    let k = -1. / 32.;
    let m = 3_f64.sqrt();
    let slopes = 1. / (m - k) + 1. / (m + k);
    let depth = 0.8 * p / slopes;
    let spec = AuthoredSpec {
        major_radius_mm: 5.,
        pitch_mm: p,
        length_mm: 8.,
        radial_taper: k,
        section: AuthoredSection::NptFlat {
            depth_mm: depth,
            flat_width_mm: (p - depth * slopes) / 2.,
        },
        crest_radius_limit_mm: None,
    };
    AuthoredThread::new(spec, 1e-7, 1e-6)
        .unwrap()
        .raw_boundary_ownership(1e-7, None)
        .unwrap()
        .ownership
        .bind_surfaces(1e-7, None)
        .unwrap()
        .certify_axial(1e-7, 1e-6, None)
        .unwrap()
        .certify_caps(None)
        .unwrap()
        .certify_ordering(None)
        .unwrap()
        .certify_material(1e-7, None)
        .unwrap()
}
#[test]
fn pipe_sleeve_step_keeps_the_exact_cavity_and_stock_and_refuses_bad_stock() {
    let cavity = internal_npt_cavity();
    let stock = SleeveStock {
        origin: [0.; 3],
        normal: [0., 0., 1.],
        reference: [1., 0., 0.],
        radius_mm: 7.,
        height_mm: 8.,
    };
    let text = write_pipe_sleeve(&cavity, stock, 1e-7).unwrap();
    let imported = crate::step_read::read_pipe_sleeve(&text, 1e-7).unwrap();
    assert_eq!(imported.stock, stock);
    assert_eq!(
        canonical(imported.cavity.input()),
        canonical(cavity.input())
    );
    // The imported cavity re-earns the full radial material certificate.
    let restored = imported
        .cavity
        .bind_surfaces(1e-7, None)
        .unwrap()
        .certify_axial(1e-7, 1e-6, None)
        .unwrap()
        .certify_caps(None)
        .unwrap()
        .certify_ordering(None)
        .unwrap()
        .certify_material(1e-7, None)
        .unwrap();
    assert_eq!(
        restored.volume_with_interrupt(1e-4, None).unwrap(),
        cavity.volume_with_interrupt(1e-4, None).unwrap()
    );
    // One cylinder, two annuli; the unchanged radial reader does not own it.
    assert_eq!(text.matches("CYLINDRICAL_SURFACE").count(), 1);
    assert_eq!(text.matches("FACE_BOUND(").count(), 2);
    assert!(crate::step_read::read_radial_boundary(&text, 1e-7).is_err());
    // Caps must be the cavity's own: a shorter stock is refused, as are
    // non-finite or degenerate stock values.
    for bad in [
        SleeveStock {
            height_mm: 7.,
            ..stock
        },
        SleeveStock {
            radius_mm: 0.,
            ..stock
        },
        SleeveStock {
            origin: [f64::NAN, 0., 0.],
            ..stock
        },
    ] {
        assert!(write_pipe_sleeve(&cavity, bad, 1e-7).is_err());
    }
    let mut writer = Writer::new();
    assert!(
        writer
            .append_pipe_sleeve(
                &cavity,
                SleeveStock {
                    height_mm: 7.,
                    ..stock
                },
                1e-7
            )
            .is_err()
    );
    writer.append_pipe_sleeve(&cavity, stock, 1e-7).unwrap();
    assert_eq!(writer.finish().unwrap(), text);
}

fn pipe_boundary(rounded: bool, external: bool) -> RadialBoundary {
    use crate::axial_pipe_thread::authored::{AuthoredSection, AuthoredSpec, AuthoredThread};
    AuthoredThread::new(
        AuthoredSpec {
            major_radius_mm: 5.,
            pitch_mm: 1.,
            length_mm: 6.,
            radial_taper: if external { 1. / 32. } else { -1. / 32. },
            section: if rounded {
                AuthoredSection::BsptRounded { depth_mm: 0.6 }
            } else {
                AuthoredSection::NptFlat {
                    depth_mm: 0.6,
                    flat_width_mm: 0.2,
                }
            },
            crest_radius_limit_mm: None,
        },
        1e-7,
        1e-6,
    )
    .unwrap()
    .raw_boundary_ownership(1e-7, None)
    .unwrap()
    .ownership
    .bind_surfaces(1e-7, None)
    .unwrap()
    .certify_axial(1e-7, 1e-6, None)
    .unwrap()
    .certify_caps(None)
    .unwrap()
    .certify_ordering(None)
    .unwrap()
    .certify_material(1e-7, None)
    .unwrap()
}
#[test]
fn pipe_hex_hardware_keeps_exact_thread_tables_and_stock_for_nuts_and_heads() {
    use crate::step_read::{HexStock, HexStockKind, read_pipe_hex_hardware};
    let tilt = 0.3_f64;
    let (s, c) = (tilt.sin(), tilt.cos());
    for rounded in [false, true] {
        for kind in [HexStockKind::Nut, HexStockKind::Head] {
            let local = pipe_boundary(rounded, kind == HexStockKind::Head);
            for placed in [false, true] {
                let (boundary, origin, normal, reference) = if placed {
                    // A rigid occurrence: the boundary and stock move together.
                    let u = [c, s, 0.];
                    let v = [0., 0., 1.];
                    let o = [3., -2., 7.];
                    let n = [
                        u[1] * v[2] - u[2] * v[1],
                        u[2] * v[0] - u[0] * v[2],
                        u[0] * v[1] - u[1] * v[0],
                    ];
                    (
                        local
                            .transformed_with_interrupt([o, u, v], 1e-7, 1e-6, None)
                            .unwrap(),
                        o,
                        n,
                        u,
                    )
                } else {
                    (local.clone(), [0.; 3], [0., 0., 1.], [1., 0., 0.])
                };
                let stock = HexStock {
                    kind,
                    origin,
                    normal,
                    reference,
                    across_flats_mm: 16.,
                    span_mm: if kind == HexStockKind::Nut {
                        [0., 6.]
                    } else {
                        [6., 10.]
                    },
                };
                let text = write_pipe_hex_hardware(&boundary, stock, 1e-7)
                    .unwrap_or_else(|e| panic!("{rounded} {kind:?} {placed}: {e:?}"));
                verify_pipe_hex_hardware(&text, &boundary, stock, 1e-7).unwrap();
                let imported = read_pipe_hex_hardware(&text, 1e-7).unwrap();
                assert_eq!(imported.kind, kind);
                assert_eq!(imported.corners.len(), 12);
                // The decoded radial tables re-earn the full material proof.
                let restored = imported
                    .boundary
                    .bind_surfaces(1e-7, None)
                    .unwrap()
                    .certify_axial(1e-7, 1e-6, None)
                    .unwrap()
                    .certify_caps(None)
                    .unwrap()
                    .certify_ordering(None)
                    .unwrap()
                    .certify_material(1e-7, None)
                    .unwrap();
                assert_eq!(
                    restored.volume_with_interrupt(1e-4, None).unwrap(),
                    boundary.volume_with_interrupt(1e-4, None).unwrap()
                );
                let annuli = if kind == HexStockKind::Nut { 2 } else { 1 };
                assert_eq!(text.matches("FACE_BOUND(").count(), annuli);
                assert_eq!(text.matches("LINE(").count(), 18);
                assert!(crate::step_read::read_radial_boundary(&text, 1e-7).is_err());
                assert!(crate::step_read::read_pipe_sleeve(&text, 1e-7).is_err());
                // Readback refuses another stock, the other arrangement and
                // a different boundary.
                for other in [
                    HexStock {
                        across_flats_mm: 16.000001,
                        ..stock
                    },
                    HexStock {
                        origin: [origin[0] + 1e-9, origin[1], origin[2]],
                        ..stock
                    },
                ] {
                    assert_eq!(
                        verify_pipe_hex_hardware(&text, &boundary, other, 1e-7),
                        Err(StepError::Unsupported)
                    );
                }
                let other = pipe_boundary(!rounded, kind == HexStockKind::Head);
                assert!(verify_pipe_hex_hardware(&text, &other, stock, 1e-7).is_err());
                // Bad stocks: wrong span ends, degenerate width, non-finite,
                // and the wrong arrangement for this boundary.
                let wrong = if kind == HexStockKind::Nut {
                    HexStockKind::Head
                } else {
                    HexStockKind::Nut
                };
                for bad in [
                    HexStock {
                        span_mm: [0., 5.],
                        ..stock
                    },
                    HexStock {
                        across_flats_mm: 0.,
                        ..stock
                    },
                    HexStock {
                        origin: [f64::NAN, 0., 0.],
                        ..stock
                    },
                    HexStock {
                        normal: [0., 0., 2.],
                        ..stock
                    },
                    HexStock {
                        kind: wrong,
                        ..stock
                    },
                ] {
                    assert!(write_pipe_hex_hardware(&boundary, bad, 1e-7).is_err());
                }
                // Appends are transactional.
                let mut writer = Writer::new();
                assert!(
                    writer
                        .append_pipe_hex_hardware(
                            &boundary,
                            HexStock {
                                span_mm: [0., 5.],
                                ..stock
                            },
                            1e-7
                        )
                        .is_err()
                );
                writer
                    .append_pipe_hex_hardware(&boundary, stock, 1e-7)
                    .unwrap();
                assert_eq!(writer.finish().unwrap(), text);
            }
        }
    }
}
