use super::*;
use crate::axial_pipe_thread::authored::{AuthoredSection, AuthoredSpec, AuthoredThread};
use crate::step::npt_stock::{CurveOptions, CurvePlan, SurfacePlan};
pub(super) fn with_source<R>(quarter: bool, f: impl FnOnce(&SurfacePlan<'_>) -> R) -> R {
    let (tpi, e1, l1, l2, od, max) = if quarter {
        (18., 0.49163, 0.2278, 0.4018, 0.540, 0.088)
    } else {
        (27., 0.37360, 0.1615, 0.2639, 0.405, 0.096)
    };
    let pitch = 25.4 / tpi;
    let k = 1. / 32.;
    let m = 3f64.sqrt();
    let slopes = 1. / (m - k) + 1. / (m + k);
    let depth = pitch / slopes - (0.033 + max) * pitch;
    let clip = AuthoredThread::new(
        AuthoredSpec {
            major_radius_mm: (e1 * 25.4 - 2. * k * l1 * 25.4 + depth) / 2.,
            pitch_mm: pitch,
            length_mm: l2 * 25.4,
            radial_taper: k,
            section: AuthoredSection::NptFlat {
                depth_mm: depth,
                flat_width_mm: (pitch - depth * slopes) / 2.,
            },
            crest_radius_limit_mm: Some(od * 25.4 / 2.),
        },
        1e-7,
        1e-6,
    )
    .unwrap()
    .stock_clip_plan(Default::default(), None)
    .unwrap();
    let curves = CurvePlan::new(&clip, CurveOptions::default(), None).unwrap();
    let source = SurfacePlan::new(curves, CurveOptions::default().into(), None).unwrap();
    f(&source)
}
fn candidate(quarter: bool, across: f64) -> String {
    let input = super::super::tests::fixture(quarter);
    let doc = CompoundDocument::parse(&input, DocumentOptions::default(), None).unwrap();
    let parsed = Parsed::parse(
        &doc,
        doc.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )
    .unwrap();
    let cap = parsed
        .faces
        .iter()
        .filter_map(|f| match f.support {
            Support::Plane { height } => Some((parsed.points[height].coordinates[2].value(), f.id)),
            _ => None,
        })
        .max_by(|a, b| a.0.total_cmp(&b.0))
        .unwrap()
        .1;
    let head = HeadInput::regular_hex(&parsed, cap, across, 3.).unwrap();
    let proof = head::certify(&parsed, head, 1e-7, 1e-6, Default::default(), None).unwrap();
    crate::step::npt_stock::full_head::propose(&proof, Default::default(), None)
        .unwrap()
        .as_str()
        .to_owned()
}
fn material(input: &str) -> Result<()> {
    let packet = DecodedFile::parse(input, Default::default(), None)?;
    let doc = CompoundDocument::parse(
        packet.derived_auxiliary_text(),
        DocumentOptions::default(),
        None,
    )
    .map_err(neutral_error)?;
    let parsed = Parsed::parse(
        &doc,
        doc.represented_members()[0].item_entity(),
        Default::default(),
        None,
    )?;
    packet.certify(&parsed, 1e-7, 1e-6, None)?;
    Ok(())
}
fn change(input: &str, id: u32, body: &str) -> String {
    let prefix = format!("#{id}=");
    let start = input.find(&prefix).unwrap();
    let end = start + input[start..].find(';').unwrap() + 1;
    format!("{}{}{};{}", &input[..start], prefix, body, &input[end..])
}
#[test]
fn actual_full_clipped_source_matches_original_law_and_directed_head_author_corners() {
    for quarter in [false, true] {
        let actual = candidate(quarter, 18.);
        with_source(quarter, |source| {
            let q = qualify_actual_file(
                &actual,
                source,
                HexHeadIntent {
                    across_flats_mm: 18.,
                    height_mm: 3.,
                },
                Default::default(),
                None,
            )
            .unwrap();
            assert_eq!(q.as_str(), actual);
            assert!(std::ptr::eq(q.source_plan(), source));
            let r = q.receipt();
            assert!(r.original_thread.position_error_mm() <= 1e-7);
            assert!(
                r.original_thread
                    .original_derivative_error_mm()
                    .iter()
                    .all(|x| *x <= 1e-6)
            );
            assert!(r.head_author_corner_error_mm > 0. && r.head_author_corner_error_mm < 1e-7);
            assert!(r.sewn_reference_position_error_mm <= 1e-7);
            let owned = OwnedActualHead::read(
                std::sync::Arc::from(actual.as_str()),
                Default::default(),
                None,
            )
            .unwrap();
            assert_eq!(
                r.thread_axial_transfer_displacement_mm,
                owned.axial_sewing_transfer().domain_displacement_mm()
            );
            assert!(
                r.sewn_reference_position_error_mm
                    >= r.original_thread.position_error_mm()
                        + r.thread_axial_transfer_displacement_mm
            );
            assert!(r.aggregate_work <= 272_000_000 && r.accounted_peak_bytes <= 64 * 1024 * 1024);
            assert_eq!(
                r.aggregate_work,
                r.material.aggregate_work
                    + r.head_matching_work
                    + r.original_thread.aggregate_work_used()
            );
            eprintln!(
                "ACTUAL_CLIPPED_SOURCE quarter={quarter} {r:?} sourcework={} sourcebytes={}",
                source.work_used(),
                source.retained_budget_bytes()
            );
        });
    }
}
#[test]
fn actual_full_clipped_source_rejects_wrong_head_and_real_material_deformation() {
    use super::super::records::{list, reference};
    use crate::step::pipe_hardware::numbers;
    let actual = candidate(false, 18.);
    let wrong_head = candidate(false, 20.);
    assert!(material(&wrong_head).is_ok());
    let doc = CompoundDocument::parse(&actual, DocumentOptions::default(), None).unwrap();
    let net = doc
        .records()
        .values()
        .find_map(|p| {
            matches!(p.as_slice(),[p] if p.name=="B_SPLINE_SURFACE_WITH_KNOTS").then_some(&p[0])
        })
        .unwrap();
    let rows = list(&net.args[3]).unwrap();
    let row = list(&rows[rows.len() / 2]).unwrap();
    let id = reference(&row[row.len() / 2]).unwrap();
    let e = &doc.records()[&id][0];
    let tokens = crate::step_read::part21::numeric_spellings(e).unwrap();
    let xyz = list(&e.args[1]).unwrap();
    let Value::Number(x) = xyz[0] else { panic!() };
    let altered = change(
        &actual,
        id,
        &format!(
            "CARTESIAN_POINT('',({},{},{}))",
            numbers::literal(x + 2e-6).unwrap(),
            tokens[1],
            tokens[2]
        ),
    );
    drop(doc);
    assert!(
        material(&altered).is_ok(),
        "independent actual material must survive small genuine interior-net deformation"
    );
    with_source(false, |source| {
        let intent = HexHeadIntent {
            across_flats_mm: 18.,
            height_mm: 3.,
        };
        assert!(matches!(
            qualify_actual_file(&wrong_head, source, intent, Default::default(), None),
            Err(Error::Unsupported)
        ));
        assert!(matches!(
            qualify_actual_file(&altered, source, intent, Default::default(), None),
            Err(Error::Unsupported)
        ));
        assert!(qualify_actual_file(&actual, source, intent, Default::default(), None).is_ok());
    });
}

#[test]
fn actual_full_clipped_source_resource_frontiers_and_cancellation_are_atomic() {
    use std::cell::Cell;
    let actual = candidate(false, 18.);
    let unchanged = actual.clone();
    with_source(false, |source| {
        let intent = HexHeadIntent {
            across_flats_mm: 18.,
            height_mm: 3.,
        };
        let q = qualify_actual_file(&actual, source, intent, Default::default(), None).unwrap();
        let good = q.receipt();
        for maximum_work in [good.aggregate_work - 1, good.aggregate_work] {
            let options = SourceOptions {
                resources: Options {
                    maximum_work,
                    ..Default::default()
                },
                ..Default::default()
            };
            let result = qualify_actual_file(&actual, source, intent, options, None);
            if maximum_work < good.aggregate_work {
                assert!(matches!(result, Err(Error::Limit)));
            } else {
                assert_eq!(
                    result.unwrap().receipt().aggregate_work,
                    good.aggregate_work
                );
            }
        }
        for maximum_bytes in [good.accounted_peak_bytes - 1, good.accounted_peak_bytes] {
            let options = SourceOptions {
                resources: Options {
                    maximum_bytes,
                    ..Default::default()
                },
                ..Default::default()
            };
            let result = qualify_actual_file(&actual, source, intent, options, None);
            if maximum_bytes < good.accounted_peak_bytes {
                assert!(matches!(result, Err(Error::Limit)));
            } else {
                assert_eq!(
                    result.unwrap().receipt().accounted_peak_bytes,
                    good.accounted_peak_bytes
                );
            }
        }
        let calls = Cell::new(0usize);
        let count = || {
            calls.set(calls.get() + 1);
            false
        };
        qualify_actual_file(&actual, source, intent, Default::default(), Some(&count)).unwrap();
        let total = calls.get();
        for at in [1, total / 2, total] {
            let calls = Cell::new(0usize);
            let cancel = || {
                calls.set(calls.get() + 1);
                calls.get() >= at
            };
            assert!(matches!(
                qualify_actual_file(&actual, source, intent, Default::default(), Some(&cancel)),
                Err(Error::Cancelled)
            ));
            assert_eq!(
                qualify_actual_file(&actual, source, intent, Default::default(), None)
                    .unwrap()
                    .receipt()
                    .aggregate_work,
                good.aggregate_work
            );
        }
        let extra = qualify_actual_file(
            &actual,
            source,
            intent,
            SourceOptions {
                resources: Options {
                    additional_resident_bytes: 4096,
                    ..Default::default()
                },
                ..Default::default()
            },
            None,
        )
        .unwrap()
        .receipt();
        assert_eq!(extra.aggregate_work, good.aggregate_work);
        assert_eq!(extra.accounted_peak_bytes, good.accounted_peak_bytes + 4096);
        for options in [
            SourceOptions {
                position_error_mm: f64::NAN,
                ..Default::default()
            },
            SourceOptions {
                material_tolerance_mm: 0.,
                ..Default::default()
            },
            SourceOptions {
                minimum_wall_mm: 0.,
                ..Default::default()
            },
        ] {
            assert!(matches!(
                qualify_actual_file(&actual, source, intent, options, None),
                Err(Error::Invalid)
            ));
        }
        assert_eq!(actual, unchanged);
    });
}
