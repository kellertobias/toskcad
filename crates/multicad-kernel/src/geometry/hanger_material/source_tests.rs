use super::*;
use multicad_rust_kernel::{
    authored_flat_thread::{AuthoredThread, SourceOptions},
    authored_hanger_bolt::Spec,
    authored_metric_thread::{AuthoredMetricThread, Spec as MetricSpec},
    helical_thread::Spec as WoodSpec,
};

pub(super) fn source(left: bool, nonbinary: bool) -> Hanger {
    let major = (2.874 + 2.980) / 4.;
    let pitch_mid = (2.580 + 2.655) / 2.;
    Hanger::new(
        Spec {
            wood: AuthoredThread::new(
                WoodSpec {
                    major_radius_mm: 3.,
                    pitch_mm: 2.,
                    length_mm: if nonbinary { 12.3 } else { 12. },
                    depth_mm: 1.,
                    crest_width_mm: 0.4,
                    lower_flank_deg: 30.,
                    upper_flank_deg: 30.,
                    starts: 1,
                    left_hand: left,
                    root_radius_mm: None,
                },
                SourceOptions {
                    position_limit_mm: 1e-7,
                    normal_angle_limit_rad: 1e-7,
                },
            )
            .unwrap(),
            metric: AuthoredMetricThread::new(MetricSpec {
                major_radius_mm: major,
                pitch_mm: 0.5,
                crest_width_mm: 0.25 - (2. * major - pitch_mid) / 3_f64.sqrt(),
                length_mm: if nonbinary { 12.7 } else { 12. },
                starts: 1,
                left_hand: left,
            })
            .unwrap(),
            plain_length_mm: if nonbinary { 2.3 } else { 2. },
            plain_radius_mm: 1.5,
        },
        None,
    )
    .unwrap()
}
pub(super) fn volume_options() -> QueryOptions {
    QueryOptions {
        volume_error_mm3: 1e-7,
        maximum_cells: 65536,
    }
}

#[test]
fn actual_authored_hanger_sources_retain_world_bounds_volume_and_frame() {
    for (left, nonbinary, frame) in [
        (false, false, PlaneFrame::default()),
        (
            true,
            true,
            PlaneFrame {
                origin: [11., -7., 5.],
                u: [0.8, 0.6, 0.],
                v: [-0.36, 0.48, 0.8],
            },
        ),
        (
            false,
            true,
            PlaneFrame {
                origin: [10., -7., 3.],
                u: [0., 1., 0.],
                v: [0., 0., 1.],
            },
        ),
    ] {
        let source = source(left, nonbinary);
        let before = source;
        let body = HangerMaterial::new(source, frame, GeometryPolicy::default(), None).unwrap();
        let bounds = body.bounds().unwrap();
        let local = source.bounds_mm();
        let axes = [frame.u, frame.v, frame.normal()];
        for mask in 0..8 {
            let point = std::array::from_fn::<_, 3, _>(|i| local[i][(mask >> i) & 1]);
            for (i, origin) in frame.origin.into_iter().enumerate() {
                let world = (0..3).fold(origin, |sum, j| sum + axes[j][i] * point[j]);
                assert!(bounds.min[i] <= world && world <= bounds.max[i]);
            }
        }
        let volume = body.volume_bounds(volume_options(), None).unwrap();
        assert!(volume[0] > 0. && volume[1] >= volume[0] && volume[1] - volume[0] <= 1e-7);
        if frame.u == [0., 1., 0.] {
            let local = source.volume(volume_options(), None).unwrap().material_mm3;
            assert!(volume[0] <= local[0] && volume[1] >= local[1]);
        }
        assert_eq!(*body.source(), before);
        assert_eq!(body.frame(), frame);
    }
}

#[test]
fn real_source_queries_and_native_policy_fail_atomically() {
    let source = source(false, true);
    let before = source;
    let frame = PlaneFrame::default();
    let body = HangerMaterial::new(source, frame, GeometryPolicy::default(), None).unwrap();
    assert_eq!(
        HangerMaterial::new(source, frame, GeometryPolicy::default(), Some(&|| true))
            .unwrap_err()
            .kind,
        GeometryErrorKind::Cancelled
    );
    assert!(
        HangerMaterial::new(
            source,
            PlaneFrame {
                u: [1.01, 0., 0.],
                ..frame
            },
            GeometryPolicy::default(),
            None
        )
        .is_err()
    );
    assert!(
        HangerMaterial::new(
            source,
            PlaneFrame {
                origin: [1e7, 0., 0.],
                ..frame
            },
            GeometryPolicy::default(),
            None
        )
        .is_err()
    );
    assert_eq!(
        body.volume_bounds(
            QueryOptions {
                volume_error_mm3: 0.,
                ..volume_options()
            },
            None
        )
        .unwrap_err()
        .kind,
        GeometryErrorKind::Input
    );
    assert_eq!(
        body.volume_bounds(
            QueryOptions {
                maximum_cells: 1,
                ..volume_options()
            },
            None
        )
        .unwrap_err()
        .kind,
        GeometryErrorKind::Limit
    );
    assert_eq!(
        body.volume_bounds(volume_options(), Some(&|| true))
            .unwrap_err()
            .kind,
        GeometryErrorKind::Cancelled
    );
    assert_eq!(*body.source(), before);
}

#[test]
fn actual_positive_shoulders_below_native_feature_size_refuse_before_meshing() {
    let source = source(false, false);
    for radius in [
        source.metric().maximum_radius_mm().next_up(),
        source.wood().minimum_radius_mm().next_down(),
    ] {
        let mut spec = source.spec();
        spec.plain_radius_mm = radius;
        let thin = Hanger::new(spec, None).unwrap();
        assert!(
            HangerMaterial::new(thin, PlaneFrame::default(), GeometryPolicy::default(), None)
                .is_err()
        );
        assert_eq!(thin.spec(), spec);
    }
}

#[test]
fn true_affine_precision_floor_and_last_native_checkpoint_refuse_cleanly() {
    let source = source(true, true);
    let frame = PlaneFrame {
        origin: [11., -7., 5.],
        u: [0.8, 0.6, 0.],
        v: [-0.36, 0.48, 0.8],
    };
    let body = HangerMaterial::new(source, frame, GeometryPolicy::default(), None).unwrap();
    let local = source.volume(volume_options(), None).unwrap().material_mm3;
    let tight = QueryOptions {
        volume_error_mm3: (local[1] - local[0]).next_up(),
        ..volume_options()
    };
    assert!(source.volume(tight, None).is_ok());
    assert_eq!(
        body.volume_bounds(tight, None).unwrap_err().kind,
        GeometryErrorKind::Unsupported
    );
    let calls = std::cell::Cell::new(0usize);
    let count = || {
        calls.set(calls.get() + 1);
        false
    };
    body.volume_bounds(volume_options(), Some(&count)).unwrap();
    let last = calls.get();
    calls.set(0);
    let cancel = || {
        calls.set(calls.get() + 1);
        calls.get() == last
    };
    assert_eq!(
        body.volume_bounds(volume_options(), Some(&cancel))
            .unwrap_err()
            .kind,
        GeometryErrorKind::Cancelled
    );
    assert_eq!(calls.get(), last);
    assert_eq!(*body.source(), source);
}
