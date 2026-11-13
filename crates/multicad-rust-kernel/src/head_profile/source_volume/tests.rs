use super::*;
fn quadrature(profile: &Profile) -> f64 {
    profile
        .segments()
        .map(|s| {
            let value = |t: f64| {
                let [r, _] = s.point(t);
                let dz = match s.curve {
                    Curve::Line => s.end[1] - s.start[1],
                    Curve::Arc {
                        center,
                        radius_mm,
                        sweep,
                    } => {
                        let a = (s.start[1] - center[1]).atan2(s.start[0] - center[0]);
                        radius_mm * (a + sweep * t).cos() * sweep
                    }
                };
                PI * r * r * dz
            };
            let n = 4096;
            let mut sum = value(0.) + value(1.);
            for i in 1..n {
                sum += if i % 2 == 0 { 2. } else { 4. } * value(i as f64 / n as f64);
            }
            sum / (3. * n as f64)
        })
        .sum()
}
#[test]
fn actual_crown_neck_cone_and_cylinder_support_volumes_enclose_independent_quadrature() {
    for head in [
        Head::Cylindrical {
            radius_mm: 3.,
            height_mm: 2.4,
        },
        Head::Crowned {
            radius_mm: 2.8,
            height_mm: 2.4,
            crown_radius_mm: 3.,
            top_land_radius_mm: 2.1,
            neck_radius_mm: 0.1,
        },
        Head::Countersunk {
            radius_mm: 2.75,
            height_mm: 1.65,
            included_angle_deg: 90.,
            edge_chamfer_mm: 0.05,
            neck_radius_mm: 0.1,
        },
    ] {
        for plain in [0., 0.2, 1.] {
            let p = Profile::new_for_thread_union(head, 1.5, plain).unwrap();
            let interval = volume_bounds(&p).unwrap();
            let numeric = quadrature(&p);
            assert!(
                interval.lo <= numeric && numeric <= interval.hi,
                "{head:?}/{plain} {interval:?}/{numeric}"
            );
            assert!(interval.lo <= p.volume_mm3() && p.volume_mm3() <= interval.hi);
            assert!(interval.width() < 1e-6);
        }
    }
}
#[test]
fn source_precision_and_true_numeric_underflow_refuse_mass_claims() {
    let p = Profile::new(
        Head::Crowned {
            radius_mm: 2.8,
            height_mm: 2.4,
            crown_radius_mm: 3.,
            top_land_radius_mm: 2.1,
            neck_radius_mm: 0.1,
        },
        1.5,
        1.,
    )
    .unwrap();
    let segments: Vec<_> = p.segments().collect();
    let loose = Profile::from_transported_segments(
        &segments,
        TransportError {
            numeric_conversion_mm: 0.005,
            arithmetic_mm: 0.,
        },
    )
    .unwrap();
    assert!(precision_envelope(&loose).is_none());
    assert!(volume_bounds(&loose).is_none());
    let tiny = Profile::new(
        Head::Cylindrical {
            radius_mm: 2e-200,
            height_mm: 1e-200,
        },
        1e-200,
        0.,
    )
    .unwrap();
    assert!(volume_bounds(&tiny).is_none());
}
