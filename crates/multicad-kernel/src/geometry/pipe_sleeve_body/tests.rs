use super::*;
use multicad_rust_kernel::axial_pipe_thread::authored::AuthoredSection;

/// 1/8-27 NPT on the NBS plain-face gauge: external pipe end and internal
/// receiver with `allowance` taken radially from each side.
fn npt(internal: bool, allowance: f64) -> AuthoredSpec {
    let (e1, l1, l2) = (0.37360_f64, 0.1615, 0.2639);
    let p = 25.4 / 27.;
    let k: f64 = if internal { -1. / 32. } else { 1. / 32. };
    let m = 3f64.sqrt();
    let s = 1. / (m - k) + 1. / (m + k);
    let d = p / s - (0.033 + 0.096) * p;
    let gauge = if internal { 0. } else { l1 * 25.4 };
    let pitch_diameter = e1 * 25.4 - 2. * k * gauge;
    let sign = if internal { 1. } else { -1. };
    AuthoredSpec {
        major_radius_mm: (pitch_diameter + d) / 2. + sign * allowance,
        pitch_mm: p,
        length_mm: if internal {
            l1 * 25.4 + 4. * p
        } else {
            l2 * 25.4
        },
        radial_taper: k,
        section: AuthoredSection::NptFlat {
            depth_mm: d,
            flat_width_mm: (p - d * s) / 2.,
        },
        crest_radius_limit_mm: (!internal).then_some(0.405 * 25.4 / 2.),
    }
}
/// The receiving chart on the opposing axis, with the gauge planes and the
/// helical phase of the external chart aligned.
fn receiving_frame() -> PlaneFrame {
    let origin = 0.1615 * 25.4;
    let phi = std::f64::consts::TAU * origin / (25.4 / 27.);
    PlaneFrame {
        origin: [0., 0., origin],
        u: [phi.portable_cos(), phi.portable_sin(), 0.],
        v: [phi.portable_sin(), -phi.portable_cos(), 0.],
    }
}
fn sleeve(allowance: f64, radius: f64) -> Result<Solid, GeometryError> {
    Solid::rust_tapered_pipe_sleeve_with_interrupt(
        npt(true, allowance),
        radius,
        EntityId::new(),
        receiving_frame(),
        GeometryPolicy::default(),
        None,
    )
}

#[test]
fn sleeve_is_exact_stock_minus_the_retained_cavity_with_annular_caps_and_wall() {
    let spec = npt(true, 0.005);
    let radius = spec.major_radius_mm + 2.;
    let solid = sleeve(0.005, radius).unwrap();
    let cavity = Solid::rust_authored_axial_pipe_thread_with_interrupt(
        spec,
        EntityId::new(),
        PlaneFrame::default(),
        GeometryPolicy::default(),
        None,
    )
    .unwrap();
    let stock = PI * radius * radius * spec.length_mm;
    let volume = solid.volume_mm3().unwrap();
    assert!((volume - (stock - cavity.volume_mm3().unwrap())).abs() < 1e-6 * stock);
    assert_eq!(solid.solid_count().unwrap(), 1);
    let faces = solid.faces().unwrap();
    let inner = cavity.faces().unwrap();
    assert_eq!(faces.len(), inner.len() + 1);
    let wall = faces.last().unwrap();
    assert!((wall.area_mm2 - 2. * PI * radius * spec.length_mm).abs() < 1e-9);
    assert!(wall.plane.is_none());
    let (body, frame) = solid.pipe_sleeve_body().unwrap();
    let n = frame.normal();
    for (end, &cap) in body.caps.iter().enumerate() {
        let annulus = &faces[cap];
        let disc = PI * radius * radius;
        assert!((annulus.area_mm2 - (disc - inner[cap].area_mm2)).abs() < 1e-6 * disc);
        let (origin, normal) = annulus.plane.unwrap();
        let sign = if end == 0 { -1. } else { 1. };
        for k in 0..3 {
            assert!((normal[k] - sign * n[k]).abs() < 1e-15);
        }
        let z = [0., spec.length_mm][end];
        assert!((origin[2] - (frame.origin[2] - z)).abs() < 1e-12);
    }
    // World box of the opposed stock: [origin - L, origin] along z.
    let b = solid.bounds().unwrap();
    assert!(b.min[2] < frame.origin[2] - spec.length_mm + 1e-9);
    assert!(b.max[2] > frame.origin[2] && b.max[2] < frame.origin[2] + 1e-9);
    let top = solid.topology_snapshot(body.feature, 1).unwrap();
    top.validate().unwrap();
    for role in ["stock.wall", "stock.rim/start", "stock.rim/end"] {
        assert!(
            top.elements
                .iter()
                .any(|e| e.origins.iter().any(|o| o.role == role))
        );
    }
    assert!(top.elements.iter().any(|e| {
        e.origins
            .iter()
            .any(|o| o.role.starts_with("receiving.inward/thread/"))
    }));
    let moved = solid
        .transformed(RigidTransform {
            rotation: [[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]],
            translation_mm: [3., 4., 5.],
        })
        .unwrap();
    assert!((moved.volume_mm3().unwrap() - volume).abs() < 1e-9 * volume);
}

#[test]
fn sleeve_step_reads_back_independently_with_the_same_material() {
    let spec = npt(true, 0.005);
    let radius = spec.major_radius_mm + 2.;
    let solid = sleeve(0.005, radius).unwrap();
    assert!(solid.supports_step_export());
    let file = tempfile::NamedTempFile::new().unwrap();
    solid.write_step(file.path()).unwrap();
    let text = std::fs::read_to_string(file.path()).unwrap();
    let back = Solid::import_rust_step(&text, EntityId::new(), GeometryPolicy::default()).unwrap();
    let (body, _) = back.pipe_sleeve_body().unwrap();
    // An import has material but never an inferred authored law.
    assert!(body.law().is_none());
    assert_eq!(body.radius_mm(), radius);
    assert!((body.height_mm() - spec.length_mm).abs() < 1e-12);
    let area: f64 = solid.faces().unwrap().iter().map(|f| f.area_mm2).sum();
    let v = solid.volume_mm3().unwrap();
    assert!((back.volume_mm3().unwrap() - v).abs() < 1e-7 * area + 1e-7 * v);
    let a = solid.bounds().unwrap();
    let b = back.bounds().unwrap();
    for k in 0..3 {
        assert!((a.min[k] - b.min[k]).abs() < 1e-6 && (a.max[k] - b.max[k]).abs() < 1e-6);
    }
}

#[test]
fn sleeve_refuses_thin_stock_and_external_or_clipped_sources() {
    let spec = npt(true, 0.005);
    let error = sleeve(0.005, spec.major_radius_mm).err().unwrap();
    assert_eq!(error.kind, GeometryErrorKind::Input);
    for source in [
        npt(false, 0.005),
        AuthoredSpec {
            crest_radius_limit_mm: Some(spec.major_radius_mm + 1.),
            ..spec
        },
    ] {
        let error = Solid::rust_tapered_pipe_sleeve_with_interrupt(
            source,
            spec.major_radius_mm + 2.,
            EntityId::new(),
            receiving_frame(),
            GeometryPolicy::default(),
            None,
        )
        .err()
        .unwrap();
        assert_eq!(error.kind, GeometryErrorKind::Unsupported);
    }
}

#[test]
fn gauge_aligned_external_pipe_thread_clears_the_sleeve_in_either_operand_order() {
    let spec = npt(true, 0.005);
    let female = sleeve(0.005, spec.major_radius_mm + 2.).unwrap();
    let male = Solid::rust_authored_axial_pipe_thread_with_interrupt(
        npt(false, 0.005),
        EntityId::new(),
        PlaneFrame::default(),
        GeometryPolicy::default(),
        None,
    )
    .unwrap();
    assert!(male.clipped_npt_body().is_some());
    let policy = GeometryPolicy::default();
    for (a, b) in [(&male, &female), (&female, &male)] {
        assert!(
            a.boolean(b, BooleanOperation::Intersect, policy)
                .unwrap()
                .is_none()
        );
        let kept = a
            .boolean(b, BooleanOperation::Subtract, policy)
            .unwrap()
            .unwrap();
        assert_eq!(kept.volume_mm3().unwrap(), a.volume_mm3().unwrap());
    }
    // Without allowance the flanks touch: no clearance is admitted, so this
    // provider never reports an empty intersection for a contact pair.
    let tight = sleeve(0., spec.major_radius_mm + 2.).unwrap();
    let flush = Solid::rust_authored_axial_pipe_thread_with_interrupt(
        npt(false, 0.),
        EntityId::new(),
        PlaneFrame::default(),
        policy,
        None,
    )
    .unwrap();
    assert!(
        flush
            .pipe_sleeve_clearance_boolean(&tight, BooleanOperation::Intersect, policy)
            .is_none()
    );
    // A half-pitch phase error puts crests into crests: refused as well.
    let (body, frame) = female.pipe_sleeve_body().unwrap();
    let half = PI;
    let shifted = Solid {
        body: SolidBody::PipeSleeve(
            Rc::new(PipeSleeveBody {
                feature: body.feature,
                cavity: body.cavity.share().unwrap(),
                radius_mm: body.radius_mm,
                height_mm: body.height_mm,
                caps: body.caps,
                law: body.law,
            }),
            PlaneFrame {
                u: rotate(frame, [half.portable_cos(), half.portable_sin(), 0.]),
                v: rotate(frame, [-half.portable_sin(), half.portable_cos(), 0.]),
                ..frame
            },
        ),
        policy,
        provenance: vec![],
        rust_origins: vec![],
        retained_placement: None,
    };
    assert!(
        male.pipe_sleeve_clearance_boolean(&shifted, BooleanOperation::Intersect, policy)
            .is_none()
    );
}
