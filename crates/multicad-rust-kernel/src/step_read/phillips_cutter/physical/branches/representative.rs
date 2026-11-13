//! Test-only data-driven closed-family feasibility. No import admission.
use super::*;
use crate::phillips_drive::{Cutter, Domain, Inner, Spec};
fn dot(a: Point, b: Point) -> f64 {
    a.into_iter().zip(b).map(|(a, b)| a * b).sum()
}
fn sub(a: Point, b: Point) -> Point {
    std::array::from_fn(|i| a[i] - b[i])
}
fn cross(a: Point, b: Point) -> Point {
    std::array::from_fn(|i| a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3])
}
fn unit(a: Point) -> Point {
    let n = dot(a, a).sqrt();
    a.map(|a| a / n)
}
struct Candidate {
    source: Cutter,
    pose: Pose,
}
// Scalar fitting is an untrusted proposal, not exact decoding of an authoring
// recipe. The original raw Body remains unchanged. Only compare() earns the
// separate finite-graph and whole-source measured correspondence. Selection
// thresholds below are hints/refusals; they never accept geometric incidence.
fn candidate(b: &Boundary) -> std::result::Result<Candidate, PhillipsReadError> {
    let fail = PhillipsReadError::Unsupported;
    if (b.faces().len(), b.edges().len(), b.vertices().len()) != (29, 76, 49) {
        return Err(fail);
    }
    let (axis, major, minor) = b
        .faces()
        .iter()
        .find_map(|f| match f.support {
            Surface::Torus {
                axis,
                major_radius_mm,
                minor_radius_mm,
            } => Some((axis, major_radius_mm, minor_radius_mm)),
            _ => None,
        })
        .ok_or(fail)?;
    let n = unit(axis.normal);
    let v = unit(cross(n, axis.reference));
    let u = cross(v, n);
    let origin = std::array::from_fn(|i| axis.origin[i] - minor * n[i]);
    let local = |p: Point| {
        [
            dot(sub(p, origin), u),
            dot(sub(p, origin), v),
            dot(sub(p, origin), n),
        ]
    };
    let vector = |p: Point| [dot(p, u), dot(p, v), dot(p, n)];
    let mut cones = b
        .faces()
        .iter()
        .filter_map(|f| match f.support {
            Surface::Cone {
                axis,
                radius_mm,
                semi_angle_rad,
            } => Some((axis, radius_mm, semi_angle_rad)),
            _ => None,
        })
        .collect::<Vec<_>>();
    cones.sort_by(|a, b| a.2.total_cmp(&b.2));
    let (_, mhalf, alpha) = *cones.first().ok_or(fail)?;
    let (_, floor0, phi) = *cones.last().ok_or(fail)?;
    let k = alpha.tan();
    let kc = phi.tan();
    // Intersect the two actual meridian cone laws r=a-k*z.
    let floor = (mhalf - floor0) / (k - kc);
    let g = 2. * (mhalf - k * floor);
    let tip = floor0 / kc;
    let (caxis, inner) = b
        .faces()
        .iter()
        .find_map(|f| match f.support {
            Surface::Cylinder { axis, radius_mm }
                if local(axis.origin)[0] > 0. && local(axis.origin)[1].abs() < 1e-8 =>
            {
                Some((axis, radius_mm))
            }
            _ => None,
        })
        .ok_or(fail)?;
    let cn = vector(unit(caxis.normal));
    let beta = (-cn[0]).atan2(cn[2]);
    let theta = b
        .faces()
        .iter()
        .find_map(|f| match f.support {
            Surface::Plane { axis } => {
                let x = vector(unit(axis.normal));
                (x[0] > 0. && x[1] < 0. && x[2] > 0.).then(|| x[0].hypot(x[2]).atan2(-x[1]))
            }
            _ => None,
        })
        .ok_or(fail)?;
    // Cylinder center x = b/2 + tan(beta)*(tip-b/(2*kc))
    //                    + inner_radius/cos(beta).
    let bwidth = 2. * (local(caxis.origin)[0] - tip * beta.tan() - inner / beta.cos())
        / (1. - beta.tan() / kc);
    let spec = Spec {
        mouth_diameter_mm: 2. * mhalf,
        mouth_radius_mm: minor,
        b_width_mm: bwidth,
        g_diameter_mm: g,
        tilt_angle_deg: beta.to_degrees(),
        outer_flank_angle_deg: alpha.to_degrees(),
        floor_angle_deg: (1. / kc).atan().to_degrees(),
        outer_corner_angle_deg: 2. * theta.to_degrees(),
        inner: Inner::Rounded { radius_mm: inner },
    };
    let expected_major = mhalf + minor * (1. - alpha.sin()) / alpha.cos();
    if (expected_major - major).abs() > 1e-8 {
        return Err(fail);
    }
    let reach = b
        .vertices()
        .iter()
        .flat_map(|v| local(v.point))
        .map(f64::abs)
        .fold(0., f64::max)
        .max(1.);
    // These are deliberately new inactive padding bounds, not recovered
    // original prism limits. Full topology/trim matching below must reject
    // any padding that changes the actual visible finite component.
    let domain = Domain {
        corner_start_depth_mm: -4. * reach,
        corner_end_depth_mm: 4. * reach,
        corner_reach_mm: 16. * reach,
    };
    eprintln!("derived={spec:?} domain={domain:?} origin={origin:?}");
    let source = Cutter::new(spec, domain).map_err(|_| fail)?;
    let pose = Pose::new(origin, u, v).map_err(|_| fail)?;
    Ok(Candidate { source, pose })
}
fn compare(
    body: &Body<'_>,
    c: &Candidate,
) -> std::result::Result<PhillipsSourceKnownProof, PhillipsReadError> {
    let b = body.cavity();
    // Reuses the qualified source-known comparator only after constructing a
    // candidate from data. This is still NOT exact raw sewing, face filling,
    // independent mass, containment, embedding or manufacturing admission.
    let entities = part21::parse(body.literal_input()).map_err(PhillipsReadError::from)?;
    let mut shell = Shell::new(&entities, 1e-7, true);
    let views = b
        .faces()
        .iter()
        .map(|f| SourceFaceView {
            surface: f.surface_entity,
            forward: f.cavity_forward,
            uses: f
                .uses
                .iter()
                .map(|&(e, s)| (b.edges()[e].entity, s))
                .collect(),
        })
        .collect::<Vec<_>>();
    verify_component(&mut shell, &views, &c.source, &c.pose, 1e-7).map_err(PhillipsReadError::from)
}
// Fixture construction is intentionally outside candidate(): none of these
// original source dimensions are available to the reconstruction function.
fn scaled_text(scale: f64, posed: bool) -> String {
    use crate::{
        head_profile::{Head, Profile},
        helical_thread::{Spec as ThreadSpec, Thread},
        phillips_drive::{Cavity, PhillipsScrew},
        profiled_screw::SourceSpec,
    };
    let (mut s, mut d) = crate::phillips_drive::test_fixture(0);
    s.mouth_diameter_mm *= scale;
    s.mouth_radius_mm *= scale;
    s.b_width_mm *= scale;
    s.g_diameter_mm *= scale;
    let Inner::Rounded { radius_mm } = s.inner else {
        unreachable!()
    };
    s.inner = Inner::Rounded {
        radius_mm: radius_mm * scale,
    };
    d.corner_start_depth_mm *= scale;
    d.corner_end_depth_mm *= scale;
    d.corner_reach_mm *= scale;
    let cavity = Cavity::new(Cutter::new(s, d).unwrap(), 1e-5, 100_000).unwrap();
    let head = Head::Crowned {
        radius_mm: 2. * scale,
        height_mm: 1.6 * scale,
        crown_radius_mm: 3.2 * scale,
        top_land_radius_mm: 1.25 * scale,
        neck_radius_mm: 0.1 * scale,
    };
    let thread = Thread::new(ThreadSpec {
        major_radius_mm: 0.97 * scale,
        pitch_mm: 0.4 * scale,
        length_mm: 2. * scale,
        depth_mm: 0.24 * scale,
        crest_width_mm: 0.05 * scale,
        lower_flank_deg: 30.,
        upper_flank_deg: 30.,
        starts: 1,
        left_hand: posed,
        root_radius_mm: None,
    })
    .unwrap();
    let screw = PhillipsScrew::new(
        thread,
        Profile::new_for_thread_union(head, scale, 0.).unwrap(),
        SourceSpec {
            head: Some(head),
            shank_length_mm: 0.,
            minimum_wall_mm: 1e-5,
        },
        cavity,
        100_000,
        &|| false,
    )
    .unwrap();
    let (o, u, v) = if posed {
        ([105., -20., 60.], [0.8, 0.6, 0.], [0.36, -0.48, -0.8])
    } else {
        ([0.; 3], [1., 0., 0.], [0., 1., 0.])
    };
    crate::step::write_phillips_screw(&screw, o, u, v, 1e-7).unwrap()
}
#[test]
fn actual_h0_candidate_whole_correspondence_probe() {
    for scale in [1., 1.137] {
        for posed in [false, true] {
            let text = scaled_text(scale, posed);
            let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
            let c = candidate(body.cavity()).unwrap();
            let r = compare(&body, &c);
            eprintln!("posed{posed} comparison {r:?}");
            let r = r.unwrap();
            assert_eq!(
                (r.face_count(), r.edge_count(), r.vertex_count()),
                (29, 76, 49)
            );
        }
    }
}

#[test]
fn altered_actual_support_or_internal_trim_control_refuses_whole_candidate() {
    let text = scaled_text(1., false);
    let body = read_phillips_body(&text, 1e-7, 1_000_000).unwrap();
    let entities = part21::parse(&text).unwrap();
    let torus = body
        .cavity()
        .faces()
        .iter()
        .find(|f| matches!(f.support, Surface::Torus { .. }))
        .unwrap();
    let words = part21::numeric_spellings(&entities[&torus.surface_entity][0]).unwrap();
    let record = text
        .lines()
        .find(|l| l.starts_with(&format!("#{}=", torus.surface_entity)))
        .unwrap();
    let changed_record = record.replacen(
        words[1],
        &format!("{:.17E}", words[1].parse::<f64>().unwrap() + 0.001),
        1,
    );
    assert_ne!(record, changed_record);
    let changed = text.replacen(record, &changed_record, 1);
    let b = read_phillips_body(&changed, 1e-7, 1_000_000).unwrap();
    assert!(candidate(b.cavity()).and_then(|c| compare(&b, &c)).is_err());
    // A side-plane offset is not used in dimension fitting. Its perturbation
    // must therefore be rejected by the complete support comparator itself.
    let plane = body.cavity().faces().iter().find(|f|matches!(f.support,Surface::Plane{axis} if axis.normal[0]>0. && axis.normal[1]<0. && axis.normal[2]>0.)).unwrap();
    let Value::Reference(axis) = entities[&plane.surface_entity][0].args[1] else {
        panic!()
    };
    let Value::Reference(point) = entities[&axis][0].args[1] else {
        panic!()
    };
    let words = part21::numeric_spellings(&entities[&point][0]).unwrap();
    let record = text
        .lines()
        .find(|l| l.starts_with(&format!("#{point}=")))
        .unwrap();
    let changed_record = record.replacen(
        words[0],
        &format!("{:.17E}", words[0].parse::<f64>().unwrap() + 0.001),
        1,
    );
    assert_ne!(record, changed_record);
    let changed = text.replacen(record, &changed_record, 1);
    let b = read_phillips_body(&changed, 1e-7, 1_000_000).unwrap();
    let c = candidate(b.cavity()).unwrap();
    assert!(compare(&b, &c).is_err());
    let point = body
        .cavity()
        .edges()
        .iter()
        .find_map(|e| match &e.curve {
            Curve::CubicIntersection { pole_entities, .. } => Some(pole_entities[1]),
            _ => None,
        })
        .unwrap();
    let words = part21::numeric_spellings(&entities[&point][0]).unwrap();
    let record = text
        .lines()
        .find(|l| l.starts_with(&format!("#{point}=")))
        .unwrap();
    let changed_record = record.replacen(
        words[0],
        &format!("{:.17E}", words[0].parse::<f64>().unwrap() + 0.001),
        1,
    );
    assert_ne!(record, changed_record);
    let changed = text.replacen(record, &changed_record, 1);
    let b = read_phillips_body(&changed, 1e-7, 1_000_000).unwrap();
    let c = candidate(b.cavity()).unwrap();
    assert!(compare(&b, &c).is_err());
    assert_eq!(body.literal_input(), text);
}
