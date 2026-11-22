use super::*;
use crate::helical_tube_cut::{Cut, Spec};
use crate::step_read::brep::pcurve::{Binding, bindings, edge, trim_surface, uv};

const TOLERANCE: f64 = 1e-7;

/// The kernel's own helical tube cut export: its cavity support is a
/// degree-11 tensor surface whose interior knots are only C0.
fn exported() -> String {
    let cut = Cut::new(
        Spec {
            stock_radius_mm: 8.,
            stock_height_mm: 30.,
            tube_radius_mm: 0.6,
            pitch_mm: 2.,
            path_height_mm: 15.,
            start_height_mm: 16.2,
            start_angle_rad: 0.,
        },
        TOLERANCE,
    )
    .unwrap();
    crate::step::write_helical_tube_cut(
        cut,
        [17., 4., -12.],
        [0.8, 0., -0.6],
        [0., 1., 0.],
        TOLERANCE,
    )
    .unwrap()
}
type Owned = (Binding, super::super::super::uv::Definition, Surface, Edge);
/// Every association on a C0 tensor surface, with its decoded inputs.
fn c0_pairs(source: &Source<'_>) -> Vec<Owned> {
    let mut b = Budget::new(512_000_000, None);
    let mut out = Vec::new();
    for binding in bindings(source, &mut b).unwrap() {
        let definition = uv::definition(source, binding.pcurve, &mut b).unwrap();
        let surface = trim_surface::read(source, &definition, &mut b).unwrap();
        let Surface::Nurbs { surface: net, .. } = &surface else {
            continue;
        };
        if !has_c0_joint(net) {
            continue;
        }
        let edge = edge(source, binding, &mut b).unwrap();
        out.push((binding, definition, surface, edge));
    }
    out
}
fn run(source: &Source<'_>, pair: &Owned, tol: f64, work: usize) -> Result<Option<Proof>> {
    let (binding, d, s, e) = pair;
    certify(
        source,
        binding.curve,
        d,
        s,
        e,
        tol,
        &mut Budget::new(work, None),
    )
}

#[test]
fn straight_iso_traces_and_the_circle_rim_certify_on_the_exported_c0_support() {
    let text = exported();
    let source = Source::parse(&text, TOLERANCE).unwrap();
    let pairs = c0_pairs(&source);
    let mut certified = 0;
    let mut circles = 0;
    for pair in &pairs {
        let straight =
            matches!(&pair.1.geometry, Curve::Nurbs { curve, .. } if curve.degree() == 1);
        let result = run(&source, pair, TOLERANCE, 512_000_000).unwrap();
        if straight {
            let Some(Proof::Certified(error)) = result else {
                panic!("straight trace {} not certified", pair.0.pcurve)
            };
            assert!(error <= TOLERANCE);
            certified += 1;
            circles += usize::from(matches!(pair.3.curve, Curve::Circle { .. }));
        } else {
            // Curved traces are not this route's: they fall through.
            assert!(result.is_none());
        }
    }
    // Two rail iso lines and the floor rim circle.
    assert_eq!((certified, circles), (3, 1));
}

#[test]
fn a_tighter_tolerance_or_a_moved_rim_is_never_certified() {
    let text = exported();
    let source = Source::parse(&text, TOLERANCE).unwrap();
    for pair in c0_pairs(&source) {
        if !matches!(&pair.1.geometry, Curve::Nurbs { curve, .. } if curve.degree() == 1) {
            continue;
        }
        let result = run(&source, &pair, 1e-12, 512_000_000).unwrap();
        assert!(matches!(result, Some(Proof::Unsupported)));
    }
    // Move ONE interior control of the support's v=0 column (row 5, which
    // shapes the rim) by 1e-5 mm. Vertices and the 3D circle are unchanged.
    let start = text.find("B_SPLINE_SURFACE_WITH_KNOTS('',11,11,(").unwrap();
    let rows = &text[start..];
    let row = rows.split("),(").nth(5).unwrap();
    let id = row
        .trim_start_matches(|c| c != '#')
        .split(',')
        .next()
        .unwrap();
    let line = text
        .lines()
        .find(|l| l.starts_with(&format!("{id}=CARTESIAN_POINT")))
        .unwrap();
    let open = line.find("',(").unwrap() + 3;
    let comma = open + line[open..].find(',').unwrap();
    let x: f64 = line[open..comma].parse().unwrap();
    let changed = format!("{}{:.16E}{}", &line[..open], x + 1e-5, &line[comma..]);
    let moved = text.replacen(line, &changed, 1);
    assert_ne!(moved, text);
    let source = Source::parse(&moved, TOLERANCE).unwrap();
    let mut refused = 0;
    for pair in c0_pairs(&source) {
        if matches!(pair.3.curve, Curve::Circle { .. }) {
            assert!(matches!(
                run(&source, &pair, TOLERANCE, 512_000_000),
                Ok(Some(Proof::Unsupported)) | Err(ReadError::Unsupported)
            ));
            refused += 1;
        }
    }
    assert_eq!(refused, 1);
}

#[test]
fn every_charge_is_metered_and_exhaustion_is_a_limit() {
    let text = exported();
    let source = Source::parse(&text, TOLERANCE).unwrap();
    let pairs = c0_pairs(&source);
    let pair = pairs
        .iter()
        .find(|p| matches!(p.3.curve, Curve::Circle { .. }))
        .unwrap();
    let mut full = Budget::new(512_000_000, None);
    let (binding, d, s, e) = pair;
    assert!(matches!(
        certify(&source, binding.curve, d, s, e, TOLERANCE, &mut full).unwrap(),
        Some(Proof::Certified(_))
    ));
    let used = full.used();
    assert_eq!(
        run(&source, pair, TOLERANCE, used - 1).err(),
        Some(ReadError::Limit)
    );
    assert!(matches!(
        run(&source, pair, TOLERANCE, used).unwrap(),
        Some(Proof::Certified(_))
    ));
}

#[test]
fn knot_merging_requires_one_exact_domain_and_keeps_order() {
    let q = |n: i64| BigRational::from_integer(n.into());
    let mut b = Budget::new(1_000_000, None);
    let merged = merge(&[q(0), q(2), q(5)], &[q(0), q(1), q(2), q(5)], &mut b).unwrap();
    assert_eq!(merged, vec![q(0), q(1), q(2), q(5)]);
    assert_eq!(
        merge(&[q(0), q(5)], &[q(0), q(4)], &mut b).err(),
        Some(ReadError::Unsupported)
    );
    assert_eq!(
        merge(&[], &[q(0)], &mut b).err(),
        Some(ReadError::Unsupported)
    );
}
