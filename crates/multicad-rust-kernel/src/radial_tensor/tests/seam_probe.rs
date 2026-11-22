//! Actual finite source reproduction only; no seam admission or certifier bypass.
use super::*;
use crate::Rational;

fn radial_step() -> BoundaryInput {
    let mut input = super::cam::wedge(&[1]);
    input.vertices[1] = [30., -40., 0.];
    input.vertices[6] = [30., -40., 4.];
    for edge in &input.edges {
        if edge.vertices.contains(&1) || edge.vertices.contains(&6) {
            let Curve::Spline { spans, poles } = &mut input.curves[edge.curve] else {
                panic!("actual linear edge");
            };
            let [a, b] = edge.vertices.map(|v| input.vertices[v]);
            *poles = (0..=9 * *spans)
                .map(|i| {
                    std::array::from_fn(|k| a[k] + (b[k] - a[k]) * i as f64 / (9 * *spans) as f64)
                })
                .collect();
        }
    }
    for face in &input.faces {
        if let Surface::Tensor(t) = &mut input.surfaces[face.surface] {
            let edge = input.edges[face.outer[0].0];
            let Curve::Spline { poles: base, .. } = &input.curves[edge.curve] else {
                panic!("actual bottom rim");
            };
            let poles = base
                .iter()
                .flat_map(|p| (0..=9).map(move |j| [p[0], p[1], 4. * j as f64 / 9.]))
                .collect();
            *t = RetainedTensor::new([1, 1], poles).unwrap();
        }
    }
    input
}
fn trimmed(input: BoundaryInput) -> TrimmedBoundary {
    BoundaryOwnership::new(input, 1e-7, None)
        .unwrap()
        .bind_surfaces(1e-7, None)
        .unwrap()
}
fn exact_cross(a: Point, b: Point) -> Rational {
    let q = crate::exact_binary::binary_fraction;
    q(a[0])
        .unwrap()
        .checked_mul(q(b[1]).unwrap())
        .unwrap()
        .checked_sub(q(a[1]).unwrap().checked_mul(q(b[0]).unwrap()).unwrap())
        .unwrap()
}

#[test]
fn actual_fixed_phase_radial_step_has_positive_area_but_strict_chart_refuses() {
    let input = radial_step();
    let original = input.clone();
    assert_eq!(input.faces.len(), 6);
    let Surface::Tensor(seam) = &input.surfaces[2] else {
        panic!("actual seam");
    };
    // Every retained coefficient lies on the same physical radial half-plane;
    // exact binary-rational products prove equality, without a tolerance test.
    for pole in seam.poles() {
        assert_eq!(
            exact_cross([3., -4., 0.], *pole),
            Rational::new(0, 1).unwrap()
        );
    }
    let derivatives = seam.derivative_bounds([0, 0], [[0., 1.]; 2]).unwrap();
    let du = derivatives.du.map(|v| I::new(v[0], v[1]));
    let dv = derivatives.dv.map(|v| I::new(v[0], v[1]));
    let normal = cross(du, dv);
    let area = (normal[0].square() + normal[1].square() + normal[2].square()).sqrt();
    assert!(area.lo > 179. && area.hi < 181.);
    assert_eq!(
        seam.certify_radial_chart(true, None),
        Err(Error::Unresolved)
    );
    for surface in &input.surfaces[3..5] {
        let Surface::Tensor(t) = surface else {
            panic!("strict adjacent face");
        };
        t.certify_radial_chart(true, None).unwrap();
    }
    // This is a real owned and support-bound six-face source; refusal occurs
    // at strict radial-chart admission, not a malformed/unsupported raw graph.
    let bound = trimmed(input.clone());
    assert_eq!(bound.input(), &input);
    assert!(matches!(
        bound.certify_axial(1e-7, 1e-6, None),
        Err(Error::Unresolved)
    ));
    assert_eq!(input, original);
}

fn pose(input: &BoundaryInput) -> BoundaryInput {
    let vector = |p: Point| [0.6 * p[0] - 0.8 * p[1], 0.8 * p[0] + 0.6 * p[1], p[2]];
    let axis = |a: Axis| Axis {
        origin: vector(a.origin),
        normal: vector(a.normal),
        reference: vector(a.reference),
    };
    let mut out = input.clone();
    out.vertices.iter_mut().for_each(|p| *p = vector(*p));
    for curve in &mut out.curves {
        match curve {
            Curve::Spline { poles, .. } => poles.iter_mut().for_each(|p| *p = vector(*p)),
            Curve::Circle { axis: a, .. } => *a = axis(*a),
        }
    }
    for surface in &mut out.surfaces {
        match surface {
            Surface::Tensor(t) => {
                *t = RetainedTensor::new(t.spans(), t.poles().iter().copied().map(vector).collect())
                    .unwrap()
            }
            Surface::Plane(a) | Surface::Cone { axis: a, .. } => *a = axis(*a),
        }
    }
    out
}

#[test]
fn actual_proper_pose_can_break_literal_rank_and_must_not_be_welded() {
    let original = radial_step();
    let posed = pose(&original);
    let Surface::Tensor(seam) = &posed.surfaces[2] else {
        panic!("actual posed seam");
    };
    assert_ne!(
        exact_cross(seam.poles()[0], seam.poles()[20]),
        Rational::new(0, 1).unwrap()
    );
    assert_eq!(
        seam.certify_radial_chart(true, None),
        Err(Error::Unresolved)
    );
    let bound = trimmed(posed.clone());
    assert_eq!(bound.input(), &posed);
    assert!(matches!(
        bound.certify_axial(1e-7, 1e-6, None),
        Err(Error::Unresolved)
    ));
    assert_eq!(original, radial_step());
}
