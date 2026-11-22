//! Raw finite-pipe boundary preflight, without manufacturing admission.
#[cfg(test)]
mod tests {
    use super::super::*;
    fn add(a: f64, b: f64) -> (f64, f64) {
        let sum = a + b;
        let split = sum - a;
        (sum, (a - (sum - split)) + (b - split))
    }
    #[test]
    fn nearest_f64_pipe_caps_cannot_represent_a_proved_literal_radial_closure() {
        let p = 25.4 / 27.;
        let k = 1. / 32.;
        let m = 3_f64.sqrt();
        let h = p / (1. / (m - k) + 1. / (m + k));
        let depth = h - 2. * 0.0645 * p;
        let width = (p - depth * (1. / (m - k) + 1. / (m + k))) / 2.;
        let spec = Spec {
            major_radius_mm: (0.37360 * 25.4 - 0.1615 * 25.4 / 16. + depth) / 2.,
            pitch_mm: p,
            length_mm: 0.2639 * 25.4,
            radial_taper: k,
            section: Section::NptFlat {
                depth_mm: depth,
                flat_width_mm: width,
            },
            crest_radius_limit_mm: None,
        };
        let t = Thread::new(spec, 1e-7, 1e-6).unwrap();
        let q = t.bands()[0].phase_mm[1];
        let (anchor, residual) = add(spec.major_radius_mm, k * q);
        assert_eq!(
            t.bands()[1].curve,
            Curve::Line {
                phase_mm: q,
                radius_mm: anchor,
                slope: -m,
            }
        );
        assert_ne!(
            residual, 0.,
            "This fixture must contain actual distinct literal incident supports"
        );
        let angle = -std::f64::consts::TAU * q / p;
        let (sin, cos) = angle.sin_cos();
        let right_radius = anchor - k * q;
        assert_eq!(right_radius, spec.major_radius_mm);
        let left = [spec.major_radius_mm * cos, spec.major_radius_mm * sin, 0.];
        let right = [right_radius * cos, right_radius * sin, 0.];
        assert_eq!(
            left, right,
            "Nearest f64 endpoints collapse the nonzero source radial closure"
        );
        // At fixed phase q the literal radial closure has P_s=delta*e_r
        // and P_z=k*e_r+omega*r*e_theta+e_z. Thus |P_s cross P_z|
        // is at least |delta| everywhere, despite zero angular projection.
        let area_lower = (I::point(residual.abs()) * I::point(spec.length_mm)).lo;
        assert!(area_lower > 0.);
        use crate::radial_tensor::{
            BoundaryInput, BoundaryOwnership, Curve as BoundaryCurve, Edge,
        };
        let raw = BoundaryInput {
            vertices: vec![left, right],
            curves: vec![BoundaryCurve::Spline {
                spans: 1,
                poles: vec![left; 10],
            }],
            edges: vec![Edge {
                vertices: [0, 1],
                curve: 0,
                forward: true,
            }],
            surfaces: vec![],
            faces: vec![],
        };
        // The actual endpoint preflight runs before whole-manifold checks;
        // this minimal raw edge is deliberately not a fabricated body proof.
        assert!(matches!(
            BoundaryOwnership::new(raw, 1e-7, None),
            Err(Error::Invalid)
        ));
        println!(
            "NPT {spec:?}; q={q:?}; left-minus-right exact radial residual={residual:?}; positive closure area lower={area_lower:?}; proposed nearest cap vertices={left:?} and {right:?}"
        );
        let a = 27.5_f64.to_radians();
        let p = 25.4 / 28.;
        let m = 1. / a.tan();
        let depth = p / (3. * a.tan());
        let r =
            (p * (m * m - k * k) / (2. * m) - depth) / (2. * (1. / a.sin() - (1. + k * k).sqrt()));
        let spec = Spec {
            major_radius_mm: 9.728 / 2. - k * 4.,
            pitch_mm: p,
            length_mm: 6.5,
            radial_taper: k,
            section: Section::BsptRounded {
                depth_mm: depth,
                rounding_radius_mm: r,
            },
            crest_radius_limit_mm: None,
        };
        let t = Thread::new(spec, 1e-7, 1e-6).unwrap();
        for i in 0..5 {
            let q = t.bands()[i].phase_mm[1];
            let next = (i + 1) % 5;
            let a = t.unclipped(i, I::point(q), I::point(0.));
            let b = t.unclipped(next, I::point(if next == 0 { 0. } else { q }), I::point(0.));
            let angle = -std::f64::consts::TAU * q / spec.pitch_mm;
            let (sin, cos) = angle.sin_cos();
            let a = [a.midpoint() * cos, a.midpoint() * sin, 0.];
            let b = [b.midpoint() * cos, b.midpoint() * sin, 0.];
            println!(
                "BSPT junction {i}: q={q:?}; proposed midpoint cap vertices={a:?},{b:?}; equal={}",
                a == b
            );
        }
    }
}
