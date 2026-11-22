//! Exact axial trims of a rounded thread by a coaxial slab.
//!
//! A rectangular prism whose faces are perpendicular or parallel to the
//! thread axis, and which covers the thread's whole radial support, meets
//! the thread exactly in the same helix over the overlapping axial
//! interval. The result is rebuilt from the unchanged profile with the
//! trimmed length; a moved start keeps the helix by turning the chart about
//! the axis by `start / helical_rate`. Nothing is sampled or meshed.
use super::*;

/// Prism corners in the thread chart, if every prism axis is (to rounding)
/// one of the thread chart's axes.
fn aligned_box(
    prism: &Prism,
    prism_frame: PlaneFrame,
    thread: PlaneFrame,
) -> Option<([[f64; 2]; 3], f64)> {
    let size = prism.dimensions().map(|d| d.as_f64());
    let prism_axes = [prism_frame.u, prism_frame.v, prism_frame.normal()];
    let thread_axes = [thread.u, thread.v, thread.normal()];
    for axis in prism_axes {
        let aligned = thread_axes
            .iter()
            .any(|t| (dot(axis, *t).abs() - 1.).abs() <= 16. * f64::EPSILON);
        if !aligned {
            return None;
        }
    }
    let mut local = [[f64::INFINITY, f64::NEG_INFINITY]; 3];
    let mut magnitude = 1_f64;
    for corner in 0..8 {
        let offset: [f64; 3] =
            std::array::from_fn(|k| if corner >> k & 1 == 1 { size[k] } else { 0. });
        let world: [f64; 3] = std::array::from_fn(|i| {
            prism_frame.origin[i] + (0..3).map(|k| prism_axes[k][i] * offset[k]).sum::<f64>()
        });
        let delta: [f64; 3] = std::array::from_fn(|i| world[i] - thread.origin[i]);
        for (j, axis) in thread_axes.iter().enumerate() {
            let value = dot(delta, *axis);
            local[j][0] = local[j][0].min(value);
            local[j][1] = local[j][1].max(value);
        }
        magnitude = world
            .into_iter()
            .chain(thread.origin)
            .map(f64::abs)
            .fold(magnitude, f64::max);
    }
    // Axis misalignment of 16 ulps over the box, plus the evaluation.
    let pad = 128. * f64::EPSILON * magnitude;
    local
        .iter()
        .flatten()
        .all(|v| v.is_finite())
        .then_some((local, pad))
}

impl Solid {
    pub(in crate::geometry) fn rounded_thread_axial_clip(
        &self,
        other: &Self,
        operation: BooleanOperation,
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        if operation != BooleanOperation::Intersect
            || policy.boolean_fuzzy_mm != 0.
            || self.policy != policy
            || other.policy != policy
        {
            return None;
        }
        let ((body, frame), probe) = match self.rounded_thread_body() {
            Some(found) => (found, other),
            None => (other.rounded_thread_body()?, self),
        };
        if body.import_identity.is_some() {
            return None;
        }
        let (prism, prism_frame) = probe.rust_prism_frame()?;
        if prism.rounded_end().is_some() || prism.inner_loops().len() != 0 {
            return None;
        }
        let ([x, y, z], pad) = aligned_box(prism, prism_frame, frame)?;
        let spec = body.geometry.spec();
        let radius = (spec.major_radius_mm * (1. + 16. * f64::EPSILON)).next_up() + pad;
        if !(x[0] <= -radius && x[1] >= radius && y[0] <= -radius && y[1] >= radius) {
            return None;
        }
        let length = spec.length_mm;
        let start = z[0].max(0.);
        let end = z[1].min(length);
        if end <= start + pad {
            // Disjoint or only touching along the axis: no common interior.
            return Some(Ok(None));
        }
        let keeps_start = z[0] <= pad;
        let keeps_end = z[1] >= length - pad;
        if keeps_start && keeps_end {
            return Some(self.share_operand(other, probe).map(Some));
        }
        let (start, end) = (
            if keeps_start { 0. } else { start },
            if keeps_end { length } else { end },
        );
        let angle = start / body.geometry.helical_rate();
        let (sin, cos) = angle.portable_sin_cos();
        let normal = frame.normal();
        let trimmed = PlaneFrame {
            origin: std::array::from_fn(|i| frame.origin[i] + normal[i] * start),
            u: std::array::from_fn(|i| frame.u[i] * cos + frame.v[i] * sin),
            v: std::array::from_fn(|i| frame.v[i] * cos - frame.u[i] * sin),
        };
        Some(
            Self::rust_rounded_thread(
                multicad_rust_kernel::rounded_thread::Spec {
                    length_mm: end - start,
                    ..spec
                },
                body.feature,
                trimmed,
                policy,
            )
            .map(Some),
        )
    }

    /// The thread operand itself, unchanged, whichever side it is on.
    fn share_operand(&self, other: &Self, probe: &Self) -> Result<Self, GeometryError> {
        if std::ptr::eq(probe, other) {
            self.share()
        } else {
            other.share()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use multicad_rust_kernel::rounded_thread::Spec;

    fn spec() -> Spec {
        let angle = 27.5_f64.to_radians();
        Spec {
            major_radius_mm: 5.,
            pitch_mm: 1.,
            length_mm: 4.,
            depth_mm: 1. / (2. * angle.portable_tan()) - 0.274 * (1. / angle.portable_sin() - 1.),
            crest_radius_mm: 0.137,
            root_radius_mm: 0.137,
            flank_angle_deg: 27.5,
            starts: 1,
            left_hand: false,
        }
    }

    fn slab(z: [f64; 2], half: f64) -> Solid {
        Solid::box_solid(
            [2. * half, 2. * half, z[1] - z[0]],
            GeometryPolicy::default(),
        )
        .unwrap()
        .transformed(RigidTransform {
            translation_mm: [-half, -half, z[0]],
            ..Default::default()
        })
        .unwrap()
    }

    #[test]
    fn slabs_trim_the_same_helix_at_either_end() {
        let p = GeometryPolicy::default();
        let thread =
            Solid::rust_rounded_thread(spec(), EntityId::new(), PlaneFrame::default(), p).unwrap();
        let (geometry, _) = thread.rounded_thread_body().unwrap();
        let full = &geometry.geometry;
        for z in [[-1., 2.], [1.25, 9.], [0.5, 3.]] {
            let clipped = thread
                .boolean(&slab(z, 8.), BooleanOperation::Intersect, p)
                .unwrap()
                .unwrap();
            let bounds = clipped.bounds().unwrap();
            let (trimmed, trimmed_frame) = clipped.rounded_thread_body().unwrap();
            assert!(
                (bounds.min[2] - z[0].max(0.)).abs() < 1e-9,
                "{z:?} {bounds:?}"
            );
            assert!(
                (bounds.max[2] - z[1].min(4.)).abs() < 1e-9,
                "{z:?} {bounds:?}"
            );
            // The trimmed material agrees with the original helix pointwise.
            for k in 0..64 {
                let theta = k as f64 * 0.37;
                let zz =
                    z[0].max(0.) + 0.05 + (k as f64 / 64.) * (z[1].min(4.) - z[0].max(0.) - 0.1);
                let r = full.radius(theta, zz);
                for (offset, inside) in [(-1e-3, true), (1e-3, false)] {
                    let point = [
                        (r + offset) * theta.portable_cos(),
                        (r + offset) * theta.portable_sin(),
                        zz,
                    ];
                    assert_eq!(
                        trimmed.contains(trimmed_frame, point),
                        inside,
                        "{z:?} {point:?}"
                    );
                }
            }
            // Volume is proportional to the kept length of whole pitches.
            let kept = z[1].min(4.) - z[0].max(0.);
            assert!(
                (clipped.volume_mm3().unwrap() / thread.volume_mm3().unwrap() - kept / 4.).abs()
                    < 1e-6
            );
        }
        // A slab enclosing the whole thread returns it unchanged; a narrow
        // slab or one beyond the caps is not a trim.
        let whole = slab([-1., 5.], 8.)
            .boolean(&thread, BooleanOperation::Intersect, p)
            .unwrap()
            .unwrap();
        assert!((whole.volume_mm3().unwrap() - thread.volume_mm3().unwrap()).abs() < 1e-9);
        assert!(
            thread
                .rounded_thread_axial_clip(&slab([1., 2.], 4.), BooleanOperation::Intersect, p)
                .is_none()
        );
        assert!(
            thread
                .rounded_thread_axial_clip(&slab([5., 6.], 8.), BooleanOperation::Intersect, p)
                .unwrap()
                .unwrap()
                .is_none()
        );
        assert!(
            thread
                .rounded_thread_axial_clip(&slab([1., 2.], 8.), BooleanOperation::Subtract, p)
                .is_none()
        );
    }
}
