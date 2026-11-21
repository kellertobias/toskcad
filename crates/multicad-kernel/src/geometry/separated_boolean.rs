//! Booleans of operands whose enclosing world boxes are strictly separated.
//!
//! No surface is intersected. A positive gap between two boxes that each
//! enclose an operand proves that the operands share no point, so the
//! regularized result is one of the unchanged operands, nothing, or both
//! operands as independent members of an exact compound. Boxes that only
//! touch still prove that the operands share no interior, which settles an
//! intersection or a difference but not a union (the contact would glue).
use super::*;

/// Bodies whose `bounds()` is the exact support function of their own
/// analytic source geometry (or the hull of exact planar corners), and
/// therefore encloses all material up to floating-point rounding.
fn has_enclosing_bounds(solid: &Solid) -> bool {
    match &solid.body {
        SolidBody::Prism(..)
        | SolidBody::PolygonPrism(..)
        | SolidBody::AxialCylinder(..)
        | SolidBody::Orthogonal(..)
        | SolidBody::Linear(..)
        | SolidBody::NormalSweep(..)
        | SolidBody::SpurGear(..)
        | SolidBody::InternalRing(..)
        | SolidBody::Sphere(..)
        | SolidBody::Torus(..) => true,
        SolidBody::RustCompound(compound) => compound.members().all(has_enclosing_bounds),
        _ => false,
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Separation {
    /// A gap wider than the tolerance: no common point.
    Apart,
    /// At most rounding-level overlap: no common interior.
    Touching,
}

/// Rounding margin proportional to the coordinates involved.
fn rounding_margin(a: Bounds, b: Bounds) -> f64 {
    64. * f64::EPSILON
        * a.min
            .into_iter()
            .chain(a.max)
            .chain(b.min)
            .chain(b.max)
            .map(f64::abs)
            .fold(1., f64::max)
}

fn separation(a: Bounds, b: Bounds, policy: GeometryPolicy) -> Option<Separation> {
    let rounding = rounding_margin(a, b);
    let gap = (0..3)
        .flat_map(|i| [b.min[i] - a.max[i], a.min[i] - b.max[i]])
        .filter(|gap| gap.is_finite())
        .fold(f64::NEG_INFINITY, f64::max);
    if gap > policy.linear_tolerance_mm + rounding {
        Some(Separation::Apart)
    } else if gap >= -rounding {
        Some(Separation::Touching)
    } else {
        None
    }
}

impl Solid {
    /// `None` leaves touching, overlapping or unenclosed operands to the
    /// remaining providers.
    pub(in crate::geometry) fn separated_boolean(
        &self,
        other: &Self,
        operation: BooleanOperation,
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        if policy.boolean_fuzzy_mm != 0.
            || self.policy != policy
            || other.policy != policy
            || !has_enclosing_bounds(self)
            || !has_enclosing_bounds(other)
        {
            return None;
        }
        let (a, b) = match (self.bounds(), other.bounds()) {
            (Ok(a), Ok(b)) => (a, b),
            (Err(error), _) | (_, Err(error)) => return Some(Err(error)),
        };
        let separation = separation(a, b, policy)?;
        if operation == BooleanOperation::Union && separation != Separation::Apart {
            return None;
        }
        Some(match operation {
            BooleanOperation::Intersect => Ok(None),
            BooleanOperation::Subtract => self.share().map(Some),
            BooleanOperation::Union => {
                let feature = match self.rust_topology() {
                    Some(Ok(topology)) => topology.owner,
                    Some(Err(error)) => return Some(Err(error)),
                    None => EntityId::new(),
                };
                RustCompound::assemble(&[self, other], feature, policy, None).map(Some)
            }
        })
    }
}

impl Solid {
    /// Orthogonal stock whose every stock cell box is separated from, or
    /// only touches, the enclosing box of the other operand. Cuts only
    /// remove stock, so the cells, in any placement, enclose the material;
    /// no common interior remains for an intersection or a difference.
    pub(in crate::geometry) fn cellwise_separated_boolean(
        &self,
        other: &Self,
        operation: BooleanOperation,
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        if operation == BooleanOperation::Union
            || policy.boolean_fuzzy_mm != 0.
            || self.policy != policy
            || other.policy != policy
        {
            return None;
        }
        let (stock, probe) = match (self.orthogonal_body(), other.orthogonal_body()) {
            (Some(_), _) => (self, other),
            (None, Some(_)) if operation == BooleanOperation::Intersect => (other, self),
            _ => return None,
        };
        if !has_enclosing_bounds(probe) {
            return None;
        }
        let (body, placement) = stock.orthogonal_body()?;
        let probe_bounds = match probe.bounds() {
            Ok(bounds) => bounds,
            Err(error) => return Some(Err(error)),
        };
        let cells = body.geometry.cells();
        if cells.is_empty() {
            return None;
        }
        for cell in cells.iter() {
            let corners = (0..8).map(|n: usize| {
                crate::rust_prism_topology::position(
                    placement,
                    std::array::from_fn(|i| {
                        if n >> i & 1 == 0 {
                            cell.min[i]
                        } else {
                            cell.max[i]
                        }
                    }),
                )
            });
            let mut bounds = Bounds {
                min: [f64::INFINITY; 3],
                max: [f64::NEG_INFINITY; 3],
            };
            for corner in corners {
                for (i, value) in corner.into_iter().enumerate() {
                    bounds.min[i] = bounds.min[i].min(value);
                    bounds.max[i] = bounds.max[i].max(value);
                }
            }
            separation(bounds, probe_bounds, policy)?;
        }
        Some(match operation {
            BooleanOperation::Intersect => Ok(None),
            _ => self.share().map(Some),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p() -> GeometryPolicy {
        GeometryPolicy::default()
    }

    fn bounds(min: [f64; 3], max: [f64; 3]) -> Bounds {
        Bounds { min, max }
    }

    fn cylinder(x: f64) -> Solid {
        Solid::cylinder(
            PlaneFrame {
                origin: [x, 0., 0.],
                ..Default::default()
            },
            3.,
            5.,
            p(),
        )
        .unwrap()
    }

    #[test]
    fn a_gap_requires_more_than_tolerance_on_some_axis() {
        use Separation::*;
        let a = bounds([0.; 3], [1.; 3]);
        let b = bounds([3., 0., 0.], [4., 1., 1.]);
        assert_eq!(separation(a, b, p()), Some(Apart));
        assert_eq!(separation(b, a, p()), Some(Apart));
        assert_eq!(
            separation(a, bounds([0., 0., -5.], [1., 1., -1.]), p()),
            Some(Apart)
        );
        assert_eq!(separation(a, bounds([0.5; 3], [2.; 3]), p()), None);
        // A sub-tolerance gap is only contact.
        let near = bounds([1. + 1.05e-7, 0., 0.], [2., 1., 1.]);
        assert_eq!(separation(a, near, p()), Some(Apart));
        // The rounding margin grows with the coordinates involved.
        let far = bounds([1e6; 3], [1e6 + 1.; 3]);
        let shifted = bounds(
            [1e6 + 1. + 1.05e-7, 1e6, 1e6],
            [1e6 + 2., 1e6 + 1., 1e6 + 1.],
        );
        assert_eq!(separation(far, shifted, p()), Some(Touching));
    }

    #[test]
    fn contact_is_distinguished_from_a_gap_and_from_overlap() {
        let a = bounds([0.; 3], [1.; 3]);
        let separated = |min, max| separation(a, bounds(min, max), p());
        assert_eq!(
            separated([2., 0., 0.], [3., 1., 1.]),
            Some(Separation::Apart)
        );
        assert_eq!(
            separated([0., 0., 1.], [1., 1., 2.]),
            Some(Separation::Touching)
        );
        assert_eq!(
            separated([0., 0., 1. + 5e-8], [1., 1., 2.]),
            Some(Separation::Touching)
        );
        assert_eq!(
            separated([0., 0., 1. - 1e-15], [1., 1., 2.]),
            Some(Separation::Touching)
        );
        assert_eq!(separated([0., 0., 1. - 1e-9], [1., 1., 2.]), None);
        assert_eq!(separated([f64::NAN; 3], [f64::NAN; 3]), None);
    }

    #[test]
    fn touching_blocks_intersect_to_nothing_but_do_not_glue() {
        let block = |z: f64| {
            Solid::box_solid([2., 3., 1.], p())
                .unwrap()
                .transformed(RigidTransform {
                    translation_mm: [0.5, 0., z],
                    ..Default::default()
                })
                .unwrap()
        };
        let disc = cylinder(1.);
        let resting = block(5.);
        // The block is rectilinear but the cylinder is not, so no other
        // provider answers before this one.
        assert!(
            disc.boolean(&resting, BooleanOperation::Intersect, p())
                .unwrap()
                .is_none()
        );
        let kept = disc
            .boolean(&resting, BooleanOperation::Subtract, p())
            .unwrap()
            .unwrap();
        assert!((kept.volume_mm3().unwrap() - disc.volume_mm3().unwrap()).abs() < 1e-9);
        assert!(
            disc.separated_boolean(&resting, BooleanOperation::Union, p())
                .is_none()
        );
        // Real overlap is never claimed.
        assert!(
            disc.separated_boolean(&block(4.5), BooleanOperation::Intersect, p())
                .is_none()
        );
    }

    #[test]
    fn overlapping_curved_operands_get_a_typed_refusal_and_far_spheres_do_not() {
        let sphere = |x: f64| {
            Solid::rust_sphere_named(
                2.,
                PlaneFrame {
                    origin: [x, 0., 0.],
                    ..Default::default()
                },
                EntityId::new(),
                p(),
            )
            .unwrap()
        };
        let error = sphere(0.)
            .boolean(&sphere(1.), BooleanOperation::Intersect, p())
            .err()
            .unwrap();
        assert_eq!(error.kind, GeometryErrorKind::Unsupported);
        assert_eq!(error.operation, "Boolean");
        assert!(
            error.message.contains("Sphere Intersect Sphere"),
            "{}",
            error.message
        );
        assert!(!error.message.contains("rectilinear"));
        // Separated spheres are settled exactly.
        let pair = sphere(0.)
            .boolean(&sphere(5.), BooleanOperation::Union, p())
            .unwrap()
            .unwrap();
        assert_eq!(pair.solid_count().unwrap(), 2);
    }

    #[test]
    fn separated_cylinders_form_each_exact_operation() {
        let [left, right] = [cylinder(0.), cylinder(10.)];
        let union = left
            .boolean(&right, BooleanOperation::Union, p())
            .unwrap()
            .unwrap();
        assert_eq!(union.solid_count().unwrap(), 2);
        let expected = 2. * std::f64::consts::PI * 9. * 5.;
        assert!((union.volume_mm3().unwrap() - expected).abs() < 1e-9);
        assert_eq!(union.bounds().unwrap().max[0], 13.);
        assert!(
            left.boolean(&right, BooleanOperation::Intersect, p())
                .unwrap()
                .is_none()
        );
        let difference = left
            .boolean(&right, BooleanOperation::Subtract, p())
            .unwrap()
            .unwrap();
        assert_eq!(difference.solid_count().unwrap(), 1);
        assert!((difference.volume_mm3().unwrap() - expected / 2.).abs() < 1e-9);
        // A compound of enclosed members is itself enclosed.
        let third = cylinder(-20.);
        let triple = union
            .boolean(&third, BooleanOperation::Union, p())
            .unwrap()
            .unwrap();
        assert_eq!(triple.solid_count().unwrap(), 3);
    }

    #[test]
    fn touching_unions_mismatched_or_fuzzy_operands_are_left_to_other_providers() {
        let left = cylinder(0.);
        assert!(
            left.separated_boolean(&cylinder(6.), BooleanOperation::Union, p())
                .is_none()
        );
        let fuzzy = GeometryPolicy {
            boolean_fuzzy_mm: 1e-6,
            ..p()
        };
        assert!(
            left.separated_boolean(&cylinder(10.), BooleanOperation::Union, fuzzy)
                .is_none()
        );
        let coarse = GeometryPolicy {
            linear_tolerance_mm: 1e-6,
            ..p()
        };
        assert!(
            left.separated_boolean(&cylinder(10.), BooleanOperation::Union, coarse)
                .is_none()
        );
    }
}
