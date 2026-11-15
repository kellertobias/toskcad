//! Actual counterbored stock; no nominal radius/box replaces the retained body.
use super::*;
impl Solid {
    pub fn rust_wood_insert_clamp_contact(
        &self,
        screw: &Self,
        insert: &Self,
        minimum_wall_mm: f64,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<source::InsertClampContact, GeometryError> {
        self.policy.validate()?;
        if self.policy != screw.policy
            || self.policy != insert.policy
            || self.policy.boolean_fuzzy_mm != 0.
        {
            return Err(input("Clamp contact requires equal actual source policies"));
        }
        let unsupported = || {
            error(
                "Insert clamp contact",
                GeometryErrorKind::Unsupported,
                "Clamp contact requires an actual one-bore counterbored polygon, original regular screw and paired insert",
            )
        };
        let stock = self
            .polygon_prism()
            .filter(|p| !p.reverse_sweep && p.pockets.is_empty())
            .ok_or_else(unsupported)?;
        let cut = stock.bore_treatments.as_ref().ok_or_else(unsupported)?;
        let (screw, screw_frame) = screw
            .rust_machine_screw_geometry()
            .ok_or_else(unsupported)?;
        let insert = insert.rust_wood_insert_data().ok_or_else(unsupported)?;
        source::certify_insert_clamp(
            source::InsertClampSources {
                prism: &stock.prism,
                cut: &cut.geometry,
                clamp_frame: frame(stock.frame),
                screw,
                screw_frame: frame(screw_frame),
                insert: &insert.source,
                insert_frame: frame(insert.frame),
            },
            source::InsertClampRequirements {
                minimum_wall_mm,
                minimum_bearing_land_mm: self
                    .policy
                    .minimum_feature_mm
                    .max(16. * self.policy.linear_tolerance_mm),
                maximum_bearing_gap_mm: self.policy.linear_tolerance_mm,
                linear_tolerance_mm: self.policy.linear_tolerance_mm,
                minimum_feature_mm: self.policy.minimum_feature_mm,
            },
            &|| interrupt.is_some_and(|f| f()),
        )
        .map_err(source_error)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{CircularPocketTool, RigidTransform};
    use multicad_core::{document::EntityId, references::Origin};
    fn clamp() -> Solid {
        let policy = GeometryPolicy::default();
        Solid::box_solid([20., 20., 4.1], policy)
            .unwrap()
            .rust_connection_counterbore(
                &CircularPocketTool {
                    frame: PlaneFrame {
                        origin: [10., 10., 4.10001],
                        u: [1., 0., 0.],
                        v: [0., -1., 0.],
                    },
                    center_uv: [0.; 2],
                    radius_mm: 1.7,
                    displacement_mm: [0., 0., -4.10002],
                    origin: Origin {
                        source: EntityId::new(),
                        role: "actual.clamp".into(),
                    },
                },
                3.5,
                1.,
                1.,
                Origin {
                    source: EntityId::new(),
                    role: "actual.seat".into(),
                },
                policy,
            )
            .unwrap()
            .0
    }
    #[test]
    fn actual_native_clamp_source_pose_and_policy_guards_remain_independent() {
        for left in [false, true] {
            let (s, i) = super::super::mating::tests::fixture(left, 8., 0., PlaneFrame::default());
            let place = |solid: &Solid, z| {
                solid
                    .transformed(RigidTransform {
                        rotation: [[1., 0., 0.], [0., -1., 0.], [0., 0., -1.]],
                        translation_mm: [10., 10., z],
                    })
                    .unwrap()
            };
            let screw = place(&s, 3.1 + 1e-10);
            let insert = place(&i, -1e-10);
            let clamp = clamp();
            let proof = clamp
                .rust_wood_insert_clamp_contact(&screw, &insert, 1., None)
                .unwrap();
            assert!(proof.bearing_face_gap_mm[0] > 0. && proof.insert_to_clamp_gap_mm[0] > 0.);
            assert_eq!(
                clamp
                    .rust_wood_insert_clamp_contact(&screw, &insert, 1., Some(&|| true))
                    .err()
                    .unwrap()
                    .kind,
                GeometryErrorKind::Cancelled
            );
            let plain = Solid::box_solid([20.; 3], GeometryPolicy::default()).unwrap();
            for (stock, screw, insert) in [
                (&plain, &screw, &insert),
                (&clamp, &plain, &insert),
                (&clamp, &screw, &plain),
                (&screw, &screw, &insert),
            ] {
                assert_eq!(
                    stock
                        .rust_wood_insert_clamp_contact(screw, insert, 1., None)
                        .err()
                        .unwrap()
                        .kind,
                    GeometryErrorKind::Unsupported
                );
            }
            let mut policy = GeometryPolicy::default();
            policy.linear_tolerance_mm *= 2.;
            let other = Solid::box_solid([20.; 3], policy).unwrap();
            assert_eq!(
                clamp
                    .rust_wood_insert_clamp_contact(&screw, &other, 1., None)
                    .err()
                    .unwrap()
                    .kind,
                GeometryErrorKind::Input
            );
            assert!(
                clamp
                    .rust_wood_insert_clamp_contact(&place(&s, 3.1 - 1e-10), &insert, 1., None)
                    .is_err()
            );
            assert!(
                clamp
                    .rust_wood_insert_clamp_contact(&screw, &place(&i, 0.), 1., None)
                    .is_err()
            );
        }
    }
}
