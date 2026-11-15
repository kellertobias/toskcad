//! Measured final insert pose; no hardware or stock is moved by the query.
use super::*;
impl Solid {
    pub fn rust_inset_wood_insert_contact(
        &self,
        insert: &Self,
        requirements: Requirements,
        maximum_inset_mm: f64,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<source::InsetInsertContact, GeometryError> {
        self.policy.validate()?;
        if self.policy != insert.policy
            || self.policy.boolean_fuzzy_mm != 0.
            || requirements.linear_tolerance_mm != self.policy.linear_tolerance_mm
            || requirements.minimum_feature_mm != self.policy.minimum_feature_mm
        {
            return Err(input(
                "Inset contact requires equal actual source geometry policies",
            ));
        }
        let unsupported = || {
            error(
                "Inset wood contact",
                GeometryErrorKind::Unsupported,
                "Inset contact requires original paired insert and one actual blind pilot",
            )
        };
        let (body, placed) = self.orthogonal_body().ok_or_else(unsupported)?;
        let (material, pilot) = body.wood_contact_source().ok_or_else(unsupported)?;
        let insert = insert.rust_wood_insert_data().ok_or_else(unsupported)?;
        source::certify_inset_insert(
            Stock {
                material,
                pilot,
                frame: frame(placed),
            },
            &insert.source,
            frame(insert.frame),
            requirements,
            maximum_inset_mm,
            &|| interrupt.is_some_and(|f| f()),
        )
        .map_err(source_error)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::RigidTransform;
    #[test]
    fn actual_inset_source_policy_and_cancellation_are_not_metadata_approvals() {
        let (stock, _) = super::super::clearance::tests::fixture(4.);
        let (_, insert) =
            super::super::mating::tests::fixture(false, 8., 0., PlaneFrame::default());
        let insert = insert
            .transformed(RigidTransform {
                rotation: [[1., 0., 0.], [0., -1., 0.], [0., 0., -1.]],
                translation_mm: [10., 10., 20. - 1e-9],
            })
            .unwrap();
        let p = GeometryPolicy::default();
        let req = Requirements {
            minimum_wall_mm: 1.,
            minimum_floor_mm: 1.,
            minimum_engagement_mm: 2.,
            tip_clearance_mm: 1.,
            linear_tolerance_mm: p.linear_tolerance_mm,
            minimum_feature_mm: p.minimum_feature_mm,
        };
        let proof = stock
            .rust_inset_wood_insert_contact(&insert, req, 1e-7, None)
            .unwrap();
        assert!(proof.inset_mm[0] > 0.);
        assert_eq!(
            stock
                .rust_inset_wood_insert_contact(&insert, req, 1e-7, Some(&|| true))
                .err()
                .unwrap()
                .kind,
            GeometryErrorKind::Cancelled
        );
        let plain = Solid::box_solid([20.; 3], p).unwrap();
        for (wood, insert) in [(&plain, &insert), (&stock, &plain), (&insert, &insert)] {
            assert_eq!(
                wood.rust_inset_wood_insert_contact(insert, req, 1e-7, None)
                    .err()
                    .unwrap()
                    .kind,
                GeometryErrorKind::Unsupported
            );
        }
        let mut bad = req;
        bad.minimum_feature_mm *= 2.;
        assert_eq!(
            stock
                .rust_inset_wood_insert_contact(&insert, bad, 1e-7, None)
                .err()
                .unwrap()
                .kind,
            GeometryErrorKind::Input
        );
        assert!(
            stock
                .rust_inset_wood_insert_contact(&insert, req, 1e-10, None)
                .is_err()
        );
    }
}
