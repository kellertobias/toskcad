//! Actual retained screw/pilot contact; no generic Boolean exception.
mod insert;
use super::*;
use multicad_rust_kernel::wood_thread_contact::{
    self as source, Contact, Frame, Requirements, Stock,
};

fn frame(value: PlaneFrame) -> Frame {
    Frame {
        origin_mm: value.origin,
        axes: [value.u, value.v, value.normal()],
    }
}
fn source_error(cause: source::Error) -> GeometryError {
    let kind = match cause {
        source::Error::Invalid => GeometryErrorKind::Input,
        source::Error::Unsupported => GeometryErrorKind::Unsupported,
        source::Error::Limit => GeometryErrorKind::Limit,
        source::Error::Cancelled => GeometryErrorKind::Cancelled,
    };
    error(
        "Wood material proof",
        kind,
        format!("Actual source material proof failed: {cause:?}"),
    )
}
impl Solid {
    /// A strict outward source-envelope separation proof. Unknown source
    /// families and overlapping envelopes grant no disjointness authority.
    pub fn rust_wood_separated_from(
        &self,
        other: &Self,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Option<bool>, GeometryError> {
        let Some((screw, screw_frame)) = self.rust_machine_screw_geometry() else {
            return Ok(None);
        };
        let interrupted = || interrupt.is_some_and(|f| f());
        if let Some((body, other_frame)) = other.orthogonal_body() {
            let Some(material) = body.wood_envelope_source() else {
                return Ok(None);
            };
            source::separated_from_stock(
                material,
                frame(other_frame),
                screw,
                frame(screw_frame),
                &interrupted,
            )
            .map(Some)
            .map_err(source_error)
        } else if let Some((prism, other_frame)) = other.rust_prism_frame() {
            source::separated_from_rectangular_prism(
                prism,
                frame(other_frame),
                screw,
                frame(screw_frame),
                &interrupted,
            )
            .map_err(source_error)
        } else {
            Ok(None)
        }
    }

    /// A positive retained-cell intersection witness; None is inconclusive.
    pub fn rust_wood_obstruction_witness(
        &self,
        other: &Self,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Option<f64>, GeometryError> {
        let Some((screw, screw_frame)) = self.rust_machine_screw_geometry() else {
            return Ok(None);
        };
        let interrupted = || interrupt.is_some_and(|f| f());
        let result = if let Some((body, other_frame)) = other.orthogonal_body() {
            let Some(material) = body.wood_obstruction_source() else {
                return Ok(None);
            };
            source::obstruction_witness(
                material,
                frame(other_frame),
                screw,
                frame(screw_frame),
                &interrupted,
            )
        } else if let Some((prism, other_frame)) = other.rust_prism_frame() {
            source::rectangular_obstruction_witness(
                prism,
                frame(other_frame),
                screw,
                frame(screw_frame),
                &interrupted,
            )
        } else {
            return Ok(None);
        };
        result.map_err(source_error)
    }
    /// Prove the only intentional material overlap against the actual current
    /// drilled stock. Other component collisions remain the caller's concern.
    pub fn rust_direct_wood_contact(
        &self,
        screw: &Self,
        requirements: Requirements,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Contact, GeometryError> {
        let fail = |kind, message| error("Direct wood contact", kind, message);
        self.policy.validate()?;
        if self.policy != screw.policy
            || self.policy.boolean_fuzzy_mm != 0.
            || requirements.linear_tolerance_mm != self.policy.linear_tolerance_mm
            || requirements.minimum_feature_mm != self.policy.minimum_feature_mm
        {
            return Err(fail(
                GeometryErrorKind::Input,
                "Contact requires the actual source geometry policies",
            ));
        }
        let unsupported = || {
            fail(
                GeometryErrorKind::Unsupported,
                "Contact requires an actual regular headed screw and one retained orthogonal blind pilot",
            )
        };
        let (stock, stock_frame) = self.orthogonal_body().ok_or_else(unsupported)?;
        let (material, pilot) = stock.wood_contact_source().ok_or_else(unsupported)?;
        let (screw, screw_frame) = screw
            .rust_machine_screw_geometry()
            .ok_or_else(unsupported)?;
        stock_frame.checked(self.policy)?;
        screw_frame.checked(self.policy)?;
        source::certify(
            Stock {
                material,
                pilot,
                frame: frame(stock_frame),
            },
            screw,
            frame(screw_frame),
            requirements,
            &|| interrupt.is_some_and(|f| f()),
        )
        .map_err(source_error)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::HelicalThread;
    use multicad_core::references::Origin;
    use multicad_rust_kernel::headed_screw::{ScrewHead, ScrewSpec};
    pub(super) fn fixture() -> (Solid, Solid, Requirements) {
        let policy = GeometryPolicy::default();
        let stock = Solid::box_solid([20.; 3], policy).unwrap();
        let tool = CircularPocketTool {
            frame: PlaneFrame {
                origin: [10., 10., 20.001],
                u: [1., 0., 0.],
                v: [0., -1., 0.],
            },
            center_uv: [0.; 2],
            radius_mm: 2.,
            displacement_mm: [0., 0., -16.001],
            origin: Origin {
                source: EntityId::new(),
                role: "actual.pilot".into(),
            },
        };
        let recipient = stock
            .rust_connection_blind_bore(&tool, 1., 1., policy)
            .unwrap()
            .0;
        let screw = Solid::rust_machine_screw_at_thread_frame(
            HelicalThread {
                major_radius_mm: 3.,
                pitch_mm: 2.,
                length_mm: 12.,
                depth_mm: 1.,
                crest_width_mm: 0.4,
                lower_flank_deg: 30.,
                upper_flank_deg: 30.,
                radial_taper: 0.,
                starts: 1,
                left_hand: false,
            },
            None,
            ScrewSpec {
                head: ScrewHead::ExternalHex {
                    across_flats_mm: 10.,
                    height_mm: 4.,
                },
                shank_length_mm: 0.,
                shank_radius_mm: 3.,
            },
            EntityId::new(),
            PlaneFrame {
                origin: [10., 10., 20.],
                u: [1., 0., 0.],
                v: [0., -1., 0.],
            },
            policy,
        )
        .unwrap();
        (
            recipient,
            screw,
            Requirements {
                minimum_wall_mm: 1.,
                minimum_floor_mm: 1.,
                minimum_engagement_mm: 8.,
                tip_clearance_mm: 1.,
                linear_tolerance_mm: policy.linear_tolerance_mm,
                minimum_feature_mm: policy.minimum_feature_mm,
            },
        )
    }
    #[test]
    fn retained_sources_are_required_and_cancelled_queries_publish_no_claims() {
        let (recipient, screw, req) = fixture();
        let proof = recipient
            .rust_direct_wood_contact(&screw, req, None)
            .unwrap();
        assert!(proof.witness_volume_lower_mm3 > 0.);
        assert_eq!(
            recipient
                .rust_direct_wood_contact(&screw, req, Some(&|| true))
                .unwrap_err()
                .kind,
            GeometryErrorKind::Cancelled
        );
        let box_ = Solid::box_solid([2.; 3], Default::default()).unwrap();
        assert!(box_.rust_direct_wood_contact(&screw, req, None).is_err());
        assert!(
            recipient
                .rust_direct_wood_contact(&box_, req, None)
                .is_err()
        );
        let uncut = Solid::box_solid([20.; 3], Default::default())
            .unwrap()
            .orthogonal_stock()
            .unwrap();
        assert!(uncut.rust_direct_wood_contact(&screw, req, None).is_err());
        let stock = Solid::box_solid([20.; 3], Default::default()).unwrap();
        let through = stock
            .rust_connection_through_bore(
                &CircularPocketTool {
                    frame: PlaneFrame {
                        origin: [10., 10., 20.001],
                        u: [1., 0., 0.],
                        v: [0., -1., 0.],
                    },
                    center_uv: [0.; 2],
                    radius_mm: 2.,
                    displacement_mm: [0., 0., -20.002],
                    origin: Origin {
                        source: EntityId::new(),
                        role: "through-not-blind".into(),
                    },
                },
                1.,
                Default::default(),
            )
            .unwrap()
            .0;
        assert!(through.rust_direct_wood_contact(&screw, req, None).is_err());

        let mut bad = req;
        bad.minimum_wall_mm = f64::NAN;
        assert_eq!(
            recipient
                .rust_direct_wood_contact(&screw, bad, None)
                .unwrap_err()
                .kind,
            GeometryErrorKind::Input
        );
        bad = req;
        bad.minimum_engagement_mm = 13.;
        assert_eq!(
            recipient
                .rust_direct_wood_contact(&screw, bad, None)
                .unwrap_err()
                .kind,
            GeometryErrorKind::Unsupported
        );
        bad = req;
        bad.linear_tolerance_mm *= 2.;
        assert_eq!(
            recipient
                .rust_direct_wood_contact(&screw, bad, None)
                .unwrap_err()
                .kind,
            GeometryErrorKind::Input
        );
        assert!(
            recipient
                .rust_direct_wood_contact(&screw, req, None)
                .is_ok()
        );
        for (source, expected) in [
            (source::Error::Invalid, GeometryErrorKind::Input),
            (source::Error::Unsupported, GeometryErrorKind::Unsupported),
            (source::Error::Limit, GeometryErrorKind::Limit),
            (source::Error::Cancelled, GeometryErrorKind::Cancelled),
        ] {
            assert_eq!(source_error(source).kind, expected);
        }
    }
    #[test]
    fn actual_obstructions_have_positive_material_and_unknown_is_not_disjointness() {
        let (recipient, screw, _) = fixture();
        let obstacle = Solid::box_solid([2.; 3], Default::default())
            .unwrap()
            .transformed(RigidTransform {
                rotation: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
                translation_mm: [9., 9., 14.],
            })
            .unwrap();
        assert!(
            screw
                .rust_wood_obstruction_witness(&obstacle, None)
                .unwrap()
                .unwrap()
                > 0.
        );
        assert!(
            screw
                .rust_wood_obstruction_witness(&obstacle.orthogonal_stock().unwrap(), None)
                .unwrap()
                .is_some()
        );
        assert!(
            screw
                .rust_wood_obstruction_witness(&recipient, None)
                .unwrap()
                .is_none()
        );
        assert!(
            obstacle
                .rust_wood_obstruction_witness(&screw, None)
                .unwrap()
                .is_none()
        );
        assert_eq!(
            screw
                .rust_wood_obstruction_witness(&obstacle, Some(&|| true))
                .unwrap_err()
                .kind,
            GeometryErrorKind::Cancelled
        );
        let round = Solid::cylinder(PlaneFrame::default(), 2., 2., Default::default()).unwrap();
        assert!(
            screw
                .rust_wood_obstruction_witness(&round, None)
                .unwrap()
                .is_none()
        );
    }
    #[test]
    fn only_outward_owned_envelopes_can_publish_disjointness() {
        let (recipient, screw, _) = fixture();
        let far = Solid::box_solid([2.; 3], Default::default())
            .unwrap()
            .transformed(RigidTransform {
                rotation: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
                translation_mm: [30.; 3],
            })
            .unwrap();
        assert_eq!(
            screw.rust_wood_separated_from(&far, None).unwrap(),
            Some(true)
        );
        assert_eq!(
            screw
                .rust_wood_separated_from(&far.orthogonal_stock().unwrap(), None)
                .unwrap(),
            Some(true)
        );
        assert_eq!(
            screw.rust_wood_separated_from(&recipient, None).unwrap(),
            Some(false)
        );
        assert_eq!(far.rust_wood_separated_from(&screw, None).unwrap(), None);
        assert_eq!(
            screw
                .rust_wood_separated_from(&far, Some(&|| true))
                .unwrap_err()
                .kind,
            GeometryErrorKind::Cancelled
        );
        let round = Solid::cylinder(
            PlaneFrame {
                origin: [100.; 3],
                ..Default::default()
            },
            1.,
            2.,
            Default::default(),
        )
        .unwrap();
        assert_eq!(screw.rust_wood_separated_from(&round, None).unwrap(), None);
        let base = Solid::box_solid([4., 4., 1.], Default::default()).unwrap();
        let top = Solid::box_solid([2., 2., 1.], Default::default())
            .unwrap()
            .transformed(RigidTransform {
                rotation: [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
                translation_mm: [0., 0., 1.],
            })
            .unwrap();
        let stepped = base
            .rust_union_orthogonal_tools(&[top], Default::default())
            .unwrap()
            .unwrap();
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("actual-imported-stock.step");
        stepped.write_step(&path).unwrap();
        let text = std::fs::read_to_string(path).unwrap();
        let imported = Solid::import_rust_step(&text, EntityId::new(), Default::default()).unwrap();
        assert!(imported.orthogonal_body().is_some());
        assert_eq!(
            screw.rust_wood_separated_from(&imported, None).unwrap(),
            None
        );
    }
}
