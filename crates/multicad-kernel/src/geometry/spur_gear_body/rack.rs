//! Placed fixed-pose proofs for a Rust spur pinion and a plain Rust rack sweep.
use super::*;
use crate::rust_kernel::spur_gear::clearance::{self, ResultKind};

const OPERATION: &str = "rack pair interference";

/// A placed spur pinion and plain normal sweep, in that order.
pub(super) struct RackPair<'a> {
    pub pinion: (&'a SpurGearBody, PlaneFrame),
    pub rack: (&'a [normal_sweep::Edge], f64, PlaneFrame),
    pub roles: Box<dyn Fn(usize) -> &'a str + 'a>,
    pub rack_support: Support<'a>,
    pub swapped: bool,
}

impl Solid {
    pub(super) fn rack_pair<'a>(&'a self, other: &'a Self) -> Option<RackPair<'a>> {
        let split = |gear: &'a Self, bar: &'a Self, swapped| {
            let pinion = gear.spur_gear_body()?;
            let (body, frame) = bar.normal_sweep_body()?;
            let (edges, height, roles) = body.plain_prism()?;
            Some(RackPair {
                pinion,
                rack: (edges, height, frame),
                roles: Box::new(roles),
                rack_support: Box::new(move |direction| body.sweep_support(direction)),
                swapped,
            })
        };
        split(self, other, false).or_else(|| split(other, self, true))
    }

    /// Prove separation or a strict common-interior point for an actual
    /// placed Rust spur pinion and an unmodified Rust rack sweep sharing an
    /// axial direction, in either operand order.
    pub fn certify_rust_rack_pair(
        &self,
        other: &Self,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<ResultKind, GeometryError> {
        policy.validate()?;
        let unsupported = |message| error(OPERATION, GeometryErrorKind::Unsupported, message);
        if self.policy != policy || other.policy != policy || policy.boolean_fuzzy_mm != 0. {
            return Err(unsupported(
                "Rust rack interference requires a shared zero-fuzz policy",
            ));
        }
        let Some(pair) = self.rack_pair(other) else {
            return Err(unsupported(
                "Rust rack interference requires an actual spur pinion and an unmodified single-ring rack sweep",
            ));
        };
        let placements = pair.chart(policy)?.0;
        clearance::certify_rack(
            &pair.pinion.0.profile,
            pair.rack.0,
            placements,
            policy.linear_tolerance_mm,
            interrupt,
        )
        .map_err(|cause| {
            let kind = match cause {
                clearance::Error::Invalid => GeometryErrorKind::Unsupported,
                clearance::Error::Limit => GeometryErrorKind::Limit,
                clearance::Error::Cancelled => GeometryErrorKind::Cancelled,
            };
            error(
                OPERATION,
                kind,
                format!("Rust rack interference certificate: {cause:?}"),
            )
        })
    }
}

impl RackPair<'_> {
    /// The pinion's hub lies inside its root disk, so enclosing drift over
    /// its whole axial extent keeps a Clear band proof valid for the hub.
    pub(super) fn chart(
        &self,
        policy: GeometryPolicy,
    ) -> Result<([clearance::Placement; 2], [f64; 2]), GeometryError> {
        let profile = &self.pinion.0.profile;
        coaxial_chart(
            (OPERATION, "Rust rack interference"),
            policy,
            [self.pinion.1, self.rack.2],
            [
                Band {
                    height_mm: profile.height_mm(),
                    extent_mm: profile.axial_bounds(),
                    support: Box::new(|direction| profile.support(direction)),
                },
                Band {
                    height_mm: self.rack.1,
                    extent_mm: [0., self.rack.1],
                    support: Box::new(|direction| (self.rack_support)(direction)),
                },
            ],
        )
    }
}
