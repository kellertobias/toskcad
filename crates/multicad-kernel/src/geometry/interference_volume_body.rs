//! The common material of two placed curved Rust gears whose pair
//! certificate proved a common interior point.
//!
//! The result keeps both exact placed members and the volume of their
//! common material; its bounds are the common box of the members. The
//! common boundary itself is not built, so meshes, faces and STEP are typed
//! refusals.
//!
//! - A worm with its generated wheel takes its volume from
//!   `worm_gear::interference`: section areas by Green's theorem over the
//!   exact wheel pieces and the explicit worm flank and end graphs,
//!   integrated across the wheel face width by adaptive Gauss-Kronrod
//!   quadrature to the policy volume tolerance.
//! - Two helical or herringbone gears on parallel axes take theirs from
//!   `spur_gear::clearance::twisted_interference_volume`: the exact common
//!   region of the two turned transverse sections at every height,
//!   integrated by adaptive Gauss-Kronrod quadrature, together with a
//!   certified enclosure of the exact volume.
use super::herringbone_body::{HerringboneBody, clearance::aligned_chart};
use super::worm_body::{WormBody, WormGeometry};
use super::*;
use multicad_rust_kernel::spur_gear::clearance as spur;
use multicad_rust_kernel::worm_gear::{clearance::Placement, interference};

const OPERATION: &str = "Boolean";

type Member<T> = (Rc<T>, PlaneFrame);

enum Members {
    Worm {
        worm: Member<WormBody>,
        wheel: Member<WormBody>,
    },
    Twisted([Member<HerringboneBody>; 2]),
}

pub(crate) struct InterferenceVolumeBody {
    members: Members,
    pub(super) feature: EntityId,
    volume_mm3: f64,
    /// Certified bounds of the exact volume, where the provider has them.
    enclosure_mm3: Option<[f64; 2]>,
}

fn placement(frame: PlaneFrame) -> Placement {
    Placement {
        origin_mm: frame.origin,
        axes: [frame.u, frame.v, frame.normal()],
    }
}

impl InterferenceVolumeBody {
    pub(super) fn volume_mm3(&self) -> f64 {
        self.volume_mm3
    }

    pub(super) fn kind_name(&self) -> &'static str {
        match self.members {
            Members::Worm { .. } => "Worm interference",
            Members::Twisted(_) => "Twisted gear interference",
        }
    }

    pub(super) fn bounds(&self) -> Bounds {
        let [a, b] = match &self.members {
            Members::Worm { worm, wheel } => [worm, wheel].map(|(body, frame)| body.bounds(*frame)),
            Members::Twisted(members) => {
                [&members[0], &members[1]].map(|(body, frame)| body.bounds(*frame))
            }
        };
        Bounds {
            min: std::array::from_fn(|i| a.min[i].max(b.min[i])),
            max: std::array::from_fn(|i| a.max[i].min(b.max[i])),
        }
    }

    pub(super) fn transformed(&self, transform: RigidTransform) -> Self {
        fn moved<T>((body, frame): &Member<T>, transform: RigidTransform) -> Member<T> {
            (Rc::clone(body), transform.frame(*frame))
        }
        Self {
            members: match &self.members {
                Members::Worm { worm, wheel } => Members::Worm {
                    worm: moved(worm, transform),
                    wheel: moved(wheel, transform),
                },
                Members::Twisted([a, b]) => {
                    Members::Twisted([moved(a, transform), moved(b, transform)])
                }
            },
            feature: self.feature,
            volume_mm3: self.volume_mm3,
            enclosure_mm3: self.enclosure_mm3,
        }
    }

    pub(super) fn unsupported(what: &'static str) -> GeometryError {
        error(
            what,
            GeometryErrorKind::Unsupported,
            "The Rust gear interference solid keeps its exact members and volume but has no built boundary",
        )
    }

    fn solid(self, policy: GeometryPolicy) -> Solid {
        Solid {
            body: SolidBody::InterferenceVolume(Rc::new(self)),
            policy,
            provenance: vec![],
            rust_origins: vec![],
            retained_placement: None,
        }
    }
}

fn no_volume() -> GeometryError {
    error(
        OPERATION,
        GeometryErrorKind::InvalidResult,
        "A proved common interior point of the gear pair gave no interference volume",
    )
}

impl Solid {
    pub(crate) fn interference_volume_body(&self) -> Option<&InterferenceVolumeBody> {
        match &self.body {
            SolidBody::InterferenceVolume(body) => Some(body),
            _ => None,
        }
    }

    /// Certified lower and upper bounds of the exact volume of a gear
    /// interference solid whose provider encloses it; `None` otherwise.
    pub fn certified_volume_bounds_mm3(&self) -> Option<[f64; 2]> {
        self.interference_volume_body()?.enclosure_mm3
    }

    /// The common material of a placed Rust worm and wheel, in either
    /// order. `None` leaves other operands to the typed refusal.
    pub(in crate::geometry) fn worm_interference(
        &self,
        other: &Self,
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        let member = |solid: &Self| match &solid.body {
            SolidBody::Worm(body, frame) => Some((Rc::clone(body), *frame)),
            _ => None,
        };
        let (a, b) = (member(self)?, member(other)?);
        let (worm, wheel) = match (&a.0.geometry, &b.0.geometry) {
            (WormGeometry::Worm(_), WormGeometry::Wheel(_)) => (a, b),
            (WormGeometry::Wheel(_), WormGeometry::Worm(_)) => (b, a),
            _ => return None,
        };
        let (WormGeometry::Worm(w), WormGeometry::Wheel(g)) = (&worm.0.geometry, &wheel.0.geometry)
        else {
            return None;
        };
        let m = g.module_mm();
        let found = interference::interference_volume(
            w,
            placement(worm.1),
            g,
            placement(wheel.1),
            policy.linear_tolerance_mm * m * m,
            None,
        );
        let volume = match found {
            Ok(volume) => volume,
            Err(cause) => {
                let kind = match cause {
                    interference::Error::Invalid => GeometryErrorKind::Unsupported,
                    interference::Error::Limit => GeometryErrorKind::Limit,
                    interference::Error::Cancelled => GeometryErrorKind::Cancelled,
                };
                return Some(Err(error(
                    OPERATION,
                    kind,
                    format!("Rust worm interference volume: {cause:?}"),
                )));
            }
        };
        if volume.volume_mm3 <= 0. {
            return Some(Err(no_volume()));
        }
        let feature = self
            .worm_body()
            .map_or_else(EntityId::new, |(body, _)| body.feature);
        Some(Ok(Some(
            InterferenceVolumeBody {
                members: Members::Worm { worm, wheel },
                feature,
                volume_mm3: volume.volume_mm3,
                enclosure_mm3: None,
            }
            .solid(policy),
        )))
    }

    /// The common material of two placed helical or herringbone gears on
    /// parallel axes with aligned starts. `None` leaves gapped, hubbed or
    /// unresolved pairs to the typed refusal.
    pub(in crate::geometry) fn twisted_gear_interference(
        &self,
        other: &Self,
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        let member = |solid: &Self| match &solid.body {
            SolidBody::Herringbone(body, frame) => Some((Rc::clone(body), *frame)),
            _ => None,
        };
        let members = [member(self)?, member(other)?];
        let placements = aligned_chart(
            OPERATION,
            policy,
            [
                (&*members[0].0, members[0].1),
                (&*members[1].0, members[1].1),
            ],
        )
        .ok()?;
        let gears = [&members[0].0.geometry, &members[1].0.geometry];
        // The policy tolerance over a square module, as for a worm pair;
        // the tip diameter spans two modules more than the teeth.
        let module = 2. * gears[0].tip_radius_mm() / (f64::from(gears[0].teeth()) + 2.);
        let found = spur::twisted_interference_volume(
            gears,
            placements,
            policy.linear_tolerance_mm,
            policy.linear_tolerance_mm * module * module,
            None,
        );
        let volume = match found {
            Ok(Some(volume)) => volume,
            Ok(None) | Err(spur::Error::Invalid) => return None,
            Err(cause) => {
                let kind = match cause {
                    spur::Error::Limit => GeometryErrorKind::Limit,
                    _ => GeometryErrorKind::Cancelled,
                };
                return Some(Err(error(
                    OPERATION,
                    kind,
                    format!("Rust twisted gear interference volume: {cause:?}"),
                )));
            }
        };
        if volume.enclosure_mm3[1] <= 0. || volume.volume_mm3 <= 0. {
            return Some(Err(no_volume()));
        }
        let feature = members[0].0.feature;
        Some(Ok(Some(
            InterferenceVolumeBody {
                members: Members::Twisted(members),
                feature,
                volume_mm3: volume.volume_mm3,
                enclosure_mm3: Some(volume.enclosure_mm3),
            }
            .solid(policy),
        )))
    }
}
