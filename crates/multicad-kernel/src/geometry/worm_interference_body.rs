//! The common material of a placed Rust worm and generated wheel whose
//! pair certificate proved a common interior point.
//!
//! The result keeps both exact placed members and its volume from
//! `worm_gear::interference`: section areas by Green's theorem over the
//! exact wheel pieces and the explicit worm flank and end graphs, integrated
//! across the wheel face width by adaptive Gauss-Kronrod quadrature to the
//! policy volume tolerance. Its bounds are the common box of the members.
//! The common boundary itself is not built, so meshes, faces and STEP are
//! typed refusals.
use super::worm_body::{WormBody, WormGeometry};
use super::*;
use multicad_rust_kernel::worm_gear::{
    clearance::Placement,
    interference::{self, Volume},
};

const OPERATION: &str = "Boolean";

pub(crate) struct WormInterferenceBody {
    worm: (Rc<WormBody>, PlaneFrame),
    wheel: (Rc<WormBody>, PlaneFrame),
    pub(super) feature: EntityId,
    volume: Volume,
}

fn placement(frame: PlaneFrame) -> Placement {
    Placement {
        origin_mm: frame.origin,
        axes: [frame.u, frame.v, frame.normal()],
    }
}

impl WormInterferenceBody {
    pub(super) fn volume_mm3(&self) -> f64 {
        self.volume.volume_mm3
    }

    pub(super) fn bounds(&self) -> Bounds {
        let [a, b] = [&self.worm, &self.wheel].map(|(body, frame)| body.bounds(*frame));
        Bounds {
            min: std::array::from_fn(|i| a.min[i].max(b.min[i])),
            max: std::array::from_fn(|i| a.max[i].min(b.max[i])),
        }
    }

    pub(super) fn transformed(&self, transform: RigidTransform) -> Self {
        Self {
            worm: (Rc::clone(&self.worm.0), transform.frame(self.worm.1)),
            wheel: (Rc::clone(&self.wheel.0), transform.frame(self.wheel.1)),
            feature: self.feature,
            volume: self.volume,
        }
    }

    pub(super) fn unsupported(what: &'static str) -> GeometryError {
        error(
            what,
            GeometryErrorKind::Unsupported,
            "The Rust worm interference solid keeps its exact members and volume but has no built boundary",
        )
    }
}

impl Solid {
    pub(crate) fn worm_interference_body(&self) -> Option<&WormInterferenceBody> {
        match &self.body {
            SolidBody::WormInterference(body) => Some(body),
            _ => None,
        }
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
            return Some(Err(error(
                OPERATION,
                GeometryErrorKind::InvalidResult,
                "A proved common interior point of the worm pair gave no interference volume",
            )));
        }
        let feature = self
            .worm_body()
            .map_or_else(EntityId::new, |(body, _)| body.feature);
        Some(Ok(Some(Self {
            body: SolidBody::WormInterference(Rc::new(WormInterferenceBody {
                worm,
                wheel,
                feature,
                volume,
            })),
            policy,
            provenance: vec![],
            rust_origins: vec![],
            retained_placement: None,
        })))
    }
}
