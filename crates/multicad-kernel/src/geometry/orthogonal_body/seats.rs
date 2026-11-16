//! Retained orthogonal stock with independently qualified finite cone seats.
use super::*;
use multicad_rust_kernel::orthogonal_bore_seats::{self as source, Seat, Seats};
use std::f64::consts::{PI, TAU};
mod receiving;
mod step;
#[cfg(test)]
mod tests;
mod topology;
mod transport;
pub use transport::OrthogonalCountersinkMeasures;
#[derive(Clone, Debug)]
pub struct OrthogonalCountersink {
    pub geometry: Seat,
    pub origin: Origin,
}
pub(crate) struct OrthogonalSeatBody {
    pub(in crate::geometry) feature: EntityId,
    original: Rc<OrthogonalBody>,
    pub(in crate::geometry) source: Seats,
    origins: Vec<Origin>,
    qualified_step: std::cell::RefCell<
        Option<(
            GeometryPolicy,
            PlaneFrame,
            Rc<multicad_rust_kernel::step::orthogonal_seats::QualifiedOrthogonalSeatsStep>,
        )>,
    >,
}
fn failure(cause: source::Error) -> GeometryError {
    error(
        "Orthogonal countersinks",
        match cause {
            source::Error::Invalid => GeometryErrorKind::Input,
            source::Error::Unsupported => GeometryErrorKind::Unsupported,
            source::Error::Limit => GeometryErrorKind::Limit,
            source::Error::Cancelled => GeometryErrorKind::Cancelled,
        },
        format!("Original finite stock/seat proof failed: {cause:?}"),
    )
}
impl Solid {
    /// Author an already admitted local cell source with stable feature names.
    /// No old placed body is reconstructed, projected or snapped by this API.
    pub fn rust_orthogonal_stock_named(
        feature: EntityId,
        geometry: OrthogonalSolid,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        policy.validate()?;
        let frame = PlaneFrame::default();
        if feature.is_nil() {
            return Err(input("Orthogonal stock requires a feature owner"));
        }
        let local = geometry.bounds();
        frame.checked_metric_bounds(
            Bounds {
                min: local.min,
                max: local.max,
            },
            policy,
        )?;
        let body = OrthogonalBody::build(feature, geometry, None)?;
        let candidate = Self {
            body: SolidBody::Orthogonal(Rc::new(body), frame),
            policy,
            provenance: vec![],
            rust_origins: vec![],
            retained_placement: None,
        };
        let b = candidate.bounds()?;
        policy.point(b.min)?;
        policy.point(b.max)?;
        Ok(candidate)
    }
    pub fn rust_countersink_orthogonal(
        &self,
        seats: &[OrthogonalCountersink],
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, GeometryError> {
        policy.validate()?;
        if policy != self.policy || policy.boolean_fuzzy_mm != 0. {
            return Err(input(
                "Countersinks require the original exact stock policy",
            ));
        }
        let unsupported = || {
            unsupported(
                "Countersinks require original orthogonal through-bores without other machining or import approximation",
            )
        };
        let SolidBody::Orthogonal(original, frame) = &self.body else {
            return Err(unsupported());
        };
        if original.imported_roundoff_mm.is_some()
            || !original.pockets.is_empty()
            || !original.slots.is_empty()
        {
            return Err(unsupported());
        }
        if seats.len() > 32 || seats.iter().any(|s| s.origin.role.len() > 1024) {
            return Err(limit("Countersinks exceed the finite owner/cut budget"));
        }
        if seats
            .iter()
            .map(|s| &s.origin)
            .collect::<BTreeSet<_>>()
            .len()
            != seats.len()
        {
            return Err(input(
                "Each finite countersink needs a distinct source owner",
            ));
        }
        let geometry: Vec<_> = seats.iter().map(|s| s.geometry).collect();
        if seats
            .iter()
            .any(|s| s.origin.source.is_nil() || s.origin.role.trim().is_empty())
        {
            return Err(input("Each countersink needs a physical source origin"));
        }
        for s in seats {
            policy.length(s.geometry.radius_mm)?;
            policy.length(s.geometry.depth_mm)?;
        }
        let source = Seats::new(
            Rc::clone(&original.geometry),
            Rc::clone(&original.through_bores),
            &geometry,
            source::Options {
                clearance_mm: policy.linear_tolerance_mm,
                minimum_feature_mm: policy.minimum_feature_mm,
                maximum_work: 16_000_000,
                maximum_bytes: 64 * 1024 * 1024,
            },
            &|| interrupt.is_some_and(|f| f()),
        )
        .map_err(failure)?;
        // Source sorting must preserve each caller's owner; index by the same
        // unique original bore/end identity instead of input list position.
        let origins = source
            .seats()
            .iter()
            .map(|s| {
                seats
                    .iter()
                    .find(|t| t.geometry.bore == s.bore && t.geometry.end == s.end)
                    .expect("qualified unique seat")
                    .origin
                    .clone()
            })
            .collect();
        let body = OrthogonalSeatBody {
            feature: original.feature,
            original: Rc::clone(original),
            source,
            origins,
            qualified_step: std::cell::RefCell::new(None),
        };
        body.validate_frame(*frame, policy)?;
        let candidate = Self {
            body: SolidBody::OrthogonalSeats(Rc::new(body), *frame),
            policy,
            provenance: vec![],
            rust_origins: vec![],
            retained_placement: None,
        };
        candidate.topology_snapshot(original.feature, 0)?;
        let b = candidate.bounds()?;
        policy.point(b.min)?;
        policy.point(b.max)?;
        Ok(candidate)
    }
    pub fn rust_orthogonal_countersink_measures(
        &self,
    ) -> Result<OrthogonalCountersinkMeasures, GeometryError> {
        let (body, frame) = self
            .orthogonal_seat_body()
            .ok_or_else(|| unsupported("Expected retained orthogonal countersinks"))?;
        body.measures(frame)
    }
    pub fn rust_orthogonal_countersink_source(&self) -> Option<(&Seats, PlaneFrame)> {
        self.orthogonal_seat_body().map(|(b, f)| (&b.source, f))
    }
}
impl OrthogonalSeatBody {
    pub(in crate::geometry) fn bounds(&self, frame: PlaneFrame) -> Result<Bounds, GeometryError> {
        self.world_bounds(frame)
    }
    pub(in crate::geometry) fn volume_mm3(&self, frame: PlaneFrame) -> Result<f64, GeometryError> {
        let v = self.world_volume_bounds(frame)?;
        Ok(v[0].midpoint(v[1]))
    }
    pub(in crate::geometry) fn mesh(
        &self,
        frame: PlaneFrame,
        options: multicad_rust_kernel::display::Tessellation,
    ) -> Result<Mesh, GeometryError> {
        self.world_mesh(frame, options)
    }
    pub(in crate::geometry) fn faces(&self, frame: PlaneFrame) -> Result<Vec<Face>, GeometryError> {
        let mut faces = self.original.faces(PlaneFrame::default());
        for (seat_index, seat) in self.source.seats().iter().enumerate() {
            let bore = self.source.bores()[seat.bore];
            let end = usize::from(seat.end);
            let center = bore.centers_mm[end];
            let normal = self.source.boundary().faces[bore.faces[end]].normal;
            let radial = seat.radius_mm - bore.radius_mm;
            let loss = PI * (seat.radius_mm - bore.radius_mm) * (seat.radius_mm + bore.radius_mm);
            let cap = &mut faces[bore.faces[end]];
            let area = cap.area_mm2 - loss;
            cap.center_mm = std::array::from_fn(|i| {
                (cap.center_mm[i] * cap.area_mm2 - center[i] * loss) / area
            });
            cap.area_mm2 = area;
            let wall = &mut faces[self.source.boundary().faces.len() + seat.bore];
            let removed = TAU * bore.radius_mm * seat.depth_mm;
            let removed_center =
                std::array::from_fn::<_, 3, _>(|i| center[i] - normal[i] * seat.depth_mm / 2.);
            let area = wall.area_mm2 - removed;
            wall.center_mm = std::array::from_fn(|i| {
                (wall.center_mm[i] * wall.area_mm2 - removed_center[i] * removed) / area
            });
            wall.area_mm2 = area;
            let cone_center = seat.depth_mm * (seat.radius_mm + 2. * bore.radius_mm)
                / (3. * (seat.radius_mm + bore.radius_mm));
            faces.push(Face {
                index: (self.source.boundary().faces.len() + self.source.bores().len() + seat_index)
                    as u32,
                plane: None,
                center_mm: std::array::from_fn(|i| center[i] - normal[i] * cone_center),
                area_mm2: PI * (seat.radius_mm + bore.radius_mm) * seat.depth_mm.hypot(radial),
            });
        }
        self.world_faces(faces, frame)
    }
}
