//! A real rectangular-stock receiver, with no bounding-box substitution.
use super::*;
use axial_profile_stack::RectangularBlindCounterbore;

fn cardinal(axis: Point) -> bool {
    axis.iter().filter(|x| x.abs() == 1.).count() == 1
        && axis.iter().all(|x| *x == 0. || x.abs() == 1.)
}
fn exact_add(a: f64, b: f64) -> Option<f64> {
    let sum = a + b;
    let bb = sum - a;
    ((a - (sum - bb)) + (b - bb) == 0.).then_some(sum)
}
impl Solid {
    /// Cut one two-radius blind recess in an unmachined rectangular prism.
    /// Only exact cardinal source charts and exactly transferable dimensions
    /// are admitted here; other stock is never replaced by its bounds.
    pub fn rust_cut_rectangular_blind_counterbore(
        &self,
        spec: CylindricalBlindPocketSpec,
        entry: PlaneFrame,
        feature: EntityId,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, GeometryError> {
        let check = || {
            if interrupt.is_some_and(|f| f()) {
                Err(error(
                    "blind counterbore",
                    GeometryErrorKind::Cancelled,
                    "Receiving cut superseded",
                ))
            } else {
                Ok(())
            }
        };
        let unsupported = || {
            error(
                "blind counterbore",
                GeometryErrorKind::Unsupported,
                "A rectangular blind counterbore requires one unmachined box with an exact cardinal entry chart",
            )
        };
        check()?;
        policy.validate()?;
        entry.checked(policy)?;
        if self.policy != policy || policy.boolean_fuzzy_mm != 0. || feature.is_nil() {
            return Err(input(
                "Receiving cut requires a distinct identity and the stock's exact policy",
            ));
        }
        let (prism, placement) = self.rust_prism_frame().ok_or_else(unsupported)?;
        if prism.feature() == feature.uuid_bits() {
            return Err(input("Receiving and stock identities must differ"));
        }
        let axes = [entry.u, entry.v, entry.normal()];
        if prism.rounded_end().is_some()
            || prism.inner_loops().len() != 0
            || !axes
                .into_iter()
                .chain([placement.u, placement.v, placement.normal()])
                .all(cardinal)
        {
            return Err(unsupported());
        }
        for length in prism.dimensions() {
            let (n, d) = length.ratio();
            if !(d as u64).is_power_of_two() || n.unsigned_abs() > (1u64 << 53) {
                return Err(unsupported());
            }
        }
        let stock = rectangular_operand(self, entry).ok_or_else(unsupported)?;
        // Recheck every coordinate transfer, rather than accepting rounded
        // offsets from the general profile-union recognition helper.
        let delta: Point = std::array::from_fn(|i| placement.origin[i] - entry.origin[i]);
        for (i, value) in delta.iter().enumerate() {
            if exact_add(placement.origin[i], -entry.origin[i]) != Some(*value) {
                return Err(unsupported());
            }
        }
        for (length, axis) in
            prism
                .dimensions()
                .into_iter()
                .zip([placement.u, placement.v, placement.normal()])
        {
            let mapped = axes.map(|a| dot(axis, a));
            let coordinate = mapped
                .iter()
                .position(|v| v.abs() == 1.)
                .ok_or_else(unsupported)?;
            let start = dot(delta, axes[coordinate]);
            if exact_add(start, mapped[coordinate] * length.as_f64()).is_none() {
                return Err(unsupported());
            }
        }
        let Profile::Rectangle { min, max } = stock.profile else {
            return Err(unsupported());
        };
        if stock.span_mm[0] != 0. {
            return Err(unsupported());
        }
        let seat = spec
            .seat
            .ok_or_else(|| input("A counterbore requires its actual seat dimensions"))?;
        let geometry = ProfileStack::blind_counterbore(
            RectangularBlindCounterbore {
                min_mm: min,
                max_mm: max,
                height_mm: stock.span_mm[1],
                bore_radius_mm: spec.radius_mm,
                bore_depth_mm: spec.depth_mm,
                seat_radius_mm: seat.radius_mm,
                seat_depth_mm: seat.depth_mm,
                minimum_wall_mm: spec.minimum_wall_mm,
                minimum_floor_mm: spec.minimum_floor_mm,
            },
            policy.linear_tolerance_mm,
            policy.minimum_feature_mm,
        )
        .map_err(|cause| {
            error(
                "blind counterbore",
                if cause == axial_profile_stack::Error::Invalid {
                    GeometryErrorKind::Input
                } else {
                    GeometryErrorKind::Unsupported
                },
                "Counterbore dimensions must preserve the actual stock wall, shoulder and floor",
            )
        })?;
        let original = self.rust_topology().ok_or_else(unsupported)??;
        let mut inherited_frames = Vec::with_capacity(geometry.faces().len());
        let origins = geometry
            .faces()
            .iter()
            .enumerate()
            .map(|(i, face)| {
                let mut ancestry = BTreeSet::new();
                let mut anchors = Vec::new();
                if let Support::Plane { normal } = face.support
                    && face.sources[0] <= 5
                {
                    let center = position(entry, face.center_mm);
                    let normal = vector(entry, normal);
                    for original in &original.elements {
                        if let Some(AnalyticGeometry::Plane { frame }) = original.geometry
                            && dot(normal, frame.normal()) == 1.
                            && dot(
                                std::array::from_fn(|j| center[j] - frame.origin_mm[j]),
                                normal,
                            ) == 0.
                        {
                            ancestry.extend(original.origins.iter().cloned());
                            anchors.extend(original.frames.iter().map(|a| AnchoredFrame {
                                origin: a.origin.clone(),
                                frame: ReferenceFrame {
                                    origin_mm: local(entry, a.frame.origin_mm),
                                    u: axes.map(|axis| dot(axis, a.frame.u)),
                                    v: axes.map(|axis| dot(axis, a.frame.v)),
                                },
                            }));
                        }
                    }
                }
                inherited_frames.push(anchors);
                std::iter::once(Origin {
                    source: feature,
                    role: format!("receiving.face.{i}"),
                })
                .chain(ancestry)
                .collect()
            })
            .collect();
        check()?;
        AxialStackBody::place_with_frames(
            geometry,
            entry,
            feature,
            origins,
            inherited_frames,
            policy,
        )
    }
}

#[cfg(test)]
mod tests;
