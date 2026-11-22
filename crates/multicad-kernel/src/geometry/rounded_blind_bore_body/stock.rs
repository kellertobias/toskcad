//! Cap-entry placement of the original full cutter into exact stock.
use super::*;
impl Solid {
    pub(in crate::geometry) fn cut_rounded_blind_stock(
        &self,
        spec: Spec,
        minimum_wall_mm: f64,
        feature: EntityId,
        requested: PlaneFrame,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        let SolidBody::AxialCylinder(cylinder, stock, _) = &self.body else {
            return Err(failure(
                GeometryErrorKind::Unsupported,
                "Blind rounded receiver needs one actual cylindrical stock",
            ));
        };
        if self.policy != policy || policy.boolean_fuzzy_mm != 0. {
            return Err(failure(
                GeometryErrorKind::Unsupported,
                "Blind rounded stock policy differs",
            ));
        }
        if feature.uuid_bits() == cylinder.feature()
            || !minimum_wall_mm.is_finite()
            || minimum_wall_mm < 0.
        {
            return Err(failure(
                GeometryErrorKind::Input,
                "Blind stock identity or minimum wall is invalid",
            ));
        }
        // This source preserves the complete original cutter and exact stock;
        // neither height is substituted into the other source's definition.
        let thread = multicad_rust_kernel::rounded_thread::Thread::new(spec).map_err(|e| {
            failure(
                GeometryErrorKind::Unsupported,
                format!("Original rounded cutter refused: {e:?}"),
            )
        })?;
        let bore =
            BlindBore::new(thread, cylinder.circle().radius(), cylinder.height()).map_err(|e| {
                failure(
                    GeometryErrorKind::Unsupported,
                    format!("Actual finite blind closure refused: {e:?}"),
                )
            })?;
        let radius = cylinder.circle().radius().as_f64();
        let height = cylinder.height().as_f64();
        let wall =
            (bore.stock_radius_bounds_mm()[0] - bore.cutter_radial_bounds_mm()[1]).next_down();
        if wall < minimum_wall_mm.max(policy.minimum_feature_mm) {
            return Err(failure(
                GeometryErrorKind::Input,
                "Actual rounded cutter leaves insufficient stock wall",
            ));
        }
        let alignment = dot(requested.normal(), stock.normal());
        if (alignment.abs() - 1.).abs() > 1e-10 {
            return Err(failure(
                GeometryErrorKind::Unsupported,
                "Blind rounded cutter must enter a coaxial stock cap",
            ));
        }
        let reverse = alignment < 0.;
        let origin = position(*stock, [0., 0., if reverse { height } else { 0. }]);
        let chart = helical_body::stock_frame::Chart::new(*stock);
        let [x, y, _] = chart.vector(requested.u);
        let norm = x.portable_hypot(y);
        let [c, s] = [x / norm, y / norm];
        let sign = if reverse { -1. } else { 1. };
        let placed = PlaneFrame {
            origin,
            u: std::array::from_fn(|i| c * stock.u[i] + s * stock.v[i]),
            v: std::array::from_fn(|i| sign * (-s * stock.u[i] + c * stock.v[i])),
        };
        helical_body::stock_frame::certify_tool_domain(
            requested,
            placed,
            bore.cutter_radial_bounds_mm()[1],
            spec.length_mm,
            policy,
        )?;
        helical_body::stock_frame::certify_circular_stock(
            *stock,
            placed,
            [c, s],
            radius,
            height,
            reverse,
            policy,
        )?;
        helical_body::stock_frame::store(
            placed,
            *stock,
            Bounds {
                min: [-radius, -radius, 0.],
                max: [radius, radius, height],
            },
            policy,
        )?;
        let mut result = Self::rust_rounded_blind_bore(
            spec,
            cylinder.circle().radius(),
            cylinder.height(),
            feature,
            placed,
            policy,
        )?;
        let inherited = self.topology_snapshot(feature, 0)?;
        helical_body::stock_frame::store_anchors(placed, inherited.clone(), policy)?;
        let SolidBody::RoundedBlindBore(body, _) = &mut result.body else {
            unreachable!("blind source constructor")
        };
        Rc::get_mut(body)
            .expect("unshared blind stock result")
            .inherited = Some((inherited, placed, reverse));
        Ok(result)
    }
}
