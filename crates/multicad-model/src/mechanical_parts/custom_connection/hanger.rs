//! Canonical dimensions create one retained continuous two-ended source.
use super::*;
use multicad_kernel::HelicalThread;
use multicad_kernel::rust_kernel::{
    authored_flat_thread::{AuthoredThread, SourceOptions},
    authored_hanger_bolt::{Hanger, Spec as HangerSpec},
    authored_metric_thread::{AuthoredMetricThread, Spec as MetricSpec, TransferOptions},
    helical_thread::{Spec, Thread},
};
fn spec(p: HelicalThread, root_radius_mm: Option<f64>) -> Spec {
    Spec {
        major_radius_mm: p.major_radius_mm,
        pitch_mm: p.pitch_mm,
        length_mm: p.length_mm,
        depth_mm: p.depth_mm,
        crest_width_mm: p.crest_width_mm,
        lower_flank_deg: p.lower_flank_deg,
        upper_flank_deg: p.upper_flank_deg,
        starts: p.starts,
        left_hand: p.left_hand,
        root_radius_mm,
    }
}
pub(super) fn build(
    feature: EntityId,
    wood: &ResolvedThread,
    metric: &ResolvedThread,
    plain_length_mm: f64,
    plain_radius_mm: f64,
    frame: Frame,
    policy: GeometryPolicy,
) -> Result<Solid, ModelError> {
    let (wood, wood_root, wood_frame) = crate::threads::rust_profile(feature, wood)?;
    let (metric, metric_root, metric_frame) = crate::threads::rust_profile(feature, metric)?;
    if wood_root.is_some() || metric_root.is_none() {
        return Err(error(
            Some(feature),
            "Rust hanger currently requires flat wood and canonical rounded metric thread sources",
        ));
    }
    if wood_frame
        != (PlaneFrame {
            origin: frame.origin,
            u: frame.u,
            v: frame.v,
        })
        || metric_frame.u != wood_frame.u
        || metric_frame.v != wood_frame.v
    {
        return Err(error(
            Some(feature),
            "Hanger threads must retain their common original phase frame",
        ));
    }
    let wood = AuthoredThread::new(
        spec(wood, None),
        SourceOptions {
            position_limit_mm: policy.linear_tolerance_mm / 32.,
            normal_angle_limit_rad: policy.linear_tolerance_mm
                / (32. * wood.major_radius_mm.max(1.)),
        },
    )
    .map_err(|e| {
        error(
            Some(feature),
            format!("Actual hanger wood source cannot be certified: {e:?}"),
        )
    })?;
    let new_metric = AuthoredMetricThread::new(MetricSpec {
        major_radius_mm: metric.major_radius_mm,
        pitch_mm: metric.pitch_mm,
        length_mm: metric.length_mm,
        crest_width_mm: metric.crest_width_mm,
        starts: metric.starts,
        left_hand: metric.left_hand,
    })
    .map_err(|e| {
        error(
            Some(feature),
            format!("Actual hanger metric source cannot be certified: {e:?}"),
        )
    })?;
    // The independently rounded legacy law is only an input compatibility
    // check. The new retained law owns its shared exact tangent expressions.
    let legacy = Thread::new(spec(metric, metric_root)).map_err(|e| {
        error(
            Some(feature),
            format!("Original hanger metric dimensions are invalid: {e:?}"),
        )
    })?;
    new_metric
        .certify_legacy(
            &legacy,
            TransferOptions {
                position_limit_mm: policy.linear_tolerance_mm / 32.,
                normal_angle_limit_rad: policy.linear_tolerance_mm
                    / (32. * metric.major_radius_mm.max(1.)),
            },
            None,
        )
        .map_err(|e| {
            error(
                Some(feature),
                format!("Hanger canonical metric transfer exceeds its precision: {e:?}"),
            )
        })?;
    let source = Hanger::new(
        HangerSpec {
            wood,
            metric: new_metric,
            plain_length_mm,
            plain_radius_mm,
        },
        None,
    )
    .map_err(|e| {
        error(
            Some(feature),
            format!("Actual connected hanger source cannot be certified: {e:?}"),
        )
    })?;
    // The canonical resolved start has floating construction arithmetic.
    // Bound it against the retained body's exact symbolic station; never
    // refit or separately place a component to that rounded point.
    let z = source.sections()[0].start_mm;
    let mut error_bound = 0_f64;
    for ((origin, resolved), normal) in frame
        .origin
        .into_iter()
        .zip(metric_frame.origin)
        .zip(wood_frame.normal())
    {
        let products = z.map(|v| normal * v);
        let lower = (origin + products[0].min(products[1]).next_down()).next_down();
        let upper = (origin + products[0].max(products[1]).next_up()).next_up();
        let delta = (resolved - lower)
            .abs()
            .max((resolved - upper).abs())
            .next_up();
        error_bound = (error_bound + delta).next_up();
    }
    if error_bound > policy.linear_tolerance_mm / 16. {
        return Err(error(
            Some(feature),
            "Hanger metric start exceeds the original source-placement precision",
        ));
    }
    Solid::rust_hanger(source, feature, wood_frame, policy, None)
        .map_err(|e| error(Some(feature), e))
}
