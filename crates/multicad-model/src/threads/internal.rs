//! Finite cap-normal internal cuts retain the complete authored thread tool.
use super::*;
use multicad_core::references::{AnalyticGeometry, TopologyKind};
use multicad_kernel::{BlindThreadCut, GeometryErrorKind};

pub(super) fn build(
    feature: EntityId,
    family: ThreadFamily,
    thread: &ResolvedThread,
    target: &Solid,
    (profile, root, requested): (HelicalThread, Option<f64>, PlaneFrame),
    policy: GeometryPolicy,
) -> Result<Solid, ModelError> {
    let original = match target.rust_cut_internal_thread(
        profile,
        thread.minimum_wall_mm,
        feature,
        requested,
        policy,
    ) {
        Ok(solid) => return Ok(solid),
        Err(cause) if cause.kind != GeometryErrorKind::Unsupported => {
            return Err(error(Some(feature), cause));
        }
        Err(cause) => cause,
    };
    if family != ThreadFamily::IsoMetric || root.is_some() {
        return Err(error(Some(feature), original));
    }
    // There is no authored floor dimension in a ThreadDefinition. These are
    // intrinsic positive-material policy reserves, not new canonical values.
    let result = target
        .rust_cut_blind_internal_thread(
            BlindThreadCut {
                profile,
                entry: requested,
                receiving_radius_mm: profile.major_radius_mm,
                minimum_wall_mm: thread.minimum_wall_mm.max(policy.minimum_feature_mm),
                minimum_floor_mm: policy.minimum_feature_mm,
            },
            feature,
            policy,
        )
        .map_err(|cause| {
            error(
                Some(feature),
                if cause.kind == GeometryErrorKind::Unsupported {
                    original
                } else {
                    cause
                },
            )
        })?;
    // The cap's chart is emitted from the actual immutable thread placement.
    // Its outward V is reversed; restore it to compare the full original tool
    // domain, including phase, rather than trusting a centerline or a recipe.
    let topology = result
        .topology_snapshot(feature, 0)
        .map_err(|cause| error(Some(feature), cause))?;
    let cap = topology
        .elements
        .iter()
        .find(|e| e.key.kind == TopologyKind::Face && e.key.ordinal == 0)
        .and_then(|e| match &e.geometry {
            Some(AnalyticGeometry::Plane { frame }) => Some(*frame),
            _ => None,
        })
        .ok_or_else(|| {
            error(
                Some(feature),
                "Blind internal thread has no actual entry chart",
            )
        })?;
    let actual = PlaneFrame {
        origin: cap.origin_mm,
        u: cap.u,
        v: cap.v.map(|x| -x),
    };
    let distance = |a: [f64; 3], b: [f64; 3]| {
        (a[0] - b[0])
            .hypot(a[1] - b[1])
            .hypot(a[2] - b[2])
            .next_up()
    };
    let mut residual = distance(requested.origin, actual.origin);
    for ((a, b), extent) in [requested.u, requested.v, requested.normal()]
        .into_iter()
        .zip([actual.u, actual.v, actual.normal()])
        .zip([
            profile.major_radius_mm,
            profile.major_radius_mm,
            profile.length_mm,
        ])
    {
        residual = (residual + (distance(a, b) * extent).next_up()).next_up();
    }
    if !residual.is_finite() || residual > policy.linear_tolerance_mm / 8. {
        return Err(error(
            Some(feature),
            "Blind internal thread cannot preserve the full original tool frame and phase",
        ));
    }
    Ok(result)
}
