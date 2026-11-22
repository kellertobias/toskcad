//! Positive-clearance certificates for external spur and twisted pairs,
//! including profile-shifted members and spur hubs, internal ring/spur
//! pinion pairs and finite straight racks with a spur pinion. Shifted and
//! rack engagements also need actual flank proximity in every tooth band.
use super::*;
use multicad_core::model::GearFamily;
use multicad_kernel::rust_kernel::spur_gear::clearance::ResultKind;

#[derive(Clone, Copy, PartialEq)]
enum Mesh {
    External(GearFamily),
    Internal,
    Rack,
}

/// The engagement witness a pose must carry, with the authored total
/// thinning it must come within.
#[derive(Clone, Copy)]
struct Engagement {
    kind: &'static str,
    chevron: bool,
    clearance_mm: f64,
}

struct Pose<'a> {
    pair: EntityId,
    mesh: Mesh,
    engagement: Option<Engagement>,
    policy: GeometryPolicy,
    interrupt: Option<&'a dyn Fn() -> bool>,
}

impl Pose<'_> {
    fn check(
        &self,
        [a, b]: [&Solid; 2],
        drive: f64,
    ) -> Result<Option<GearEngagementWitness>, ModelError> {
        let Self {
            pair,
            mesh,
            engagement,
            policy,
            interrupt,
        } = *self;
        let at = drive.to_degrees();
        let result = match mesh {
            Mesh::External(GearFamily::Spur) => a.certify_rust_spur_pair(b, policy, interrupt),
            Mesh::External(_) => a.certify_rust_herringbone_pair(b, policy, interrupt),
            Mesh::Internal => a.certify_rust_internal_pair(b, policy, interrupt),
            Mesh::Rack => a.certify_rust_rack_pair(b, policy, interrupt),
        }
        .map_err(|cause| error(Some(pair), cause))?;
        let reason = match result {
            ResultKind::Clear => {
                let Some(engagement) = engagement else {
                    return Ok(None);
                };
                return witness(pair, [a, b], at, engagement, policy, interrupt).map(Some);
            }
            ResultKind::CommonInterior => "overlaps with a certified common-interior point",
            ResultKind::Unresolved => "does not have certified positive clearance",
        };
        Err(error(
            Some(pair),
            format!(
                "Rust gear pair {reason} at tested drive {at:.6} degrees; repair its phase or clearance"
            ),
        ))
    }
}

/// A Clear result also accepts separated bodies. Require actual flank
/// points within the authored total thinning in every tooth band.
fn witness(
    pair: EntityId,
    [a, b]: [&Solid; 2],
    at: f64,
    engagement: Engagement,
    policy: GeometryPolicy,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<GearEngagementWitness, ModelError> {
    let clearance = engagement.clearance_mm;
    let reach = clearance + policy.linear_tolerance_mm * 10.;
    let bands = a
        .rust_gear_flank_witnesses(b, policy, reach, interrupt)
        .map_err(|cause| error(Some(pair), cause))?;
    if bands.len() != if engagement.chevron { 2 } else { 1 } || bands.iter().any(Option::is_none) {
        return Err(error(
            Some(pair),
            format!(
                "Gear members have no Rust flank proximity in every tooth band within authored total thinning {clearance} mm at tested drive {at:.6} degrees"
            ),
        ));
    }
    let widest = bands
        .into_iter()
        .flatten()
        .max_by(|p, q| p.distance_mm.total_cmp(&q.distance_mm))
        .expect("checked bands");
    Ok(GearEngagementWitness {
        distance_mm: widest.distance_mm,
        kind: engagement.kind,
        faces: widest.faces,
        points_mm: widest.points_mm,
    })
}

fn refuse(pair: EntityId) -> ModelError {
    error(
        Some(pair),
        "Rust gear-pair interference requires external spur, helical, continuous herringbone or equal-gap double-helical members of one family, an internal ring with a hubless spur pinion, or a straight rack with a spur pinion",
    )
}

/// Classify the pair and the engagement witness each pose must carry.
fn classify(
    document: &Document,
    pair: EntityId,
    source: &GearPairDefinition,
    resolved: &model::ResolvedGearPair,
) -> Result<(Mesh, Option<Engagement>), ModelError> {
    let rotary = resolved.gears.each_ref().map(|gear| gear.rotary());
    let thinning = |gear: &model::ResolvedGear| {
        gear.nominal_tooth_thickness_mm - gear.physical_tooth_thickness_mm
    };
    if let Some(travel) = &resolved.rack_travel {
        let rack = resolved.gears[travel.rack_member.index()]
            .rack()
            .ok_or_else(|| refuse(pair))?;
        let pinion = rotary[travel.pinion_member.index()]
            .filter(|pinion| pinion.family == GearFamily::Spur)
            .ok_or_else(|| refuse(pair))?;
        return Ok((
            Mesh::Rack,
            Some(Engagement {
                kind: "rack_pinion_flanks",
                chevron: false,
                clearance_mm: rack.nominal_tooth_thickness_mm - rack.physical_tooth_thickness_mm
                    + thinning(pinion),
            }),
        ));
    }
    let [Some(first), Some(second)] = rotary else {
        return Err(refuse(pair));
    };
    let families = [first.family, second.family];
    if families.contains(&GearFamily::InternalRing) {
        let pinion = if families[0] == GearFamily::InternalRing {
            second
        } else {
            first
        };
        if pinion.family != GearFamily::Spur || pinion.hub.is_some() {
            return Err(refuse(pair));
        }
        return Ok((Mesh::Internal, None));
    }
    let family = families[0];
    if family != families[1]
        || !matches!(
            family,
            GearFamily::Spur
                | GearFamily::Helical
                | GearFamily::Herringbone
                | GearFamily::DoubleHelical
        )
    {
        return Err(refuse(pair));
    }
    let engagement = shifted_engagement(document, source.parameters.gears).then(|| Engagement {
        kind: match family {
            GearFamily::Spur => "shifted_spur_flank_faces",
            GearFamily::Helical => "shifted_helical_flank_faces",
            _ => "shifted_chevron_both_bands",
        },
        chevron: matches!(family, GearFamily::Herringbone | GearFamily::DoubleHelical),
        clearance_mm: thinning(first) + thinning(second),
    });
    Ok((Mesh::External(family), engagement))
}

fn scope(mesh: Mesh, shifted: bool) -> &'static str {
    match (mesh, shifted) {
        (Mesh::Rack, _) => {
            "Rust positive-clearance certificates for the current pose, nine sampled phases over one pinion tooth period and both conservative finite travel endpoints of a straight rack and spur pinion, each with actual flank points within the authored total thinning. Continuous motion, strength and manufacturing fit are outside this check."
        }
        (Mesh::Internal, _) => {
            "Rust positive-clearance certificates for the current pose and nine sampled phases over one tooth cycle of an unshifted internal ring and spur pinion, including nesting of the pinion inside the ring cavity. Continuous motion, strength and manufacturing fit are outside this check."
        }
        (Mesh::External(GearFamily::Spur), true) => {
            "Rust positive-clearance certificates for the current pose and nine sampled phases over one matched profile-shifted external spur tooth cycle, including any hubs, each with actual named flank faces within the authored total thinning. Continuous motion, strength and manufacturing fit are outside this check."
        }
        (Mesh::External(GearFamily::Spur), false) => {
            "Rust positive-clearance certificates for the current pose and nine sampled phases over one unshifted external spur tooth cycle, including any hubs. Continuous motion, strength and manufacturing fit are outside this check."
        }
        (Mesh::External(GearFamily::Helical), true) => {
            "Rust full-axial positive-clearance certificates for the current pose and nine sampled phases of one matched profile-shifted helical tooth cycle, each with actual named flank faces within the authored total thinning. Continuous motion, strength and manufacturing fit are outside this check."
        }
        (Mesh::External(GearFamily::Helical), false) => {
            "Rust full-axial positive-clearance certificates for the current pose and nine sampled phases of one unshifted helical tooth cycle. Continuous motion, strength and manufacturing fit are outside this check."
        }
        (Mesh::External(_), true) => {
            "Rust full-axial positive-clearance certificates for the current pose and nine sampled phases of one matched profile-shifted chevron tooth cycle, each with actual named flank faces within the authored total thinning in both bands. Continuous motion, strength and manufacturing fit are outside this check."
        }
        (Mesh::External(_), false) => {
            "Rust full-axial positive-clearance certificates for the current pose and nine sampled phases of one unshifted continuous herringbone or equal-gap double-helical tooth cycle. Continuous motion, strength and manufacturing fit are outside this check."
        }
    }
}

fn rotary(
    inputs: &ValidationInputs<'_>,
    pair: EntityId,
    source: &GearPairDefinition,
    resolved: &model::ResolvedGearPair,
    (phases, reused): (Vec<f64>, bool),
) -> Result<GearPairCheck, ModelError> {
    let ValidationInputs {
        document,
        evaluated,
        solids,
        policy,
        interrupt,
        ..
    } = *inputs;
    let (mesh, engagement) = classify(document, pair, source, resolved)?;
    // A ring/pinion pair is certified with the ring as the first operand.
    let ring_second =
        mesh == Mesh::Internal && resolved.gears[1].family() == GearFamily::InternalRing;
    let order = if ring_second { [1, 0] } else { [0, 1] };
    let pose = Pose {
        pair,
        mesh,
        engagement,
        policy,
        interrupt,
    };
    let start = std::time::Instant::now();
    if !reused {
        let local = order.map(|index| solids.get(&source.parameters.gears[index]));
        let [Some(a), Some(b)] = local else {
            return Err(error(
                Some(pair),
                "Rust gear-pair interference requires both exact member solids",
            ));
        };
        for &drive in &phases {
            checkpoint(interrupt)?;
            let sample = sample_document(document, pair, source, drive)?;
            let [pa, pb] = order.map(|index| {
                model::resolve_instance_pose(&sample, source.outputs.instances[index], evaluated)
                    .map_err(|e| error(Some(pair), e))
            });
            pose.check([&placed(a, pa?, pair)?, &placed(b, pb?, pair)?], drive)?;
        }
    }
    let sampled_cycle_ms = if reused {
        0.
    } else {
        start.elapsed().as_secs_f64() * 1000.
    };
    checkpoint(interrupt)?;
    let current = source
        .parameters
        .drive_angle
        .value_in(multicad_core::parameters::Unit::Radian, Some(evaluated))
        .map_err(|e| error(Some(pair), e))?;
    let start = std::time::Instant::now();
    let current_engagement = pose.check(
        [
            inputs.current_member(pair, source, order[0])?,
            inputs.current_member(pair, source, order[1])?,
        ],
        current,
    )?;
    Ok(GearPairCheck {
        pair,
        center_distance_mm: resolved.center_distance_mm,
        bevel: None,
        worm: None,
        current_engagement,
        cycle_drive_angles_rad: phases,
        cycle_reused: reused,
        sampled_cycle_ms,
        current_pose_ms: start.elapsed().as_secs_f64() * 1000.,
        current_drive_angle_rad: current,
        current_pose_checked: true,
        maximum_overlap_volume_mm3: 0.,
        scope: scope(mesh, engagement.is_some() && mesh != Mesh::Rack),
    })
}

pub(crate) fn validate(
    inputs: ValidationInputs<'_>,
    mut cache: Option<&mut crate::pair_cycle_proofs::CycleTransaction>,
) -> Result<BTreeMap<EntityId, GearPairCheck>, ModelError> {
    let ValidationInputs {
        document,
        evaluated,
        interrupt,
        ..
    } = inputs;
    let pairs = document
        .entities()
        .filter(|entity| matches!(&entity.data, EntityData::Extension { definition, .. } if definition.type_id == model::GEAR_PAIR_TYPE))
        .collect::<Vec<_>>();
    if pairs.len() > MAX_PAIRS {
        return Err(error(
            None,
            "Rust gear-pair validation exceeds the 64-pair budget",
        ));
    }
    let mut reports = BTreeMap::new();
    for entity in pairs {
        checkpoint(interrupt)?;
        let pair = entity.id;
        let source = model::read_gear_pair(document, pair).map_err(|e| error(Some(pair), e))?;
        let resolved = model::resolve_gear_pair_parameters(document, evaluated, &source.parameters)
            .map_err(|e| error(Some(pair), e))?;
        // A completed cycle proof is reused only for identical constructor
        // inputs; the current pose is always certified again.
        let key = cycle_key("rust-gear-pair-cycle/1", &inputs, pair, &source, &resolved)?;
        let previous = cache.as_deref().and_then(|c| c.lookup(&key));
        let reused = previous.is_some();
        let phases = previous.unwrap_or_else(|| cycle_phases(&resolved));
        let cycle = (phases.clone(), reused);
        let report = match resolved.gears[0].family() {
            GearFamily::Worm => {
                super::rust_worm::validate(&inputs, pair, &source, &resolved, cycle)?
            }
            GearFamily::Bevel => {
                super::rust_bevel::validate(&inputs, pair, &source, &resolved, cycle)?
            }
            _ => rotary(&inputs, pair, &source, &resolved, cycle)?,
        };
        if !reused && let Some(cache) = cache.as_deref_mut() {
            cache.stage(key, phases);
        }
        reports.insert(pair, report);
    }
    Ok(reports)
}

#[cfg(test)]
mod tests;
