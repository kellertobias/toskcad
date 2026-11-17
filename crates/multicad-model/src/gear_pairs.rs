//! Bounded interference evidence for periodic, unmodified generator pairs.
use crate::{ModelError, OccurrenceSnapshot, compiler::checkpoint, error};
use multicad_core::{
    constraint_sketch::Dimension,
    document::{Document, EntityData, EntityId, Mutation},
    model::{self, GearPairDefinition},
    parameters::EvaluatedParameters,
    references::TopologyKind,
};
use multicad_kernel::{BRepSelection, BooleanOperation, GeometryPolicy, RigidTransform, Solid};
use serde::Serialize;
use std::{
    collections::BTreeMap,
    f64::consts::{PI, TAU},
};
mod rust_bevel;
mod rust_spur;
mod rust_worm;

pub(crate) use rust_spur::validate as validate_rust;

const SAMPLES: usize = 9;
const MAX_PAIRS: usize = 64;
pub(crate) fn validate_parameters(
    document: &Document,
    evaluated: &EvaluatedParameters,
) -> Result<(), ModelError> {
    let mut count = 0;
    for entity in document.entities() {
        if !matches!(&entity.data,EntityData::Extension{definition,..} if definition.type_id==model::GEAR_PAIR_TYPE)
        {
            continue;
        }
        count += 1;
        if count > MAX_PAIRS {
            return Err(error(
                None,
                "Native gear-pair validation exceeds the 64-pair budget",
            ));
        }
        let source =
            model::read_gear_pair(document, entity.id).map_err(|e| error(Some(entity.id), e))?;
        model::resolve_gear_pair_parameters(document, evaluated, &source.parameters)
            .map_err(|e| error(Some(entity.id), e))?;
    }
    Ok(())
}
#[derive(Clone, Debug, Serialize)]
pub struct GearEngagementWitness {
    pub distance_mm: f64,
    pub kind: &'static str,
    /// Evaluation-local native face ordinals, never persistent source identity.
    pub faces: Option<[u32; 2]>,
    pub points_mm: [[f64; 3]; 2],
}
#[derive(Clone, Debug, Serialize)]
pub struct GearPairCheck {
    pub pair: EntityId,
    pub center_distance_mm: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub bevel: Option<model::ResolvedBevelPair>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub worm: Option<model::ResolvedWormPair>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub current_engagement: Option<GearEngagementWitness>,
    pub cycle_drive_angles_rad: Vec<f64>,
    pub cycle_reused: bool,
    /// Wall time in sampled interference checks; zero when a completed cycle is reused.
    pub sampled_cycle_ms: f64,
    pub current_pose_ms: f64,
    pub current_drive_angle_rad: f64,
    pub current_pose_checked: bool,
    pub maximum_overlap_volume_mm3: f64,
    pub scope: &'static str,
}
fn placed(solid: &Solid, pose: model::RigidPose, pair: EntityId) -> Result<Solid, ModelError> {
    solid
        .transformed(RigidTransform {
            rotation: pose.rotation,
            translation_mm: pose.translation_mm,
        })
        .map_err(|e| error(Some(pair), e))
}

fn overlap_checked(
    pair: EntityId,
    a: &Solid,
    b: &Solid,
    drive: f64,
    policy: GeometryPolicy,
    parallel: bool,
) -> Result<(), ModelError> {
    if let Some(common) = (if parallel {
        a.intersect_parallel(b, policy)
    } else {
        a.boolean(b, BooleanOperation::Intersect, policy)
    })
    .map_err(|e| {
        error(
            Some(pair),
            format!("Native gear-pair interference check: {e}"),
        )
    })? {
        let volume = common.volume_mm3().map_err(|e| error(Some(pair), e))?;
        return Err(error(
            Some(pair),
            format!(
                "Gear pair overlaps by {volume:.8} mm³ at tested drive {:.6} degrees; repair its geometry or clearance",
                drive.to_degrees()
            ),
        ));
    }
    Ok(())
}
fn engagement_flanks(
    solid: &Solid,
    feature: EntityId,
    teeth: u32,
    phase: f64,
    meridian: f64,
    pair: EntityId,
    names: (&str, &str),
) -> Result<Vec<u32>, ModelError> {
    let center = ((meridian - phase).rem_euclid(TAU) * f64::from(teeth) / TAU).round() as i64;
    let topology = solid
        .topology_snapshot(feature, 0)
        .map_err(|e| error(Some(pair), e))?;
    let mut faces = vec![];
    for offset in [0, -1, 1, -2, 2] {
        let tooth = (center + offset).rem_euclid(i64::from(teeth));
        let prefix = format!("{}/{tooth}/{}.", names.0, names.1);
        let mut found = [false; 2];
        for face in &topology.elements {
            if face.key.kind != TopologyKind::Face {
                continue;
            }
            for (i, suffix) in ["left", "right"].iter().enumerate() {
                if face.origins.iter().any(|origin| {
                    origin.source == feature && origin.role == format!("{prefix}{suffix}")
                }) {
                    found[i] = true;
                    if !faces.contains(&face.key.ordinal) {
                        faces.push(face.key.ordinal);
                    }
                }
            }
        }
        if !found.into_iter().all(|value| value) || faces.len() > 32 {
            return Err(error(
                Some(pair),
                "Gear engagement needs bounded native tooth-flank ancestry",
            ));
        }
    }
    Ok(faces)
}
fn rotary_engagement(
    pair: EntityId,
    features: [EntityId; 2],
    resolved: &model::ResolvedGearPair,
    solids: [&Solid; 2],
    policy: GeometryPolicy,
    bevel: bool,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<GearEngagementWitness, ModelError> {
    if bevel {
        return rotary_engagement_band(
            pair,
            features,
            resolved,
            solids,
            policy,
            ("bevel/tooth", "spherical_involute", true),
            interrupt,
        );
    }
    let family = resolved.gears[0].rotary().expect("external pair").family;
    if matches!(
        family,
        model::GearFamily::Herringbone | model::GearFamily::DoubleHelical
    ) {
        // Opposite bands each need their own actual native flank witness.
        // Whole-body proximity cannot establish engagement in both bands.
        let lower = rotary_engagement_band(
            pair,
            features,
            resolved,
            solids,
            policy,
            ("gear/lower/tooth", "involute", false),
            interrupt,
        )?;
        let upper = rotary_engagement_band(
            pair,
            features,
            resolved,
            solids,
            policy,
            ("gear/upper/tooth", "involute", false),
            interrupt,
        )?;
        let mut result = if lower.distance_mm >= upper.distance_mm {
            lower
        } else {
            upper
        };
        result.kind = "shifted_chevron_both_bands";
        return Ok(result);
    }
    let mut result = rotary_engagement_band(
        pair,
        features,
        resolved,
        solids,
        policy,
        ("gear/whole/tooth", "involute", true),
        interrupt,
    )?;
    if family == model::GearFamily::Helical && result.kind != "whole_solid_fallback" {
        result.kind = "shifted_helical_flank_faces";
    }
    Ok(result)
}
fn rotary_engagement_band(
    pair: EntityId,
    features: [EntityId; 2],
    resolved: &model::ResolvedGearPair,
    solids: [&Solid; 2],
    policy: GeometryPolicy,
    band: (&str, &str, bool),
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<GearEngagementWitness, ModelError> {
    let bevel = band.0 == "bevel/tooth";
    let names = (band.0, band.1);
    let allow_whole_fallback = band.2;
    let members = resolved.gears.each_ref().map(|gear| {
        if bevel {
            let g = gear.bevel().expect("resolved bevel pair");
            (
                g.teeth,
                g.nominal_tooth_thickness_mm - g.physical_tooth_thickness_mm,
            )
        } else {
            let g = gear.rotary().expect("resolved shifted external pair");
            (
                g.teeth,
                g.nominal_tooth_thickness_mm - g.physical_tooth_thickness_mm,
            )
        }
    });
    let faces = [
        engagement_flanks(
            solids[0],
            features[0],
            members[0].0,
            resolved.phases_rad[0],
            0.,
            pair,
            names,
        )?,
        engagement_flanks(
            solids[1],
            features[1],
            members[1].0,
            resolved.phases_rad[1],
            PI,
            pair,
            names,
        )?,
    ];
    let clearance = members.iter().map(|gear| gear.1).sum::<f64>();
    let mut nearest = f64::INFINITY;
    // Whole-solid intersection was checked before entering this function.
    // We need an existence witness of close native flanks, not the global
    // filled-solid minimum (whose containment work would duplicate that check).
    // Every fast-path witness comes from actual named BRep faces. A miss falls
    // back to the original whole-solid query, preserving acceptance outside
    // this small search window and recording the different witness scope.
    for &a in &faces[0] {
        for &b in &faces[1] {
            checkpoint(interrupt)?;
            let distance = solids[0]
                .minimum_distance(BRepSelection::Face(a), solids[1], BRepSelection::Face(b))
                .map_err(|e| error(Some(pair), e))?;
            nearest = nearest.min(distance.distance_mm);
            if distance.distance_mm <= clearance + policy.linear_tolerance_mm * 10. {
                let witness = &distance.witnesses[0];
                return Ok(GearEngagementWitness {
                    distance_mm: distance.distance_mm,
                    kind: if bevel {
                        "bevel_flank_faces"
                    } else {
                        "shifted_spur_flank_faces"
                    },
                    faces: Some([a, b]),
                    points_mm: [witness.a_mm, witness.b_mm],
                });
            }
        }
    }
    if !allow_whole_fallback {
        return Err(error(
            Some(pair),
            format!(
                "Gear band {} has no named native engagement proximity: nearest candidate flank gap {nearest} mm exceeds authored total thinning {clearance} mm",
                names.0,
            ),
        ));
    }
    checkpoint(interrupt)?;
    let distance = solids[0]
        .minimum_distance(BRepSelection::Whole, solids[1], BRepSelection::Whole)
        .map_err(|e| error(Some(pair), e))?;
    if distance.distance_mm <= clearance + policy.linear_tolerance_mm * 10. {
        let witness = &distance.witnesses[0];
        return Ok(GearEngagementWitness {
            distance_mm: distance.distance_mm,
            kind: "whole_solid_fallback",
            faces: None,
            points_mm: [witness.a_mm, witness.b_mm],
        });
    }
    Err(error(
        Some(pair),
        format!(
            "Gear members have no native engagement proximity: whole-solid gap {} mm and nearest candidate flank gap {nearest} mm exceed authored total thinning {clearance} mm",
            distance.distance_mm
        ),
    ))
}
fn worm_engagement(
    pair: EntityId,
    features: [EntityId; 2],
    resolved: &model::ResolvedGearPair,
    solids: [&Solid; 2],
    policy: GeometryPolicy,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<GearEngagementWitness, ModelError> {
    let context = resolved.worm.as_ref().expect("resolved worm pair");
    let wi = context.worm_member.index();
    let gi = context.wheel_member.index();
    let clearance = resolved
        .gears
        .iter()
        .map(|g| g.worm().expect("worm member").thinning_mm)
        .sum::<f64>();
    let mut selected = [Vec::new(), Vec::new()];
    let worm_topology = solids[wi]
        .topology_snapshot(features[wi], 0)
        .map_err(|e| error(Some(pair), e))?;
    for face in &worm_topology.elements {
        if face.key.kind == TopologyKind::Face
            && face.origins.iter().any(|o| {
                o.source == features[wi]
                    && matches!(o.role.as_str(), "worm/flank.lower" | "worm/flank.upper")
            })
        {
            selected[wi].push(face.key.ordinal);
        }
    }
    let wheel_topology = solids[gi]
        .topology_snapshot(features[gi], 0)
        .map_err(|e| error(Some(pair), e))?;
    let teeth = i64::from(context.wheel_teeth);
    let center = (-resolved.phases_rad[gi] * teeth as f64 / TAU).round() as i64;
    for offset in [0, -1, 1, -2, 2] {
        let sector = (center + offset).rem_euclid(teeth);
        for face in &wheel_topology.elements {
            if face.key.kind == TopologyKind::Face
                && face.origins.iter().any(|o| {
                    o.source == features[gi]
                        && (o.role == format!("worm_wheel/sector/{sector}/flank.minus")
                            || o.role == format!("worm_wheel/sector/{sector}/flank.plus"))
                })
                && !selected[gi].contains(&face.key.ordinal)
            {
                selected[gi].push(face.key.ordinal);
            }
        }
    }
    // Prioritize screw faces near the central contact section. Numeric face
    // ordinals and distances remain evaluation-local, never persistent IDs.
    let centers = solids[wi]
        .faces()
        .map_err(|e| error(Some(pair), e))?
        .into_iter()
        .map(|f| (f.index, f.center_mm))
        .collect::<BTreeMap<_, _>>();
    let origin = resolved.mating[wi].origin_mm;
    selected[wi].sort_by(|a, b| {
        let score = |id| {
            centers.get(id).map_or(f64::INFINITY, |p| {
                (0..3).map(|j| (p[j] - origin[j]).powi(2)).sum::<f64>()
            })
        };
        score(a).total_cmp(&score(b))
    });
    // The whole-solid overlap check has already passed. These actual named
    // faces provide only an existence witness within the declared clearance;
    // a candidate miss always runs the unchanged whole-solid distance query.
    if selected[wi].len() <= 256 && selected[gi].len() <= 32 {
        for &a in &selected[0] {
            for &b in &selected[1] {
                checkpoint(interrupt)?;
                let distance = solids[0]
                    .minimum_distance(BRepSelection::Face(a), solids[1], BRepSelection::Face(b))
                    .map_err(|e| error(Some(pair), e))?;
                if distance.distance_mm <= clearance + policy.linear_tolerance_mm * 10. {
                    let witness = &distance.witnesses[0];
                    return Ok(GearEngagementWitness {
                        distance_mm: distance.distance_mm,
                        kind: "worm_flank_faces",
                        faces: Some([a, b]),
                        points_mm: [witness.a_mm, witness.b_mm],
                    });
                }
            }
        }
    }
    checkpoint(interrupt)?;
    let distance = solids[0]
        .minimum_distance(BRepSelection::Whole, solids[1], BRepSelection::Whole)
        .map_err(|e| error(Some(pair), e))?;
    if distance.distance_mm > clearance + policy.linear_tolerance_mm * 10. {
        return Err(error(
            Some(pair),
            format!(
                "Gear members have no native engagement proximity: gap {} mm exceeds authored total thinning {clearance} mm",
                distance.distance_mm
            ),
        ));
    }
    let witness = &distance.witnesses[0];
    Ok(GearEngagementWitness {
        distance_mm: distance.distance_mm,
        kind: "whole_solid_fallback",
        faces: None,
        points_mm: [witness.a_mm, witness.b_mm],
    })
}
fn chevron_pair(resolved: &model::ResolvedGearPair) -> bool {
    resolved.gears.iter().all(|gear| {
        matches!(
            gear.family(),
            model::GearFamily::Herringbone | model::GearFamily::DoubleHelical
        )
    })
}

struct PoseCheck<'a> {
    pair: EntityId,
    features: [EntityId; 2],
    policy: GeometryPolicy,
    shifted: bool,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl PoseCheck<'_> {
    fn with_placement(
        &self,
        document: &Document,
        evaluated: &EvaluatedParameters,
        source: &GearPairDefinition,
        resolved: &model::ResolvedGearPair,
        solids: [&Solid; 2],
        drive: f64,
    ) -> Result<Option<GearEngagementWitness>, ModelError> {
        let mut placed_solids = Vec::with_capacity(2);
        for (instance, solid) in source.outputs.instances.into_iter().zip(solids) {
            checkpoint(self.interrupt)?;
            let pose = model::resolve_instance_pose(document, instance, evaluated)
                .map_err(|e| error(Some(self.pair), e))?;
            placed_solids.push(placed(solid, pose, self.pair)?);
        }
        self.placed(resolved, [&placed_solids[0], &placed_solids[1]], drive)
    }

    fn placed(
        &self,
        resolved: &model::ResolvedGearPair,
        solids: [&Solid; 2],
        drive: f64,
    ) -> Result<Option<GearEngagementWitness>, ModelError> {
        let Self {
            pair,
            features,
            policy,
            shifted,
            interrupt,
        } = *self;
        checkpoint(interrupt)?;
        // One exact COMMON at a time on this owner. OCCT owns any workers
        // within that query; Rust native handles stay on their owner thread.
        overlap_checked(
            pair,
            solids[0],
            solids[1],
            drive,
            policy,
            resolved.worm.is_some() || chevron_pair(resolved),
        )?;
        checkpoint(interrupt)?;
        let engagement = if resolved.worm.is_some() {
            Some(worm_engagement(
                pair, features, resolved, solids, policy, interrupt,
            )?)
        } else if resolved.bevel.is_some() || shifted {
            Some(rotary_engagement(
                pair,
                features,
                resolved,
                solids,
                policy,
                resolved.bevel.is_some(),
                interrupt,
            )?)
        } else {
            if let Some(travel) = &resolved.rack_travel {
                let rack = resolved.gears[travel.rack_member.index()]
                    .rack()
                    .expect("rack member");
                let pinion = resolved.gears[travel.pinion_member.index()]
                    .rotary()
                    .expect("pinion member");
                // A no-overlap result also accepts separated bodies. Require
                // actual native proximity within the total authored thinning.
                let clearance = rack.nominal_tooth_thickness_mm - rack.physical_tooth_thickness_mm
                    + pinion.nominal_tooth_thickness_mm
                    - pinion.physical_tooth_thickness_mm;
                let distance = solids[0]
                    .minimum_distance(BRepSelection::Whole, solids[1], BRepSelection::Whole)
                    .map_err(|e| error(Some(pair), e))?;
                if distance.distance_mm > clearance + policy.linear_tolerance_mm * 10. {
                    return Err(error(
                        Some(pair),
                        format!(
                            "Gear members have no native engagement proximity: gap {} mm exceeds authored total thinning {} mm",
                            distance.distance_mm, clearance
                        ),
                    ));
                }
            }
            None
        };
        checkpoint(interrupt)?;
        Ok(engagement)
    }
}

fn shifted_engagement(document: &Document, features: [EntityId; 2]) -> bool {
    features.into_iter().any(|feature| {
        model::read_gear(document, feature).is_ok_and(|g| g.parameters.shift_mate.is_some())
    })
}

/// Process-local signature of every actual constructor input. Debug is
/// derived on these closed structs and never persisted as a file format.
/// Drive is excluded only from the periodic cycle proof; its actual current
/// pose is checked independently on every compile. Each backend has its own
/// prefix, so one backend's proof never stands in for the other's.
fn cycle_key(
    prefix: &str,
    inputs: &ValidationInputs<'_>,
    pair: EntityId,
    source: &GearPairDefinition,
    resolved: &model::ResolvedGearPair,
) -> Result<String, ModelError> {
    let document = inputs.document;
    let specifications = [
        crate::gears::specification(source.parameters.gears[0], &resolved.gears[0])?,
        crate::gears::specification(source.parameters.gears[1], &resolved.gears[1])?,
    ];
    let base = source
        .parameters
        .frame
        .resolve(inputs.evaluated)
        .map_err(|e| error(Some(pair), e))?;
    let shifted_engagement = shifted_engagement(document, source.parameters.gears);
    let key = format!(
        "{prefix}/{SAMPLES}|{:?}|{pair:?}|{:?}|{:?}|{specifications:?}|{:?}|{base:?}|{:?}|shifted-engagement:{shifted_engagement}",
        document.id(),
        source.parameters.gears,
        source.outputs,
        resolved.gears.each_ref().map(|g| g.frame()),
        inputs.policy,
    );
    Ok(if let Some(worm) = &resolved.worm {
        format!("{key}|{worm:?}")
    } else if let Some(travel) = &resolved.rack_travel {
        format!(
            "{key}|{:?}|{:?}|{}|{}|{}",
            travel.rack_member,
            travel.pinion_member,
            travel.min_drive_angle_rad,
            travel.max_drive_angle_rad,
            travel.pinion_period_rad
        )
    } else {
        key
    })
}

/// Nine phases over one tooth (or worm) period; a finite rack adds both
/// conservative travel endpoints.
fn cycle_phases(resolved: &model::ResolvedGearPair) -> Vec<f64> {
    if let Some(travel) = &resolved.rack_travel {
        let mut phases = (0..SAMPLES)
            .map(|i| travel.pinion_period_rad * (i as f64 / (SAMPLES - 1) as f64 - 0.5))
            .collect::<Vec<_>>();
        phases.extend([travel.min_drive_angle_rad, travel.max_drive_angle_rad]);
        phases
    } else if let Some(worm) = &resolved.worm {
        (0..SAMPLES)
            .map(|i| worm.worm_period_rad * i as f64 / (SAMPLES - 1) as f64)
            .collect()
    } else {
        (0..SAMPLES)
            .map(|i| TAU / f64::from(resolved.gears[0].teeth()) * i as f64 / (SAMPLES - 1) as f64)
            .collect()
    }
}

pub(crate) struct ValidationInputs<'a> {
    pub document: &'a Document,
    pub evaluated: &'a EvaluatedParameters,
    pub solids: &'a BTreeMap<EntityId, Solid>,
    pub occurrences: &'a BTreeMap<EntityId, std::sync::Arc<OccurrenceSnapshot>>,
    pub revision: u64,
    pub policy: GeometryPolicy,
    pub interrupt: Option<&'a dyn Fn() -> bool>,
}

impl ValidationInputs<'_> {
    fn current_member(
        &self,
        pair: EntityId,
        source: &GearPairDefinition,
        index: usize,
    ) -> Result<&Solid, ModelError> {
        // The Full compiler created this native operand and occurrence receipt
        // together. Confirm the pair's fixed source still names precisely that
        // instance, Part, producer and evaluated world pose.
        let instance = source.outputs.instances[index];
        let producer = source.parameters.gears[index];
        let component = model::read_gear(self.document, producer)
            .map_err(|e| error(Some(pair), e))?
            .output
            .body;
        let expected_pose = model::resolve_instance_pose(self.document, instance, self.evaluated)
            .map_err(|e| error(Some(pair), e))?;
        let invalid = || {
            error(
                Some(pair),
                "Current gear occurrence does not match its checked Part, producer, pose or revision",
            )
        };
        let occurrence = self.occurrences.get(&instance).ok_or_else(invalid)?;
        let solid = self.solids.get(&instance).ok_or_else(invalid)?;
        if occurrence.instance != instance
            || occurrence.component != component
            || occurrence.producer != producer
            || occurrence.pose != expected_pose
            || occurrence.topology.owner != instance
            || occurrence.topology.revision != self.revision
            || solid.geometry_policy() != self.policy
            || solid.bounds().map_err(|e| error(Some(pair), e))? != occurrence.bounds
        {
            return Err(invalid());
        }
        Ok(solid)
    }
}

fn sample_document(
    document: &Document,
    pair: EntityId,
    source: &GearPairDefinition,
    drive: f64,
) -> Result<Document, ModelError> {
    let mut source = source.clone();
    source.parameters.drive_angle = Dimension::Literal {
        value: drive,
        unit: "rad".into(),
    };
    let mut entity = document
        .get(pair)
        .ok_or_else(|| error(Some(pair), "Gear pair is missing"))?
        .clone();
    entity.data = EntityData::Extension {
        definition: source.definition().map_err(|e| error(Some(pair), e))?,
        inputs: source
            .parameters
            .inputs()
            .map_err(|e| error(Some(pair), e))?,
    };
    document
        .edited(&[Mutation::Replace(entity)])
        .map_err(|e| error(Some(pair), e))
}
pub(crate) fn validate(
    inputs: ValidationInputs<'_>,
    mut cache: Option<&mut crate::pair_cycle_proofs::CycleTransaction>,
) -> Result<BTreeMap<EntityId, GearPairCheck>, ModelError> {
    let ValidationInputs {
        document,
        evaluated,
        solids,
        policy,
        interrupt,
        ..
    } = inputs;
    let pairs=document.entities().filter(|e|matches!(&e.data,EntityData::Extension{definition,..} if definition.type_id==model::GEAR_PAIR_TYPE)).collect::<Vec<_>>();
    if pairs.len() > MAX_PAIRS {
        return Err(error(
            None,
            "Native gear-pair validation exceeds the 64-pair budget",
        ));
    }
    let mut reports = BTreeMap::new();
    for entity in pairs {
        checkpoint(interrupt)?;
        let pair = entity.id;
        let source = model::read_gear_pair(document, pair).map_err(|e| error(Some(pair), e))?;
        let resolved = model::resolve_gear_pair_parameters(document, evaluated, &source.parameters)
            .map_err(|e| error(Some(pair), e))?;
        let a = solids
            .get(&source.parameters.gears[0])
            .ok_or_else(|| error(Some(pair), "First paired native gear is unavailable"))?;
        let b = solids
            .get(&source.parameters.gears[1])
            .ok_or_else(|| error(Some(pair), "Second paired native gear is unavailable"))?;
        let shifted_engagement = shifted_engagement(document, source.parameters.gears);
        let checker = PoseCheck {
            pair,
            features: source.parameters.gears,
            policy,
            shifted: shifted_engagement,
            interrupt,
        };
        let key = cycle_key("gear-pair-cycle/6", &inputs, pair, &source, &resolved)?;
        let previous = cache.as_deref().and_then(|c| c.lookup(&key));
        let reused = previous.is_some();
        let cycle_started = std::time::Instant::now();
        let phases = if let Some(phases) = previous {
            phases
        } else {
            let phases = cycle_phases(&resolved);
            for drive in &phases {
                checkpoint(interrupt)?;
                let sample = sample_document(document, pair, &source, *drive)?;
                // Use this sample's drive for both native placement and its
                // named flank witnesses, including wrapped and rack endpoints.
                let sampled =
                    model::read_gear_pair(&sample, pair).map_err(|e| error(Some(pair), e))?;
                let sampled_resolved =
                    model::resolve_gear_pair_parameters(&sample, evaluated, &sampled.parameters)
                        .map_err(|e| error(Some(pair), e))?;
                checker.with_placement(
                    &sample,
                    evaluated,
                    &sampled,
                    &sampled_resolved,
                    [a, b],
                    *drive,
                )?;
            }
            if let Some(cache) = cache.as_deref_mut() {
                cache.stage(key, phases.clone());
            }
            phases
        };
        let sampled_cycle_ms = if reused {
            0.
        } else {
            cycle_started.elapsed().as_secs_f64() * 1000.
        };
        checkpoint(interrupt)?;
        let current = source
            .parameters
            .drive_angle
            .value_in(multicad_core::parameters::Unit::Radian, Some(evaluated))
            .map_err(|e| error(Some(pair), e))?;
        let current_started = std::time::Instant::now();
        let current_engagement = if chevron_pair(&resolved) {
            let placed = [
                inputs.current_member(pair, &source, 0)?,
                inputs.current_member(pair, &source, 1)?,
            ];
            checker.placed(&resolved, placed, current)?
        } else {
            checker.with_placement(document, evaluated, &source, &resolved, [a, b], current)?
        };
        let current_pose_ms = current_started.elapsed().as_secs_f64() * 1000.;
        reports.insert(pair,GearPairCheck{pair,current_engagement,center_distance_mm:resolved.center_distance_mm,bevel:resolved.bevel.clone(),worm:resolved.worm.clone(),cycle_drive_angles_rad:phases,cycle_reused:reused,sampled_cycle_ms,current_pose_ms,current_drive_angle_rad:current,current_pose_checked:true,maximum_overlap_volume_mm3:0.,scope:if resolved.worm.is_some() { "Native current pose and nine sampled phases over one full ZA worm turn with native no-overlap and engagement proximity witnesses for its generated wheel. No continuous collision-free motion, strength or manufacturing-fit certification." } else if resolved.rack_travel.is_some() { "Native current pose, nine sampled phases over one pinion tooth period and both conservative finite travel endpoints. No continuous collision-free motion, strength or manufacturing-fit certification." } else if resolved.bevel.is_some() { "Native current pose and nine sampled phases over one spherical-involute bevel tooth cycle with native no-overlap and engagement proximity witnesses. No continuous collision-free motion, strength or manufacturing-fit certification." } else if shifted_engagement { "Native current pose and nine sampled phases over one matched shifted external tooth cycle, with whole-solid no-overlap and native named-flank proximity in each chevron band, or a whole-solid fallback for single-band gears. No continuous collision-free motion, strength or manufacturing-fit certification." } else { "Native current pose and nine sampled phases over one pinion tooth cycle, under the unmodified periodic generator assumption. No continuous collision-free motion, strength or manufacturing-fit certification." }});
    }
    Ok(reports)
}
