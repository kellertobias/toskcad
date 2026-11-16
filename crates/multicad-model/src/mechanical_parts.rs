//! Native mounting geometry and nut threads from retained canonical parameters.
mod custom_connection;
mod profiled_heads;
mod screw_profiles;
mod wood_screws;
use crate::{ModelError, error};
use multicad_core::{
    document::EntityId,
    model::{
        MechanicalPartDefinition, MechanicalPartParameters, ResolvedMechanicalPart,
        ResolvedScrewDrive, ResolvedScrewHead, ResolvedThread, ThreadMode, ThreadRepresentation,
    },
    parameters::EvaluatedParameters,
    references::Origin,
    sketch::Frame,
};
use multicad_kernel::{
    BooleanOperation, CurvedEdge2d, CurvedProfile, GeometryPolicy, PlaneFrame, RigidTransform,
    Solid,
    rust_kernel::{
        cosmetic_fastener::Countersunk,
        headed_screw::{ScrewHead, ScrewSpec},
    },
};

use std::{cell::RefCell, collections::BTreeMap};

const SCREW_FACTORY_ENTRIES: usize = 16;
const SCREW_FACTORY_BYTES: usize = 32 * 1024 * 1024;
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord)]
struct ScrewKey {
    numeric: Vec<u64>,
    // Derived Debug preserves every resolved field, including fit semantics;
    // retaining it also makes future resolved-field additions invalidate hits.
    resolved_guard: String,
}
fn screw_key(resolved: &ResolvedMechanicalPart, policy: GeometryPolicy) -> Option<ScrewKey> {
    let ResolvedMechanicalPart::MachineScrew {
        underhead_length_mm,
        shank_radius_mm,
        head,
        drive,
        thread,
        ..
    } = resolved
    else {
        return None;
    };
    let mut numbers = vec![
        *underhead_length_mm,
        *shank_radius_mm,
        head.height_mm(),
        thread.nominal_major_radius_mm,
        thread.major_radius_mm,
        thread.pitch_mm,
        thread.lead_mm,
        thread.length_mm,
        thread.minimum_wall_mm,
        thread.profile.depth_mm,
        thread.profile.crest_width_mm,
        thread.profile.lower_flank_deg,
        thread.profile.upper_flank_deg,
        thread.profile.radial_taper,
        policy.linear_tolerance_mm,
        policy.boolean_fuzzy_mm,
        policy.minimum_feature_mm,
        policy.maximum_coordinate_mm,
    ];
    numbers.extend(thread.frame.origin);
    numbers.extend(thread.frame.u);
    numbers.extend(thread.frame.v);
    numbers.extend(thread.frame.n);
    match head {
        ResolvedScrewHead::SocketCap { radius_mm, .. } => numbers.push(*radius_mm),
        ResolvedScrewHead::ExternalHex {
            across_flats_mm, ..
        } => numbers.push(*across_flats_mm),
        ResolvedScrewHead::Pan {
            radius_mm,
            crown_radius_mm,
            top_land_radius_mm,
            neck_radius_mm,
            ..
        }
        | ResolvedScrewHead::Button {
            radius_mm,
            crown_radius_mm,
            top_land_radius_mm,
            neck_radius_mm,
            ..
        } => numbers.extend([
            *radius_mm,
            *crown_radius_mm,
            *top_land_radius_mm,
            *neck_radius_mm,
        ]),
        ResolvedScrewHead::Countersunk {
            radius_mm,
            included_angle_deg,
            edge_chamfer_mm,
            neck_radius_mm,
            ..
        } => numbers.extend([
            *radius_mm,
            *included_angle_deg,
            *edge_chamfer_mm,
            *neck_radius_mm,
        ]),
    }
    match drive {
        ResolvedScrewDrive::Hexalobular {
            major_radius_mm,
            minor_radius_mm,
            outer_lobe_radius_mm,
            depth_mm,
            bottom_angle_deg,
        } => numbers.extend([
            *major_radius_mm,
            *minor_radius_mm,
            *outer_lobe_radius_mm,
            *depth_mm,
            *bottom_angle_deg,
        ]),
        ResolvedScrewDrive::PhillipsH { profile: p } => {
            numbers.extend([
                p.mouth_diameter_mm,
                p.mouth_radius_mm,
                p.b_width_mm,
                p.g_diameter_mm,
                p.tilt_angle_deg,
                p.outer_flank_angle_deg,
                p.floor_angle_deg,
                p.outer_corner_angle_deg,
            ]);
            match p.inner {
                multicad_core::model::ResolvedPhillipsInnerForm::Rounded { radius_mm } => {
                    numbers.push(radius_mm)
                }
                multicad_core::model::ResolvedPhillipsInnerForm::ShallowV {
                    angle_deg,
                    e_width_mm,
                } => numbers.extend([angle_deg, e_width_mm]),
            }
        }
        _ => {}
    }
    if let ResolvedScrewDrive::HexSocket {
        across_flats_mm,
        depth_mm,
    } = drive
    {
        numbers.extend([*across_flats_mm, *depth_mm]);
    }
    let mut numeric = numbers.into_iter().map(f64::to_bits).collect::<Vec<_>>();
    for value in [
        thread.nominal_pitch_diameter_mm,
        thread.actual_pitch_diameter_mm,
        thread.profile.rounding_radius_mm,
    ] {
        numeric.push(u64::from(value.is_some()));
        numeric.push(value.unwrap_or(0.).to_bits());
    }
    numeric.push(u64::from(thread.starts));
    Some(ScrewKey {
        numeric,
        resolved_guard: format!("machine-screw-factory/1|{resolved:?}"),
    })
}
struct ScrewEntry {
    solid: Solid,
    owner: EntityId,
    bytes: usize,
    used: u64,
}
#[derive(Default)]
struct ScrewFactory {
    entries: BTreeMap<ScrewKey, ScrewEntry>,
    clock: u64,
    bytes: usize,
    requests: u64,
    hits: u64,
    builds: u64,
}
thread_local! {static SCREWS: RefCell<ScrewFactory> = RefCell::new(ScrewFactory::default());}
fn cached_screw(key: &ScrewKey, feature: EntityId) -> Option<Solid> {
    SCREWS.with(|slot| {
        let mut cache = slot.borrow_mut();
        cache.clock += 1;
        cache.requests += 1;
        let clock = cache.clock;
        let entry = cache.entries.get_mut(key)?;
        let result = entry
            .solid
            .rebind_exclusive_owner(entry.owner, feature)
            .ok()?;
        entry.used = clock;
        cache.hits += 1;
        Some(result)
    })
}
fn remember_screw(key: ScrewKey, feature: EntityId, solid: &Solid) {
    // Retain only an already completed, fully owned native result. Optional
    // bookkeeping errors cannot turn successful geometry into a model error.
    let Ok(owned) = solid.rebind_exclusive_owner(feature, feature) else {
        return;
    };
    let Ok(brep) = owned.to_brep() else {
        return;
    };
    let Ok(topology) = owned.topology_snapshot(feature, 0) else {
        return;
    };
    let Ok(receipt) = serde_json::to_vec(&topology) else {
        return;
    };
    let Some(bytes) = brep.len().checked_add(receipt.len()) else {
        return;
    };
    if bytes > SCREW_FACTORY_BYTES {
        return;
    }
    SCREWS.with(|slot| {
        let mut cache = slot.borrow_mut();
        cache.builds += 1;
        if let Some(old) = cache.entries.remove(&key) {
            cache.bytes -= old.bytes;
        }
        while cache.entries.len() >= SCREW_FACTORY_ENTRIES
            || cache.bytes > SCREW_FACTORY_BYTES - bytes
        {
            let Some(oldest) = cache
                .entries
                .iter()
                .min_by_key(|(_, entry)| entry.used)
                .map(|(key, _)| key.clone())
            else {
                break;
            };
            if let Some(old) = cache.entries.remove(&oldest) {
                cache.bytes -= old.bytes;
            }
        }
        cache.bytes += bytes;
        let used = cache.clock;
        cache.entries.insert(
            key,
            ScrewEntry {
                solid: owned,
                owner: feature,
                bytes,
                used,
            },
        );
    });
}

fn hexagon(across_flats: f64) -> Vec<CurvedEdge2d> {
    let radius = across_flats / 3_f64.sqrt();
    let point = |i: usize| {
        let angle = std::f64::consts::PI / 6. + i as f64 * std::f64::consts::TAU / 6.;
        [radius * angle.cos(), radius * angle.sin()]
    };
    (0..6)
        .map(|i| CurvedEdge2d::Line {
            a: point(i),
            b: point((i + 1) % 6),
        })
        .collect()
}

/// Screw components are built in the underhead coordinate system. Each owns
/// distinct profile and cap roles before native ancestry is combined.
fn screw_prism(
    feature: EntityId,
    component: &str,
    z: f64,
    height: f64,
    outer: Vec<CurvedEdge2d>,
    policy: GeometryPolicy,
) -> Result<Solid, ModelError> {
    let origins = (0..outer.len())
        .map(|i| Origin {
            source: feature,
            role: format!("screw.{component}/edge/{i}"),
        })
        .collect::<Vec<_>>();
    let profile = CurvedProfile {
        frame: PlaneFrame {
            origin: [0., 0., z],
            u: [1., 0., 0.],
            v: [0., 1., 0.],
        },
        outer,
        holes: vec![],
    };
    Solid::extrude_curved_profile_vector_named(
        &profile,
        [0., 0., height],
        feature,
        &origins,
        policy,
    )
    .and_then(|s| {
        s.with_feature_cap_roles(
            feature,
            &format!("screw.{component}/start"),
            &format!("screw.{component}/end"),
        )
    })
    .map_err(|e| error(Some(feature), e))
}
fn screw_boolean(
    feature: EntityId,
    a: &Solid,
    b: &Solid,
    operation: BooleanOperation,
    policy: GeometryPolicy,
) -> Result<Solid, ModelError> {
    let result = a
        .boolean(b, operation, policy)
        .map_err(|e| error(Some(feature), e))?
        .ok_or_else(|| error(Some(feature), "Machine screw operation produced no solid"))?;
    if result.solid_count().map_err(|e| error(Some(feature), e))? != 1 {
        return Err(error(
            Some(feature),
            "Machine screw must remain one connected solid",
        ));
    }
    Ok(result)
}
fn extrude(
    feature: EntityId,
    frame: Frame,
    height: f64,
    outer: Vec<CurvedEdge2d>,
    holes: Vec<Vec<CurvedEdge2d>>,
    policy: GeometryPolicy,
) -> Result<Solid, ModelError> {
    let mut origins = (0..outer.len())
        .map(|i| Origin {
            source: feature,
            role: format!("hardware.outer/{i}"),
        })
        .collect::<Vec<_>>();
    for (loop_index, edges) in holes.iter().enumerate() {
        origins.extend((0..edges.len()).map(|i| Origin {
            source: feature,
            role: format!("hardware.bore/{loop_index}/{i}"),
        }));
    }
    let profile = CurvedProfile {
        frame: PlaneFrame {
            origin: frame.origin,
            u: frame.u,
            v: frame.v,
        },
        outer,
        holes,
    };
    Solid::extrude_curved_profile_vector_named(
        &profile,
        frame.n.map(|n| n * height),
        feature,
        &origins,
        policy,
    )
    .map_err(|e| error(Some(feature), e))
}
fn ring(
    feature: EntityId,
    frame: Frame,
    bore: f64,
    outside: f64,
    height: f64,
    policy: GeometryPolicy,
) -> Result<Solid, ModelError> {
    let circle = |radius_mm| CurvedEdge2d::Circle {
        center: [0., 0.],
        radius_mm,
    };
    extrude(
        feature,
        frame,
        height,
        vec![circle(outside)],
        vec![vec![circle(bore)]],
        policy,
    )
}
pub(crate) fn build(
    feature: EntityId,
    source: &MechanicalPartDefinition,
    evaluated: &EvaluatedParameters,
    policy: GeometryPolicy,
) -> Result<Solid, ModelError> {
    let resolved = source
        .parameters
        .resolve(evaluated)
        .map_err(|e| error(Some(feature), e))?;
    build_resolved(feature, resolved, policy)
}

pub(crate) fn build_rust(
    feature: EntityId,
    source: &MechanicalPartDefinition,
    evaluated: &EvaluatedParameters,
    policy: GeometryPolicy,
) -> Result<Solid, ModelError> {
    if let MechanicalPartParameters::MachineScrew { thread, .. }
    | MechanicalPartParameters::HexNut { thread, .. } = &source.parameters
    {
        crate::threads::validate_rust_family(feature, thread.family)?;
    }
    let resolved = source
        .parameters
        .resolve(evaluated)
        .map_err(|cause| error(Some(feature), cause))?;
    build_rust_resolved(feature, resolved, policy)
}

pub(crate) fn build_rust_resolved(
    feature: EntityId,
    resolved: ResolvedMechanicalPart,
    policy: GeometryPolicy,
) -> Result<Solid, ModelError> {
    match resolved {
        ResolvedMechanicalPart::CustomConnection { parameters } => {
            custom_connection::build_rust(feature, *parameters, policy)
        }
        ResolvedMechanicalPart::WoodScrew { parameters } => {
            wood_screws::build_rust(feature, &parameters, policy)
        }
        ResolvedMechanicalPart::MachineScrew {
            underhead_length_mm,
            shank_radius_mm,
            head,
            drive,
            thread,
            ..
        } => rust_machine_screw(
            feature,
            underhead_length_mm,
            shank_radius_mm,
            head,
            drive,
            &thread,
            policy,
        ),
        ResolvedMechanicalPart::HexNut {
            across_flats_mm,
            thread,
        } => {
            if thread.representation == ThreadRepresentation::CosmeticEnvelope {
                let frame = rust_cosmetic_frame(feature, &thread)?;
                return Solid::rust_cosmetic_hex_nut(
                    across_flats_mm,
                    thread.major_radius_mm - thread.profile.depth_mm,
                    thread.length_mm,
                    thread.minimum_wall_mm + thread.profile.depth_mm,
                    feature,
                    frame,
                    policy,
                )
                .map_err(|cause| error(Some(feature), cause));
            }
            let (profile, _, frame) = crate::threads::rust_profile(feature, &thread)?;
            Solid::rust_hex_nut(
                profile,
                across_flats_mm,
                thread.minimum_wall_mm,
                feature,
                frame,
                policy,
            )
            .map_err(|cause| error(Some(feature), cause))
        }
        ResolvedMechanicalPart::Washer {
            bore_radius_mm,
            outside_radius_mm,
            thickness_mm,
            frame,
        }
        | ResolvedMechanicalPart::BearingEnvelope {
            bore_radius_mm,
            outside_radius_mm,
            width_mm: thickness_mm,
            frame,
        } => Solid::rust_annulus(
            bore_radius_mm,
            outside_radius_mm,
            thickness_mm,
            feature,
            PlaneFrame {
                origin: frame.origin,
                u: frame.u,
                v: frame.v,
            },
            policy,
        )
        .map_err(|cause| error(Some(feature), cause)),
    }
}

fn rust_cosmetic_frame(
    feature: EntityId,
    thread: &ResolvedThread,
) -> Result<PlaneFrame, ModelError> {
    if thread.profile.radial_taper != 0. || thread.end_treatments.is_some() {
        return Err(error(
            Some(feature),
            "Rust cosmetic fasteners require a cylindrical thread envelope without end treatments",
        ));
    }
    Ok(PlaneFrame {
        origin: thread.frame.origin,
        u: thread.frame.u,
        v: thread.frame.v,
    })
}

fn rust_machine_screw(
    feature: EntityId,
    underhead_length_mm: f64,
    shank_radius_mm: f64,
    head: ResolvedScrewHead,
    drive: ResolvedScrewDrive,
    thread: &ResolvedThread,
    policy: GeometryPolicy,
) -> Result<Solid, ModelError> {
    let legacy_head = rust_legacy_head(&head, &drive);
    if thread.representation == ThreadRepresentation::CosmeticEnvelope {
        let frame = rust_cosmetic_frame(feature, thread)?;
        if shank_radius_mm != thread.major_radius_mm {
            return Err(error(
                Some(feature),
                "Rust cosmetic screws require matching shank and thread envelope diameters",
            ));
        }
        if let Some(head) = legacy_head {
            return Solid::rust_cosmetic_headed_screw(
                head,
                shank_radius_mm,
                underhead_length_mm,
                feature,
                frame,
                policy,
            )
            .map_err(|cause| error(Some(feature), cause));
        }
        if matches!(drive, ResolvedScrewDrive::Hexalobular { .. }) {
            let mut physical = thread.clone();
            physical.representation = ThreadRepresentation::Modeled;
            let source = profiled_heads::build(
                feature,
                underhead_length_mm,
                shank_radius_mm,
                &head,
                &drive,
                &physical,
                policy,
            )?
            .ok_or_else(|| error(Some(feature), "Unsupported cosmetic physical head"))?;
            return source
                .cosmetic_profiled_screw_envelope()
                .map_err(|cause| error(Some(feature), cause));
        }
        let (
            ResolvedScrewHead::Countersunk {
                radius_mm,
                height_mm,
                included_angle_deg,
                edge_chamfer_mm,
                neck_radius_mm: 0.,
            },
            ResolvedScrewDrive::HexSocket {
                across_flats_mm,
                depth_mm,
            },
        ) = (head, drive)
        else {
            return Err(error(
                Some(feature),
                "Rust cosmetic screws support external hex heads or socket-cap and straight countersunk heads with flat hex sockets",
            ));
        };
        return Solid::rust_cosmetic_countersunk_screw(
            Countersunk {
                head_radius_mm: radius_mm,
                head_height_mm: height_mm,
                included_angle_deg,
                edge_chamfer_mm,
                shaft_radius_mm: shank_radius_mm,
                shaft_length_mm: underhead_length_mm,
                socket_across_flats_mm: across_flats_mm,
                socket_depth_mm: depth_mm,
            },
            feature,
            frame,
            policy,
        )
        .map_err(|cause| error(Some(feature), cause));
    }
    if let Some(solid) = profiled_heads::build(
        feature,
        underhead_length_mm,
        shank_radius_mm,
        &head,
        &drive,
        thread,
        policy,
    )? {
        return Ok(solid);
    }
    let (profile, root_radius, frame) = crate::threads::rust_profile(feature, thread)?;
    let head = legacy_head.ok_or_else(|| {
        error(
            Some(feature),
            "Rust machine screws support a socket cap with a flat hex socket or an external hex head",
        )
    })?;
    Solid::rust_machine_screw(
        profile,
        root_radius,
        ScrewSpec {
            head,
            shank_length_mm: underhead_length_mm - thread.length_mm,
            shank_radius_mm,
        },
        feature,
        frame,
        policy,
    )
    .map_err(|cause| error(Some(feature), cause))
}

fn rust_legacy_head(head: &ResolvedScrewHead, drive: &ResolvedScrewDrive) -> Option<ScrewHead> {
    Some(match (head, drive) {
        (
            ResolvedScrewHead::SocketCap {
                radius_mm,
                height_mm,
            },
            ResolvedScrewDrive::HexSocket {
                across_flats_mm,
                depth_mm,
            },
        ) => ScrewHead::SocketCap {
            radius_mm: *radius_mm,
            height_mm: *height_mm,
            socket_across_flats_mm: *across_flats_mm,
            socket_depth_mm: *depth_mm,
        },
        (
            ResolvedScrewHead::ExternalHex {
                across_flats_mm,
                height_mm,
            },
            ResolvedScrewDrive::ExternalHex,
        ) => ScrewHead::ExternalHex {
            across_flats_mm: *across_flats_mm,
            height_mm: *height_mm,
        },
        _ => return None,
    })
}

pub(super) fn build_resolved(
    feature: EntityId,
    resolved: ResolvedMechanicalPart,
    policy: GeometryPolicy,
) -> Result<Solid, ModelError> {
    let factory_key = screw_key(&resolved, policy);
    match resolved {
        ResolvedMechanicalPart::CustomConnection { parameters } => {
            custom_connection::build(feature, *parameters, policy)
        }
        ResolvedMechanicalPart::WoodScrew { parameters } => {
            wood_screws::build(feature, &parameters, policy)
        }
        ResolvedMechanicalPart::MachineScrew {
            underhead_length_mm,
            shank_radius_mm,
            head,
            drive,
            mut thread,
            ..
        } => {
            let frame = thread.frame;
            let shank_length = underhead_length_mm - thread.length_mm;
            if let ResolvedScrewDrive::HexSocket {
                across_flats_mm,
                depth_mm,
            } = &drive
                && matches!(head, ResolvedScrewHead::SocketCap { .. })
                && (head.height_mm() - depth_mm < policy.minimum_feature_mm
                    || head.minimum_radius_mm() - across_flats_mm / 3_f64.sqrt()
                        < policy.minimum_feature_mm)
            {
                return Err(error(
                    Some(feature),
                    "Hex socket floor or corner wall is below native minimum feature size",
                ));
            }
            policy.validate().map_err(|e| error(Some(feature), e))?;
            screw_profiles::validate_native_clearance(
                feature,
                &head,
                &drive,
                shank_radius_mm,
                underhead_length_mm,
                &thread,
                policy,
            )?;
            let factory_key = factory_key
                .ok_or_else(|| error(Some(feature), "Missing resolved screw factory key"))?;
            if let Some(cached) = cached_screw(&factory_key, feature) {
                return Ok(cached);
            }
            let circle = |radius_mm| {
                vec![CurvedEdge2d::Circle {
                    center: [0., 0.],
                    radius_mm,
                }]
            };
            let new_profile = screw_profiles::new_head(&head);
            let height = head.height_mm();
            let early_socket = matches!(head, ResolvedScrewHead::SocketCap { .. })
                && matches!(drive, ResolvedScrewDrive::HexSocket { .. });
            let mut result = if new_profile {
                screw_profiles::blank(feature, &head, shank_radius_mm, shank_length, policy)?
            } else {
                let head_outer = match head {
                    ResolvedScrewHead::SocketCap { radius_mm, .. } => circle(radius_mm),
                    ResolvedScrewHead::ExternalHex {
                        across_flats_mm, ..
                    } => hexagon(across_flats_mm),
                    _ => unreachable!("New head profile dispatched above"),
                };
                screw_prism(feature, "head", -height, height, head_outer, policy)?
            };
            // Preserve the old flat-cap/flat-hex-socket construction exactly.
            if early_socket
                && let ResolvedScrewDrive::HexSocket {
                    across_flats_mm,
                    depth_mm,
                } = &drive
            {
                let socket = screw_prism(
                    feature,
                    "socket",
                    -height,
                    *depth_mm,
                    hexagon(*across_flats_mm),
                    policy,
                )?;
                result = screw_boolean(
                    feature,
                    &result,
                    &socket,
                    BooleanOperation::Subtract,
                    policy,
                )?;
            }
            if !new_profile && shank_length > 0. {
                let shank = screw_prism(
                    feature,
                    "shank",
                    0.,
                    shank_length,
                    circle(shank_radius_mm),
                    policy,
                )?;
                result = screw_boolean(feature, &result, &shank, BooleanOperation::Union, policy)?;
            }
            thread.frame = Frame {
                origin: [0., 0., shank_length],
                u: [1., 0., 0.],
                v: [0., 1., 0.],
                n: [0., 0., 1.],
            };
            let threaded = crate::threads::build_resolved(
                feature,
                *thread,
                ThreadMode::NewExternal,
                None,
                policy,
            )
            .map_err(|e| error(Some(feature), format!("Screw standalone thread: {e}")))?
            .with_feature_cap_roles(feature, "screw.thread/start", "screw.thread/end")
            .map_err(|e| error(Some(feature), e))?;
            result = screw_boolean(feature, &result, &threaded, BooleanOperation::Union, policy)
                .map_err(|e| error(Some(feature), format!("Screw thread fuse: {e}")))?;
            // New recesses may legitimately reach through a countersunk head
            // into its neck. Cut the complete head/shank/thread body so a later
            // union cannot refill the retained cavity.
            if !early_socket
                && let Some(cutter) = screw_profiles::cutter(feature, &head, &drive, policy)?
            {
                result = screw_boolean(
                    feature,
                    &result,
                    &cutter,
                    BooleanOperation::Subtract,
                    policy,
                )?;
            }
            let topology = result
                .topology_snapshot(feature, 0)
                .map_err(|e| error(Some(feature), e))?;
            if topology.elements.iter().any(|e| {
                e.key.kind == multicad_core::references::TopologyKind::Face && e.origins.is_empty()
            }) {
                return Err(error(
                    Some(feature),
                    "Machine screw has incomplete native face ancestry",
                ));
            }
            let completed = result
                .transformed(RigidTransform {
                    rotation: std::array::from_fn(|i| [frame.u[i], frame.v[i], frame.n[i]]),
                    translation_mm: frame.origin,
                })
                .map_err(|e| error(Some(feature), format!("Screw world transform: {e}")))?;
            remember_screw(factory_key, feature, &completed);
            Ok(completed)
        }
        ResolvedMechanicalPart::Washer {
            bore_radius_mm,
            outside_radius_mm,
            thickness_mm,
            frame,
        } => ring(
            feature,
            frame,
            bore_radius_mm,
            outside_radius_mm,
            thickness_mm,
            policy,
        ),
        ResolvedMechanicalPart::BearingEnvelope {
            bore_radius_mm,
            outside_radius_mm,
            width_mm,
            frame,
        } => ring(
            feature,
            frame,
            bore_radius_mm,
            outside_radius_mm,
            width_mm,
            policy,
        ),
        ResolvedMechanicalPart::HexNut {
            across_flats_mm,
            thread,
        } => {
            let radius = across_flats_mm / 3_f64.sqrt();
            let point = |i: usize| {
                let angle = std::f64::consts::PI / 6. + i as f64 * std::f64::consts::TAU / 6.;
                [radius * angle.cos(), radius * angle.sin()]
            };
            let outer = (0..6)
                .map(|i| CurvedEdge2d::Line {
                    a: point(i),
                    b: point((i + 1) % 6),
                })
                .collect();
            let blank = extrude(
                feature,
                thread.frame,
                thread.length_mm,
                outer,
                vec![],
                policy,
            )?;
            crate::threads::build_resolved(
                feature,
                *thread,
                ThreadMode::CutInternal,
                Some(&blank),
                policy,
            )
        }
    }
}

#[cfg(test)]
mod screw_factory_tests {
    use super::*;

    use multicad_core::{
        document::Document,
        model::{self, BodyOutput, ThreadRepresentation},
        parameters,
    };

    fn source() -> MechanicalPartDefinition {
        MechanicalPartDefinition {
            output: BodyOutput {
                body: EntityId::new(),
            },
            parameters: model::default_machine_screw(),
            provenance: None,
        }
    }

    #[test]
    fn exact_screw_key_covers_resolved_fields_pose_representation_and_policy() {
        let document = Document::new("Exact factory keys").unwrap();
        let values = parameters::evaluate(&document).unwrap();
        let resolved = source().parameters.resolve(&values).unwrap();
        let policy = GeometryPolicy::default();
        let base = screw_key(&resolved, policy).unwrap();
        let mut changed = resolved.clone();
        if let ResolvedMechanicalPart::MachineScrew { thread, .. } = &mut changed {
            thread.frame.origin[0] = -0.;
        }
        assert!(
            screw_key(&changed, policy).unwrap() != base,
            "Exact zero sign must not alias"
        );
        for field in 0..10 {
            let mut changed = resolved.clone();
            let ResolvedMechanicalPart::MachineScrew {
                underhead_length_mm,
                shank_radius_mm,
                head,
                drive,
                thread,
                ..
            } = &mut changed
            else {
                unreachable!()
            };
            let value = match field {
                0 => underhead_length_mm,
                1 => shank_radius_mm,
                2 => &mut thread.length_mm,
                3 => &mut thread.pitch_mm,
                4 => &mut thread.profile.depth_mm,
                5 => &mut thread.frame.origin[2],
                6 => &mut thread.frame.u[0],
                7 => match head {
                    ResolvedScrewHead::SocketCap { radius_mm, .. } => radius_mm,
                    _ => unreachable!(),
                },
                8 => match drive {
                    ResolvedScrewDrive::HexSocket { depth_mm, .. } => depth_mm,
                    _ => unreachable!(),
                },
                _ => &mut thread.minimum_wall_mm,
            };
            *value = f64::from_bits(value.to_bits() + 1);
            assert!(screw_key(&changed, policy).unwrap() != base);
        }
        let mut cosmetic = resolved.clone();
        if let ResolvedMechanicalPart::MachineScrew { thread, .. } = &mut cosmetic {
            thread.representation = ThreadRepresentation::CosmeticEnvelope;
        }
        assert!(screw_key(&cosmetic, policy).unwrap() != base);
        for field in 0..4 {
            let mut p = policy;
            let value = match field {
                0 => &mut p.linear_tolerance_mm,
                1 => &mut p.boolean_fuzzy_mm,
                2 => &mut p.minimum_feature_mm,
                _ => &mut p.maximum_coordinate_mm,
            };
            *value = f64::from_bits(value.to_bits() + 1);
            assert!(screw_key(&resolved, p).unwrap() != base);
        }
    }
}
