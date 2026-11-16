//! A nominal wood screw envelope does not recover thread or point geometry.
use super::*;
use multicad_core::model::{ResolvedWoodScrew, WoodPointIntent, WoodScrewRepresentation};
use multicad_kernel::rust_kernel::Rational;

const MAX_ENTRIES: usize = 16;
const MAX_BYTES: usize = 8 * 1024 * 1024;

struct Entry {
    solid: Solid,
    owner: EntityId,
    bytes: usize,
    used: u64,
}
#[derive(Default)]
struct Cache {
    entries: BTreeMap<String, Entry>,
    bytes: usize,
    clock: u64,
}
thread_local! {static ENVELOPES: RefCell<Cache> = RefCell::new(Cache::default());}

fn neutral_owner() -> EntityId {
    EntityId::parse("d067b697-a074-4bc7-8ae8-4c05d18b3a90")
        .expect("fixed private wood envelope owner")
}
fn authored_role(role: &str) -> bool {
    ["head", "stock"].into_iter().any(|component| {
        role == format!("wood.envelope.{component}/start")
            || role == format!("wood.envelope.{component}/end")
            || (0..if component == "head" { 6 } else { 1 })
                .any(|i| role == format!("wood.envelope.{component}/side/{i}"))
    })
}
fn valid_face_origins(origins: &std::collections::BTreeSet<Origin>, owner: EntityId) -> bool {
    let side_prefix = format!("side/{owner}/");
    !origins.is_empty()
        && origins.iter().any(|origin| authored_role(&origin.role))
        && origins.iter().all(|origin| {
            origin.source == owner
                && (authored_role(&origin.role)
                    || origin
                        .role
                        .strip_prefix(&side_prefix)
                        .is_some_and(|authored| {
                            authored.contains("/side/")
                                && authored_role(authored)
                                && origins.contains(&Origin {
                                    source: owner,
                                    role: authored.into(),
                                })
                        }))
        })
}
fn key(p: &ResolvedWoodScrew, policy: GeometryPolicy) -> String {
    // Include source-intent receipts even when they do not change the stock
    // envelope. They must never inherit another profile's completed proof.
    format!("wood-screw-envelope/1|{p:?}|{policy:?}")
}
fn cached(key: &str, feature: EntityId) -> Option<Solid> {
    ENVELOPES.with(|slot| {
        let mut cache = slot.borrow_mut();
        cache.clock = cache.clock.wrapping_add(1);
        let clock = cache.clock;
        let entry = cache.entries.get_mut(key)?;
        let solid = entry
            .solid
            .rebind_exclusive_owner(entry.owner, feature)
            .ok()?;
        entry.used = clock;
        Some(solid)
    })
}
fn remember(key: String, feature: EntityId, solid: &Solid) {
    let Ok(owned) = solid.rebind_exclusive_owner(feature, feature) else {
        return;
    };
    let Ok(brep) = owned.to_brep() else { return };
    let Ok(topology) = owned.topology_snapshot(feature, 0) else {
        return;
    };
    let Ok(receipt) = serde_json::to_vec(&topology) else {
        return;
    };
    let Some(bytes) = brep.len().checked_add(receipt.len()) else {
        return;
    };
    if bytes > MAX_BYTES {
        return;
    }
    ENVELOPES.with(|slot| {
        let mut cache = slot.borrow_mut();
        if let Some(old) = cache.entries.remove(&key) {
            cache.bytes -= old.bytes;
        }
        while cache.entries.len() >= MAX_ENTRIES || cache.bytes + bytes > MAX_BYTES {
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
        cache.clock = cache.clock.wrapping_add(1);
        let used = cache.clock;
        cache.entries.insert(
            key,
            Entry {
                solid: owned,
                owner: feature,
                bytes,
                used,
            },
        );
        cache.bytes += bytes;
    });
}
fn prism(
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
            role: format!("wood.envelope.{component}/side/{i}"),
        })
        .collect::<Vec<_>>();
    Solid::extrude_curved_profile_vector_named(
        &CurvedProfile {
            frame: PlaneFrame {
                origin: [0., 0., z],
                u: [1., 0., 0.],
                v: [0., 1., 0.],
            },
            outer,
            holes: vec![],
        },
        [0., 0., height],
        feature,
        &origins,
        policy,
    )
    .and_then(|solid| {
        solid.with_feature_cap_roles(
            feature,
            &format!("wood.envelope.{component}/start"),
            &format!("wood.envelope.{component}/end"),
        )
    })
    .map_err(|e| error(Some(feature), e))
}
// New authored stock retains the evaluator's actual binary64 dimensions.
// This is construction from parameters, never normalization of existing stock.
fn binary_dimension(feature: EntityId, value: f64) -> Result<Rational, ModelError> {
    let failed = || {
        error(
            Some(feature),
            "Wood stock dimensions exceed the exact arithmetic budget",
        )
    };
    if !value.is_finite() || value <= 0. {
        return Err(failed());
    }
    let bits = value.to_bits();
    let exponent = ((bits >> 52) & 0x7ff) as i32;
    let mut numerator = bits & ((1_u64 << 52) - 1);
    let mut shift = if exponent == 0 {
        -1074
    } else {
        exponent - 1075
    };
    if exponent != 0 {
        numerator |= 1_u64 << 52;
    }
    let zeros = numerator.trailing_zeros();
    numerator >>= zeros;
    shift += zeros as i32;
    let (n, d) = if shift >= 0 {
        (
            u128::from(numerator)
                .checked_mul(1_u128.checked_shl(shift as u32).ok_or_else(failed)?)
                .ok_or_else(failed)?,
            1_u128,
        )
    } else {
        (
            u128::from(numerator),
            1_u128.checked_shl((-shift) as u32).ok_or_else(failed)?,
        )
    };
    let n = i64::try_from(n).map_err(|_| failed())?;
    let d = i64::try_from(d).map_err(|_| failed())?;
    Rational::new(n, d).map_err(|_| failed())
}
fn rust_stock(
    feature: EntityId,
    radius: f64,
    height: f64,
    policy: GeometryPolicy,
) -> Result<Solid, ModelError> {
    let zero = Rational::new(0, 1).expect("zero");
    let circle =
        multicad_kernel::rust_kernel::Circle::new([zero; 2], binary_dimension(feature, radius)?)
            .map_err(|e| error(Some(feature), format!("Invalid wood stock circle: {e:?}")))?;
    let cylinder = multicad_kernel::rust_kernel::cylinder::AxialCylinder::new(
        feature.uuid_bits(),
        circle,
        zero,
        binary_dimension(feature, height)?,
    )
    .map_err(|e| error(Some(feature), format!("Invalid wood stock: {e:?}")))?;
    Solid::from_axial_cylinder_sweep(
        cylinder,
        PlaneFrame::default(),
        [0.; 2],
        false,
        Origin {
            source: feature,
            role: "wood.envelope.stock/side/0".into(),
        },
        policy,
    )
    .and_then(|solid| {
        solid.with_feature_cap_roles(
            feature,
            "wood.envelope.stock/start",
            "wood.envelope.stock/end",
        )
    })
    .map_err(|e| error(Some(feature), e))
}
pub(super) fn build(
    feature: EntityId,
    p: &ResolvedWoodScrew,
    policy: GeometryPolicy,
) -> Result<Solid, ModelError> {
    build_inner(feature, p, policy, false)
}
pub(super) fn build_rust(
    feature: EntityId,
    p: &ResolvedWoodScrew,
    policy: GeometryPolicy,
) -> Result<Solid, ModelError> {
    build_inner(feature, p, policy, true)
}
fn build_inner(
    feature: EntityId,
    p: &ResolvedWoodScrew,
    policy: GeometryPolicy,
    rust: bool,
) -> Result<Solid, ModelError> {
    policy.validate().map_err(|e| error(Some(feature), e))?;
    let ResolvedScrewHead::ExternalHex {
        across_flats_mm,
        height_mm,
    } = &p.head
    else {
        return Err(error(
            Some(feature),
            "Wood envelope supports only an external hex head",
        ));
    };
    let (across_flats_mm, height_mm) = (*across_flats_mm, *height_mm);
    if p.drive != ResolvedScrewDrive::ExternalHex
        || p.representation != WoodScrewRepresentation::SimplifiedEnvelope
        || p.point != WoodPointIntent::ManufacturerUnspecified
    {
        return Err(error(
            Some(feature),
            "Unsupported wood screw envelope intent",
        ));
    }
    let radius = p.nominal_diameter_mm / 2.;
    if [
        radius,
        p.underhead_length_mm,
        height_mm,
        across_flats_mm / 2. - radius,
    ]
    .into_iter()
    .any(|value| !value.is_finite() || value < policy.minimum_feature_mm)
    {
        return Err(error(
            Some(feature),
            "Wood envelope stock, head or bearing land is below native minimum feature size",
        ));
    }
    let cache_key = key(p, policy);
    if !rust && let Some(hit) = cached(&cache_key, feature) {
        return Ok(hit);
    }
    let neutral = if rust { feature } else { neutral_owner() };
    let head = prism(
        neutral,
        "head",
        -height_mm,
        height_mm,
        super::hexagon(across_flats_mm),
        policy,
    )?;
    let stock = if rust {
        rust_stock(neutral, radius, p.underhead_length_mm, policy)?
    } else {
        prism(
            neutral,
            "stock",
            0.,
            p.underhead_length_mm,
            vec![CurvedEdge2d::Circle {
                center: [0., 0.],
                radius_mm: radius,
            }],
            policy,
        )?
    };
    let solid = if head.uses_rust_geometry() && stock.uses_rust_geometry() {
        Solid::rust_join_planar_head_and_round_stock(&head, &stock, neutral, policy)
            .map_err(|e| error(Some(feature), e))?
    } else {
        head.boolean(&stock, BooleanOperation::Union, policy)
            .map_err(|e| error(Some(feature), e))?
            .ok_or_else(|| error(Some(feature), "Wood envelope union produced no solid"))?
    };
    if solid.solid_count().map_err(|e| error(Some(feature), e))? != 1 {
        return Err(error(
            Some(feature),
            "Wood envelope must remain one connected solid",
        ));
    }
    let topology = solid
        .topology_snapshot(neutral, 0)
        .map_err(|e| error(Some(feature), e))?;
    if topology.elements.iter().any(|element| {
        element.key.kind == multicad_core::references::TopologyKind::Face
            && !valid_face_origins(&element.origins, neutral)
    }) {
        return Err(error(
            Some(feature),
            "Wood envelope has incomplete native face ancestry",
        ));
    }
    let frame = p.frame;
    let completed = solid
        .transformed(RigidTransform {
            rotation: std::array::from_fn(|i| [frame.u[i], frame.v[i], frame.n[i]]),
            translation_mm: frame.origin,
        })
        .map_err(|e| error(Some(feature), e))?;
    // Rust's immutable Model evaluation cache already retains completed
    // feature geometry. Keep its actual owner instead of native-handle
    // serialization and neutral-owner rebinding.
    if rust {
        return Ok(completed);
    }
    let result = completed
        .rebind_exclusive_owner(neutral, feature)
        .map_err(|e| error(Some(feature), e))?;
    remember(cache_key, neutral, &completed);
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn resolved() -> ResolvedWoodScrew {
        let session = multicad_core::commands::Session::new(
            multicad_core::document::Document::new("Wood cache test").unwrap(),
        )
        .unwrap();
        let ResolvedMechanicalPart::WoodScrew { parameters } =
            multicad_core::model::default_wood_screw()
                .resolve(session.parameters())
                .unwrap()
        else {
            unreachable!()
        };
        *parameters
    }
    #[test]
    fn authored_binary_dimensions_are_exact_and_refuse_out_of_budget_values() {
        let feature = EntityId::new();
        for value in [0.125, 3.3, 8.7, 12345.6789] {
            assert_eq!(
                binary_dimension(feature, value).unwrap().as_f64().to_bits(),
                value.to_bits()
            );
        }
        assert_eq!(
            binary_dimension(feature, 3.3).unwrap().ratio(),
            (3715469692580659, 1125899906842624)
        );
        for value in [
            0.,
            -1.,
            f64::NAN,
            f64::INFINITY,
            1e-300,
            f64::from_bits(1),
            0.001,
            1e-20,
            1e19,
            f64::MAX,
        ] {
            assert!(binary_dimension(feature, value).is_err());
        }
    }

    #[test]
    fn envelope_receipts_key_reference_intent_pose_and_all_policy_fields() {
        let p = resolved();
        let policy = GeometryPolicy::default();
        let original = key(&p, policy);
        let mut changes = vec![];
        let mut q = p.clone();
        q.reference_profile.pitch_mm += 0.01;
        changes.push(q);
        let mut q = p.clone();
        q.reference_profile.core_diameter_mm += 0.01;
        changes.push(q);
        let mut q = p.clone();
        q.reference_profile.included_angle_deg += 0.1;
        changes.push(q);
        let mut q = p.clone();
        q.minimum_tip_inclusive_thread_extent_mm += 1.;
        q.maximum_possible_plain_length_mm -= 1.;
        changes.push(q);
        let mut q = p.clone();
        q.frame.origin[0] += 1.;
        changes.push(q);
        let mut q = p.clone();
        q.authored_length_mm += 1.;
        changes.push(q);
        let mut q = p.clone();
        q.underhead_length_mm += 1.;
        changes.push(q);
        let mut q = p.clone();
        q.nominal_diameter_mm += 0.01;
        changes.push(q);
        assert!(changes.iter().all(|q| key(q, policy) != original));
        for q in [
            GeometryPolicy {
                linear_tolerance_mm: 2e-7,
                ..policy
            },
            GeometryPolicy {
                boolean_fuzzy_mm: 1e-8,
                ..policy
            },
            GeometryPolicy {
                minimum_feature_mm: 2e-5,
                ..policy
            },
            GeometryPolicy {
                maximum_coordinate_mm: 1e6,
                ..policy
            },
        ] {
            assert_ne!(key(&p, q), original);
        }
    }
}
