use super::architecture;
use crate::{
    commands::Prepared,
    document::{Document, EntityId, Placement},
    project_file::ContentHash,
};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SymbolKind {
    Static,
    Dynamic,
}
/// The author-defined simplified plan of a reusable symbol: a captured
/// sketch in the definition plus its explicit alignment to the symbol's local
/// coordinates. Plan views draw this fixed linework instead of projecting the
/// detailed 3D content; 3D views keep the detailed model.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SymbolPlan {
    pub sketch: EntityId,
    /// Identity alignment is omitted from storage, so definitions authored
    /// before alignment existed keep their pins.
    #[serde(default, skip_serializing_if = "PlanAlignment::is_identity")]
    pub alignment: PlanAlignment,
}
impl SymbolPlan {
    pub fn new(sketch: EntityId) -> Self {
        Self {
            sketch,
            alignment: PlanAlignment::default(),
        }
    }
    pub fn aligned(sketch: EntityId, alignment: PlanAlignment) -> Self {
        Self { sketch, alignment }
    }
}
/// A rigid in-plane placement of the plan sketch in the symbol's local
/// coordinates: rotate the sketch about local Z by `rotation_deg`
/// (counter-clockwise, seen from +Z), then shift it by `offset_mm` along
/// local X and Y. True size; plan symbols are never scaled.
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PlanAlignment {
    pub offset_mm: [f64; 2],
    pub rotation_deg: f64,
}
/// The largest authored plan offset from the symbol origin.
pub const MAX_PLAN_OFFSET_MM: f64 = 100_000.;
impl PlanAlignment {
    pub fn is_identity(&self) -> bool {
        self.offset_mm == [0., 0.] && self.rotation_deg == 0.
    }
    pub fn validate(&self) -> Result<(), crate::document::DocumentError> {
        if !self
            .offset_mm
            .iter()
            .all(|v| v.is_finite() && v.abs() <= MAX_PLAN_OFFSET_MM)
            || !self.rotation_deg.is_finite()
            || self.rotation_deg.abs() > 360.
        {
            return Err(crate::document::DocumentError::Invalid(format!(
                "Plan alignment needs finite offsets within ±{MAX_PLAN_OFFSET_MM} mm and a rotation within ±360°"
            )));
        }
        Ok(())
    }
    /// The local rigid transform from plan-sketch coordinates to symbol
    /// coordinates; compose the occurrence pose after it.
    pub fn pose(&self) -> crate::model::RigidPose {
        use multicad_precision::PortableMath as _;
        let (s, c) = self.rotation_deg.to_radians().portable_sin_cos();
        crate::model::RigidPose {
            rotation: [[c, -s, 0.], [s, c, 0.], [0., 0., 1.]],
            translation_mm: [self.offset_mm[0], self.offset_mm[1], 0.],
        }
    }
}
/// A placed instance's authored plan, resolved to canonical identities.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SymbolPlanRepresentation {
    /// The placed canonical object (the outer Instance) both views pick.
    pub instance: EntityId,
    /// The instance-owned copy of the definition's plan sketch.
    pub sketch: EntityId,
    pub alignment: PlanAlignment,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SymbolDefinition {
    pub symbol: EntityId,
    pub revision: String,
    pub name: String,
    pub kind: SymbolKind,
    pub template: Document,
    pub root: EntityId,
    pub exposed: Vec<EntityId>,
    pub plan: Option<SymbolPlan>,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SymbolValue {
    pub value: f64,
    pub unit: String,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SymbolInstance {
    pub definition: EntityId,
    pub pin: ContentHash,
    pub instance: EntityId,
    pub wrapper: EntityId,
    pub wrapper_metadata: EntityId,
    pub member: EntityId,
    pub ids: BTreeMap<EntityId, EntityId>,
    pub overrides: BTreeMap<EntityId, SymbolValue>,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SymbolAuthorRequest {
    pub definition: EntityId,
    pub symbol: EntityId,
    pub revision: String,
    pub name: String,
    pub kind: SymbolKind,
    pub root: EntityId,
    pub exposed: Vec<EntityId>,
    pub plan: Option<SymbolPlan>,
}
impl SymbolAuthorRequest {
    pub fn new(
        name: impl Into<String>,
        root: EntityId,
        kind: SymbolKind,
        exposed: Vec<EntityId>,
        plan: Option<SymbolPlan>,
    ) -> Self {
        Self {
            definition: EntityId::new(),
            symbol: EntityId::new(),
            revision: "1".into(),
            name: name.into(),
            kind,
            root,
            exposed,
            plan,
        }
    }
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SymbolPlaceRequest {
    pub record: EntityId,
    pub definition: EntityId,
    pub instance: EntityId,
    pub wrapper: EntityId,
    pub wrapper_metadata: EntityId,
    pub member: EntityId,
    pub ids: BTreeMap<EntityId, EntityId>,
    pub name: String,
    pub overrides: BTreeMap<EntityId, SymbolValue>,
    pub placement: Placement,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SymbolEditRequest {
    pub record: EntityId,
    pub overrides: BTreeMap<EntityId, SymbolValue>,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SymbolUpdateRequest {
    pub definition: EntityId,
    pub replacement: SymbolDefinition,
    pub maps: BTreeMap<EntityId, BTreeMap<EntityId, EntityId>>,
    pub instance_overrides: BTreeMap<EntityId, BTreeMap<EntityId, SymbolValue>>,
    #[serde(default)]
    pub architecture_inputs: BTreeMap<EntityId, architecture::ArchitectureInputs>,
    #[serde(default)]
    pub architecture_contexts: BTreeMap<EntityId, architecture::ArchitectureContext>,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SymbolOverrideConflict {
    pub instance: EntityId,
    pub parameter: EntityId,
    pub reason: String,
}
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SymbolUpgradeReport {
    pub old_pin: ContentHash,
    pub new_pin: ContentHash,
    pub affected_instances: Vec<EntityId>,
    pub added_parameters: Vec<EntityId>,
    pub removed_parameters: Vec<EntityId>,
    pub changed_parameters: Vec<EntityId>,
    pub changed_entities: Vec<EntityId>,
    pub conflicts: Vec<SymbolOverrideConflict>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub architecture_changes: Vec<architecture::ArchitectureUpgradeChange>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub architecture_error: Option<String>,
}
pub struct SymbolPrepared {
    pub prepared: Prepared,
    pub definition: EntityId,
    pub instances: Vec<EntityId>,
    pub plans: BTreeMap<EntityId, EntityId>,
}
#[cfg(test)]
mod plan_alignment_tests {
    use super::*;
    #[test]
    fn identity_alignment_keeps_stored_plans_and_pins_of_older_definitions() {
        let sketch = EntityId::new();
        let identity = serde_json::to_value(SymbolPlan::new(sketch)).unwrap();
        assert_eq!(identity, serde_json::json!({"sketch": sketch}));
        let older: SymbolPlan =
            serde_json::from_value(serde_json::json!({"sketch": sketch})).unwrap();
        assert_eq!(older, SymbolPlan::new(sketch));
        let aligned = SymbolPlan::aligned(
            sketch,
            PlanAlignment {
                offset_mm: [12.5, -3.],
                rotation_deg: 30.,
            },
        );
        let value = serde_json::to_value(&aligned).unwrap();
        assert_eq!(
            value["alignment"],
            serde_json::json!({"offset_mm":[12.5,-3.0],"rotation_deg":30.0})
        );
        assert_eq!(
            serde_json::from_value::<SymbolPlan>(value).unwrap(),
            aligned
        );
        assert!(
            serde_json::from_value::<SymbolPlan>(
                serde_json::json!({"sketch": sketch, "alignment": {"offset_mm":[0,0],"rotation_deg":0,"scale":2}})
            )
            .is_err(),
            "plan symbols are true size; unknown fields are refused"
        );
    }
    #[test]
    fn alignment_rotates_about_local_z_before_offsetting_and_rejects_bad_values() {
        let alignment = PlanAlignment {
            offset_mm: [100., -20.],
            rotation_deg: 90.,
        };
        let p = alignment.pose().point([30., 0., 5.]);
        for (a, b) in p.into_iter().zip([100., 10., 5.]) {
            assert!((a - b).abs() < 1e-12, "{p:?}");
        }
        assert!(alignment.pose().validate().is_ok());
        for bad in [
            PlanAlignment {
                offset_mm: [f64::NAN, 0.],
                rotation_deg: 0.,
            },
            PlanAlignment {
                offset_mm: [MAX_PLAN_OFFSET_MM * 2., 0.],
                rotation_deg: 0.,
            },
            PlanAlignment {
                offset_mm: [0., 0.],
                rotation_deg: 361.,
            },
            PlanAlignment {
                offset_mm: [0., 0.],
                rotation_deg: f64::INFINITY,
            },
        ] {
            assert!(bad.validate().is_err(), "{bad:?}");
        }
    }
}
