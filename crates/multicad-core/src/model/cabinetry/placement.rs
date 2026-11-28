//! Panel position and orientation as the inspector shows them (CAD-430):
//! the world position of the panel's origin corner, the direction its face
//! looks and its turn about that direction, all measured where the Part is
//! now, after every Move / rotate in its history.
//!
//! An edit is one canonical transaction. A panel without later placement
//! features whose new orientation is one of its own principal planes keeps
//! that placement in its panel feature (plane and origin). Any other edit
//! keeps the panel feature and carries the panel by a Move / rotate at the
//! end of its history: the trailing literal Move / rotate is edited in
//! place, otherwise one is appended. Joints follow either way (CAD-431).
use super::*;
use crate::{
    commands::{Actor, CommandError, Request, Session},
    model::{
        BodyInput, ModelPrepared, NumericTransform, RigidPose, TransformRequest,
        common::{error, revision},
        length,
    },
    parameters::EvaluatedParameters,
};
use multicad_precision::PortableMath as _;

/// The six directions a panel face can look along, as presets.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum PanelFacing {
    PosZ,
    NegZ,
    PosY,
    NegY,
    PosX,
    NegX,
}
impl PanelFacing {
    pub const ALL: [Self; 6] = [
        Self::PosZ,
        Self::NegZ,
        Self::PosY,
        Self::NegY,
        Self::PosX,
        Self::NegX,
    ];
    /// The panel's length (u) and width (v) directions at rotation 0; the
    /// thickness runs along u × v, which is this direction. +Z, −Y and +X
    /// are the XY, XZ and YZ panel planes.
    pub fn frame(self) -> ([f64; 3], [f64; 3]) {
        match self {
            Self::PosZ => ([1., 0., 0.], [0., 1., 0.]),
            Self::NegZ => ([1., 0., 0.], [0., -1., 0.]),
            Self::NegY => ([1., 0., 0.], [0., 0., 1.]),
            Self::PosY => ([1., 0., 0.], [0., 0., -1.]),
            Self::PosX => ([0., 1., 0.], [0., 0., 1.]),
            Self::NegX => ([0., 1., 0.], [0., 0., -1.]),
        }
    }
    pub fn normal(self) -> [f64; 3] {
        let (u, v) = self.frame();
        cross(u, v)
    }
    pub fn label(self) -> &'static str {
        match self {
            Self::PosZ => "+Z · face up",
            Self::NegZ => "−Z · face down",
            Self::PosY => "+Y · face back",
            Self::NegY => "−Y · face front",
            Self::PosX => "+X · face right",
            Self::NegX => "−X · face left",
        }
    }
}
/// Where a rectangle panel is: its origin corner in world millimetres, the
/// preset its face looks along (None when it is turned off every principal
/// direction) and its turn in degrees about that direction, counter-
/// clockwise looking against it, from the preset's length direction.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PanelPlacement {
    pub origin_mm: [f64; 3],
    pub facing: Option<PanelFacing>,
    pub rotation_deg: f64,
}
/// An inspector placement edit: the origin corner's world position, the
/// face direction and the turn about it. Dimensions may bind parameters.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PanelPlacementRequest {
    pub position: [Dimension; 3],
    pub facing: PanelFacing,
    pub rotation: Dimension,
}
impl PanelPlacementRequest {
    pub fn resolve(
        &self,
        parameters: &EvaluatedParameters,
    ) -> Result<PanelPlacement, DocumentError> {
        let mm = |d: &Dimension| crate::model::common::number(d, Unit::Millimeter, parameters);
        Ok(PanelPlacement {
            origin_mm: [
                mm(&self.position[0])?,
                mm(&self.position[1])?,
                mm(&self.position[2])?,
            ],
            facing: Some(self.facing),
            rotation_deg: crate::model::common::number(&self.rotation, Unit::Radian, parameters)?
                .to_degrees(),
        })
    }
}
const TOL: f64 = 1e-9;

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|i| a[i] * b[i]).sum()
}
fn near(a: [f64; 3], b: [f64; 3]) -> bool {
    (0..3).all(|i| (a[i] - b[i]).abs() <= TOL)
}
/// Round values within 1e-9 of a whole multiple of `step` onto it, so exact
/// presets read and write back exactly.
fn snap(x: f64, step: f64) -> f64 {
    let k = (x / step).round();
    if (x - k * step).abs() <= TOL {
        k * step
    } else {
        x
    }
}
fn rotate(v: [f64; 3], axis: [f64; 3], deg: f64) -> [f64; 3] {
    let (sin, cos) = crate::model::common::cardinal_sin_cos(deg.to_radians());
    let c = cross(axis, v);
    let d = dot(axis, v);
    std::array::from_fn(|i| v[i] * cos + c[i] * sin + axis[i] * d * (1. - cos))
}

/// The placement of a panel frame: origin corner, length and width directions.
pub fn placement_of(origin: [f64; 3], u: [f64; 3], v: [f64; 3]) -> PanelPlacement {
    let n = cross(u, v);
    let facing = PanelFacing::ALL.into_iter().find(|f| near(f.normal(), n));
    let rotation_deg = facing.map_or(0., |f| {
        let (ru, rv) = f.frame();
        let deg = dot(u, rv).portable_atan2(dot(u, ru)).to_degrees();
        let deg = snap(deg, 90.);
        if deg <= -180. + TOL { 180. } else { deg }
    });
    PanelPlacement {
        origin_mm: origin.map(|x| snap(x, 1.)),
        facing,
        rotation_deg,
    }
}
/// The length and width directions a placement asks for.
pub fn placement_frame(placement: &PanelPlacement) -> Result<([f64; 3], [f64; 3]), DocumentError> {
    let facing = placement
        .facing
        .ok_or_else(|| invalid("Choose the direction the panel face looks along"))?;
    if !placement.rotation_deg.is_finite()
        || placement.rotation_deg.abs() > 1e6
        || placement
            .origin_mm
            .iter()
            .any(|x| !x.is_finite() || x.abs() > 1e9)
    {
        return Err(invalid("Panel position or rotation is out of range"));
    }
    let (u, v) = facing.frame();
    let n = facing.normal();
    let clean = |w: [f64; 3]| w.map(|c| snap(c, 1.));
    Ok((
        clean(rotate(u, n, placement.rotation_deg)),
        clean(rotate(v, n, placement.rotation_deg)),
    ))
}

/// Where a rectangle panel Part is now.
pub fn panel_placement(
    document: &Document,
    parameters: &EvaluatedParameters,
    part: EntityId,
) -> Result<PanelPlacement, DocumentError> {
    let b = super::joint::panel_box(document, parameters, part)?;
    Ok(placement_of(b.origin, b.u, b.v))
}

/// The rigid pose with rotation `[u v n] · [u0 v0 n0]ᵀ` carrying `o0` to `o`.
fn pose_between(
    from: ([f64; 3], [f64; 3], [f64; 3]),
    to: ([f64; 3], [f64; 3], [f64; 3]),
) -> Result<RigidPose, DocumentError> {
    let (o0, u0, v0) = from;
    let (o, u, v) = to;
    let (n0, n) = (cross(u0, v0), cross(u, v));
    let rotation: [[f64; 3]; 3] = std::array::from_fn(|i| {
        std::array::from_fn(|j| snap(u[i] * u0[j] + v[i] * v0[j] + n[i] * n0[j], 1.))
    });
    let mut pose = RigidPose {
        rotation,
        translation_mm: [0.; 3],
    };
    let r = pose.vector(o0);
    pose.translation_mm = std::array::from_fn(|i| o[i] - r[i]);
    pose.validate()?;
    Ok(pose)
}
/// Literal Move / rotate parameters for a rigid pose: axis-angle about the
/// world origin, then a translation.
pub fn numeric_transform(pose: RigidPose) -> Result<NumericTransform, DocumentError> {
    pose.validate()?;
    let r = pose.rotation;
    let trace = r[0][0] + r[1][1] + r[2][2];
    let cos = ((trace - 1.) / 2.).clamp(-1., 1.);
    let angle = snap(cos.portable_acos().to_degrees(), 90.);
    let skew = [r[2][1] - r[1][2], r[0][2] - r[2][0], r[1][0] - r[0][1]];
    let len = dot(skew, skew).sqrt();
    let axis = if angle.abs() <= TOL {
        [0., 0., 1.]
    } else if len > 1e-6 {
        skew.map(|c| snap(c / len, 1.))
    } else {
        // A half turn: the axis is the column of R + I with the largest norm.
        let k = (0..3)
            .max_by(|&a, &b| r[a][a].total_cmp(&r[b][b]))
            .expect("three axes");
        let col: [f64; 3] = std::array::from_fn(|i| r[i][k] + if i == k { 1. } else { 0. });
        let len = dot(col, col).sqrt();
        col.map(|c| snap(c / len, 1.))
    };
    Ok(NumericTransform {
        translation: pose.translation_mm.map(|x| length(snap(x, 1.))),
        pivot: std::array::from_fn(|_| length(0.)),
        axis,
        angle: angle_dimension(angle),
    })
}
fn angle_dimension(deg: f64) -> Dimension {
    crate::model::angle(deg)
}
fn literal(transform: &NumericTransform) -> bool {
    transform
        .dimensions()
        .into_iter()
        .all(|d| matches!(d, Dimension::Literal { .. }))
}

/// Edit a rectangle panel and place it, in one transaction: `panel` is the
/// edited panel definition (its own plane and origin are replaced only when
/// the placement stays in the panel feature) and `placement` the target.
pub fn prepare_panel_placement(
    session: &Session,
    expected_revision: u64,
    actor: Actor,
    feature: EntityId,
    panel: &PanelDefinition,
    request: &PanelPlacementRequest,
) -> Result<ModelPrepared, CommandError> {
    revision(session, expected_revision)?;
    let document = session.document();
    let parameters = session.parameters();
    let part = panel.output.body;
    let placement = &request.resolve(parameters)?;
    let PanelSource::Rectangle { .. } = &panel.source else {
        return Err(error(
            "Only rectangle panels are positioned here; move other panels with Move / rotate",
        ));
    };
    let (u, v) = placement_frame(placement)?;
    let target = (placement.origin_mm, u, v);
    let placements = super::joint::panel_placements(document, part);
    let principal = [PanelPlane::Xy, PanelPlane::Xz, PanelPlane::Yz]
        .into_iter()
        .find(|p| {
            let (pu, pv) = p.axes();
            near(pu, u) && near(pv, v)
        });
    let mut edited = panel.clone();
    let mut in_feature = false;
    if placements.is_empty()
        && let Some(p) = principal
    {
        // The whole placement stays in the panel feature.
        let PanelSource::Rectangle {
            plane: new_plane,
            origin: new_origin,
            ..
        } = &mut edited.source
        else {
            unreachable!()
        };
        *new_plane = p;
        *new_origin = request.position.clone();
        in_feature = true;
    }
    // Stage the panel edit, then the Move / rotate, on derived sessions.
    let stage = Session::new(document.clone())?;
    let edit = prepare_panel_edit(&stage, stage.revision(), actor.clone(), feature, &edited)?;
    let mut candidate = edit.prepared.document().clone();
    // The panel feature is the edit's one output, so a feature-edit scope
    // keeps its producer identity; a trailing Move / rotate rides along.
    let outputs = edit.outputs;
    if !in_feature {
        let stage = Session::new(candidate.clone())?;
        let b = super::joint::unplaced_panel_box(stage.document(), stage.parameters(), part)?;
        let want = pose_between((b.origin, b.u, b.v), target)?;
        let current =
            super::joint::panel_pose_after(stage.document(), stage.parameters(), part, feature)?;
        let EntityData::Part { features, .. } = &stage.document().get(part).expect("Part").data
        else {
            return Err(error("The panel output must be a Part"));
        };
        let last = *features.last().expect("panel feature");
        let trailing = crate::model::read_transform(stage.document(), last)
            .ok()
            .filter(|t| literal(&t.parameters))
            .filter(|_| {
                matches!(
                    stage.document().get(last).map(|e| &e.data),
                    Some(EntityData::Feature {
                        suppressed: false,
                        ..
                    })
                )
            });
        let result = if let Some(t) = trailing {
            // Everything before the trailing Move / rotate stays; only it changes.
            let own = t.parameters.resolve(stage.parameters())?;
            let before = current.then(own.inverse()?)?;
            let parameters = numeric_transform(before.inverse()?.then(want)?)?;
            crate::model::prepare_transform_edit(
                &stage,
                stage.revision(),
                actor.clone(),
                last,
                &parameters,
            )?
        } else {
            let parameters = numeric_transform(current.inverse()?.then(want)?)?;
            let mut request = TransformRequest::new(
                BodyInput {
                    body: part,
                    producer: last,
                },
                parameters,
            );
            request.name = "Panel placement".into();
            crate::model::prepare_transform(&stage, stage.revision(), actor.clone(), &request)?
        };
        candidate = result.prepared.document().clone();
    }
    let prepared = session.prepare(Request {
        expected_revision,
        actor,
        label: "Place panel".into(),
        edits: crate::model::extrude_boolean::delta(document, &candidate)?,
    })?;
    Ok(ModelPrepared { prepared, outputs })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_preset_reads_back_with_its_turns() {
        for facing in PanelFacing::ALL {
            for rotation_deg in [0., 90., 180., -90., 30.] {
                let p = PanelPlacement {
                    origin_mm: [1., 2., 3.],
                    facing: Some(facing),
                    rotation_deg,
                };
                let (u, v) = placement_frame(&p).unwrap();
                let back = placement_of(p.origin_mm, u, v);
                assert_eq!(back.facing, Some(facing));
                let want = if rotation_deg == -90. {
                    -90.
                } else {
                    rotation_deg
                };
                assert!(
                    (back.rotation_deg - want).abs() < 1e-9,
                    "{facing:?} {back:?}"
                );
            }
        }
        // The three principal planes are the presets at rotation 0.
        for (plane, facing) in [
            (PanelPlane::Xy, PanelFacing::PosZ),
            (PanelPlane::Xz, PanelFacing::NegY),
            (PanelPlane::Yz, PanelFacing::PosX),
        ] {
            assert_eq!(plane.axes(), facing.frame());
        }
    }

    #[test]
    fn literal_transforms_round_trip_poses() {
        let parameters =
            crate::parameters::evaluate(&Document::new("Transforms").unwrap()).unwrap();
        for (axis, deg) in [
            ([0., 0., 1.], 90.),
            ([1., 0., 0.], 180.),
            ([0., 1., 0.], -90.),
            ([1., 1., 1.], 120.),
            ([0., 0., 1.], 0.),
            ([0., 1., 0.], 30.),
        ] {
            let t = NumericTransform {
                translation: [length(5.), length(-7.), length(11.)],
                axis,
                angle: crate::model::angle(deg),
                ..Default::default()
            };
            let pose = t.resolve(&parameters).unwrap();
            let again = numeric_transform(pose)
                .unwrap()
                .resolve(&parameters)
                .unwrap();
            for i in 0..3 {
                assert!(near(pose.rotation[i], again.rotation[i]), "{axis:?} {deg}");
            }
            assert!(near(pose.translation_mm, again.translation_mm));
        }
    }
}
