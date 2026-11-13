//! Exact finite surfaces for a retained head, sharing every thread start rim.
use super::*;
use crate::{
    head_profile::{Curve, Role},
    profiled_screw::Screw,
};
use std::f64::consts::TAU;
mod outer;
mod source;
pub(crate) use source::HeadSource;
pub(super) mod phillips;
#[cfg(test)]
mod tests;
mod torx;

#[derive(Clone, Copy)]
struct Station {
    center: Entity,
    axis: Option<Entity>,
    rim: Option<Entity>,
    vertex: Vertex,
    radius: f64,
}
impl From<CircularRim> for Station {
    fn from(rim: CircularRim) -> Self {
        Self {
            center: rim.center,
            axis: Some(rim.axis),
            rim: Some(rim.edge),
            vertex: rim.vertex,
            radius: rim.radius,
        }
    }
}

pub(super) struct Junction {
    pub(super) crests: Vec<Entity>,
    pub(super) complements: Vec<Entity>,
    vertex: Vertex,
    reference: Point,
    positive: Vec<(Entity, bool)>,
    pub(super) outer: Option<CircularRim>,
}

pub(super) fn admission(screw: &Screw, origin: Point, tolerance: f64) -> Result<(), StepError> {
    admission_source(HeadSource::from(screw), origin, tolerance)
}
pub(super) fn admission_source(
    source: HeadSource<'_>,
    origin: Point,
    tolerance: f64,
) -> Result<(), StepError> {
    outer::admission(source.data, origin, tolerance)?;
    match source.kind {
        source::Kind::Existing(screw) => match screw.recess() {
            crate::screw_recess::Recess::FlatHex(hex) => {
                let (drive, socket_depth) = hex.geometry();
                if socket_depth <= 10. * tolerance
                    || drive.across_flats_mm() <= 10. * tolerance
                    || !screw.admits_recess_envelope(
                        drive.across_flats_mm() / 3_f64.sqrt(),
                        socket_depth,
                        10. * tolerance,
                    )
                {
                    return Err(StepError::Limit);
                }
            }
            crate::screw_recess::Recess::ConicalHexalobular(recess) => {
                torx::admission(screw, recess, tolerance)?
            }
        },
        source::Kind::Phillips(screw) => {
            screw
                .cavity()
                .certify_host_envelope(
                    screw.profile(),
                    screw.thread(),
                    screw.spec().shank_length_mm,
                    10. * tolerance.min(MAX_ERROR_MM),
                    100_000,
                    &|| false,
                )
                .map_err(phillips::source_error)?;
        }
    }
    Ok(())
}

pub(super) fn placement(screw: &Screw, p: &Placement, tolerance: f64) -> Result<(), StepError> {
    placement_source(HeadSource::from(screw), p, tolerance)
}
pub(super) fn placement_source(
    source: HeadSource<'_>,
    p: &Placement,
    tolerance: f64,
) -> Result<(), StepError> {
    outer::placement(source.data, p, tolerance)?;
    match source.kind {
        source::Kind::Existing(screw) => torx::placement(screw, p),
        source::Kind::Phillips(_) => Ok(()),
    }
}

pub(super) fn junction(
    file: &mut StepFile,
    p: &Placement,
    thread: &Thread,
    shank_radius: f64,
    z: f64,
    corners: &[[Vertex; 2]],
) -> Result<Junction, StepError> {
    let spec = thread.spec();
    let [_, start, _, _, crest, _] = thread.profile_breaks();
    let rate = thread.helical_rate();
    let radial = |phase: f64| {
        let (sin, cos) = ((z - phase) / rate).sin_cos();
        p.vector([cos, sin, 0.])
    };
    let center = p.point_f64([0.; 2], z);
    let axis = p.normal.map(|x| if spec.left_hand { x } else { -x });
    let mut crests = Vec::new();
    let mut complements = Vec::new();
    for (k, &[first, last]) in corners.iter().enumerate() {
        let shift = k as f64 * spec.pitch_mm;
        let next = corners[(k + 1) % corners.len()][0];
        crests.push(
            circular_arc_edge(
                file,
                [last, next],
                center,
                axis,
                radial(crest + shift),
                spec.major_radius_mm,
                (start + spec.pitch_mm - crest) / rate.abs(),
            )?
            .0,
        );
        if shank_radius == spec.major_radius_mm {
            complements.push(
                circular_arc_edge(
                    file,
                    [first, last],
                    center,
                    axis,
                    radial(start + shift),
                    spec.major_radius_mm,
                    (crest - start) / rate.abs(),
                )?
                .0,
            );
        }
    }
    let mut positive: Vec<_> = complements
        .iter()
        .zip(&crests)
        .flat_map(|(&a, &b)| [(a, true), (b, true)])
        .collect();
    if !spec.left_hand {
        positive = positive.into_iter().rev().map(|(e, s)| (e, !s)).collect();
    }
    let outer = if shank_radius > spec.major_radius_mm {
        Some(CircularRim::new(
            file,
            center,
            p.normal,
            radial(start),
            shank_radius,
        )?)
    } else {
        None
    };
    if let Some(outer) = &outer {
        positive = vec![(outer.edge, true)];
    }
    Ok(Junction {
        crests,
        complements,
        vertex: corners[0][0],
        reference: radial(start),
        positive,
        outer,
    })
}

pub(super) fn head(
    file: &mut StepFile,
    p: &Placement,
    source: HeadSource<'_>,
    join: Junction,
    tolerance_mm: f64,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<Vec<Entity>, StepError> {
    match source.kind {
        source::Kind::Existing(screw) => existing_head(file, p, screw, join, tolerance_mm),
        source::Kind::Phillips(screw) => {
            let shell = outer::shell(file, p, source.data, join)?;
            let pose = crate::step::phillips_cutter::support::Pose::new(
                p.point_f64([0.; 2], screw.head_top_z()),
                p.u,
                p.v,
            )
            .map_err(phillips::plan_error)?;
            let plan = crate::step::phillips_cutter::plan::Plan::new(
                screw.cavity().source(),
                &pose,
                tolerance_mm,
                interrupt,
            )
            .map_err(phillips::plan_error)?;
            let cavity = plan
                .emit_reversed_cavity(file, interrupt)
                .map_err(phillips::plan_error)?;
            let outer = circular_bound(file, shell.top_edge, false, true);
            let top = file.face_on(
                &[outer, cavity.mouth],
                cavity.mouth_surface,
                cavity.mouth_forward,
            );
            let mut faces = vec![top];
            faces.extend(shell.walls);
            faces.extend(cavity.faces);
            Ok(faces)
        }
    }
}
fn existing_head(
    file: &mut StepFile,
    p: &Placement,
    screw: &Screw,
    join: Junction,
    tolerance_mm: f64,
) -> Result<Vec<Entity>, StepError> {
    let data = outer::Data::from(screw);
    let outer::Shell {
        top_edge,
        top_center,
        walls,
    } = outer::shell(file, p, data, join)?;
    let top = data.top_z();
    let (mouth, mut cavity) = match screw.recess() {
        crate::screw_recess::Recess::FlatHex(recess) => {
            let (drive, socket_depth) = recess.geometry();
            let ([mouth, floor], sides) =
                encode_drive_prism(file, p, drive.corners_mm(), [top, top + socket_depth], true)?;
            let floor_center = file.point(p.point_f64([0.; 2], top + socket_depth));
            let mut faces = vec![file.face(&[floor], floor_center, p.normal.map(|x| -x), p.u)];
            faces.extend(sides);
            (mouth, faces)
        }
        crate::screw_recess::Recess::ConicalHexalobular(recess) => {
            torx::encode(file, p, recess, top, tolerance_mm)?
        }
    };
    let outer = circular_bound(file, top_edge, false, true);
    let down = p.normal.map(|x| -x);
    let mut faces = vec![file.face(&[outer, mouth], top_center, down, p.u)];
    faces.extend(walls);
    faces.append(&mut cavity);
    Ok(faces)
}

fn seam(
    file: &mut StepFile,
    surface: Entity,
    circle: Entity,
    [start, end]: [Entity; 2],
    v: f64,
    sign: f64,
) -> Entity {
    let context=file.entity(format_args!("(GEOMETRIC_REPRESENTATION_CONTEXT(2)PARAMETRIC_REPRESENTATION_CONTEXT()REPRESENTATION_CONTEXT('2D SPACE',''))"));
    let curves = [0., TAU].map(|u| {
        let point = file.entity(format_args!("CARTESIAN_POINT('',({u:.16E},{v:.16E}))"));
        let direction = file.entity(format_args!("DIRECTION('',(0.,{sign:.16E}))"));
        let vector = file.entity(format_args!("VECTOR('',{direction},1.)"));
        let line = file.entity(format_args!("LINE('',{point},{vector})"));
        let representation = file.entity(format_args!(
            "DEFINITIONAL_REPRESENTATION('',({line}),{context})"
        ));
        file.entity(format_args!("PCURVE('',{surface},{representation})"))
    });
    let [a, b] = curves;
    let seam = file.entity(format_args!(
        "SEAM_CURVE('',{circle},({a},{b}),.PCURVE_S1.)"
    ));
    file.entity(format_args!("EDGE_CURVE('',{start},{end},{seam},.T.)"))
}
