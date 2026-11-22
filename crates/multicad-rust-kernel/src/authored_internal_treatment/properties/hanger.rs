//! Composition of original component measures, with original source phase and
//! symbolic axial stations. No Cavity, legacy Thread or nominal centroid proxy.
use super::*;
use crate::authored_hanger_bolt::{
    Error as HangerError, FaceKey as HangerFace, FaceProperties as HangerProperties, Hanger,
};
fn failure(e: Error) -> HangerError {
    match e {
        Error::Invalid => HangerError::Invalid,
        Error::Limit => HangerError::Limit,
        Error::Unresolved => HangerError::Unresolved,
        Error::Cancelled => HangerError::Cancelled,
    }
}
impl Hanger {
    pub fn local_face_properties(
        &self,
        face: HangerFace,
        options: FacePropertyOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<HangerProperties, HangerError> {
        hanger_query(self, face, options, interrupt).map_err(failure)
    }
}
fn hanger_query(
    source: &Hanger,
    face: HangerFace,
    options: FacePropertyOptions,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<HangerProperties, Error> {
    let mut work = integrate::Work::new(options, interrupt)?;
    if source.face_index(face).is_none() {
        return Err(Error::Invalid);
    }
    let spec = source.spec();
    let radius = I::point(spec.plain_radius_mm);
    let plain = I::point(spec.plain_length_mm);
    let shift = -plain - I::point(source.metric().full_length_mm());
    let metric_chart = metric::Chart {
        shift,
        cap_body_z: None,
        outer: None,
    };
    let values = match face {
        HangerFace::MetricTip => metric::evaluate(
            source.metric(),
            FaceKey::StartCap,
            metric::Chart {
                cap_body_z: Some(shift),
                ..metric_chart
            },
            options,
            &mut work,
        )?,
        HangerFace::MetricShoulder => metric::evaluate(
            source.metric(),
            FaceKey::EndCap,
            metric::Chart {
                cap_body_z: Some(-plain),
                outer: Some(radius),
                ..metric_chart
            },
            options,
            &mut work,
        )?,
        HangerFace::MetricThread { start, band } => metric::evaluate(
            source.metric(),
            FaceKey::Thread { start, band },
            metric_chart,
            options,
            &mut work,
        )?,
        HangerFace::WoodTip => source::evaluate(
            source::Domain::bare(source.wood()),
            None,
            FaceKey::EndCap,
            options,
            &mut work,
        )?,
        HangerFace::WoodShoulder => {
            let mut values = source::evaluate(
                source::Domain::bare(source.wood()),
                None,
                FaceKey::StartCap,
                options,
                &mut work,
            )?;
            let circle = source::circle(radius, None, 0.);
            for i in 0..4 {
                values[i] = values[i] - circle[i];
            }
            values
        }
        HangerFace::WoodThread { start, band } => source::evaluate(
            source::Domain::bare(source.wood()),
            None,
            FaceKey::Thread { start, band },
            options,
            &mut work,
        )?,
        HangerFace::PlainWall => source::cylinder(radius, -spec.plain_length_mm, 0.),
    };
    let values = finish_values(values, options)?;
    let result = HangerProperties {
        face,
        area_mm2: values.area_mm2,
        first_moment_mm3: values.first_moment_mm3,
        centroid_mm: values.centroid_mm,
    };
    work.checkpoint()?;
    Ok(result)
}
