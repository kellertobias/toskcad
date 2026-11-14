//! Paired insert source contact with one actual retained blind pilot.
use super::*;
use source::InsertContact;
mod clearance;
mod mating;
impl Solid {
    pub fn rust_wood_insert_contact(
        &self,
        insert: &Self,
        requirements: Requirements,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<InsertContact, GeometryError> {
        let fail = |kind, message| error("Wood insert contact", kind, message);
        self.policy.validate()?;
        if self.policy != insert.policy
            || self.policy.boolean_fuzzy_mm != 0.
            || requirements.linear_tolerance_mm != self.policy.linear_tolerance_mm
            || requirements.minimum_feature_mm != self.policy.minimum_feature_mm
        {
            return Err(fail(
                GeometryErrorKind::Input,
                "Contact requires the actual source geometry policies",
            ));
        }
        let unsupported = || {
            fail(
                GeometryErrorKind::Unsupported,
                "Insert contact requires an actual paired insert and one retained orthogonal blind pilot",
            )
        };
        let (stock, stock_frame) = self.orthogonal_body().ok_or_else(unsupported)?;
        let (material, pilot) = stock.wood_contact_source().ok_or_else(unsupported)?;
        let insert = insert.rust_wood_insert_data().ok_or_else(unsupported)?;
        stock_frame.checked(self.policy)?;
        insert.frame.checked(self.policy)?;
        source::certify_insert(
            Stock {
                material,
                pilot,
                frame: frame(stock_frame),
            },
            &insert.source,
            frame(insert.frame),
            requirements,
            &|| interrupt.is_some_and(|f| f()),
        )
        .map_err(source_error)
    }
}
#[cfg(test)]
mod tests;
