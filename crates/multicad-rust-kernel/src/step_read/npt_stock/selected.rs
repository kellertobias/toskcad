//! Selection retains neutral represented ownership independently of geometry.
use super::*;
use crate::step_read::owned_insert::{MappingLink, Member, Options as MappingOptions, mapping};
use crate::tapered_thread::interval::I;
mod placement;
pub use placement::Placement;
pub struct Selected<'document, 'input> {
    parsed: Parsed<'document, 'input>,
    member: Member,
}
impl<'document, 'input> Selected<'document, 'input> {
    pub fn parse(
        document: &'document CompoundDocument<'input>,
        item: u32,
        options: Options,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self> {
        let (parsed, member) = parse::selected_body(document, item, options, interrupt)?;
        Ok(Self { parsed, member })
    }
    pub fn parsed(&self) -> &Parsed<'document, 'input> {
        &self.parsed
    }
    /// Sorted original dependency IDs; these retain the selected source's
    /// provenance without copying or reparsing the neutral AST.
    pub fn entity_ids(&self) -> &[u32] {
        &self.parsed.dependencies
    }
    pub fn member(&self) -> Member {
        self.member
    }
}
#[cfg(test)]
mod tests;
