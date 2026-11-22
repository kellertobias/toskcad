//! Conservative exact-decimal authority for complete analytic support laws.
//! Equality is syntactic, before binary64 conversion. It does not normalize
//! frames or merge merely close supports, and retains every raw owner separately.
use super::*;

#[derive(Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub struct ExactSupportKey(Vec<String>);
impl ExactSupportKey {
    pub fn retained_capacity_bytes(&self) -> usize {
        std::mem::size_of::<Self>()
            + self.0.capacity() * std::mem::size_of::<String>()
            + self.0.iter().map(String::capacity).sum::<usize>()
    }
}
impl Boundary {
    /// None means exact numeric provenance is unavailable; binary64 equality
    /// must never be substituted. Equal keys authorize only identical literal
    /// parameter realizations, not arbitrary independent uncertainty choices.
    pub fn exact_support_key(&self, face_index: usize) -> Option<&ExactSupportKey> {
        self.exact_support_keys.get(face_index)?.as_ref()
    }
    /// Dynamic key payload plus its vector allocation (excluding this Boundary).
    pub fn exact_support_keys_capacity_bytes(&self) -> usize {
        self.exact_support_keys.capacity() * std::mem::size_of::<Option<ExactSupportKey>>()
            + self
                .exact_support_keys
                .iter()
                .flatten()
                .map(|k| k.retained_capacity_bytes() - std::mem::size_of::<ExactSupportKey>())
                .sum::<usize>()
    }
}
fn append(key: &mut Vec<String>, entity: &Entity<'_>) -> Option<()> {
    let tokens = crate::step_read::part21::numeric_spellings(entity)?;
    key.push(entity.name.to_owned());
    key.extend(tokens.into_iter().map(str::to_owned));
    Some(())
}
pub(super) fn read(
    shell: &Shell<'_, '_>,
    id: u32,
    work: &mut Work<'_>,
) -> Result<Option<ExactSupportKey>> {
    let support = shell.simple(id)?;
    let axis = reference(support.args.get(1).ok_or(ReadError::Invalid)?)?;
    let axis = shell.simple(axis)?;
    let [Value::String(_), origin, normal, reference_axis] = axis.args.as_slice() else {
        return Err(ReadError::Invalid);
    };
    let entities = [
        support,
        axis,
        shell.simple(reference(origin)?)?,
        shell.simple(reference(normal)?)?,
        shell.simple(reference(reference_axis)?)?,
    ];
    // Charge all original text before token scans or owned key allocations.
    work.step(entities.iter().map(|e| e.source.map_or(0, str::len)).sum())?;
    let mut key = Vec::new();
    for entity in entities {
        if append(&mut key, entity).is_none() {
            return Ok(None);
        }
    }
    Ok(Some(ExactSupportKey(key)))
}
#[cfg(test)]
mod tests;
