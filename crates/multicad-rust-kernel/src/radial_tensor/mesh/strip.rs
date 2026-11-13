//! Bounded source-only U restrictions within one Tensor and check pass.
use super::*;
#[cfg(test)]
pub(super) mod tests;
pub(super) struct Cache<'a> {
    owner: &'a Tensor,
    entries: BTreeMap<[u64; 6], certificate::SourceStrip>,
    // Separate source-position and normal tables each obey the same fixed
    // entry limit; neither carries a facet or a successful acceptance result.
    normals: BTreeMap<[u64; 6], certificate::SourceStrip>,
    limit: usize,
}
impl<'a> Cache<'a> {
    pub fn new(owner: &'a Tensor, limit: usize) -> Self {
        Self {
            owner,
            entries: BTreeMap::new(),
            normals: BTreeMap::new(),
            limit,
        }
    }
    pub fn advance(&mut self, lower_v: f64, work: &mut Work<'_>) -> Result<(), Error> {
        work.charge(self.entries.len() + self.normals.len())?;
        // Closed patch endpoints are retained until strictly behind the row.
        self.entries
            .retain(|key, _| f64::from_bits(key[3]) >= lower_v);
        self.normals
            .retain(|key, _| f64::from_bits(key[3]) >= lower_v);
        Ok(())
    }
    /// Retain every coefficient interval of the original normal polynomial.
    /// The U-then-V restriction arithmetic matches the uncached certificate;
    /// each actual facet must still check every coefficient independently.
    pub fn normal_coefficients(
        &mut self,
        patch: &tensor::Patch,
        range: [[f64; 2]; 2],
        work: &mut Work<'_>,
    ) -> Result<Vec<Vec<[I; 3]>>, Error> {
        work.step()?;
        if !self.owner.owns(patch) {
            return Err(Error::Invalid);
        }
        let mut key = [0; 6];
        for (slot, x) in key
            .iter_mut()
            .zip(patch.domain.into_iter().flatten().chain(range[0]))
        {
            *slot = x.to_bits();
        }
        if let Some(source) = self.normals.get(&key) {
            return source.normal_coefficients(range[1], work);
        }
        let source = certificate::SourceStrip::new(&patch.normal, patch.domain, range[0], work)?;
        let coefficients = source.normal_coefficients(range[1], work)?;
        if self.normals.len() < self.limit {
            work.charge(16)?;
            self.normals.insert(key, source);
        }
        Ok(coefficients)
    }
    pub fn cell(
        &mut self,
        patch: &tensor::Patch,
        range: [[f64; 2]; 2],
        work: &mut Work<'_>,
    ) -> Result<certificate::SourceCell, Error> {
        work.step()?;
        if !self.owner.owns(patch) {
            return Err(Error::Invalid);
        }
        let mut key = [0; 6];
        for (slot, x) in key
            .iter_mut()
            .zip(patch.domain.into_iter().flatten().chain(range[0]))
        {
            *slot = x.to_bits();
        }
        if let Some(source) = self.entries.get(&key) {
            return source.cell(range[1], work);
        }
        if self.entries.len() >= self.limit {
            return certificate::SourceCell::new(&patch.controls, patch.domain, range, work);
        }
        let source = certificate::SourceStrip::new(&patch.controls, patch.domain, range[0], work)?;
        let cell = source.cell(range[1], work)?;
        work.charge(16)?;
        self.entries.insert(key, source);
        Ok(cell)
    }
}
