//! Opaque preparation of ONE actual literal curved face. This is deliberately
//! not a partial LiteralGeometryPlan and cannot be supplied to whole meshes.
use super::*;
pub struct LiteralFaceGeometryPlan {
    source: OwnedMaterial,
    face: usize,
    tensor: usize,
    normals: Vec<normal::Normal>,
    preparation: usize,
    retained: usize,
    peak: usize,
}
impl LiteralFaceGeometryPlan {
    pub fn source(&self) -> &OwnedMaterial {
        &self.source
    }
    pub fn face_index(&self) -> usize {
        self.face
    }
    pub fn preparation_work_used(&self) -> usize {
        self.preparation
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.retained
    }
    pub fn construction_peak_bytes(&self) -> usize {
        self.peak
    }
    pub(in crate::step_read::npt_stock::selected::owned) fn original_normal_coefficients(
        &self,
        index: [usize; 2],
    ) -> Result<&[[I; 3]]> {
        let tensor = &self.source.tensors()[self.tensor];
        if index[0] >= tensor.spans[0] || index[1] >= tensor.spans[1] {
            return Err(Error::Invalid);
        }
        Ok(self.normals[index[0] * tensor.spans[1] + index[1]].coefficients())
    }
}
impl OwnedMaterial {
    pub fn face_geometry_plan(
        &self,
        face: usize,
        options: LiteralGeometryOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<LiteralFaceGeometryPlan> {
        if options.maximum_work == 0
            || options.maximum_work > 512_000_000
            || options.maximum_bytes == 0
            || options.maximum_bytes > BYTES
        {
            return Err(Error::Invalid);
        }
        let owner = self.faces().get(face).ok_or(Error::Invalid)?;
        let Support::Tensor(tensor) = owner.support() else {
            return Err(Error::Unsupported);
        };
        let bytes = self
            .retained_budget_bytes()
            .checked_add(SCRATCH)
            .and_then(|n| n.checked_add(size_of::<LiteralFaceGeometryPlan>()))
            .ok_or(Error::Limit)?;
        if bytes > options.maximum_bytes {
            return Err(Error::Limit);
        }
        let mut b = Budget {
            used: 0,
            remaining: options.maximum_work,
            bytes,
            peak: bytes,
            limit: options.maximum_bytes,
            interrupt,
        };
        b.step(256)?;
        let t = &self.tensors()[tensor];
        let count = t.spans[0].checked_mul(t.spans[1]).ok_or(Error::Limit)?;
        let mut normals = b.allocate(count)?;
        for u in 0..t.spans[0] {
            for v in 0..t.spans[1] {
                normals.push(normal::Normal::new(self, t, u, v, &mut b)?);
            }
        }
        b.step(0)?;
        Ok(LiteralFaceGeometryPlan {
            source: self.clone(),
            face,
            tensor,
            normals,
            preparation: b.used,
            retained: b.bytes - SCRATCH,
            peak: b.peak,
        })
    }
}
