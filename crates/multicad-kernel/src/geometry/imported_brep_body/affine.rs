//! Scale of curved Rust solids through their exact retained boundary.
//!
//! The source is read back from its own exact STEP encoding, or taken from
//! an already retained boundary. The boundary's affine image keeps every
//! face, loop and sense: circles become principal-axis ellipses, cylinders
//! linear extrusions of them, and other quadrics keep their family only
//! under a uniform scale. The image is revalidated and integrated afresh, so
//! its volume is an independent check of the determinant.
use super::*;
use crate::AffineTransform;
use crate::geometry::authored_sweep::{AuthoredNaming, FaceName};

impl Solid {
    /// The source boundary, then its affine image, as a retained solid.
    /// `None` leaves solids without a Rust STEP encoding to their caller.
    pub(crate) fn rust_scale_brep(
        &self,
        transform: AffineTransform,
        feature: EntityId,
        policy: GeometryPolicy,
    ) -> Option<Result<Self, GeometryError>> {
        let source = match self.imported_brep_body() {
            Some(body) => body.solid.clone(),
            None => {
                if !self.supports_step_export() {
                    return None;
                }
                match self.exact_brep(policy) {
                    Ok(Some(solid)) => solid,
                    Ok(None) => return None,
                    Err(cause) => return Some(Err(cause)),
                }
            }
        };
        Some(self.scaled_brep(&source, transform, feature, policy))
    }

    fn exact_brep(&self, policy: GeometryPolicy) -> Result<Option<brep::Solid>, GeometryError> {
        let text = self.rust_step_text("scale")?;
        let verified = super::super::step_import::pcurve_validation::validate(&text, policy)?;
        let mut assembly =
            step_read::brep::read_assembly_validated(&verified, None).map_err(failure)?;
        if assembly.nodes.len() != 1
            || !assembly.nodes[0].children.is_empty()
            || assembly.nodes[0].solids.len() != 1
        {
            return Ok(None);
        }
        Ok(assembly.nodes.remove(0).solids.pop())
    }

    fn scaled_brep(
        &self,
        source: &brep::Solid,
        transform: AffineTransform,
        feature: EntityId,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        let map = brep::Affine {
            matrix: transform.matrix,
            translation: transform.translation_mm,
        };
        let image = source.boundary().affine_image(map).map_err(failure)?;
        let scaled = brep::Solid::new(image, source.tolerance(), None).map_err(failure)?;
        let determinant = transform.determinant();
        let expected = source.volume_mm3() * determinant;
        let allowed = scaled.volume_error_mm3() + source.volume_error_mm3() * determinant;
        if !(scaled.volume_mm3() - expected)
            .abs()
            .le(&allowed.max(1e-9 * expected))
        {
            return Err(error(
                OPERATION,
                GeometryErrorKind::InvalidResult,
                "Rust Scale did not conserve the affine material volume",
            ));
        }
        match self.scaled_names(source, transform)? {
            Some(naming) => Self::from_retained_brep(scaled, feature, Rc::new(naming), policy),
            None => {
                let roles = Rc::new(roles(&scaled));
                Self::from_imported_brep(scaled, feature, roles, policy)
            }
        }
    }

    /// Source face names carried one-to-one onto the image. Each retained
    /// face must match exactly one published source face by area and
    /// centroid; otherwise the image is named from its own geometry.
    fn scaled_names(
        &self,
        source: &brep::Solid,
        transform: AffineTransform,
    ) -> Result<Option<AuthoredNaming>, GeometryError> {
        if self.imported_brep_body().is_some() {
            return Ok(None);
        }
        let Some(topology) = self.rust_topology().transpose()? else {
            return Ok(None);
        };
        let published = self.faces()?;
        let [low, high] = source.bounds_mm();
        let size = (0..3).map(|i| high[i] - low[i]).fold(1_f64, f64::max);
        let mut used = vec![false; published.len()];
        let mut faces = Vec::with_capacity(source.faces().len());
        for face in source.faces() {
            let mut matches = published.iter().enumerate().filter(|(_, candidate)| {
                (candidate.area_mm2 - face.area_mm2).abs() <= 1e-6 * face.area_mm2.max(1.)
                    && (0..3).all(|i| {
                        (candidate.center_mm[i] - face.centroid_mm[i]).abs() <= 1e-6 * size
                    })
            });
            let (Some((index, candidate)), None) = (matches.next(), matches.next()) else {
                return Ok(None);
            };
            if std::mem::replace(&mut used[index], true) {
                return Ok(None);
            }
            let Some(element) = topology.elements.iter().find(|element| {
                element.key.kind == TopologyKind::Face && element.key.ordinal == candidate.index
            }) else {
                return Ok(None);
            };
            // The anchor is published under the first origin of the face.
            let anchor = element.frames.first();
            let mut origins: Vec<_> = element.origins.iter().cloned().collect();
            if let Some(anchor) = anchor
                && let Some(position) = origins.iter().position(|o| *o == anchor.origin)
            {
                origins.swap(0, position);
            }
            let anchor = anchor
                .filter(|anchor| origins.first() == Some(&anchor.origin))
                .map(|anchor| transform.frame(anchor.frame))
                .transpose()?;
            faces.push(FaceName {
                origins,
                support: None,
                anchor,
            });
        }
        Ok(Some(AuthoredNaming { faces }))
    }
}
