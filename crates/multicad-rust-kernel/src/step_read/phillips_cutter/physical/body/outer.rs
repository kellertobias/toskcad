//! Measured matching of actual retained outer carriers. This does not certify
//! raw exact support sewing, cavity containment, closed embedding or net mass.
use super::*;
mod caps;
mod integers;
use crate::{step::helical::prepare_retained_reference, tapered_thread::interval::I};
#[derive(Clone, Copy, Debug)]
pub struct OuterCarrierOptions {
    pub maximum_distance_mm: f64,
    pub work_limit: usize,
}
#[derive(Clone, Copy, Debug)]
pub struct TensorCarrierBound {
    pub face_entity: u32,
    pub surface_entity: u32,
    pub control_rows: usize,
    pub control_columns: usize,
    /// Uniform tensor carrier distance from the reconstructed reference carrier
    /// under per-cell Bernstein reparameterization, including numeric conversion.
    pub control_deviation_mm: f64,
    /// Measured C0 support residual of the literal end-rim spline controls.
    pub tip_plane_distance_mm: f64,
    /// The real entry crest Circle and retained polynomial rim are independently
    /// retained; this measured C0 bound does not assert exact common incidence.
    pub entry_circle_distance_mm: Option<f64>,
}
#[derive(Debug)]
pub struct OuterCarrierBounds<'a, 'input> {
    body: &'a Body<'input>,
    tensors: Vec<TensorCarrierBound>,
    work_used: usize,
}
impl<'a, 'input> OuterCarrierBounds<'a, 'input> {
    pub fn body(&self) -> &'a Body<'input> {
        self.body
    }
    pub fn tensors(&self) -> &[TensorCarrierBound] {
        &self.tensors
    }
    pub fn meridian_reconstruction_error_mm(&self) -> f64 {
        self.body.profile.rounding_error_mm()
    }
    pub fn work_used(&self) -> usize {
        self.work_used
    }
    pub fn retained_capacity_bytes(&self) -> usize {
        std::mem::size_of::<Self>()
            + self.tensors.capacity() * std::mem::size_of::<TensorCarrierBound>()
    }
}
fn norm(v: [I; 3]) -> Result<f64> {
    let squared = v.into_iter().fold(I::point(0.), |s, v| s + v.square());
    require(squared.valid())?;
    let result = squared.sqrt();
    require(result.valid() && result.hi.is_finite())?;
    Ok(result.hi)
}
fn deviation(
    actual: Point,
    radial: [f64; 2],
    axial: [f64; 2],
    z: f64,
    [o, u, v]: [Point; 3],
) -> Result<f64> {
    let r = radial.map(I::point);
    let a = axial.map(I::point);
    let xy = [r[0] * a[0] - r[1] * a[1], r[0] * a[1] + r[1] * a[0]];
    let ui = u.map(I::point);
    let vi = v.map(I::point);
    let n = [
        ui[1] * vi[2] - ui[2] * vi[1],
        ui[2] * vi[0] - ui[0] * vi[2],
        ui[0] * vi[1] - ui[1] * vi[0],
    ];
    norm(std::array::from_fn(|i| {
        let literal = I::new(actual[i].next_down(), actual[i].next_up());
        literal - (I::point(o[i]) + ui[i] * xy[0] + vi[i] * xy[1] + n[i] * I::point(z))
    }))
}
impl<'input> Body<'input> {
    pub fn outer_carrier_bounds<'a>(
        &'a self,
        options: OuterCarrierOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> std::result::Result<OuterCarrierBounds<'a, 'input>, PhillipsReadError> {
        let mut work = Work {
            remaining: options.work_limit,
            interrupt,
            cancelled: Cell::new(false),
        };
        let result = (|| {
            require(options.maximum_distance_mm.is_finite() && options.maximum_distance_mm > 0.)?;
            require(self.profile.rounding_error_mm() <= options.maximum_distance_mm)?;
            work.step(1)?;
            let mut shell = Shell::new(&self.entities, self.tolerance, true);
            // Bound work from the literal arena before constructing comparison
            // nets; the file controls, knots and trims remain authoritative.
            for &id in self.tensor_face_entities() {
                let (surface, _, _) = helical::face(&shell, id)?;
                let args = shell.args(surface, "B_SPLINE_SURFACE_WITH_KNOTS")?;
                let Value::List(rows) = args.get(3).ok_or(ReadError::Invalid)? else {
                    return Err(ReadError::Invalid);
                };
                for row in rows {
                    let Value::List(values) = row else {
                        return Err(ReadError::Invalid);
                    };
                    work.step(values.len())?;
                }
            }
            integers::validate(self, &shell, &mut work)?;
            let basis = prepare_retained_reference(&self.thread, self.visible[0], self.tolerance)
                .map_err(helical::step_error)?;
            require(basis.bands.len() == self.tensor_face_count)?;
            let caps = caps::Caps::new(self, &shell)?;
            let mut tensors = Vec::with_capacity(self.tensor_face_count);
            let proof = helical::Proof {
                error: Some(options.maximum_distance_mm),
                extent: self
                    .thread
                    .spec()
                    .length_mm
                    .max(self.thread.spec().major_radius_mm),
            };
            let prefix = helical::verify_prefix(
                &mut shell,
                self.tensor_face_entities(),
                &basis,
                self.frame,
                helical::PrefixContract {
                    proof,
                    forward: self.thread.spec().left_hand,
                    legacy_shank_entry: false,
                    all_entry_crests: true,
                    exit_crests: false,
                },
                |index, controls| {
                    let mut maximum = 0_f64;
                    let radial = &basis.bands[index];
                    for (row, &r) in controls.iter().zip(&radial.poles) {
                        for ((&p, &a), &z) in row.iter().zip(&basis.axial.poles).zip(&basis.heights)
                        {
                            work.step(1)?;
                            maximum = maximum.max(deviation(p, r, a, z, self.frame)?);
                        }
                    }
                    let (tip_plane_distance_mm, entry_circle_distance_mm) =
                        caps.measure(index, controls, basis.visible_start, &mut work)?;
                    require(
                        maximum <= options.maximum_distance_mm
                            && tip_plane_distance_mm <= options.maximum_distance_mm
                            && entry_circle_distance_mm
                                .is_none_or(|e| e <= options.maximum_distance_mm),
                    )?;
                    let id = self.tensor_face_entities()[index];
                    let (surface, _, _) =
                        helical::face(&Shell::new(&self.entities, self.tolerance, true), id)?;
                    tensors.push(TensorCarrierBound {
                        face_entity: id,
                        surface_entity: surface,
                        control_rows: controls.len(),
                        control_columns: basis.axial.poles.len(),
                        control_deviation_mm: maximum,
                        tip_plane_distance_mm,
                        entry_circle_distance_mm,
                    });
                    Ok(())
                },
            )?;
            work.step(1)?;
            // The real tip cap uses the exact retained rim/vertex identities.
            helical::verify_plain_tip(
                proof,
                &shell,
                self.material_faces[self.tensor_face_count + 1],
                self.frame,
                (self.thread.spec().length_mm, self.thread.spec().left_hand),
                &prefix.rings[1],
            )?;
            helical::profiled::match_phillips_outer(
                &mut shell,
                &self.material_faces,
                &self.thread,
                (self.frame, self.visible[0]),
                &self.profile,
                &prefix,
            )?;
            Ok(OuterCarrierBounds {
                body: self,
                tensors,
                work_used: options.work_limit - work.remaining,
            })
        })();
        if work.cancelled.get() {
            Err(PhillipsReadError::Cancelled)
        } else {
            result.map_err(PhillipsReadError::from)
        }
    }
}
#[cfg(test)]
mod tests;
