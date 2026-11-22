//! Literal blind closure of the original rounded wall and its finite tip rim.
use super::*;
use crate::rounded_thread::BlindBore;
use crate::tapered_thread::interval::I;

/// Produce a literal approximation; this is not an exact boundary or STEP
/// admission. Native consumers must still earn fresh global embedding.
pub fn tessellate_blind_bore_with_interrupt<'a>(
    bore: &'a BlindBore,
    options: Tessellation,
    interrupt: Option<&'a dyn Fn() -> bool>,
) -> Result<DisplayMesh, Error> {
    options.check()?;
    // The original wall/floor approximation receives half the requested C0
    // allowance. Actual rational stock coordinates are transferred separately.
    let radius = bore.stock_radius_mm().as_f64();
    let height = bore.stock_height_mm().as_f64();
    let [rlo, rhi] = bore.stock_radius_bounds_mm();
    let [hlo, hhi] = bore.stock_height_bounds_mm();
    let transfer = (I::point(2.) * I::point((I::new(rlo, rhi) - I::point(radius)).magnitude())
        + I::point((I::new(hlo, hhi) - I::point(height)).magnitude()))
    .hi;
    if !transfer.is_finite() || transfer > options.deflection / 4. {
        return Err(MeshError::UnresolvedBoundary.into());
    }
    let local = Tessellation {
        deflection: options.deflection / 2.,
        ..options
    };
    tessellate_shape(bore.cutter(), Some(Receiver::Blind(bore)), local, interrupt)
}

pub(in crate::rounded_thread_mesh) fn max_angle(
    bore: &BlindBore,
    options: Tessellation,
    current: f64,
) -> Result<f64, Error> {
    let radius = bore.stock_radius_mm().as_f64();
    let height = bore.stock_height_mm().as_f64();
    let error = point_error(radius, height);
    let gap = (radius - bore.cutter_radial_bounds_mm()[1]).next_down();
    if error > options.deflection / 16. || gap <= 8. * error {
        return Err(MeshError::UnresolvedBoundary.into());
    }
    Ok(current
        .min((options.deflection / (2. * radius)).sqrt())
        .min((gap / (2. * radius)).sqrt()))
}

impl Builder<'_, '_> {
    pub(in crate::rounded_thread_mesh) fn add_blind_bore(
        &mut self,
        bore: &BlindBore,
    ) -> Result<(), AttemptError> {
        let radius = bore.stock_radius_mm().as_f64();
        let height = bore.stock_height_mm().as_f64();
        let mut caps = std::mem::take(&mut self.caps);
        for ring in &mut caps {
            Self::sort_ring(ring)?;
        }
        let needed = caps[0].len() + caps[1].len() + 4 * self.steps;
        if needed > self.limit.saturating_sub(self.mesh.face_indices.len()) {
            return Err(MeshError::BudgetExceeded.into());
        }
        self.work.charge(needed)?;
        // The actual tip rim bounds the floor, not a circle of nominal radius.
        // Existing cap fan proof uses exactly the original source rim points.
        let tip = std::mem::take(&mut caps[1]);
        self.add_cap(1, tip)?;
        for triangles in self.mesh.indices.as_chunks_mut::<3>().0.chunks_mut(128) {
            self.work.charge(triangles.len())?;
            for triangle in triangles {
                triangle.swap(1, 2);
            }
        }
        // Every inward wall/floor facet retains the original finite source
        // certificate. Only its occupied-side orientation has reversed.
        prove_outer_polygon(
            radius,
            height,
            bore.cutter_radial_bounds_mm()[1],
            self.steps,
        )?;
        let mut outer: [Vec<u32>; 2] = std::array::from_fn(|_| Vec::with_capacity(self.steps));
        for (cap, ring) in outer.iter_mut().enumerate() {
            for column in 0..self.steps {
                self.work.charge(1)?;
                let (s, c) = (TAU * self.fraction(column)).portable_sin_cos();
                let z = if cap == 0 { 0. } else { height };
                let id = self.mesh.vertices.len() as u32 / 3;
                self.mesh.vertices.extend([radius * c, radius * s, z]);
                self.parameters.push([TAU * self.fraction(column), z]);
                ring.push(id);
            }
        }
        self.bore_cap(0, &caps[0], &outer[0])?;
        for column in 0..self.steps {
            self.work.checkpoint()?;
            let next = (column + 1) % self.steps;
            let (s, c) = (TAU * self.fraction(column)).portable_sin_cos();
            for triangle in [
                [outer[0][column], outer[0][next], outer[1][next]],
                [outer[0][column], outer[1][next], outer[1][column]],
            ] {
                self.bore_triangle(triangle, 7, [c, s, 0.], true)?;
            }
        }
        self.work.charge(1)?;
        let center = self.mesh.vertices.len() as u32 / 3;
        self.mesh.vertices.extend([0., 0., height]);
        self.parameters.push([0., height]);
        for column in 0..self.steps {
            self.work.checkpoint()?;
            self.bore_triangle(
                [
                    center,
                    outer[1][column],
                    outer[1][(column + 1) % self.steps],
                ],
                8,
                [0., 0., 1.],
                false,
            )?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests;
