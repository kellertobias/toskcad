//! Hex stock closures of the literal rounded thread wall. Every hex vertex
//! lies on its flat (the six corners are emitted exactly), so the planar walls
//! and caps add no deflection; only the thread wall's own receipt remains.
use super::*;
use crate::rounded_thread::{HexForm, HexStock};
use std::f64::consts::{FRAC_PI_3, FRAC_PI_6};

pub fn tessellate_hex(stock: &HexStock, options: Tessellation) -> Result<DisplayMesh, Error> {
    tessellate_hex_with_interrupt(stock, options, None)
}
pub fn tessellate_hex_with_interrupt<'a>(
    stock: &'a HexStock,
    options: Tessellation,
    interrupt: Option<&'a dyn Fn() -> bool>,
) -> Result<DisplayMesh, Error> {
    tessellate_shape(
        stock.thread(),
        Some(Receiver::Hex(stock)),
        options,
        interrupt,
    )
}
/// The hex prism and its planar caps are meshed exactly; the signed-volume
/// error is the thread wall's radial-graph bound alone.
pub fn hex_volume_error_bound(stock: &HexStock, options: Tessellation) -> Result<f64, Error> {
    volume_error_bound(stock.thread(), options)
}
fn point_error(stock: &HexStock) -> f64 {
    let scale = ((stock.across_flats_mm() + stock.length_mm() + stock.hex_height_mm()).next_up()
        + 1.)
        .next_up();
    (64. * f64::EPSILON * scale).next_up()
}
pub(super) fn max_angle(
    stock: &HexStock,
    options: Tessellation,
    current: f64,
) -> Result<f64, Error> {
    let error = point_error(stock);
    let outer = stock.across_flats_mm() / 3_f64.sqrt();
    let gap = (stock.across_flats_mm() / 2. - stock.thread().spec().major_radius_mm).next_down();
    if error > options.deflection / 16. || gap <= 8. * error {
        return Err(MeshError::UnresolvedBoundary.into());
    }
    Ok(current.min((gap / (2. * outer)).sqrt()))
}
/// Corner `k` of the hexagon is at `30° + 60° k`; side `k` joins corners `k`
/// and `k + 1`, with outward normal at `60° + 60° k`.
fn side_of(angle: f64) -> usize {
    (((angle - FRAC_PI_6).rem_euclid(TAU)) / FRAC_PI_3).floor() as usize % 6
}
fn side_normal(side: usize) -> [f64; 3] {
    let angle = FRAC_PI_3 * (side as f64 + 1.);
    let (s, c) = angle.sin_cos();
    [c, s, 0.]
}
impl Builder<'_, '_> {
    /// Hex ring at height `z`: every column angle on its flat plus the six
    /// exact corners, in increasing angle from zero. Returns `(angle, id)`.
    fn hex_ring(&mut self, stock: &HexStock, z: f64) -> Result<Vec<(f64, u32)>, AttemptError> {
        let corners = stock.hex_vertices();
        let half = stock.across_flats_mm() / 2.;
        let mut angles: Vec<(f64, Option<usize>)> = (0..self.steps)
            .map(|column| (TAU * self.fraction(column), None))
            .collect();
        angles.extend((0..6).map(|k| (FRAC_PI_6 + FRAC_PI_3 * k as f64, Some(k))));
        angles.sort_by(|a, b| a.0.total_cmp(&b.0));
        // A column within a few rounding units of a corner is that corner.
        let merge = 1e-9;
        let mut ring: Vec<(f64, Option<usize>)> = Vec::with_capacity(angles.len());
        for (angle, corner) in angles {
            if let Some(last) = ring.last_mut()
                && angle - last.0 <= merge
            {
                if corner.is_some() {
                    *last = (angle, corner);
                }
                continue;
            }
            ring.push((angle, corner));
        }
        if ring.len() > 1 && ring[0].0 + TAU - ring.last().expect("ring").0 <= merge {
            let last = ring.pop().expect("ring");
            if last.1.is_some() {
                ring[0] = (last.0 - TAU, last.1);
            }
        }
        let mut out = Vec::with_capacity(ring.len());
        for (angle, corner) in ring {
            self.work.charge(1)?;
            let [x, y] = match corner {
                Some(k) => corners[k],
                None => {
                    let side = side_of(angle);
                    let mid = FRAC_PI_3 * (side as f64 + 1.);
                    let r = half / (angle - mid).cos();
                    let (s, c) = angle.sin_cos();
                    [r * c, r * s]
                }
            };
            let id = self.mesh.vertices.len() as u32 / 3;
            self.mesh.vertices.extend([x, y, z]);
            self.parameters.push([angle, z]);
            out.push((angle, id));
        }
        if out.len() < 7 || out[0].0 < 0. || out.windows(2).any(|w| w[0].0 >= w[1].0) {
            return Err(MeshError::UnresolvedBoundary.into());
        }
        Ok(out)
    }
    pub(super) fn add_hex(&mut self, stock: &HexStock) -> Result<(), AttemptError> {
        let mut caps = std::mem::take(&mut self.caps);
        for cap in &mut caps {
            Self::sort_ring(cap)?;
        }
        let ring = self.steps + 6;
        let needed = 4 * ring + 2 * (caps[0].len() + caps[1].len() + ring);
        if needed > self.limit.saturating_sub(self.mesh.face_indices.len()) {
            return Err(MeshError::BudgetExceeded.into());
        }
        self.work.charge(needed)?;
        let error = point_error(stock);
        let inner = (stock.thread().spec().major_radius_mm + self.receipt.point_error).next_up();
        if (stock.across_flats_mm() / 2. - error).next_down() <= inner {
            return Err(MeshError::UnresolvedBoundary.into());
        }
        let [low, high] = stock.hex_interval_mm();
        let rings = [self.hex_ring(stock, low)?, self.hex_ring(stock, high)?];
        if let HexForm::Nut = stock.form() {
            // Reversal keeps every checked thread triangle and its receipt;
            // only the material side changes, exactly as in a round bore.
            for chunk in self.mesh.indices.as_chunks_mut::<3>().0.chunks_mut(128) {
                self.work.charge(chunk.len())?;
                for triangle in chunk {
                    triangle.swap(1, 2);
                }
            }
        }
        let n = rings[0].len();
        for i in 0..n {
            self.work.checkpoint()?;
            let j = (i + 1) % n;
            let next = if j == 0 { TAU } else { rings[0][j].0 };
            let side = side_of(0.5 * (rings[0][i].0 + next));
            let [a, b] = [&rings[0], &rings[1]];
            for triangle in [[a[i].1, a[j].1, b[j].1], [a[i].1, b[j].1, b[i].1]] {
                self.bore_triangle(triangle, HexStock::wall_face(side), side_normal(side), true)?;
            }
        }
        let [start, end] = caps;
        match stock.form() {
            HexForm::Nut => {
                self.ring_cap(0, false, &start, &rings[0])?;
                self.ring_cap(1, true, &end, &rings[1])?;
            }
            HexForm::Head { .. } => {
                // The seat faces the shaft; the hex top closes the head.
                self.ring_cap(0, true, &start, &rings[1])?;
                self.add_cap(1, end)?;
                self.hex_top(&rings[0], low)?;
            }
        }
        Ok(())
    }
    fn hex_top(&mut self, ring: &[(f64, u32)], z: f64) -> Result<(), AttemptError> {
        let center = self.mesh.vertices.len() as u32 / 3;
        self.work.charge(1)?;
        self.mesh.vertices.extend([0., 0., z]);
        self.parameters.push([0., z]);
        for (&(_, a), &(_, b)) in ring.iter().zip(ring.iter().cycle().skip(1)) {
            self.bore_triangle([center, b, a], 13, [0., 0., -1.], false)?;
        }
        Ok(())
    }
}
impl Builder<'_, '_> {
    /// Angular sweep between an inner thread ring and an outer ring with
    /// explicit angles, both starting at angle zero (the outer ring's first
    /// vertex may be its first corner). `up` selects the +Z normal side. The
    /// edge-pairing proof checks a closed annulus with both rings as boundary.
    fn ring_cap(
        &mut self,
        face: u32,
        up: bool,
        inner: &[(f64, u32)],
        outer: &[(f64, u32)],
    ) -> Result<(), AttemptError> {
        if inner[0].0 != 0. || inner.last().expect("checked ring").0 >= TAU || outer[0].0 < 0. {
            return Err(MeshError::UnresolvedBoundary.into());
        }
        // Rotate the outer ring to begin at angle zero when its first vertex
        // is a later corner: the sweep below needs a common starting ray.
        let start = outer.iter().position(|v| v.0 == 0.);
        let outer: Vec<(f64, u32)> = match start {
            Some(s) => outer[s..].iter().chain(&outer[..s]).copied().collect(),
            None => return Err(MeshError::UnresolvedBoundary.into()),
        };
        let mut a = 0;
        let mut b = 0;
        let mut proof = bore::CapProof::default();
        while a < inner.len() || b < outer.len() {
            self.work.checkpoint()?;
            let ia = inner[a % inner.len()].1;
            let ob = outer[b % outer.len()].1;
            let next_inner = if a + 1 < inner.len() {
                inner[a + 1].0
            } else {
                TAU
            };
            let next_outer = if b + 1 < outer.len() {
                outer[b + 1].0
            } else {
                TAU
            };
            let mut triangles = [[0; 3]; 2];
            let mut count = 0;
            if next_inner <= next_outer {
                triangles[count] = [ia, ob, inner[(a + 1) % inner.len()].1];
                count += 1;
                a += 1;
            }
            if next_outer <= next_inner {
                let ia = inner[a % inner.len()].1;
                triangles[count] = [ia, ob, outer[(b + 1) % outer.len()].1];
                count += 1;
                b += 1;
            }
            for mut triangle in triangles[..count].iter().copied() {
                if !up {
                    triangle.swap(1, 2);
                }
                self.bore_triangle(triangle, face, [0., 0., if up { 1. } else { -1. }], false)?;
                proof.add(triangle, self.work)?;
            }
        }
        bore::closed_sweep([a, b], [inner.len(), outer.len()])?;
        let outer_edges = outer
            .iter()
            .zip(outer.iter().cycle().skip(1))
            .map(|(&(_, a), &(_, b))| [a, b]);
        let inner_edges = inner
            .iter()
            .zip(inner.iter().cycle().skip(1))
            .map(|(&(_, a), &(_, b))| [b, a]);
        proof.finish(
            outer_edges.chain(inner_edges).map(|mut edge| {
                if !up {
                    edge.swap(0, 1);
                }
                edge
            }),
            self.work,
        )?;
        Ok(())
    }
}
