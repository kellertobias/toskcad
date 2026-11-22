//! Indexed tessellation of the helical radial envelope in screw coordinates.
//! A cell stays in one analytic flank, crest or root band. Cap intersections
//! and the periodic seam share vertex identities with the side triangles.

use crate::{
    display::{DisplayMesh, MeshError, Tessellation},
    helical_thread::Thread,
    hex_nut::Nut,
};
use multicad_precision::PortableMath as _;
use std::{
    cell::Cell,
    collections::HashMap,
    f64::consts::{PI, TAU},
};
/// Counted logical geometry operations; allocator/CDT workspaces keep their
/// separately bounded contracts. Every proposal/refinement loop pays before
/// constructing its cell, including cells discarded by finite clipping.
struct Counted<'a> {
    remaining: Cell<usize>,
    cancelled: Cell<bool>,
    interrupt: &'a dyn Fn() -> bool,
}
impl Counted<'_> {
    fn step(&self, n: usize) -> Result<(), MeshError> {
        if (self.interrupt)() {
            self.cancelled.set(true);
            return Err(MeshError::BudgetExceeded);
        }
        self.remaining.set(
            self.remaining
                .get()
                .checked_sub(n)
                .ok_or(MeshError::BudgetExceeded)?,
        );
        Ok(())
    }
}

mod blind_counterbore;
pub use blind_counterbore::counterbore_volume_error_bound;
pub(crate) use blind_counterbore::open_counterbored_cavity;
mod receiver;
pub use receiver::{receiver_volume_error_bound, tessellate_receiver};
mod external;
pub use external::{shaft_volume_error_bound, tessellate_shaft};
mod thread_end_relief;
pub use thread_end_relief::{end_relief_volume_error_bound, tessellate_end_relief};
mod screw;
pub use screw::tessellate_profiled_cosmetic_screw;
pub use screw::{
    PhillipsCavityVertex, PhillipsHeadVertex, PhillipsThreadVertex, PhillipsThreadVertexKind,
};
pub use screw::{
    PhillipsMeshError, PhillipsScrewMesh, tessellate_cosmetic_phillips_screw,
    tessellate_phillips_screw, tessellate_phillips_screw_with_work,
};
pub use screw::{profiled_screw_volume_error_bound, tessellate_profiled_screw};
pub use screw::{screw_volume_error_bound, tessellate_screw};

/// The positional bound is for linear interpolation in each screw-coordinate
/// triangle. Angular deflection bounds the variation of the exact normal
/// within each smooth triangle, without smoothing across profile creases.
pub fn tessellate(thread: &Thread, options: Tessellation) -> Result<DisplayMesh, MeshError> {
    tessellate_envelope(thread, None, options)
}

/// A cylindrical solid with a coaxial threaded through bore. The inner
/// wall uses the thread face ordinals, followed by the outer cylinder.
pub fn tessellate_bore(
    thread: &Thread,
    outer_radius_mm: f64,
    options: Tessellation,
) -> Result<DisplayMesh, MeshError> {
    if !outer_radius_mm.is_finite()
        || outer_radius_mm > 1e7
        || outer_radius_mm - thread.spec().major_radius_mm <= 64. * f64::EPSILON * outer_radius_mm
    {
        return Err(MeshError::InvalidOptions);
    }
    tessellate_envelope(
        thread,
        Some(BoreExterior::Cylinder(outer_radius_mm)),
        options,
    )
}

/// A hexagonal solid with a coaxial threaded through bore. Six exact planar
/// wrench faces follow the cap and helical wall ordinals in corner order.
pub fn tessellate_nut(nut: &Nut, options: Tessellation) -> Result<DisplayMesh, MeshError> {
    tessellate_envelope(nut.thread(), Some(BoreExterior::Hex(nut)), options)
}

#[derive(Clone, Copy)]
enum BoreExterior<'a> {
    Cylinder(f64),
    Hex(&'a Nut),
}

fn tessellate_envelope(
    thread: &Thread,
    exterior: Option<BoreExterior<'_>>,
    options: Tessellation,
) -> Result<DisplayMesh, MeshError> {
    envelope_builder(thread, exterior, options)?.build(exterior)
}

pub(crate) fn open_blind_cavity(
    thread: &Thread,
    options: Tessellation,
) -> Result<(DisplayMesh, Vec<u32>), MeshError> {
    let mut builder = envelope_builder(thread, None, options)?;
    builder.build_wall()?;
    let [mut mouth, floor] = std::mem::take(&mut builder.caps);
    Builder::sort_ring(&mut mouth)?;
    builder.add_cap(1, floor)?;
    for triangle in builder.mesh.indices.as_chunks_mut::<3>().0 {
        triangle.swap(1, 2);
    }
    Ok((
        builder.mesh,
        mouth.into_iter().map(|(_, vertex)| vertex).collect(),
    ))
}

/// Outward external wall with both cap rims as counterclockwise vertex rings,
/// for callers that close the planar ends themselves.
pub(crate) fn open_wall(
    thread: &Thread,
    options: Tessellation,
) -> Result<(DisplayMesh, [Vec<u32>; 2]), MeshError> {
    let mut builder = envelope_builder(thread, None, options)?;
    builder.build_wall()?;
    let mut rings = std::mem::take(&mut builder.caps);
    for ring in &mut rings {
        Builder::sort_ring(ring)?;
    }
    Ok((
        builder.mesh,
        rings.map(|ring| ring.into_iter().map(|(_, vertex)| vertex).collect()),
    ))
}

fn envelope_builder<'a>(
    thread: &'a Thread,
    exterior: Option<BoreExterior<'_>>,
    options: Tessellation,
) -> Result<Builder<'a>, MeshError> {
    check_options(options)?;
    let major = thread.spec().major_radius_mm;
    let (radius, deflection, stride) = match exterior {
        // Each inscribed cylindrical chord must clear the bore crest,
        // including cap rail intersections between angular columns.
        Some(BoreExterior::Cylinder(radius)) => {
            (radius, options.deflection.min(2. * (radius - major)), 1.)
        }
        // A multiple of twelve angular steps includes every hex corner.
        // The planar stock has no curvature to refine.
        Some(BoreExterior::Hex(_)) => (major, options.deflection, 12.),
        None => (major, options.deflection, 1.),
    };
    let max_angle = (options.angle * 0.5).min((2. * deflection / radius).sqrt());
    let angular_steps = (TAU / (stride * max_angle)).ceil() * stride;
    if angular_steps * 2. > options.triangle_limit as f64 {
        return Err(MeshError::BudgetExceeded);
    }
    let angular_steps = angular_steps as usize;
    let strips = phase_strips(thread, TAU / angular_steps as f64, options)?;
    Ok(Builder::new(thread, angular_steps, strips, options))
}

fn check_options(options: Tessellation) -> Result<(), MeshError> {
    if !options.deflection.is_finite()
        || options.deflection <= 0.
        || !(0.01..=1.).contains(&options.angle)
        || !(1..=5_000_000).contains(&options.triangle_limit)
    {
        return Err(MeshError::InvalidOptions);
    }
    Ok(())
}

fn hex_side(angle: f64) -> u32 {
    ((angle - PI / 6.).rem_euclid(TAU) / (TAU / 6.)).floor() as u32
}

/// Bound volume deviation for the generated mesh. At equal depth, each
/// angular sector remains a positive, single-valued radial graph. A chord
/// perturbation of `deflection` changes radius directly and changes the exact
/// radial function through its angular phase. Integrate that radial error
/// over the full angle and length.
pub fn volume_error_bound(thread: &Thread, options: Tessellation) -> f64 {
    let spec = thread.spec();
    let min_radius = thread.root_radius_mm() * (options.angle * 0.5).portable_cos();
    let angular_error = 2.
        * (options.deflection / (2. * min_radius))
            .min(1.)
            .portable_asin();
    let max_slope = thread.maximum_profile_slope();
    let radial_error = options.deflection + thread.helical_rate().abs() * max_slope * angular_error;
    TAU * spec.major_radius_mm * spec.length_mm * radial_error
}

/// The hole keeps the external thread's error bound. The inscribed outer
/// polygon's radial deficit is at most the specified chord deflection.
pub fn bore_volume_error_bound(
    thread: &Thread,
    outer_radius_mm: f64,
    options: Tessellation,
) -> f64 {
    volume_error_bound(thread, options)
        + TAU * outer_radius_mm * thread.spec().length_mm * options.deflection
}

#[derive(Clone, Copy)]
struct Strip {
    phase: f64,
    band: usize,
}

fn phase_strips(
    thread: &Thread,
    theta: f64,
    options: Tessellation,
) -> Result<Vec<Strip>, MeshError> {
    phase_strips_counted(thread, theta, options, None)
}
fn phase_strips_counted(
    thread: &Thread,
    theta: f64,
    options: Tessellation,
    work: Option<&Counted<'_>>,
) -> Result<Vec<Strip>, MeshError> {
    let spec = thread.spec();
    let root = thread.root_radius_mm();
    let rate = thread.helical_rate().abs();
    let normal_scale = (rate / root).portable_hypot(1.);
    let mut strips = Vec::new();
    for (band, bounds) in thread.profile_breaks().windows(2).enumerate() {
        let mut pending = vec![[bounds[0], bounds[1]]];
        while let Some([lo, hi]) = pending.pop() {
            if let Some(w) = work {
                w.step(256)?;
            }
            let [a, b] = [lo, hi].map(|s| thread.slope_at_phase(s, band));
            let width = hi - lo;
            let slope = a.abs().max(b.abs());
            let curvature = if band == 2 {
                spec.root_radius_mm
                    .map_or(0., |radius| slope.portable_hypot(1.).powi(3) / radius)
            } else {
                0.
            };
            // For barycentric interpolation, each parameter's variance is
            // at most its span squared / 4. Bound the Taylor remainder
            // using |Xθθ|≤r, |Xθs|≤|r'| and |Xss|≤|r''|.
            let chord = (spec.major_radius_mm * theta * theta
                + 2. * slope * theta * width
                + curvature * width * width)
                * 0.125;
            // At fixed θ the unit normal is parallel to
            // (1, (lead/TAU) r'/r, -r') in the radial basis. The slope
            // is monotone within each band and |Δr|≤max|r'|Δs.
            let normal = normal_scale * (b.portable_atan() - a.portable_atan()).abs()
                + rate * slope * width / (root * root);
            if chord <= options.deflection && theta + normal <= options.angle {
                if strips.len() == options.triangle_limit {
                    return Err(MeshError::BudgetExceeded);
                }
                strips.push(Strip { phase: lo, band });
            } else {
                let midpoint = lo + width * 0.5;
                if !(lo < midpoint && midpoint < hi) {
                    return Err(MeshError::UnresolvedBoundary);
                }
                pending.extend([[midpoint, hi], [lo, midpoint]]);
            }
        }
    }
    Ok(strips)
}

#[derive(Clone, Copy, Eq, Hash, PartialEq)]
enum VertexKey {
    Grid(usize, i64),
    CapColumn(usize, usize),
    CapRail(i64, usize),
}

#[derive(Clone, Copy)]
struct Corner {
    column: usize,
    rail: i64,
    z: f64,
}

struct BoreRim {
    vertices: Vec<u32>,
    columns: Vec<usize>,
}

struct Builder<'a> {
    thread: &'a Thread,
    steps: usize,
    angular: Option<Vec<f64>>,
    strips: Vec<Strip>,
    limit: usize,
    lead: f64,
    length: f64,
    lower: f64,
    pitch: f64,
    snap: f64,
    vertex_ids: HashMap<VertexKey, u32>,
    cap_rails: HashMap<(usize, usize), i64>,
    caps: [Vec<(f64, u32)>; 2],
    mesh: DisplayMesh,
    work: Option<&'a Counted<'a>>,
}

impl<'a> Builder<'a> {
    fn charge(&self, n: usize) -> Result<(), MeshError> {
        if let Some(w) = self.work {
            w.step(n)?;
        }
        Ok(())
    }
    fn new(thread: &'a Thread, steps: usize, strips: Vec<Strip>, options: Tessellation) -> Self {
        let spec = thread.spec();
        let lead = if spec.left_hand { -1. } else { 1. } * thread.lead_mm();
        let smallest_strip = strips
            .iter()
            .map(|s| s.phase)
            .zip(
                strips
                    .iter()
                    .skip(1)
                    .map(|s| s.phase)
                    .chain([spec.pitch_mm]),
            )
            .map(|(lo, hi)| hi - lo)
            .fold(f64::INFINITY, f64::min);
        Self {
            thread,
            steps,
            angular: None,
            strips,
            limit: options.triangle_limit,
            lead,
            length: spec.length_mm,
            lower: 0.,
            pitch: spec.pitch_mm,
            snap: (64. * f64::EPSILON * spec.length_mm.max(lead.abs()))
                .min(smallest_strip / 8.)
                .min(options.deflection * 1e-3),
            vertex_ids: HashMap::new(),
            cap_rails: HashMap::new(),
            caps: [Vec::new(), Vec::new()],
            mesh: DisplayMesh {
                vertices: Vec::new(),
                indices: Vec::new(),
                face_indices: Vec::new(),
            },
            work: None,
        }
    }

    fn build(mut self, exterior: Option<BoreExterior<'_>>) -> Result<DisplayMesh, MeshError> {
        self.build_wall()?;
        let caps = std::mem::take(&mut self.caps);
        if let Some(exterior) = exterior {
            self.add_bore(caps, exterior)?;
        } else {
            for (cap, ring) in caps.into_iter().enumerate() {
                self.add_cap(cap, ring)?;
            }
        }
        Ok(self.mesh)
    }

    fn insert_angular_stations(
        &mut self,
        stations: &[f64],
        coupling_angle: f64,
    ) -> Result<(), MeshError> {
        if !coupling_angle.is_finite()
            || coupling_angle <= 0.
            || !stations.iter().all(|x| x.is_finite())
        {
            return Err(MeshError::UnresolvedBoundary);
        }
        let mut columns = (0..=self.steps)
            .map(|i| i as f64 / self.steps as f64)
            .collect::<Vec<_>>();
        for &angle in stations {
            let value = angle.rem_euclid(TAU) / TAU;
            let index = columns.partition_point(|&x| x < value);
            let near = columns[index.saturating_sub(1)..=(index.min(columns.len() - 1))]
                .iter()
                .any(|&x| (x - value).abs() <= coupling_angle / TAU);
            if !near {
                columns.insert(index, value);
            }
        }
        if (columns.len() - 1).saturating_mul(2) > self.limit {
            return Err(MeshError::BudgetExceeded);
        }
        self.steps = columns.len() - 1;
        self.angular = Some(columns);
        Ok(())
    }
    fn insert_cap_stations(&mut self) -> Result<(), MeshError> {
        let stations = self.strips.len().saturating_mul(2);
        self.charge(
            stations
                .saturating_mul(self.steps.saturating_add(stations))
                .saturating_add(128),
        )?;
        let mut angles = Vec::with_capacity(stations);
        for cap in [self.lower, self.length] {
            let cap_phase = cap.rem_euclid(self.pitch);
            for strip in &self.strips {
                // Every intersection of a finite axial cap and a phase rail
                // becomes a shared angular COLUMN. Adjacent bands therefore
                // split the same original rail with identical vertex keys.
                let fraction = ((cap_phase - strip.phase) / self.lead).rem_euclid(1.);
                angles.push(TAU * fraction);
            }
        }
        // The existing finite clipping snap covers this phase/angle conversion
        // and the inverse rail evaluation. Merge only within that same bound;
        // no geometric feature is proximity-welded or assigned a new owner.
        self.insert_angular_stations(&angles, TAU * self.snap / self.lead.abs())
    }

    fn build_wall(&mut self) -> Result<(), MeshError> {
        let per_period = self.strips.len() as i64;
        for column in 0..self.steps {
            let a = self.advance(column);
            let b = self.advance(column + 1);
            // Include the periodic cells straddling either end; the
            // actual axial planes clip them before adding any mesh data.
            let first = ((self.lower - a.max(b)) / self.pitch).floor() as i64;
            let last = ((self.length - a.min(b)) / self.pitch).floor() as i64;
            for period in first..=last {
                for strip in 0..self.strips.len() {
                    self.charge(256)?;
                    let rail = period * per_period + strip as i64;
                    let lo = self.rail_position(rail);
                    let hi = self.rail_position(rail + 1);
                    if hi + a.max(b) <= self.lower + self.snap
                        || lo + a.min(b) >= self.length - self.snap
                    {
                        continue;
                    }
                    let corners = [
                        self.corner(column, rail, lo + a),
                        self.corner(column + 1, rail, lo + b),
                        self.corner(column + 1, rail + 1, hi + b),
                        self.corner(column, rail + 1, hi + a),
                    ];
                    let face = self.thread.face_index(period, self.strips[strip].band);
                    self.side(corners, face)?;
                }
            }
        }
        Ok(())
    }

    fn sort_ring(ring: &mut [(f64, u32)]) -> Result<(), MeshError> {
        ring.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));
        if ring.len() < 3 || !ring.windows(2).all(|pair| pair[0].0 < pair[1].0) {
            return Err(MeshError::UnresolvedBoundary);
        }
        Ok(())
    }

    fn add_cap(&mut self, cap: usize, mut ring: Vec<(f64, u32)>) -> Result<(), MeshError> {
        Self::sort_ring(&mut ring)?;
        let center = self.mesh.vertices.len() as u32 / 3;
        self.mesh.vertices.extend([0., 0., self.cap_z(cap)]);
        for (&(_, a), &(_, b)) in ring.iter().zip(ring.iter().cycle().skip(1)) {
            let triangle = if cap == 0 {
                [center, b, a]
            } else {
                [center, a, b]
            };
            self.triangle(triangle, cap as u32)?;
        }
        Ok(())
    }

    fn add_bore(
        &mut self,
        [bottom, top]: [Vec<(f64, u32)>; 2],
        exterior: BoreExterior<'_>,
    ) -> Result<(), MeshError> {
        // Each rim point adds two cap triangles and one triangle on the
        // outer wall. Check the budget before allocating either outer rim.
        if self.mesh.face_indices.len() + 3 * (bottom.len() + top.len()) > self.limit {
            return Err(MeshError::BudgetExceeded);
        }
        for triangle in self.mesh.indices.as_chunks_mut::<3>().0 {
            triangle.swap(1, 2);
        }
        let bottom = self.add_annular_cap(0, bottom, exterior)?;
        let top = self.add_annular_cap(1, top, exterior)?;
        let first_face = 2 + 4 * self.thread.spec().starts;
        for column in 0..self.steps {
            let face = match exterior {
                BoreExterior::Cylinder(_) => first_face,
                BoreExterior::Hex(_) => {
                    let middle = TAU * (column as f64 + 0.5) / self.steps as f64;
                    first_face + hex_side(middle)
                }
            };
            let low = &bottom.vertices[bottom.columns[column]..=bottom.columns[column + 1]];
            let high = &top.vertices[top.columns[column]..=top.columns[column + 1]];
            // The two rim subdivisions can differ. Both lie on the same
            // flat exterior panel, so split the panel on its diagonal.
            for &[a, b] in low.array_windows() {
                self.triangle([a, b, high[high.len() - 1]], face)?;
            }
            for &[a, b] in high.array_windows() {
                self.triangle([low[0], b, a], face)?;
            }
        }
        Ok(())
    }

    fn add_annular_cap(
        &mut self,
        cap: usize,
        mut ring: Vec<(f64, u32)>,
        exterior: BoreExterior<'_>,
    ) -> Result<BoreRim, MeshError> {
        Self::sort_ring(&mut ring)?;
        let mut columns = (0..self.steps)
            .map(|column| {
                let angle = TAU * self.fraction(column);
                ring.binary_search_by(|sample| sample.0.total_cmp(&angle))
                    .map_err(|_| MeshError::UnresolvedBoundary)
            })
            .collect::<Result<Vec<_>, _>>()?;
        columns.push(ring.len());
        let mut outer = Vec::with_capacity(ring.len() + 1);
        let z = self.cap_z(cap);
        for (column, &[first, end]) in columns.array_windows().enumerate() {
            let lo = TAU * self.fraction(column);
            let hi = TAU * self.fraction(column + 1);
            let middle = (lo + hi) * 0.5;
            let half_chord_cos = ((hi - lo) * 0.5).portable_cos();
            for &(angle, _) in &ring[first..end] {
                // The inner and outer vertices share the same radial ray.
                // Hex corners divide columns, so every segment is planar.
                let radial = match exterior {
                    BoreExterior::Cylinder(radius) => {
                        radius * half_chord_cos / (angle - middle).portable_cos()
                    }
                    BoreExterior::Hex(nut) => nut.outer_radius(angle),
                };
                let (sin, cos) = angle.portable_sin_cos();
                outer.push(self.mesh.vertices.len() as u32 / 3);
                self.mesh.vertices.extend([radial * cos, radial * sin, z]);
            }
        }
        outer.push(outer[0]);
        for (i, (&(_, a), &(_, b))) in ring.iter().zip(ring.iter().cycle().skip(1)).enumerate() {
            let [c, d] = [outer[i], outer[i + 1]];
            for [a, b, c] in [[a, c, d], [a, d, b]] {
                let triangle = if cap == 0 { [a, c, b] } else { [a, b, c] };
                self.triangle(triangle, cap as u32)?;
            }
        }
        Ok(BoreRim {
            vertices: outer,
            columns,
        })
    }

    fn advance(&self, column: usize) -> f64 {
        self.lead * self.fraction(column)
    }

    fn fraction(&self, column: usize) -> f64 {
        self.angular
            .as_ref()
            .map_or_else(|| column as f64 / self.steps as f64, |v| v[column])
    }

    fn rail_phase(&self, rail: i64) -> f64 {
        let n = self.strips.len() as i64;
        self.strips[rail.rem_euclid(n) as usize].phase
    }

    fn rail_position(&self, rail: i64) -> f64 {
        let n = self.strips.len() as i64;
        (rail.div_euclid(n) as f64).mul_add(self.pitch, self.rail_phase(rail))
    }

    fn cap_z(&self, cap: usize) -> f64 {
        if cap == 0 { self.lower } else { self.length }
    }

    fn corner(&self, column: usize, rail: i64, z: f64) -> Corner {
        let z = if (z - self.lower).abs() <= self.snap {
            self.lower
        } else if (z - self.length).abs() <= self.snap {
            self.length
        } else {
            z
        };
        Corner { column, rail, z }
    }

    fn corner_key(&mut self, corner: Corner) -> Result<Option<VertexKey>, MeshError> {
        let column = corner.column % self.steps;
        let shift = if corner.column == self.steps {
            self.lead.signum() as i64
                * i64::from(self.thread.spec().starts)
                * self.strips.len() as i64
        } else {
            0
        };
        if corner.z == self.lower || corner.z == self.length {
            let cap = usize::from(corner.z == self.length);
            let rail = corner.rail + shift;
            if let Some(old) = self.cap_rails.insert((column, cap), rail)
                && old != rail
            {
                return Err(MeshError::UnresolvedBoundary);
            }
            Ok(Some(VertexKey::CapColumn(column, cap)))
        } else if (self.lower..self.length).contains(&corner.z) {
            // The same exact rail-owner shift closes one complete revolution.
            Ok(Some(VertexKey::Grid(column, corner.rail + shift)))
        } else {
            Ok(None)
        }
    }

    fn side(&mut self, corners: [Corner; 4], face: u32) -> Result<(), MeshError> {
        let mut polygon = [0; 6];
        let mut count = 0;
        for (a, b) in corners.into_iter().zip(corners.into_iter().cycle().skip(1)) {
            if let Some(key) = self.corner_key(a)? {
                polygon[count] = self.vertex(key);
                count += 1;
            }
            let cap_order = if a.z <= b.z { [0, 1] } else { [1, 0] };
            for cap in cap_order {
                let z = self.cap_z(cap);
                if (a.z < z && z < b.z) || (b.z < z && z < a.z) {
                    let key = if a.column == b.column {
                        VertexKey::CapColumn(a.column % self.steps, cap)
                    } else {
                        debug_assert_eq!(a.rail, b.rail);
                        VertexKey::CapRail(a.rail, cap)
                    };
                    polygon[count] = self.vertex(key);
                    count += 1;
                }
            }
        }
        for i in 2..count {
            self.triangle([polygon[0], polygon[i - 1], polygon[i]], face)?;
        }
        Ok(())
    }

    fn vertex(&mut self, key: VertexKey) -> u32 {
        if let Some(&id) = self.vertex_ids.get(&key) {
            return id;
        }
        let (angle, z, radius, cap) = match key {
            VertexKey::Grid(column, rail) => (
                TAU * self.fraction(column),
                self.rail_position(rail) + self.advance(column),
                self.thread.radius_at_phase(self.rail_phase(rail)),
                None,
            ),
            VertexKey::CapColumn(column, cap) => {
                let angle = TAU * self.fraction(column);
                let z = self.cap_z(cap);
                (angle, z, self.thread.radius(angle, z), Some(cap))
            }
            VertexKey::CapRail(rail, cap) => {
                let z = self.cap_z(cap);
                let angle = TAU * ((z - self.rail_position(rail)) / self.lead);
                (
                    angle,
                    z,
                    self.thread.radius_at_phase(self.rail_phase(rail)),
                    Some(cap),
                )
            }
        };
        let (sin, cos) = angle.portable_sin_cos();
        let id = self.mesh.vertices.len() as u32 / 3;
        self.mesh.vertices.extend([radius * cos, radius * sin, z]);
        self.vertex_ids.insert(key, id);
        if let Some(cap) = cap {
            self.caps[cap].push((angle, id));
        }
        id
    }

    fn triangle(&mut self, triangle: [u32; 3], face: u32) -> Result<(), MeshError> {
        self.charge(128)?;
        if self.mesh.face_indices.len() == self.limit {
            return Err(MeshError::BudgetExceeded);
        }
        if triangle[0] == triangle[1] || triangle[1] == triangle[2] || triangle[0] == triangle[2] {
            return Err(MeshError::UnresolvedBoundary);
        }
        self.mesh.indices.extend(triangle);
        self.mesh.face_indices.push(face);
        Ok(())
    }
}

#[cfg(test)]
mod tests;
