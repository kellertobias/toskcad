//! Actual connected canonical Hanger candidate. Native literal embedding remains
//! mandatory; source mass and literal signed volume are distinct receipts.
use super::{exact::Q, *};
use crate::{
    authored_hanger_bolt::{FaceKey, Hanger},
    authored_metric_thread::MetricShell,
};
use std::{cmp::Ordering, mem::size_of};
mod caps;
mod closed;
mod ring;
mod transfer;
use ring::{Ring, Wall};
type E = crate::authored_hanger_bolt::Error;
type R<T> = Result<T, E>;
const BYTES: usize = 128 * 1024 * 1024;
const SCRATCH: usize = 8 * 1024 * 1024;
pub struct HangerCandidate {
    source: Hanger,
    mesh: DisplayMesh,
    chord: f64,
    normal: f64,
    work: usize,
    source_volume: [f64; 2],
    literal_volume: [f64; 2],
}
impl HangerCandidate {
    pub fn source(&self) -> &Hanger {
        &self.source
    }
    pub fn mesh(&self) -> &DisplayMesh {
        &self.mesh
    }
    pub fn into_mesh(self) -> DisplayMesh {
        self.mesh
    }
    pub fn source_chord_bound_mm(&self) -> f64 {
        self.chord
    }
    pub fn source_normal_angle_bound_rad(&self) -> f64 {
        self.normal
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn source_volume_mm3(&self) -> [f64; 2] {
        self.source_volume
    }
    pub fn literal_volume_mm3(&self) -> [f64; 2] {
        self.literal_volume
    }
}
impl From<Error> for E {
    fn from(e: Error) -> Self {
        match e {
            Error::Invalid => E::Invalid,
            Error::Limit => E::Limit,
            Error::Unresolved => E::Unresolved,
            Error::Cancelled => E::Cancelled,
        }
    }
}
impl From<crate::authored_flat_thread::Error> for E {
    fn from(e: crate::authored_flat_thread::Error) -> Self {
        use crate::authored_flat_thread::Error as F;
        match e {
            F::Invalid => E::Invalid,
            F::Limit => E::Limit,
            F::Unresolved => E::Unresolved,
            F::Cancelled => E::Cancelled,
            F::Unsupported => E::Unsupported,
        }
    }
}
impl From<crate::authored_metric_thread::Error> for E {
    fn from(e: crate::authored_metric_thread::Error) -> Self {
        use crate::authored_metric_thread::Error as F;
        match e {
            F::Invalid => E::Invalid,
            F::Limit => E::Limit,
            F::Unresolved => E::Unresolved,
            F::Cancelled => E::Cancelled,
        }
    }
}
struct Memory {
    used: usize,
}
impl Memory {
    fn charge(&mut self, bytes: usize) -> R<()> {
        self.used = self
            .used
            .checked_add(bytes)
            .filter(|&v| v <= BYTES)
            .ok_or(E::Limit)?;
        Ok(())
    }
    fn allocate<T>(&mut self, count: usize, w: &mut Work<'_>) -> R<Vec<T>> {
        w.step()?;
        let bytes = count.checked_mul(size_of::<T>()).ok_or(E::Limit)?;
        if bytes > BYTES - self.used {
            return Err(E::Limit);
        }
        let mut result = Vec::new();
        result.try_reserve_exact(count).map_err(|_| E::Limit)?;
        self.charge(result.capacity() * size_of::<T>())?;
        Ok(result)
    }
}
fn point(mesh: &DisplayMesh, id: usize) -> [f64; 3] {
    mesh.vertices[id * 3..id * 3 + 3].try_into().unwrap()
}
fn add_triangle(
    mesh: &mut DisplayMesh,
    ids: [usize; 3],
    face: u32,
    options: Tessellation,
) -> R<()> {
    if mesh.face_indices.len() >= options.triangle_limit
        || mesh.indices.len() + 3 > mesh.indices.capacity()
        || mesh.face_indices.len() == mesh.face_indices.capacity()
    {
        return Err(E::Limit);
    }
    mesh.indices.extend(ids.map(|i| i as u32));
    mesh.face_indices.push(face);
    Ok(())
}
fn tau() -> I {
    I::point(2.)
        * I::new(
            std::f64::consts::PI.next_down(),
            std::f64::consts::PI.next_up(),
        )
}
fn cross(a: [I; 3], b: [I; 3]) -> [I; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
#[cfg(test)]
fn dot(a: [I; 3], b: [I; 3]) -> I {
    a.into_iter()
        .zip(b)
        .fold(I::point(0.), |s, (a, b)| s + a * b)
}
#[cfg(test)]
fn normal(p: [[f64; 3]; 3]) -> [I; 3] {
    let p = p.map(|x| x.map(I::point));
    cross(
        std::array::from_fn(|i| p[1][i] - p[0][i]),
        std::array::from_fn(|i| p[2][i] - p[0][i]),
    )
}
fn copy_side(
    mesh: &mut DisplayMesh,
    wall: &DisplayMesh,
    source: &Hanger,
    metric: bool,
    options: Tessellation,
    w: &mut Work<'_>,
) -> R<(usize, f64, f64)> {
    let offset = mesh.vertices.len() / 3;
    let mut delta = 0_f64;
    let mut angle = 0_f64;
    let p = source.spec().plain_length_mm;
    let l = source.metric().full_length_mm();
    w.charge(wall.vertices.len() + wall.face_indices.len())?;
    for raw in wall.vertices.as_chunks::<3>().0 {
        let mut v = *raw;
        if metric {
            v[2] = (v[2] - l) - p;
            let exact = I::point(raw[2]) - I::point(l) - I::point(p);
            delta = delta.max((exact - I::point(v[2])).magnitude());
        }
        mesh.vertices.extend(v);
    }
    for (ids, &face) in wall
        .indices
        .as_chunks::<3>()
        .0
        .iter()
        .zip(&wall.face_indices)
    {
        let ids = ids.map(|v| v as usize);
        if metric {
            let old_points = ids.map(|id| point(wall, id));
            let bound =
                transfer::checked(old_points, ids.map(|id| point(mesh, id + offset)), delta, w)?;
            angle = angle.max(bound);
        }
        let mapped = map_face(source, metric, face)?;
        add_triangle(mesh, ids.map(|i| i + offset), mapped, options)?;
    }
    Ok((offset, delta, angle))
}
pub(crate) fn hanger_candidate(
    source: &Hanger,
    options: Tessellation,
    maximum_work: usize,
    interrupt: Option<&dyn Fn() -> bool>,
) -> R<HangerCandidate> {
    options.check().map_err(|_| E::Invalid)?;
    #[cfg(test)]
    controls::STAGE.set(0);
    let mut w = Work {
        remaining: maximum_work.min(16_000_000),
        used: 0,
        interrupt,
    };
    w.step()?;
    let mut memory = Memory {
        used: SCRATCH + size_of::<HangerCandidate>(),
    };
    let metric = source
        .metric()
        .tessellate_shell(options, w.remaining, interrupt)?;
    #[cfg(test)]
    controls::STAGE.set(1);
    w.charge(metric.work_used())?;
    memory.charge(metric.retained_bytes())?;
    let wood = source
        .wood()
        .tessellate_shell(options, w.remaining, interrupt)?;
    #[cfg(test)]
    controls::STAGE.set(2);
    w.charge(wood.work_used())?;
    memory.charge(wood.retained_bytes())?;
    let circles = caps::circle_count(source.spec().plain_radius_mm, options)?;
    let side = metric.mesh().face_indices.len() + wood.mesh().face_indices.len();
    let rims = metric.rim_count(ShellEnd::Start)
        + metric.rim_count(ShellEnd::End)
        + wood.rim_count(ShellEnd::Start)
        + wood.rim_count(ShellEnd::End);
    let triangles = side + rims + 4 * circles;
    if triangles > options.triangle_limit {
        return Err(E::Limit);
    }
    let vertices = metric.mesh().vertices.len() + wood.mesh().vertices.len() + 6 * circles + 6;
    let mut mesh = DisplayMesh {
        vertices: memory.allocate(vertices, &mut w)?,
        indices: memory.allocate(3 * triangles, &mut w)?,
        face_indices: memory.allocate(triangles, &mut w)?,
    };
    let (mo, delta, transfer) = copy_side(&mut mesh, metric.mesh(), source, true, options, &mut w)?;
    let (wo, _, _) = copy_side(&mut mesh, wood.mesh(), source, false, options, &mut w)?;
    let metric_error = metric.source_chord_bound_mm();
    let metric_angle = metric.source_normal_angle_bound_rad();
    let (metric_chord, metric_normal) =
        quality(metric_error, delta, metric_angle, transfer, options)?;
    let mut chord = metric_chord.max(wood.source_chord_bound_mm());
    let normal = metric_normal.max(wood.source_normal_angle_bound_rad());
    #[cfg(test)]
    controls::STAGE.set(3);
    let metric_bottom = ring::read(
        Wall::Metric(&metric),
        ShellEnd::Start,
        mo,
        &mut memory,
        &mut w,
    )?;
    let metric_wall = Wall::Metric(&metric);
    let metric_top = ring::read(metric_wall, ShellEnd::End, mo, &mut memory, &mut w)?;
    let wood_bottom = ring::read(Wall::Wood(&wood), ShellEnd::Start, wo, &mut memory, &mut w)?;
    let wood_top = ring::read(Wall::Wood(&wood), ShellEnd::End, wo, &mut memory, &mut w)?;
    #[cfg(test)]
    controls::STAGE.set(4);
    let bottom_z = point(&mesh, metric_bottom.vertices[0].id)[2];
    let top_z = source.wood().full_length_mm();
    let (plain_bottom, plain_top, plain_chord) =
        caps::plain(&mut mesh, source, circles, options, &mut memory, &mut w)?;
    chord = chord.max(plain_chord);
    #[cfg(test)]
    controls::STAGE.set(5);
    for (outer, inner, key, z, plane_error) in [
        (
            &plain_bottom,
            Some(&metric_top),
            FaceKey::MetricShoulder,
            -source.spec().plain_length_mm,
            0.,
        ),
        (
            &wood_bottom,
            Some(&plain_top),
            FaceKey::WoodShoulder,
            0.,
            0.,
        ),
        (
            &metric_bottom,
            None,
            FaceKey::MetricTip,
            bottom_z,
            (I::point(bottom_z)
                + I::point(source.spec().plain_length_mm)
                + I::point(source.metric().full_length_mm()))
            .magnitude(),
        ),
        (&wood_top, None, FaceKey::WoodTip, top_z, 0.),
    ] {
        let cap_bound = caps::fill(
            &mut mesh,
            outer,
            inner,
            caps::Target {
                key,
                z,
                plane_error,
            },
            source,
            options,
            &mut w,
        );
        chord = chord.max(cap_bound?);
    }
    closed::certify(&mesh, &mut memory, &mut w)?;
    #[cfg(test)]
    controls::STAGE.set(6);
    #[cfg(test)]
    controls::STAGE.set(7);
    let query = crate::authored_hanger_bolt::QueryOptions {
        volume_error_mm3: 1e-7,
        maximum_cells: 11,
    };
    let volume = source.volume(query, interrupt)?.material_mm3;
    w.charge(256)?;
    let literal = closed::volume(&mesh, &mut w)?;
    #[cfg(test)]
    controls::STAGE.set(8);
    w.step()?;
    Ok(HangerCandidate {
        source: *source,
        mesh,
        chord,
        normal,
        work: w.used,
        source_volume: volume,
        literal_volume: literal,
    })
}
#[cfg(test)]
mod tests {
    use super::*;
    struct Explicit(bool);
    impl Explicit {
        fn new() -> Self {
            Self(transfer::FORCE_EXPLICIT.replace(true))
        }
    }
    impl Drop for Explicit {
        fn drop(&mut self) {
            transfer::FORCE_EXPLICIT.set(self.0);
        }
    }
    #[test]
    fn bounded_buffers_typed_errors_and_candidate_moves_preserve_their_contracts() {
        let mut w = Work {
            remaining: 1000000,
            used: 0,
            interrupt: None,
        };
        let mut m = Memory { used: BYTES - 8 };
        let actual: Vec<u64> = m.allocate(1, &mut w).unwrap();
        assert_eq!(actual.capacity(), 1);
        assert!(m.allocate::<u64>(1, &mut w).is_err());
        let mut m = Memory { used: 0 };
        assert!(m.allocate::<u64>(usize::MAX, &mut w).is_err());
        assert!(m.charge(usize::MAX).is_err());
        for (a, b) in [
            (Error::Invalid, E::Invalid),
            (Error::Limit, E::Limit),
            (Error::Unresolved, E::Unresolved),
            (Error::Cancelled, E::Cancelled),
        ] {
            assert_eq!(E::from(a), b);
        }
        use crate::authored_flat_thread::Error as F;
        for (a, b) in [
            (F::Invalid, E::Invalid),
            (F::Limit, E::Limit),
            (F::Unresolved, E::Unresolved),
            (F::Cancelled, E::Cancelled),
            (F::Unsupported, E::Unsupported),
        ] {
            assert_eq!(E::from(a), b);
        }
        use crate::authored_metric_thread::Error as M;
        for (a, b) in [
            (M::Invalid, E::Invalid),
            (M::Limit, E::Limit),
            (M::Unresolved, E::Unresolved),
            (M::Cancelled, E::Cancelled),
        ] {
            assert_eq!(E::from(a), b);
        }
        let mut mesh = DisplayMesh {
            vertices: vec![0.; 9],
            indices: Vec::new(),
            face_indices: Vec::new(),
        };
        assert!(
            add_triangle(
                &mut mesh,
                [0, 1, 2],
                0,
                Tessellation {
                    deflection: 0.1,
                    angle: 0.3,
                    triangle_limit: 1
                }
            )
            .is_err()
        );
    }
    #[test]
    fn actual_rounding_cone_and_projected_divergence_preserve_literal_geometry() {
        let source = Hanger::new(
            crate::authored_hanger_bolt::tests::fixture(false, false),
            None,
        )
        .unwrap();
        let options = Tessellation {
            deflection: 0.1,
            angle: 0.3,
            triangle_limit: 1_000_000,
        };
        let now = std::time::Instant::now();
        let explicit = {
            let _proof = Explicit::new();
            hanger_candidate(&source, options, 16_000_000, None).unwrap()
        };
        eprintln!(
            "explicitNormal {}work {:?}",
            explicit.work_used(),
            now.elapsed()
        );
        let now = std::time::Instant::now();
        let cone = hanger_candidate(&source, options, 16_000_000, None).unwrap();
        eprintln!("coneNormal {}work {:?}", cone.work_used(), now.elapsed());
        assert!(
            explicit.mesh().vertices.iter().map(|x| x.to_bits()).eq(cone
                .mesh()
                .vertices
                .iter()
                .map(|x| x.to_bits()))
        );
        assert_eq!(explicit.mesh().indices, cone.mesh().indices);
        assert_eq!(explicit.mesh().face_indices, cone.mesh().face_indices);
        assert!(cone.work_used() < explicit.work_used());
        assert!(cone.source_normal_angle_bound_rad() <= options.angle);
        let mut w = Work {
            remaining: 16_000_000,
            used: 0,
            interrupt: None,
        };
        let tetra = closed::tetra_volume(cone.mesh(), &mut w).unwrap();
        let flux = cone.literal_volume_mm3();
        assert!(tetra[0] <= flux[1] && flux[0] <= tetra[1]);
        assert_eq!(cone.source_volume_mm3(), explicit.source_volume_mm3());
        let count = cone.mesh().face_indices.len();
        assert_eq!(cone.into_mesh().face_indices.len(), count);
    }
}

fn upper(value: I) -> R<f64> {
    if value.valid() {
        Ok(value.hi)
    } else {
        Err(E::Unresolved)
    }
}
fn room(length: usize, extra: usize, capacity: usize) -> R<()> {
    if extra > capacity.saturating_sub(length) {
        Err(E::Limit)
    } else {
        Ok(())
    }
}
fn map_face(source: &Hanger, metric: bool, face: u32) -> R<u32> {
    let role = face.checked_sub(2).ok_or(E::Unresolved)?;
    let count = 4 * if metric {
        source.metric().spec().starts
    } else {
        source.wood().spec().starts
    };
    if role >= count {
        return Err(E::Unresolved);
    }
    // Disjoint ranges are precisely those exposed by Hanger::face_key.
    Ok(5 + role
        + if metric {
            0
        } else {
            4 * source.metric().spec().starts
        })
}
fn quality(
    chord: f64,
    delta: f64,
    angle: f64,
    transfer: f64,
    options: Tessellation,
) -> R<(f64, f64)> {
    let chord = upper(I::point(chord) + I::point(delta))?;
    let angle = upper(I::point(angle) + I::point(transfer))?;
    if chord > options.deflection || angle > options.angle {
        Err(E::Unresolved)
    } else {
        Ok((chord, angle))
    }
}
#[cfg(test)]
mod publication_tests {
    use super::*;
    #[test]
    fn actual_proof_budget_and_face_domain_failures_are_atomic() {
        let source = Hanger::new(
            crate::authored_hanger_bolt::tests::fixture(false, false),
            None,
        )
        .unwrap();
        assert!(map_face(&source, true, 6).is_err());
        assert!(map_face(&source, false, 0).is_err());
        let options = Tessellation {
            deflection: 0.1,
            angle: 0.3,
            triangle_limit: 1_000_000,
        };
        assert!(quality(0.099, 0.01, 0.2, 0., options).is_err());
        assert!(quality(0.01, 0., 0.29, 0.02, options).is_err());
        assert!(upper(I::new(0., f64::INFINITY)).is_err());
        assert!(upper(I::new(2., 1.)).is_err());
        assert!(upper(I::point(f64::NAN)).is_err());
        assert!(room(3, 1, 3).is_err());
        assert!(room(1, 1, 3).is_ok());
        let before = source;
        assert!(
            hanger_candidate(
                &source,
                Tessellation {
                    triangle_limit: 1,
                    ..options
                },
                16_000_000,
                None
            )
            .is_err()
        );
        assert!(hanger_candidate(&source, options, 0, None).is_err());
        assert!(matches!(
            hanger_candidate(&source, options, 16_000_000, Some(&|| true)),
            Err(E::Cancelled)
        ));
        assert_eq!(source, before);
    }
}

#[cfg(test)]
mod controls {
    use super::*;
    thread_local! { pub(super) static STAGE:std::cell::Cell<u8> = const {std::cell::Cell::new(0)}; }
    #[test]
    fn late_cancellation_at_every_actual_pipeline_stage_never_publishes_a_candidate() {
        let source = Hanger::new(
            crate::authored_hanger_bolt::tests::fixture(false, false),
            None,
        )
        .unwrap();
        let before = source;
        let options = Tessellation {
            deflection: 0.1,
            angle: 0.3,
            triangle_limit: 1_000_000,
        };
        for target in 0..=8 {
            STAGE.set(0);
            let interrupt = || STAGE.get() >= target;
            let result = hanger_candidate(&source, options, 16_000_000, Some(&interrupt));
            assert!(matches!(result, Err(E::Cancelled)), "stage{target}");
            assert_eq!(source, before);
        }
        assert!(matches!(
            hanger_candidate(
                &source,
                Tessellation {
                    triangle_limit: 41400,
                    ..options
                },
                16_000_000,
                None
            ),
            Err(E::Limit)
        ));
    }
}
