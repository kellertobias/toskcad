//! Directed literal-point and whole-cell correspondence. Producer arithmetic
//! is observed as data, never covered by an assumed epsilon multiplier.
use super::*;
use crate::{
    display::DisplayMesh,
    head_profile::{Curve, Segment},
    helical_thread_mesh::{PhillipsScrewMesh, PhillipsThreadVertexKind as Kind},
    phillips_drive::ScrewFace,
};
const LIMIT: usize = 512_000_000;
const FIXED: usize = 65536;
#[derive(Clone, Copy)]
enum Surface {
    Thread(usize),
    Head(usize),
}
#[derive(Clone, Copy)]
struct Sample {
    parameter: I,
    angle: I,
    axial: I,
    error: f64,
    period: i64,
}
pub struct PeriodicLiteralCertificate<'s, 'm> {
    source: &'s OwnedPeriodicMaterial,
    candidate: &'m PhillipsScrewMesh,
    distance: f64,
    outer_vertex: f64,
    work: usize,
    peak: usize,
}
impl PeriodicLiteralCertificate<'_, '_> {
    pub fn source(&self) -> &OwnedPeriodicMaterial {
        self.source
    }
    pub fn candidate(&self) -> &PhillipsScrewMesh {
        self.candidate
    }
    pub fn maximum_surface_error_mm(&self) -> f64 {
        self.distance
    }
    pub fn maximum_outer_vertex_error_mm(&self) -> f64 {
        self.outer_vertex
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn accounted_peak_bytes(&self) -> usize {
        self.peak
    }
}
fn point(mesh: &DisplayMesh, id: u32) -> [f64; 3] {
    std::array::from_fn(|i| mesh.vertices[3 * id as usize + i])
}
fn norm(v: [I; 3]) -> I {
    let q = v.into_iter().fold(I::point(0.), |s, x| s + x.square());
    I::new(0., q.hi).sqrt()
}
fn distance(a: [I; 3], b: [f64; 3]) -> f64 {
    norm(std::array::from_fn(|i| a[i] - I::point(b[i]))).hi
}
fn width(v: [I; 3]) -> I {
    I::point(v.iter().map(|x| x.hi).fold(f64::NEG_INFINITY, f64::max))
        - I::point(v.iter().map(|x| x.lo).fold(f64::INFINITY, f64::min))
}
// Min/max are monotone in each argument: this encloses the exact expression
// clamp(z,lower,upper), retaining both nonrepresentable cap expressions.
fn axial(z: I, caps: [I; 2]) -> I {
    I::new(
        z.lo.max(caps[0].lo).min(caps[1].lo),
        z.hi.max(caps[0].hi).min(caps[1].hi),
    )
}
fn unit(p: I, start: I, end: I) -> I {
    let u = (p - start) / (end - start);
    I::new(u.lo.clamp(0., 1.), u.hi.clamp(0., 1.))
}
struct Core<'o, 'i, 'w> {
    owner: &'o OwnedPeriodicMaterial,
    rotation: rotation::Rotation,
    pi: I,
    work: &'w mut Work<'i>,
}
impl Core<'_, '_, '_> {
    fn thread(&mut self, c: &PhillipsScrewMesh, id: u32, band: usize) -> FResult<Sample> {
        self.work.step(64)?;
        let witness = c
            .thread_vertices()
            .binary_search_by_key(&id, |v| v.vertex)
            .ok()
            .map(|i| c.thread_vertices()[i])
            .ok_or(PhillipsReadError::Unsupported)?;
        let s = &self.owner.thread;
        let pitch = I::point(s.pitch_mm());
        let hand = I::point(if s.left_hand() { -1. } else { 1. });
        let caps = self.owner.caps.map(expression);
        let (phase, z, event) = match witness.kind {
            Kind::Grid {
                fraction,
                phase,
                period,
            } => {
                let phase = I::point(phase) + I::point(period as f64) * pitch;
                let z = axial(phase + hand * pitch * I::point(fraction), caps);
                (phase, z, true)
            }
            Kind::CapRail { phase, period, cap } => (
                I::point(phase) + I::point(period as f64) * pitch,
                caps[cap],
                true,
            ),
            Kind::CapColumn {
                fraction,
                cap,
                rail,
            } => {
                let z = caps[cap];
                if let Some((period, phase)) = rail {
                    (I::point(phase) + I::point(period as f64) * pitch, z, true)
                } else {
                    (z - hand * pitch * I::point(fraction), z, false)
                }
            }
        };
        let nodes = s.phases();
        let center = (nodes[band] + nodes[band + 1]) / I::point(2.);
        // The integer is only a proposed representation of this periodic
        // field. Every resulting point and complete cell is bounded below;
        // a wrong period cannot earn a false distance/parameter receipt.
        let quotient = (phase - center) / pitch;
        if !quotient.valid() || quotient.magnitude() > 1e12 {
            return Err(PhillipsReadError::Unsupported);
        }
        let k = quotient.midpoint().round();
        let relative = phase - I::point(k) * pitch;
        let u = unit(relative, nodes[band], nodes[band + 1]);
        // A source-owned cap/rail event is mapped to that very real phase
        // rail. Other columns must stay within the actual assigned band;
        // refusing uncertain crossings preserves finite patch coverage.
        if !event && (relative.lo < nodes[band].hi || relative.hi > nodes[band + 1].lo) {
            return Err(PhillipsReadError::Unsupported);
        }
        let selected = nodes[band] + u * (nodes[band + 1] - nodes[band]) + I::point(k) * pitch;
        let angle = hand * I::point(2.) * self.pi * (z - selected) / pitch;
        let cs = self.rotation.sin_cos(angle, self.work)?;
        let r = s.radii(band, u);
        let exact = [r * cs[0], r * cs[1], z];
        let error = distance(exact, point(c.mesh(), id));
        if !u.valid() || !angle.valid() || !error.is_finite() {
            return Err(PhillipsReadError::Unsupported);
        }
        Ok(Sample {
            parameter: u,
            angle,
            axial: z,
            error,
            period: k as i64,
        })
    }
    fn head(&mut self, c: &PhillipsScrewMesh, id: u32, segment: usize) -> FResult<Sample> {
        self.work.step(4096)?;
        let witness = c
            .head_vertices()
            .binary_search_by_key(&(id, segment), |v| (v.vertex, v.segment))
            .ok()
            .map(|i| c.head_vertices()[i])
            .ok_or(PhillipsReadError::Unsupported)?;
        let t = I::point(witness.parameter);
        let proposal = self.owner.head.view();
        let s = proposal
            .profile_proposal()
            .segments()
            .nth(segment)
            .ok_or(PhillipsReadError::Unsupported)?;
        let p = meridian(s, t)?;
        let angle = I::point(witness.angle);
        let cs = self.rotation.sin_cos(angle, self.work)?;
        let exact = [
            p[0] * cs[0],
            p[0] * cs[1],
            p[1] - I::point(self.owner.head_shift),
        ];
        let error = distance(exact, point(c.mesh(), id));
        if !(0. ..=1.).contains(&witness.parameter) || !error.is_finite() {
            return Err(PhillipsReadError::Unsupported);
        }
        Ok(Sample {
            parameter: t,
            angle,
            axial: exact[2],
            error,
            period: 0,
        })
    }
}
fn meridian(s: Segment, t: I) -> FResult<[I; 2]> {
    match s.curve {
        Curve::Line => Ok(std::array::from_fn(|i| {
            I::point(s.start[i]) + t * (I::point(s.end[i]) - I::point(s.start[i]))
        })),
        Curve::Arc {
            center,
            radius_mm,
            sweep,
        } => {
            // This cached atan2 is the proposal's defining constant. Its
            // mathematical trigonometric evaluation is enclosed separately;
            // the observed producer output is compared with this enclosure.
            let start = (s.start[1] - center[1]).portable_atan2(s.start[0] - center[0]);
            let cs = crate::head_profile::proposal_rotation(I::point(start) + t * I::point(sweep))
                .map_err(
                    crate::step_read::phillips_cutter::physical::body::closed_head::repaired_error,
                )?;
            Ok([
                I::point(center[0]) + I::point(radius_mm) * cs[0],
                I::point(center[1]) + I::point(radius_mm) * cs[1],
            ])
        }
    }
}
fn unwrap(samples: &mut [Sample; 3], pi: I) -> FResult<[i64; 3]> {
    let mut shifts = [0_i64; 3];
    let first = samples[0].angle.midpoint();
    for (j, s) in samples[1..].iter_mut().enumerate() {
        let turns = ((first - s.angle.midpoint()) / (2. * std::f64::consts::PI)).round();
        s.angle = s.angle + I::point(2. * turns) * pi;
        shifts[j + 1] = turns as i64;
    }
    // Verified after the proposed integer choices; no angle/proximity guess
    // is used as a source correspondence or an admission condition.
    if width(samples.map(|s| s.angle)).hi >= 1. {
        return Err(PhillipsReadError::Unsupported);
    }
    Ok(shifts)
}
impl<'s> PeriodicAuthoredApproximation<'s> {
    /// New literal/cell audit stage, bounded by512M logical operations. The
    /// consumer debits preparation and producer work before supplying this
    /// remaining allowance. Bytes cover the borrowed owner/candidate and new
    /// cache; existing producer delegated workspaces retain their contracts.
    pub fn certify_literal<'m>(
        &self,
        c: &'m PhillipsScrewMesh,
        tolerance: f64,
        maximum_work: usize,
        maximum_bytes: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> FResult<PeriodicLiteralCertificate<'s, 'm>> {
        let mut work = Work {
            remaining: maximum_work,
            interrupt,
            cancelled: Cell::new(false),
        };
        let result = (|| {
            work.step(1)?;
            if maximum_work > LIMIT
                || maximum_bytes > BYTES
                || !tolerance.is_finite()
                || tolerance <= 0.
            {
                return Err(PhillipsReadError::Invalid);
            }
            let actual = c.proposal();
            let expected = &self.proposal;
            if actual.thread().spec() != expected.thread().spec()
                || actual.spec() != expected.spec()
                || actual.profile().top() != expected.profile().top()
                || actual.profile().end() != expected.profile().end()
                || !actual
                    .profile()
                    .segments()
                    .eq(expected.profile().segments())
                || actual.cavity().source().spec() != self.owner.cavity.spec()
                || actual.cavity().source().domain() != self.owner.cavity.domain()
            {
                return Err(PhillipsReadError::Unsupported);
            }
            let vertices = c.mesh().vertices.len() / 3;
            let peak = self
                .owner
                .retained_budget_bytes()
                .saturating_add(c.retained_budget_bytes())
                .saturating_add(vertices.saturating_mul(std::mem::size_of::<Option<Sample>>()))
                .saturating_add(FIXED);
            if peak > maximum_bytes {
                return Err(PhillipsReadError::Limit);
            }
            work.step(vertices)?;
            let mut cache = vec![None; vertices];
            let rotation = rotation::Rotation::new(&mut work)?;
            let pi = rotation.pi(&mut work)?;
            let mut core = Core {
                owner: self.owner,
                rotation,
                pi,
                work: &mut work,
            };
            let mut largest = 0_f64;
            let mut outer_vertex = 0_f64;
            let nodes = self.owner.thread.phases();
            let radii = self.owner.thread.radius_nodes_mm();
            for (face_index, &face) in c.faces().iter().enumerate() {
                let surface = match face {
                    ScrewFace::ThreadWall(i) => Surface::Thread([2, 3, 0, 1][i as usize]),
                    ScrewFace::ProfileWall(i) => Surface::Head(i as usize),
                    _ => continue,
                };
                core.work.step(vertices)?;
                cache.fill(None);
                for (tri, &owner) in c
                    .mesh()
                    .indices
                    .as_chunks::<3>()
                    .0
                    .iter()
                    .zip(&c.mesh().face_indices)
                {
                    core.work.step(1)?;
                    if owner as usize != face_index {
                        continue;
                    }
                    let mut samples = [Sample {
                        parameter: I::point(0.),
                        angle: I::point(0.),
                        axial: I::point(0.),
                        error: 0.,
                        period: 0,
                    }; 3];
                    for (sample, &id) in samples.iter_mut().zip(tri) {
                        if cache[id as usize].is_none() {
                            let v = match surface {
                                Surface::Thread(band) => core.thread(c, id, band)?,
                                Surface::Head(segment) => core.head(c, id, segment)?,
                            };
                            cache[id as usize] = Some(v);
                        }
                        *sample = cache[id as usize]
                            .expect("actual vertex cache filled before its immutable read");
                    }
                    let turns = unwrap(&mut samples, pi)?;
                    core.work.step(256)?;
                    let angle_span = width(samples.map(|s| s.angle));
                    let literal_error = samples.iter().map(|s| s.error).fold(0., f64::max);
                    let (remainder, transfer) = match surface {
                        Surface::Thread(band) => {
                            let hand = if self.owner.thread.left_hand() {
                                -1_i64
                            } else {
                                1_i64
                            };
                            let periods = std::array::from_fn::<_, 3, _>(|j| {
                                samples[j].period - hand * turns[j]
                            });
                            if periods[1] != periods[0] || periods[2] != periods[0] {
                                return Err(PhillipsReadError::Unsupported);
                            }
                            let slope = (I::point(radii[band + 1]) - I::point(radii[band]))
                                / (nodes[band + 1] - nodes[band]);
                            let alpha =
                                I::point(self.owner.thread.pitch_mm()) / (I::point(2.) * pi);
                            let dz = width(samples.map(|s| s.axial));
                            let mtt = I::point(radii[band].max(radii[band + 1]))
                                + I::point(2.) * alpha * I::point(slope.magnitude());
                            (
                                (mtt * angle_span.square()
                                    + I::point(2.) * I::point(slope.magnitude()) * angle_span * dz)
                                    / I::point(8.),
                                0.,
                            )
                        }
                        Surface::Head(i) => {
                            let segment =
                                actual.profile().segments().nth(i).expect(
                                    "associated proposal has every emitted profile segment",
                                );
                            let dt = width(samples.map(|s| s.parameter));
                            let (mtt, mt_theta, r) = match segment.curve {
                                Curve::Line => (
                                    I::point(0.),
                                    // Exact ruled cone/cylinder/planar sector:
                                    // its two triangles fill the same planar
                                    // trapezoid between the end-circle chords.
                                    // At every linear meridian station, the
                                    // radial deficit is <=R(1-cos(dtheta/2))
                                    // <=R*dtheta²/8. Radial slope creates no
                                    // extra surface-set error, even on a flat
                                    // bearing; vertex residuals remain below.
                                    I::point(0.),
                                    I::point(segment.start[0].max(segment.end[0])),
                                ),
                                Curve::Arc {
                                    center,
                                    radius_mm,
                                    sweep,
                                } => (
                                    I::point(radius_mm) * I::point(sweep).square(),
                                    I::point(radius_mm) * I::point(sweep.abs()),
                                    I::point(center[0].abs()) + I::point(radius_mm),
                                ),
                            };
                            (
                                (mtt * dt.square()
                                    + I::point(2.) * mt_theta * dt * angle_span
                                    + r * angle_span.square())
                                    / I::point(8.),
                                self.head_position,
                            )
                        }
                    };
                    let vertex = (I::point(literal_error) + I::point(transfer)).hi;
                    let total = (I::point(vertex) + remainder).hi;
                    if !total.is_finite() || total > tolerance {
                        return Err(PhillipsReadError::Unsupported);
                    }
                    largest = largest.max(total);
                    outer_vertex = outer_vertex.max(vertex);
                }
            }
            // The actual same Cutter produced the retained source-certified
            // wire/side candidate. Compare every literal translated coordinate
            // with its original literal point plus the REAL mouth expression.
            let mouth = expression(self.owner.mouth);
            let mut translated = 0_f64;
            for v in c.cavity_vertices() {
                core.work.step(128)?;
                let original = [
                    I::point(v.original[0]),
                    I::point(v.original[1]),
                    I::point(v.original[2]) + mouth,
                ];
                translated = translated.max(distance(original, point(c.mesh(), v.vertex)));
            }
            let cavity = (I::point(c.cavity_source_distance_mm()) + I::point(translated)).hi;
            largest = largest.max(cavity);
            // Planar tip/top/shoulder fills add no interior curvature. Their
            // exact ordered boundary uses are the same sewn side/wire uses:
            // source cap/rail owners above, or the actual traced Cutter mouth.
            // The literal CDT preserves those constraints and winding. The
            // boundary correspondence gives a homotopy in this error tube;
            // outside that tube the planar winding number cannot change.
            // Thus the symmetric difference of each planar fill is inside the
            // same boundary tube already bounded by the side-cell/wire tests.
            if largest > tolerance {
                return Err(PhillipsReadError::Unsupported);
            }
            core.work.step(1)?;
            Ok(PeriodicLiteralCertificate {
                source: self.owner,
                candidate: c,
                distance: largest,
                outer_vertex,
                work: maximum_work - core.work.remaining,
                peak,
            })
        })();
        result.map_err(|e| {
            if work.cancelled.get() {
                PhillipsReadError::Cancelled
            } else {
                e
            }
        })
    }
}
