//! Actual source charts and whole-cell bounds, preceding detached mesh publication.
use super::*;
mod mesh;
mod properties;
use crate::tapered_thread::interval::J;
pub use mesh::BoundaryMesh;
pub use properties::{FaceProperties, FacePropertyOptions};
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PatchKey {
    Envelope {
        envelope: EnvelopeRole,
        wall: WallKind,
        quarter: u8,
    },
    Wall {
        envelope: EnvelopeRole,
        wall: WallKind,
        quarter: u8,
        lower: bool,
    },
    MouthTrace {
        envelope: EnvelopeRole,
        wall: WallKind,
        quarter: u8,
        lower: bool,
    },
    MouthGap {
        quarter: u8,
    },
    FloorFan {
        quarter: u8,
    },
}
#[derive(Clone, Copy, Debug)]
enum Law {
    Envelope(usize, u8),
    Wall(usize, u8, bool),
    Trace(usize, u8, bool),
    Gap(u8),
    Floor(u8),
}
#[derive(Clone, Copy, Debug)]
pub struct SurfacePatch {
    key: PatchKey,
    face: usize,
    law: Law,
}
impl SurfacePatch {
    pub fn source_identity(&self) -> PatchKey {
        self.key
    }
    pub fn face_index(&self) -> usize {
        self.face
    }
    /// Exact chart orientation inherited from the typed source law.
    pub fn winding(&self) -> i8 {
        match self.law {
            Law::Gap(_) => 1,
            Law::Wall(_, _, true) | Law::Trace(_, _, true) => 1,
            _ => -1,
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PatchBounds {
    pub point_mm: [Bounds; 3],
    pub du_mm: [Bounds; 3],
    pub dv_mm: [Bounds; 3],
    pub duu_mm: [Bounds; 3],
    pub duv_mm: [Bounds; 3],
    pub dvv_mm: [Bounds; 3],
    pub normal: [Bounds; 3],
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct QuadErrorBound {
    pub source_to_bilinear_mm: f64,
    pub literal_vertex_error_mm: f64,
    pub bilinear_to_triangles_mm: f64,
    pub normal_radius_chord: f64,
}
impl QuadErrorBound {
    pub fn positional_error_mm(self) -> f64 {
        (I::point(self.source_to_bilinear_mm)
            + I::point(self.literal_vertex_error_mm)
            + I::point(self.bilinear_to_triangles_mm))
        .hi
    }
}
#[derive(Clone, Debug)]
pub struct PatchAtlas<'a> {
    topology: BoundaryTopology<'a>,
    patches: Vec<SurfacePatch>,
    radial_margin: I,
}
fn domain(uv: [[f64; 2]; 2]) -> Result<[I; 2], Error> {
    if uv.into_iter().flatten().any(|v| !v.is_finite())
        || uv.into_iter().any(|[a, b]| a < 0. || b > 1. || a > b)
    {
        return Err(Error::Invalid);
    }
    Ok(uv.map(b))
}
fn norm(p: [I; 3]) -> I {
    p.into_iter()
        .fold(I::point(0.), |s, x| s + x.square())
        .sqrt()
}
fn cj(x: I) -> J {
    J {
        v: x,
        ..J::constant(0.)
    }
}
fn rotated([x, y]: [J; 2], q: u8) -> [J; 2] {
    match q {
        0 => [x, y],
        1 => [-y, x],
        2 => [-x, -y],
        _ => [y, -x],
    }
}
fn pi() -> I {
    I::new(PI.next_down(), PI.next_up())
}
fn framed(point: [J; 3], normal: [I; 3], q: u8, lower: bool) -> ([J; 3], [I; 3]) {
    let [x, y, z] = point;
    let [x, y] = rotated([x, if lower { -y } else { y }], q);
    let [nx, ny, nz] = normal;
    let [nx, ny] = rotate([nx, if lower { -ny } else { ny }], q);
    ([x, y, z], [nx, ny, nz])
}
// Bilinear interpolation retains mixed variation. The complete finite mixed
// difference then bounds the actual two-triangle replacement, including a saddle.
fn quad_error(
    source: I,
    vertex: f64,
    literal: [[f64; 3]; 4],
    normals: [I; 3],
    center: [I; 3],
) -> Result<QuadErrorBound, Error> {
    let warp = norm(std::array::from_fn(|i| {
        I::point(literal[0][i]) - I::point(literal[1][i]) + I::point(literal[2][i])
            - I::point(literal[3][i])
    })) / I::point(4.);
    let normal = norm(std::array::from_fn(|i| normals[i] - center[i]));
    let bound = QuadErrorBound {
        source_to_bilinear_mm: source.hi,
        literal_vertex_error_mm: vertex,
        bilinear_to_triangles_mm: warp.hi,
        normal_radius_chord: normal.hi.min(2.),
    };
    if !bound.positional_error_mm().is_finite() || !normal.valid() {
        return Err(Error::Unresolved);
    }
    Ok(bound)
}
fn cylinder_normal(k: &Corner, y: I, radius: f64) -> Result<[I; 3], Error> {
    let radius = I::point(radius);
    let square = radius.square() - y.square();
    if square.lo <= 0. {
        return Err(Error::Unresolved);
    }
    let radial = square.sqrt() / radius;
    Ok([k.cos_beta * radial, -y / radius, k.sin_beta * radial])
}
impl<'a> PatchAtlas<'a> {
    pub fn topology(&self) -> &BoundaryTopology<'a> {
        &self.topology
    }
    pub fn patches(&self) -> &[SurfacePatch] {
        &self.patches
    }
    /// Positive lower reserve for x-y*f'(y) over all actual corner fronts.
    pub fn radial_margin_mm_bounds(&self) -> Bounds {
        bounds(self.radial_margin)
    }
    fn state(&self, index: usize, uv: [I; 2]) -> Result<([J; 3], [I; 3]), Error> {
        let patch = self.patches.get(index).ok_or(Error::Invalid)?;
        let u = J::d(uv[0]);
        let v = J::z(uv[1]);
        let events = self.topology.events();
        let c = events.source();
        let k = &c.corners[0];
        let e = &c.envelope;
        let cap = [I::point(0.), I::point(0.), I::point(-1.)];
        let (point, normal) = match patch.law {
            Law::Envelope(index, q) => {
                let p = evaluate::trim_state::<J>(events, index, uv[0])?;
                let r = evaluate::root(p[0] * p[0] + p[1] * p[1])?;
                let phi = (p[1] / p[0]).atan();
                // The symmetric form retains the authored angle while avoiding
                // two independent interval occurrences of phi near either rim.
                let theta = cj(pi() / I::point(4.))
                    + (v - J::constant(0.5)) * (cj(pi() / I::point(2.)) - J::constant(2.) * phi);
                let normal = match events.trims[index].envelope {
                    EnvelopeRole::EntryRoundover => {
                        let angle = events.trims[index].parameters[0]
                            + uv[0]
                                * (events.trims[index].parameters[1]
                                    - events.trims[index].parameters[0]);
                        [
                            angle.sin() * theta.v.cos(),
                            angle.sin() * theta.v.sin(),
                            angle.cos(),
                        ]
                    }
                    EnvelopeRole::OuterFlank => {
                        [e.cos_a * theta.v.cos(), e.cos_a * theta.v.sin(), e.sin_a]
                    }
                    EnvelopeRole::Floor => {
                        [e.sin_c * theta.v.cos(), e.sin_c * theta.v.sin(), e.cos_c]
                    }
                };
                framed([r * theta.cos(), r * theta.sin(), p[2]], normal, q, false)
            }
            Law::Wall(index, q, lower) | Law::Trace(index, q, lower) => {
                let p = evaluate::trim_state::<J>(events, index, uv[0])?;
                let beta = cj(events.beta_tangent);
                let top = [p[0] + beta * p[2], p[1], J::constant(0.)];
                if matches!(patch.law, Law::Trace(..)) {
                    framed([v * top[0], v * top[1], J::constant(0.)], cap, q, lower)
                } else {
                    let normal = match events.trims[index].wall {
                        WallKind::Outer => [
                            k.cos_beta * k.sin_theta,
                            -k.cos_theta,
                            k.sin_beta * k.sin_theta,
                        ],
                        WallKind::Inner => match k.inner {
                            SourceInner::ShallowV { sine, cosine } => {
                                [k.cos_beta * sine, -cosine, k.sin_beta * sine]
                            }
                            SourceInner::Rounded { radius_mm } => {
                                cylinder_normal(k, p[1].v, radius_mm)?
                            }
                        },
                    };
                    framed(
                        [p[0] + beta * p[2] * (J::constant(1.) - v), p[1], v * p[2]],
                        normal,
                        q,
                        lower,
                    )
                }
            }
            Law::Gap(q) => {
                let p = events.event(EventRole::MouthEnvelope).point;
                let phi = (p[1] / p[0]).atan();
                let theta = cj(phi) + u * cj(pi() / I::point(2.) - I::point(2.) * phi);
                let radius = cj(e.entry) * v;
                framed(
                    [radius * theta.cos(), radius * theta.sin(), J::constant(0.)],
                    cap,
                    q,
                    false,
                )
            }
            Law::Floor(q) => {
                let radius = cj(k.b_half) * (J::constant(1.) - u);
                let depth = cj(k.b_depth) + u * cj(e.tip - k.b_depth);
                let theta = v * cj(pi() / I::point(2.));
                let normal = [e.sin_c * theta.v.cos(), e.sin_c * theta.v.sin(), e.cos_c];
                framed(
                    [radius * theta.cos(), radius * theta.sin(), depth],
                    normal,
                    q,
                    false,
                )
            }
        };
        let point = [
            point[0].checked().map_err(|_| Error::Unresolved)?,
            point[1].checked().map_err(|_| Error::Unresolved)?,
            point[2].checked().map_err(|_| Error::Unresolved)?,
        ];
        Ok((point, normal))
    }
    pub fn patch_bounds(&self, index: usize, uv: [[f64; 2]; 2]) -> Result<PatchBounds, Error> {
        let (p, n) = self.state(index, domain(uv)?)?;
        Ok(PatchBounds {
            point_mm: p.map(|j| bounds(j.v)),
            du_mm: p.map(|j| bounds(j.d)),
            dv_mm: p.map(|j| bounds(j.z)),
            duu_mm: p.map(|j| bounds(j.dd)),
            duv_mm: p.map(|j| bounds(j.dz)),
            dvv_mm: p.map(|j| bounds(j.zz)),
            normal: n.map(bounds),
        })
    }
    /// Fixed source-chart orientation, never repaired from a generated triangle.
    pub fn patch_winding(&self, index: usize) -> Result<i8, Error> {
        let patch = self.patches.get(index).ok_or(Error::Invalid)?;
        Ok(patch.winding())
    }
    pub fn certify_quad(
        &self,
        index: usize,
        uv: [[f64; 2]; 2],
        literal: [[f64; 3]; 4],
    ) -> Result<QuadErrorBound, Error> {
        let [u, v] = domain(uv)?;
        if !literal.into_iter().flatten().all(f64::is_finite) {
            return Err(Error::Invalid);
        }
        if literal.into_iter().flatten().any(|x| x.abs() > 1e9) {
            return Err(Error::Limit);
        }
        let (p, n) = self.state(index, [u, v])?;
        let hu = I::point(u.hi) - I::point(u.lo);
        let hv = I::point(v.hi) - I::point(v.lo);
        let source = (norm(p.map(|j| j.dd)) * hu.square() + norm(p.map(|j| j.zz)) * hv.square())
            / I::point(8.);
        let coords = [[u.lo, v.lo], [u.hi, v.lo], [u.hi, v.hi], [u.lo, v.hi]];
        let mut vertex = 0_f64;
        for (pair, point) in coords.into_iter().zip(literal) {
            let (source, _) = self.state(index, pair.map(I::point))?;
            vertex = vertex.max(norm(std::array::from_fn(|i| I::point(point[i]) - source[i].v)).hi)
        }
        let (_, center) = self.state(index, [I::point(u.midpoint()), I::point(v.midpoint())])?;
        quad_error(source, vertex, literal, n, center)
    }
}
fn face_owner(faces: &[BoundaryFace], key: FaceKey) -> Result<usize, Error> {
    let found: Vec<_> = faces
        .iter()
        .enumerate()
        .filter(|(_, f)| f.source_identity() == key)
        .map(|(i, _)| i)
        .collect();
    if found.len() != 1 {
        return Err(Error::Unresolved);
    }
    Ok(found[0])
}
struct FaceOwners {
    mouth: usize,
    envelope: [[usize; 3]; 4],
    wall: [[[usize; 2]; 2]; 4],
}
impl FaceOwners {
    fn read(faces: &[BoundaryFace]) -> Result<Self, Error> {
        let mouth = face_owner(faces, FaceKey::Mouth)?;
        let mut envelope = [[0; 3]; 4];
        let mut wall = [[[0; 2]; 2]; 4];
        for q in 0..4 {
            for (j, role) in [
                EnvelopeRole::EntryRoundover,
                EnvelopeRole::OuterFlank,
                EnvelopeRole::Floor,
            ]
            .into_iter()
            .enumerate()
            {
                let key = FaceKey::Envelope { role, quarter: q };
                envelope[q as usize][j] = face_owner(faces, key)?;
            }
            for (j, kind) in [WallKind::Inner, WallKind::Outer].into_iter().enumerate() {
                for lower in [false, true] {
                    let key = FaceKey::Wall {
                        kind,
                        quarter: q,
                        lower,
                    };
                    wall[q as usize][j][usize::from(lower)] = face_owner(faces, key)?;
                }
            }
        }
        Ok(Self {
            mouth,
            envelope,
            wall,
        })
    }
    fn envelope(&self, q: u8, role: EnvelopeRole) -> usize {
        self.envelope[q as usize][match role {
            EnvelopeRole::EntryRoundover => 0,
            EnvelopeRole::OuterFlank => 1,
            EnvelopeRole::Floor => 2,
        }]
    }
    fn wall(&self, q: u8, kind: WallKind, lower: bool) -> usize {
        self.wall[q as usize][usize::from(kind == WallKind::Outer)][usize::from(lower)]
    }
}
pub(in crate::phillips_drive) fn construct<'a>(
    source: &'a Cutter,
    cancel: &dyn Fn() -> bool,
) -> Result<PatchAtlas<'a>, Error> {
    let topology = source.boundary_topology_with_cancel(cancel)?;
    let owners = FaceOwners::read(topology.faces())?;
    let mut patches = Vec::new();
    for (index, trim) in topology.events().trims().iter().enumerate() {
        for quarter in 0..4 {
            checkpoint(cancel)?;
            let envelope = trim.envelope_role();
            let wall = trim.wall_kind();
            patches.push(SurfacePatch {
                key: PatchKey::Envelope {
                    envelope,
                    wall,
                    quarter,
                },
                face: owners.envelope(quarter, envelope),
                law: Law::Envelope(index, quarter),
            });
            for lower in [false, true] {
                patches.push(SurfacePatch {
                    key: PatchKey::Wall {
                        envelope,
                        wall,
                        quarter,
                        lower,
                    },
                    face: owners.wall(quarter, wall, lower),
                    law: Law::Wall(index, quarter, lower),
                });
                patches.push(SurfacePatch {
                    key: PatchKey::MouthTrace {
                        envelope,
                        wall,
                        quarter,
                        lower,
                    },
                    face: owners.mouth,
                    law: Law::Trace(index, quarter, lower),
                })
            }
        }
    }
    for quarter in 0..4 {
        checkpoint(cancel)?;
        patches.push(SurfacePatch {
            key: PatchKey::MouthGap { quarter },
            face: owners.mouth,
            law: Law::Gap(quarter),
        });
        patches.push(SurfacePatch {
            key: PatchKey::FloorFan { quarter },
            face: owners.envelope(quarter, EnvelopeRole::Floor),
            law: Law::Floor(quarter),
        })
    }
    let k = &source.corners[0];
    let rb = k.b_half + (k.b_depth - source.envelope.tip) * k.sin_beta / k.cos_beta;
    let outer = rb + (k.join_s - k.join_y * k.cos_theta / k.sin_theta) / k.cos_beta;
    // V intercepts are rb/outer. The H0 inner radial intercept decreases to the
    // identical outer tangent intercept. Admission already proves both positive.
    let radial_margin = I::new(rb.lo.min(outer.lo), rb.hi.min(outer.hi));
    Ok(PatchAtlas {
        topology,
        patches,
        radial_margin,
    })
}
#[cfg(test)]
mod tests;

// Source-only cell-complex certificate; no literal vertices/facets are emitted.
pub(in crate::phillips_drive) struct ChartIncidence {
    pub topology: TopologyCertificate,
    pub sides: [[crate::phillips_drive::SourceChartSide; 4]; 88],
}
impl PatchAtlas<'_> {
    pub(in crate::phillips_drive) fn chart_incidence(
        &self,
        cancel: &dyn Fn() -> bool,
    ) -> Result<ChartIncidence, Error> {
        mesh::chart_incidence(self, cancel)
    }
}
