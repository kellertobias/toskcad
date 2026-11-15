//! Literal thread/tail facets share original rim IDs. Manufacturing admission
//! additionally requires the Native whole literal embedding phase.
use super::super::authored::{StockClipGridLine, StockClipPointOwner};
use super::*;
use crate::{
    display::{DisplayMesh, Tessellation},
    manufacturing::{self, MeshProof},
};
fn p(v: f64) -> I {
    I::point(v)
}
fn interval(v: [f64; 2]) -> I {
    I::new(v[0], v[1])
}
fn norm(v: [I; 3]) -> I {
    (v[0].square() + v[1].square() + v[2].square()).sqrt()
}
fn cross(a: [I; 3], b: [I; 3]) -> [I; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn at(mesh: &DisplayMesh, id: u32) -> [I; 3] {
    std::array::from_fn(|i| p(mesh.vertices[id as usize * 3 + i]))
}
fn error(literal: [f64; 3], source: [I; 3]) -> f64 {
    norm(std::array::from_fn(|i| p(literal[i]) - source[i])).hi
}
fn magnitude(v: I) -> I {
    p(v.magnitude())
}
pub struct MeshCandidate {
    source: ExternalStockCut,
    mesh: DisplayMesh,
    proof: MeshProof,
    chord: f64,
    normal: f64,
    work: usize,
    bytes: usize,
}
impl MeshCandidate {
    pub fn source(&self) -> &ExternalStockCut {
        &self.source
    }
    pub fn mesh(&self) -> &DisplayMesh {
        &self.mesh
    }
    pub fn into_mesh(self) -> DisplayMesh {
        self.mesh
    }
    pub fn connectivity(&self) -> MeshProof {
        self.proof
    }
    pub fn source_distance_mm(&self) -> f64 {
        self.chord
    }
    pub fn normal_chord(&self) -> f64 {
        self.normal
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.bytes
    }
}
struct Meter<'a> {
    used: usize,
    left: usize,
    limit: usize,
    peak: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl Meter<'_> {
    fn charge(&mut self, n: usize) -> Result<(), Error> {
        if self.interrupt.is_some_and(|f| f()) {
            return Err(Error::Cancelled);
        }
        self.left = self.left.checked_sub(n).ok_or(Error::Limit)?;
        self.used += n;
        Ok(())
    }
    fn memory(&mut self, n: usize) -> Result<(), Error> {
        if n > self.limit {
            return Err(Error::Limit);
        }
        self.peak = self.peak.max(n);
        Ok(())
    }
    fn options(&self, prefix: usize) -> Result<StockClipOptions, Error> {
        if self.left == 0 {
            return Err(Error::Limit);
        }
        Ok(StockClipOptions {
            maximum_work: 16_000_000.min(prefix.saturating_add(self.left)),
            maximum_bytes: 64 * 1024,
        })
    }
    fn allowance(&self, prefix: usize) -> Result<usize, Error> {
        if self.left == 0 {
            return Err(Error::Limit);
        }
        Ok(512_000_000.min(prefix.saturating_add(self.left)))
    }
}
#[derive(Clone, Copy)]
struct Rim {
    a: u32,
    b: u32,
    band: StockClipBand,
}
struct Corner {
    outer: u32,
    top: u32,
    inner_error: f64,
    outer_error: f64,
}
fn band_index(band: StockClipBand) -> usize {
    match band {
        StockClipBand::Lower => 0,
        StockClipBand::Crest => 1,
        StockClipBand::Upper => 2,
        StockClipBand::Root => 3,
    }
}
const BANDS: [StockClipBand; 4] = [
    StockClipBand::Lower,
    StockClipBand::Crest,
    StockClipBand::Upper,
    StockClipBand::Root,
];
const RAILS: [StockClipRail; 4] = [
    StockClipRail::RootLower,
    StockClipRail::LowerCrest,
    StockClipRail::CrestUpper,
    StockClipRail::UpperRoot,
];
fn owner_band(a: StockClipPointOwner, b: StockClipPointOwner) -> Result<StockClipBand, Error> {
    let band = |owner| {
        if let StockClipPointOwner::Curve {
            curve:
                StockClipCurve::Ring {
                    station: StockClipStation::End,
                    band,
                },
            ..
        } = owner
        {
            Some(band)
        } else {
            None
        }
    };
    match (band(a), band(b)) {
        (Some(a), Some(b)) if a != b => Err(Error::Unresolved),
        (Some(a), _) | (_, Some(a)) => Ok(a),
        _ => {
            let rail = |owner| {
                if let StockClipPointOwner::Vertex(StockClipVertex {
                    station: StockClipStation::End,
                    rail,
                }) = owner
                {
                    Some(rail)
                } else {
                    None
                }
            };
            let (a, b) = (
                rail(a).ok_or(Error::Unresolved)?,
                rail(b).ok_or(Error::Unresolved)?,
            );
            BANDS
                .into_iter()
                .find(|&band| {
                    let i = band_index(band);
                    (RAILS[i] == a && RAILS[(i + 1) % 4] == b)
                        || (RAILS[i] == b && RAILS[(i + 1) % 4] == a)
                })
                .ok_or(Error::Unresolved)
        }
    }
}
fn parameter(
    owner: StockClipPointOwner,
    band: StockClipBand,
    fields: [[I; 3]; 2],
    q0: I,
    pitch: I,
    length: I,
    angular: usize,
) -> Result<I, Error> {
    match owner {
        StockClipPointOwner::Vertex(StockClipVertex {
            station: StockClipStation::End,
            rail,
        }) => {
            let i = band_index(band);
            if rail == RAILS[i] {
                Ok(p(0.))
            } else if rail == RAILS[(i + 1) % 4] {
                Ok(p(1.))
            } else {
                Err(Error::Unresolved)
            }
        }
        StockClipPointOwner::Curve {
            curve:
                StockClipCurve::Ring {
                    station: StockClipStation::End,
                    band: b,
                },
            grid: StockClipGridLine::Angular(k),
        } if b == band => {
            let phase = q0 + p(k as f64) * pitch / p(angular as f64) + length;
            let t = (phase - fields[0][1]) / (fields[1][1] - fields[0][1]);
            // The original finite cell owner already proves membership in this
            // curve. Intersect its evaluation enclosure with that earned domain.
            let t = I::new(t.lo.max(0.), t.hi.min(1.));
            if t.valid() {
                Ok(t)
            } else {
                Err(Error::Unresolved)
            }
        }
        _ => Err(Error::Unresolved),
    }
}
fn confine(mut literal: [f64; 3], cut: &ExternalStockCut) -> Result<[f64; 3], Error> {
    let center = cut.stock.circle().center().map(Rational::as_f64);
    let radial = norm([
        p(literal[0]) - cut.center[0],
        p(literal[1]) - cut.center[1],
        p(0.),
    ]);
    if radial.hi > cut.radius.lo {
        let center_error = norm([
            p(center[0]) - cut.center[0],
            p(center[1]) - cut.center[1],
            p(0.),
        ]);
        let safe = cut.radius - p(center_error.hi);
        let local = [literal[0] - center[0], literal[1] - center[1]];
        let length = norm([p(local[0]), p(local[1]), p(0.)]);
        if safe.lo <= 0. || length.lo <= 0. {
            return Err(Error::Unresolved);
        }
        let factor = (p(safe.lo) / length).lo * (1. - 16. * f64::EPSILON);
        literal[0] = center[0] + local[0] * factor;
        literal[1] = center[1] + local[1] * factor;
    }
    if norm([
        p(literal[0]) - cut.center[0],
        p(literal[1]) - cut.center[1],
        p(0.),
    ])
    .hi > cut.radius.lo
    {
        return Err(Error::Unresolved);
    }
    Ok(literal)
}
fn push_point(
    mesh: &mut DisplayMesh,
    source: [I; 3],
    cut: &ExternalStockCut,
) -> Result<(u32, f64), Error> {
    let mut literal = source.map(|v| v.lo.midpoint(v.hi));
    literal[2] = cut.join.lo;
    let literal = confine(literal, cut)?;
    let index = u32::try_from(mesh.vertices.len() / 3).map_err(|_| Error::Limit)?;
    mesh.vertices.extend(literal);
    Ok((index, error(literal, source)))
}
fn mapped_normal_error(
    original: &[f64],
    mesh: &DisplayMesh,
    cut: &ExternalStockCut,
    meter: &mut Meter<'_>,
) -> Result<f64, Error> {
    let mut maximum = 0_f64;
    let inv = p(1.) / cut.radial;
    let transfer = p(4.) * magnitude(inv - p(1.)) / p(inv.lo.min(1.));
    for (face, ids) in mesh.indices.as_chunks::<3>().0.iter().enumerate() {
        if mesh.face_indices[face] == 6 {
            continue;
        }
        meter.charge(512)?;
        let old =
            ids.map(|id| std::array::from_fn::<_, 3, _>(|i| p(original[id as usize * 3 + i])));
        let ab = std::array::from_fn(|i| old[1][i] - old[0][i]);
        let ac = std::array::from_fn(|i| old[2][i] - old[0][i]);
        let n = cross(ab, ac);
        let ideal = [
            n[0] * cut.radial,
            n[1] * cut.radial,
            n[2] * cut.radial.square(),
        ];
        let ideal_length = norm(ideal);
        let [a, b, c] = ids.map(|id| at(mesh, id));
        let actual = cross(
            std::array::from_fn(|i| b[i] - a[i]),
            std::array::from_fn(|i| c[i] - a[i]),
        );
        let actual_length = norm(actual);
        if ideal_length.lo <= 0. || actual_length.lo <= 0. {
            return Err(Error::Unresolved);
        }
        let error = norm(std::array::from_fn(|i| {
            actual[i] / actual_length - ideal[i] / ideal_length
        }));
        maximum = maximum.max((error + transfer).hi);
    }
    Ok(maximum)
}
fn push_triangle(
    mesh: &mut DisplayMesh,
    ids: [u32; 3],
    face: u32,
    theta: I,
    allowed: f64,
    meter: &mut Meter<'_>,
) -> Result<f64, Error> {
    // These are only symbolically collapsed crest corners, identified by the
    // SAME original vertex ID; no positional tolerance removes a triangle.
    if ids[0] == ids[1] || ids[1] == ids[2] || ids[2] == ids[0] {
        return if face == 6 {
            Ok(0.)
        } else {
            Err(Error::Unresolved)
        };
    }
    meter.charge(512)?;
    let [a, b, c] = ids.map(|id| at(mesh, id));
    let ab = std::array::from_fn(|i| b[i] - a[i]);
    let ac = std::array::from_fn(|i| c[i] - a[i]);
    let n = cross(ab, ac);
    let length = norm(n);
    if length.lo <= 0. {
        return Err(Error::Unresolved);
    }
    let normal = n.map(|v| v / length);
    let target = match face {
        6 => [p(0.), p(0.), p(-1.)],
        8 => [p(0.), p(0.), p(1.)],
        _ => [theta.cos(), theta.sin(), p(0.)],
    };
    let error = norm(std::array::from_fn(|i| normal[i] - target[i])).hi;
    if !error.is_finite() || error > allowed {
        return Err(Error::Unresolved);
    }
    mesh.indices.extend(ids);
    mesh.face_indices.push(face);
    Ok(error)
}
fn connected(mesh: &DisplayMesh, meter: &mut Meter<'_>) -> Result<(), Error> {
    let mut parents = (0..mesh.vertices.len() / 3)
        .map(|v| v as u32)
        .collect::<Vec<_>>();
    let mut used = vec![false; parents.len()];
    fn root(p: &mut [u32], mut i: u32) -> u32 {
        while p[i as usize] != i {
            p[i as usize] = p[p[i as usize] as usize];
            i = p[i as usize]
        }
        i
    }
    for tri in mesh.indices.as_chunks::<3>().0 {
        meter.charge(64)?;
        for &v in tri {
            used[v as usize] = true;
        }
        for &v in &tri[1..] {
            let a = root(&mut parents, tri[0]);
            let b = root(&mut parents, v);
            parents[a.max(b) as usize] = a.min(b);
        }
    }
    let mut first = None;
    for (i, used) in used.into_iter().enumerate() {
        meter.charge(1)?;
        if used {
            let r = root(&mut parents, i as u32);
            if first.is_some_and(|v| v != r) {
                return Err(Error::Unresolved);
            }
            first = Some(r);
        }
    }
    if first.is_none() {
        return Err(Error::Unresolved);
    }
    Ok(())
}
impl ExternalStockCut {
    pub fn mesh_candidate(
        &self,
        quality: Tessellation,
        maximum_work: usize,
        maximum_bytes: usize,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<MeshCandidate, Error> {
        quality.check().map_err(|_| Error::Invalid)?;
        if maximum_work == 0
            || maximum_work > 512_000_000
            || maximum_bytes == 0
            || maximum_bytes > 64 * 1024 * 1024
        {
            return Err(Error::Invalid);
        }
        let mut meter = Meter {
            used: 0,
            left: maximum_work,
            limit: maximum_bytes,
            peak: 0,
            interrupt,
        };
        meter.charge(self.work)?;
        let spec = self.thread.source().spec();
        let radius = self.stock_radius_mm();
        let pitch = p(spec.pitch_mm);
        let length = p(spec.length_mm);
        let tau = I::new(
            std::f64::consts::TAU.next_down(),
            std::f64::consts::TAU.next_up(),
        );
        meter.charge(1024)?;
        // Steering selects a small dyadic grid in BOTH directions. The
        // existing source/facet/normal certificates still decide admission.
        let mut chosen = None;
        for n in [128usize, 256, 512] {
            for m in [32usize, 64, 128, 256] {
                let omega = std::f64::consts::TAU / spec.pitch_mm;
                let ds = spec.pitch_mm / n as f64;
                let dz = spec.length_mm / m as f64;
                let bound =
                    0.5 * ((omega * omega * radius + 4. * omega) * ds * ds + 4. * omega * ds * dz);
                if bound < quality.deflection * 0.9
                    && 24. / (n as f64) < quality.angle
                    && chosen.is_none_or(|(a, b)| n * m < a * b)
                {
                    chosen = Some((n, m));
                }
            }
        }
        let (angular, axial) = chosen.ok_or(Error::Limit)?;
        let plan = self.thread.cell_plan(
            angular,
            axial,
            meter.allowance(self.thread.work_used())?,
            maximum_bytes.checked_sub(4096).ok_or(Error::Limit)?,
            interrupt,
        )?;
        meter.charge(plan.work_used() - self.thread.work_used())?;
        let old = plan.mesh_candidate(
            quality,
            meter.allowance(plan.work_used())?,
            maximum_bytes.checked_sub(4096).ok_or(Error::Limit)?,
            interrupt,
        )?;
        meter.charge(old.work_used() - plan.work_used())?;
        meter.memory(old.retained_budget_bytes() + 4096)?;
        let mut chord = old.source_distance_mm();
        let mut normal_error = old.normal_chord();
        let mut mesh = old.into_mesh();
        let reserve = plan.cap_edge_count();
        let max_vertices = mesh.vertices.len() / 3 + 2 * reserve + 1;
        let max_triangles = mesh.indices.len() / 3 + 5 * reserve;
        if max_triangles > quality.triangle_limit {
            return Err(Error::Limit);
        }
        let old_bytes = mesh.vertices.capacity() * 8
            + mesh.indices.capacity() * 4
            + mesh.face_indices.capacity() * 4;
        let metadata = plan.vertex_count() * size_of::<Option<usize>>()
            + reserve * (size_of::<Rim>() + size_of::<Corner>())
            + 64 * 1024;
        let output = 24 * max_vertices + 16 * max_triangles;
        meter.memory(
            plan.retained_budget_bytes()
                + old_bytes
                + output
                + metadata
                + 52 * max_triangles
                + 5 * max_vertices
                + mesh.vertices.len() * 8,
        )?;
        mesh.vertices
            .try_reserve_exact(6 * reserve + 3)
            .map_err(|_| Error::Limit)?;
        mesh.indices
            .try_reserve_exact(15 * reserve)
            .map_err(|_| Error::Limit)?;
        mesh.face_indices
            .try_reserve_exact(5 * reserve)
            .map_err(|_| Error::Limit)?;
        let original = mesh.vertices.clone();
        let mut rounding = 0_f64;
        if self.join.lo <= self.bottom.hi {
            return Err(Error::Unresolved);
        }
        for vertex in mesh.vertices.as_chunks_mut::<3>().0 {
            meter.charge(256)?;
            let target = [
                p(vertex[0]) * self.radial + self.center[0],
                p(vertex[1]) * self.radial + self.center[1],
                p(vertex[2]) + self.bottom,
            ];
            let mut literal = target.map(|v| v.lo.midpoint(v.hi));
            literal[2] = literal[2].clamp(self.bottom.hi, self.join.lo);
            literal = confine(literal, self)?;
            rounding = rounding.max(error(literal, target));
            *vertex = literal;
        }
        let transfer = mapped_normal_error(&original, &mesh, self, &mut meter)?;
        normal_error = (p(normal_error) + p(transfer)).hi;
        chord = (p(chord) * p(self.radial.hi.max(1.)) + p(rounding)).hi;
        if chord > quality.deflection
            || normal_error > (p(2.) * (p(quality.angle) / p(2.)).sin()).lo
        {
            return Err(Error::Unresolved);
        }
        drop(original);
        let mut rim = Vec::with_capacity(reserve);
        let mut kept = 0;
        let end_center = mesh.vertices.len() / 3 - 1;
        for i in 0..mesh.face_indices.len() {
            meter.charge(16)?;
            let ids = [
                mesh.indices[3 * i],
                mesh.indices[3 * i + 1],
                mesh.indices[3 * i + 2],
            ];
            if mesh.face_indices[i] == 6 {
                if ids[0] as usize != end_center {
                    return Err(Error::Unresolved);
                }
                let a = plan
                    .vertex_owner(ids[1] as usize)
                    .ok_or(Error::Unresolved)?;
                let b = plan
                    .vertex_owner(ids[2] as usize)
                    .ok_or(Error::Unresolved)?;
                rim.push(Rim {
                    a: ids[1],
                    b: ids[2],
                    band: owner_band(a, b)?,
                });
            } else {
                mesh.indices[3 * kept..3 * kept + 3].copy_from_slice(&ids);
                mesh.face_indices[kept] = mesh.face_indices[i];
                kept += 1;
            }
        }
        mesh.indices.truncate(3 * kept);
        mesh.face_indices.truncate(kept);
        mesh.vertices.truncate(3 * end_center);
        let mut laws = [[[p(0.); 3]; 2]; 4];
        for (i, &band) in BANDS.iter().enumerate() {
            let law = self.thread.curve_law(
                StockClipCurve::Ring {
                    station: StockClipStation::End,
                    band,
                },
                meter.options(self.thread.work_used())?,
                interrupt,
            )?;
            meter.charge(law.work_used())?;
            laws[i] = law.affine_fields_mm().map(|r| r.map(interval));
        }
        let base = self.thread.vertex_bounds(
            StockClipVertex {
                station: StockClipStation::Start,
                rail: StockClipRail::RootLower,
            },
            meter.options(self.thread.work_used())?,
            interrupt,
        )?;
        meter.charge(base.work_used() - self.thread.work_used())?;
        let q0 = interval(base.phase_mm());
        let mut map = vec![None; plan.vertex_count()];
        let mut corners = Vec::<Corner>::with_capacity(reserve);
        for edge in &rim {
            for id in [edge.a, edge.b] {
                if map[id as usize].is_some() {
                    continue;
                }
                meter.charge(256)?;
                let owner = plan.vertex_owner(id as usize).ok_or(Error::Unresolved)?;
                let (sample, crest) = match owner {
                    StockClipPointOwner::Vertex(
                        key @ StockClipVertex {
                            station: StockClipStation::End,
                            rail,
                        },
                    ) => (
                        self.vertex_bounds(
                            Vertex::Thread(key),
                            meter.options(self.work)?,
                            interrupt,
                        )?,
                        matches!(rail, StockClipRail::LowerCrest | StockClipRail::CrestUpper),
                    ),
                    StockClipPointOwner::Curve {
                        curve:
                            StockClipCurve::Ring {
                                station: StockClipStation::End,
                                band,
                            },
                        ..
                    } => {
                        let t = parameter(
                            owner,
                            band,
                            laws[band_index(band)],
                            q0,
                            pitch,
                            length,
                            angular,
                        )?;
                        (
                            self.curve_bounds(
                                Curve::Thread(StockClipCurve::Ring {
                                    station: StockClipStation::End,
                                    band,
                                }),
                                [t.lo, t.hi],
                                meter.options(self.work)?,
                                interrupt,
                            )?,
                            band == StockClipBand::Crest,
                        )
                    }
                    _ => return Err(Error::Unresolved),
                };
                meter.charge(sample.work_used() - self.work)?;
                let source = sample.position_mm().map(interval);
                let r = interval(sample.radius_mm());
                if r.lo <= 0. {
                    return Err(Error::Unresolved);
                }
                let factor = self.radius / r;
                let projected = [
                    (source[0] - self.center[0]) * factor + self.center[0],
                    (source[1] - self.center[1]) * factor + self.center[1],
                    self.join,
                ];
                let old_literal = std::array::from_fn(|i| mesh.vertices[id as usize * 3 + i]);
                let inner_error = error(old_literal, source);
                let (outer, outer_error) = if crest {
                    (id, error(old_literal, projected))
                } else {
                    push_point(&mut mesh, projected, self)?
                };
                let top = u32::try_from(mesh.vertices.len() / 3).map_err(|_| Error::Limit)?;
                mesh.vertices.extend([
                    mesh.vertices[outer as usize * 3],
                    mesh.vertices[outer as usize * 3 + 1],
                    self.top.lo,
                ]);
                map[id as usize] = Some(corners.len());
                let top_literal = std::array::from_fn(|i| mesh.vertices[top as usize * 3 + i]);
                let top_error = error(top_literal, [projected[0], projected[1], self.top]);
                corners.push(Corner {
                    outer,
                    top,
                    inner_error,
                    outer_error: outer_error.max(top_error),
                });
            }
        }
        let top_center = u32::try_from(mesh.vertices.len() / 3).map_err(|_| Error::Limit)?;
        let center = self.stock.circle().center().map(Rational::as_f64);
        mesh.vertices.extend([center[0], center[1], self.top.lo]);
        let allowed = (p(2.) * (p(quality.angle) / p(2.)).sin()).lo;
        let root_min = (I::new(
            laws[3][0][0].lo.min(laws[3][1][0].lo),
            laws[3][0][0].hi.max(laws[3][1][0].hi),
        ) * self.radial)
            .lo;
        for edge in &rim {
            meter.charge(512)?;
            let a = &corners[map[edge.a as usize].ok_or(Error::Unresolved)?];
            let b = &corners[map[edge.b as usize].ok_or(Error::Unresolved)?];
            let fields = laws[band_index(edge.band)];
            let ta = parameter(
                plan.vertex_owner(edge.a as usize)
                    .ok_or(Error::Unresolved)?,
                edge.band,
                fields,
                q0,
                pitch,
                length,
                angular,
            )?;
            let tb = parameter(
                plan.vertex_owner(edge.b as usize)
                    .ok_or(Error::Unresolved)?,
                edge.band,
                fields,
                q0,
                pitch,
                length,
                angular,
            )?;
            let qa = fields[0][1] + (fields[1][1] - fields[0][1]) * ta;
            let qb = fields[0][1] + (fields[1][1] - fields[0][1]) * tb;
            let theta_a = (length - qa) * tau / pitch;
            let theta_b = (length - qb) * tau / pitch;
            let width = theta_b - theta_a;
            if width.lo <= 0. || width.hi >= std::f64::consts::FRAC_PI_4 {
                return Err(Error::Unresolved);
            }
            let ra = fields[0][0] + (fields[1][0] - fields[0][0]) * ta;
            let rb = fields[0][0] + (fields[1][0] - fields[0][0]) * tb;
            let inner_chord = (self.radius * width.square()
                + p(2.) * magnitude((ra - rb) * self.radial) * width)
                / p(8.);
            let outer_chord = self.radius * width.square() / p(8.);
            let ei = inner_chord + p(a.inner_error.max(b.inner_error));
            let eo = outer_chord + p(a.outer_error.max(b.outer_error));
            // Radial transport of the finite curve/chord error. The global
            // NPT slope |dr/dtheta|<=pitch/pi holds across all source junctions.
            if root_min <= 0. || ei.hi >= root_min / 4. {
                return Err(Error::Unresolved);
            }
            let transport = p(1.)
                + p(4.) * pitch * self.radial / (p(std::f64::consts::PI.next_down()) * p(root_min));
            let new_chord = ei * transport + eo;
            if !new_chord.valid() || new_chord.hi > quality.deflection {
                return Err(Error::Unresolved);
            }
            chord = chord.max(new_chord.hi);
            let margin = p(2.) * p(a.outer_error.max(b.outer_error)) / self.radius;
            let theta = I::new(theta_a.lo - margin.hi, theta_b.hi + margin.hi);
            if edge.band != StockClipBand::Crest {
                normal_error = normal_error.max(push_triangle(
                    &mut mesh,
                    [edge.a, edge.b, b.outer],
                    6,
                    theta,
                    allowed,
                    &mut meter,
                )?);
                normal_error = normal_error.max(push_triangle(
                    &mut mesh,
                    [edge.a, b.outer, a.outer],
                    6,
                    theta,
                    allowed,
                    &mut meter,
                )?);
            }
            for ids in [[a.outer, b.outer, b.top], [a.outer, b.top, a.top]] {
                normal_error = normal_error.max(push_triangle(
                    &mut mesh, ids, 7, theta, allowed, &mut meter,
                )?);
            }
            normal_error = normal_error.max(push_triangle(
                &mut mesh,
                [top_center, a.top, b.top],
                8,
                theta,
                allowed,
                &mut meter,
            )?);
        }
        meter.charge(
            mesh.indices.len() * (128 + 3 * (mesh.indices.len().max(1).ilog2() as usize + 1)),
        )?;
        let proof = manufacturing::certify(&mesh).map_err(|_| Error::Unresolved)?;
        connected(&mesh, &mut meter)?;
        meter.charge(1)?;
        Ok(MeshCandidate {
            source: *self,
            mesh,
            proof,
            chord,
            normal: normal_error,
            work: meter.used,
            bytes: meter.peak,
        })
    }
}
