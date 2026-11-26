//! Imported hex-stock pipe hardware: one certified radial thread (head) or
//! through cavity (nut) in the chart of a regular hex prism. No authored law
//! is inferred from the file; the stock is earned from its twelve decoded
//! corners, each within an eighth of the linear tolerance of the regular
//! prism the chart describes. Re-export writes the same retained tables.
use super::*;
use crate::geometry::internal_thread_treatment_body::interval::Interval as I;
use crate::rust_prism_topology::position;
use multicad_core::references::{
    AnalyticGeometry, AnchoredFrame, Origin, ReferenceFrame, TopologyElement, TopologyKey,
    TopologyKind,
};
use multicad_rust_kernel::{
    axial_pipe_thread::authored::hardware::{
        Kind,
        mesh::{self as joined, ThreadBoundaryMesh},
    },
    radial_tensor::Surface,
    step_read::{HexStock, HexStockKind, ImportedPipeHexHardware},
};
use std::collections::BTreeSet;
use std::f64::consts::PI;
#[cfg(test)]
mod tests;

const OPERATION: &str = "Rust imported pipe hex hardware";
type Point = [f64; 3];

fn failure(kind: GeometryErrorKind, message: impl Into<String>) -> GeometryError {
    error(OPERATION, kind, message)
}
fn unsupported(message: impl Into<String>) -> GeometryError {
    failure(GeometryErrorKind::Unsupported, message)
}
fn dot(a: Point, b: Point) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross(a: Point, b: Point) -> Point {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn sub(a: Point, b: Point) -> Point {
    std::array::from_fn(|i| a[i] - b[i])
}
fn length(a: Point) -> f64 {
    dot(a, a).sqrt()
}
fn unit(a: Point) -> Option<Point> {
    let l = length(a);
    (l.is_finite() && l > 0.).then(|| a.map(|v| v / l))
}
fn centroid(ring: &[Point; 6]) -> Point {
    std::array::from_fn(|k| ring.iter().map(|p| p[k]).sum::<f64>() / 6.)
}

pub(crate) struct PipeHexBody {
    pub(super) feature: EntityId,
    /// Certified thread or cavity material at the identity occurrence of the chart.
    thread: Solid,
    kind: HexStockKind,
    across_flats_mm: f64,
    /// Local axial extent of the hex prism.
    span_mm: [f64; 2],
    /// Thread faces that became stock planes: nut [lower, upper], head [shoulder].
    caps: Vec<usize>,
    /// Number of retained thread faces; stock faces follow.
    face_count: usize,
}

impl Solid {
    /// Independent STEP readback of hex-stock pipe hardware: the radial part
    /// re-earns its material certificate in the file's own world chart, the
    /// hex prism is earned from its corners, and both move into the chart.
    pub(super) fn import_pipe_hex(
        imported: ImportedPipeHexHardware,
        feature: EntityId,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        let world = Self::qualify_imported_radial(imported.boundary, policy)?;
        Self::import_qualified_pipe_hex(world, imported.kind, imported.rings, feature, policy)
    }
    /// The same import from an already certified world boundary, as a
    /// compound member is after its own qualification.
    pub(super) fn import_qualified_pipe_hex(
        world: multicad_rust_kernel::radial_tensor::RadialBoundary,
        kind: HexStockKind,
        rings: [multicad_rust_kernel::step_read::HexRing; 2],
        feature: EntityId,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        let imported_kind = kind;
        let tolerance = policy.linear_tolerance_mm;
        let axes = world.cap_axes();
        let o = axes[0].origin;
        let [a, b] = rings;
        let nearer = |r: &[Point; 6]| length(sub(centroid(r), o));
        let lower_first = match imported_kind {
            HexStockKind::Head => match (a.annulus, b.annulus) {
                (true, false) => true,
                (false, true) => false,
                _ => return Err(unsupported("A head needs one shoulder and one top")),
            },
            HexStockKind::Nut => {
                if !(a.annulus && b.annulus) {
                    return Err(unsupported("A nut needs two annular caps"));
                }
                nearer(&a.corners) <= nearer(&b.corners)
            }
        };
        let (lower, upper) = if lower_first {
            (a.corners, b.corners)
        } else {
            (b.corners, a.corners)
        };
        // The chart takes the thread's own exact cap axis wherever the stock
        // agrees with it, so identity and cardinal files keep exact charts.
        let rise = sub(centroid(&upper), centroid(&lower));
        let n = if dot(axes[0].normal, rise) >= 0. {
            axes[0].normal
        } else {
            axes[0].normal.map(|v| -v)
        };
        if unit(rise).is_none_or(|r| length(cross(n, r)) > 1e-9) {
            return Err(unsupported("The hex stock is not coaxial with the thread"));
        }
        let perp = |p: Point| {
            let d = sub(p, o);
            let h = dot(d, n);
            std::array::from_fn::<f64, 3, _>(|i| d[i] - h * n[i])
        };
        // The lower ring is counterclockwise about -n: reverse it, start at
        // the corner nearest the thread reference, then align the upper ring.
        let lower: [Point; 6] = std::array::from_fn(|k| lower[(6 - k) % 6]);
        let toward =
            |p: Point| unit(perp(p)).map_or(f64::NEG_INFINITY, |d| dot(d, axes[0].reference));
        let start = (0..6)
            .max_by(|&i, &j| toward(lower[i]).total_cmp(&toward(lower[j])))
            .expect("six corners");
        let lower: [Point; 6] = std::array::from_fn(|k| lower[(k + start) % 6]);
        let shift = (0..6)
            .min_by(|&i, &j| {
                length(sub(perp(upper[i]), perp(lower[0])))
                    .total_cmp(&length(sub(perp(upper[j]), perp(lower[0]))))
            })
            .expect("six corners");
        let upper: [Point; 6] = std::array::from_fn(|k| upper[(k + shift) % 6]);
        let corner =
            unit(perp(lower[0])).ok_or_else(|| unsupported("Hex corner lies on the axis"))?;
        let reference = if length(cross(corner, axes[0].reference)) <= 1e-9 {
            axes[0].reference
        } else {
            corner
        };
        let radius = lower
            .iter()
            .chain(&upper)
            .map(|&p| length(perp(p)))
            .sum::<f64>()
            / 12.;
        let across_flats_mm = radius * 3_f64.sqrt();
        let span_mm = [
            dot(sub(centroid(&lower), o), n),
            dot(sub(centroid(&upper), o), n),
        ];
        let frame = PlaneFrame {
            origin: o,
            u: reference,
            v: cross(n, reference),
        };
        frame.checked(policy)?;
        let stock = HexStock {
            kind: imported_kind,
            origin: o,
            normal: n,
            reference,
            across_flats_mm,
            span_mm,
        };
        let rings = stock
            .corners()
            .ok_or_else(|| unsupported("The decoded hex stock is not a regular prism"))?;
        for (expected, actual) in rings.iter().zip([lower, upper]) {
            for k in 0..6 {
                if length(sub(expected[k], actual[k])) > tolerance / 8. {
                    return Err(unsupported(
                        "The decoded hex corners do not form a regular prism",
                    ));
                }
            }
        }
        let local = if frame == PlaneFrame::default() {
            world
        } else {
            let rows = [frame.u, frame.v, frame.normal()];
            let column = |k: usize| rows.map(|r| r[k]);
            let inverse_origin = std::array::from_fn(|i| -dot(rows[i], o));
            world
                .transformed_with_interrupt(
                    [inverse_origin, column(0), column(1)],
                    tolerance,
                    policy.minimum_feature_mm,
                    None,
                )
                .map_err(|cause| {
                    unsupported(format!(
                        "Imported thread does not transfer into the stock chart: {cause:?}"
                    ))
                })?
        };
        let mut caps = vec![
            usize::MAX;
            if imported_kind == HexStockKind::Nut {
                2
            } else {
                1
            }
        ];
        for (index, face) in local.input().faces.iter().enumerate() {
            if let Surface::Plane(axis) = local.input().surfaces[face.surface] {
                let z = axis.origin[2];
                for (slot, level) in caps.iter_mut().zip(span_mm) {
                    if (z - level).abs() <= tolerance {
                        if *slot != usize::MAX {
                            return Err(unsupported("Two thread caps meet one stock plane"));
                        }
                        *slot = index;
                    }
                }
            }
        }
        if caps.contains(&usize::MAX) {
            return Err(unsupported("The thread caps do not meet the stock planes"));
        }
        // The complete thread is strictly inside the hex flats: a
        // circumscribed 32-gon of its actual supports.
        let mut support = 0f64;
        for k in 0..32 {
            let angle = std::f64::consts::TAU * k as f64 / 32.;
            let [_, hi] = local
                .support_bounds(
                    [angle.portable_cos(), angle.portable_sin(), 0.],
                    tolerance,
                    None,
                )
                .map_err(|cause| {
                    unsupported(format!("Thread support is not resolved: {cause:?}"))
                })?;
            support = support.max(hi);
        }
        let reach = (support / (PI / 32.).portable_cos() * (1. + 1e-12)).next_up();
        let wall = across_flats_mm / 2. - reach;
        if wall.is_nan() || wall <= policy.minimum_feature_mm {
            return Err(failure(
                GeometryErrorKind::Input,
                "The hex stock does not clear the complete thread",
            ));
        }
        let face_count = local.input().faces.len();
        let thread = Self::from_rust_radial_boundary(local, feature, policy)?;
        let body = PipeHexBody {
            feature,
            thread,
            kind: imported_kind,
            across_flats_mm,
            span_mm,
            caps,
            face_count,
        };
        body.volume_interval()?;
        let bounds = body.bounds(frame)?;
        policy.point(bounds.min)?;
        policy.point(bounds.max)?;
        Ok(Self {
            body: SolidBody::PipeHex(Rc::new(body), frame),
            policy,
            provenance: vec![],
            rust_origins: vec![],
            retained_placement: None,
        })
    }
    /// Fresh exact readback of ACTUAL bytes (a standalone file or one flat
    /// member's local view) against this imported body's retained thread
    /// tables and earned stock. Returns the body's volume and world bounds.
    pub fn verify_pipe_hex_step_source(
        &self,
        text: &str,
        tolerance: f64,
    ) -> Option<Result<([f64; 2], Bounds), GeometryError>> {
        let (body, frame) = self.pipe_hex_body()?;
        Some((|| {
            if !tolerance.is_finite() || tolerance <= 0. {
                return Err(failure(
                    GeometryErrorKind::Input,
                    "Positive finite STEP tolerance required",
                ));
            }
            let tol = tolerance.min(self.policy.linear_tolerance_mm);
            let boundary = body.radial().placed_boundary(frame)?;
            multicad_rust_kernel::step::verify_pipe_hex_hardware(
                text,
                &boundary,
                body.stock(frame),
                tol,
            )
            .map_err(|e| {
                failure(
                    if e == multicad_rust_kernel::step::StepError::Limit {
                        GeometryErrorKind::Limit
                    } else {
                        GeometryErrorKind::InvalidResult
                    },
                    format!("Imported pipe hex STEP readback failed: {e:?}"),
                )
            })?;
            let v = body.volume_interval()?;
            Ok(([v.lo, v.hi], body.bounds(frame)?))
        })())
    }
    pub(crate) fn pipe_hex_body(&self) -> Option<(&PipeHexBody, PlaneFrame)> {
        match &self.body {
            SolidBody::PipeHex(body, frame) => Some((body, *frame)),
            _ => None,
        }
    }
    /// Imported hex-stock pipe hardware: its kind, across-flats width and
    /// local hex span. Imports never carry an inferred authored law.
    pub fn rust_imported_pipe_hex_data(&self) -> Option<(HexStockKind, f64, [f64; 2])> {
        self.pipe_hex_body()
            .map(|(b, _)| (b.kind, b.across_flats_mm, b.span_mm))
    }
}

impl PipeHexBody {
    fn radial(&self) -> &super::radial_boundary_body::RadialBoundaryBody {
        self.thread.radial_boundary_body().expect("sealed thread").0
    }
    fn is_nut(&self) -> bool {
        self.kind == HexStockKind::Nut
    }
    fn prefix(&self) -> &'static str {
        if self.is_nut() { "nut" } else { "head" }
    }
    fn origin(&self, role: impl Into<String>) -> Origin {
        Origin {
            source: self.feature,
            role: role.into(),
        }
    }
    pub(super) fn extent(&self) -> [f64; 3] {
        let r = self.across_flats_mm / 3_f64.sqrt();
        [r, r, self.span_mm[1].abs().max(self.span_mm[0].abs())]
    }
    fn corner(&self, k: usize, z: f64) -> Point {
        let half = self.across_flats_mm / 2.;
        let h = half / 3_f64.sqrt();
        let polygon = [
            [2. * h, 0.],
            [h, half],
            [-h, half],
            [-2. * h, 0.],
            [-h, -half],
            [h, -half],
        ];
        let [x, y] = polygon[k % 6];
        [x, y, z]
    }
    fn hex_area(&self) -> I {
        I::point(3.).sqrt() / I::point(2.) * I::point(self.across_flats_mm).square()
    }
    fn volume_interval(&self) -> Result<I, GeometryError> {
        let [lo, hi] = self
            .radial()
            .volume_properties(PlaneFrame::default(), None)?
            .volume_mm3;
        let thread = I::new(lo, hi);
        let prism = self.hex_area() * (I::point(self.span_mm[1]) - I::point(self.span_mm[0]));
        let v = if self.is_nut() {
            prism - thread
        } else {
            prism + thread
        };
        if !v.valid() || v.lo <= 0. {
            return Err(unsupported("The pipe hex hardware volume is not resolved"));
        }
        Ok(v)
    }
    pub(super) fn volume_mm3(&self) -> Result<f64, GeometryError> {
        let v = self.volume_interval()?;
        Ok(v.lo + (v.hi - v.lo) * 0.5)
    }
    pub(super) fn bounds(&self, frame: PlaneFrame) -> Result<Bounds, GeometryError> {
        let mut min = [f64::INFINITY; 3];
        let mut max = [f64::NEG_INFINITY; 3];
        for z in self.span_mm {
            for k in 0..6 {
                let p = position(frame, self.corner(k, z));
                for i in 0..3 {
                    min[i] = min[i].min(p[i].next_down());
                    max[i] = max[i].max(p[i].next_up());
                }
            }
        }
        let thread = self.radial().bounds(frame, None)?;
        for i in 0..3 {
            min[i] = min[i].min(thread.min[i]);
            max[i] = max[i].max(thread.max[i]);
        }
        Ok(Bounds { min, max })
    }
    fn thread_face_count(&self) -> usize {
        self.face_count
    }
    fn direction(frame: PlaneFrame, d: Point) -> Point {
        let n = frame.normal();
        std::array::from_fn(|i| frame.u[i] * d[0] + frame.v[i] * d[1] + n[i] * d[2])
    }
    fn wall_frame(&self, frame: PlaneFrame, k: usize) -> ReferenceFrame {
        let a = self.corner(k, 0.);
        let b = self.corner(k + 1, 0.);
        let edge = sub(b, a);
        let l = length(edge);
        let mid = (self.span_mm[0] + self.span_mm[1]) / 2.;
        ReferenceFrame {
            origin_mm: position(frame, [(a[0] + b[0]) / 2., (a[1] + b[1]) / 2., mid]),
            u: Self::direction(frame, edge.map(|v| v / l)),
            v: Self::direction(frame, [0., 0., 1.]),
        }
    }
    fn plane_frame(frame: PlaneFrame, z: f64, up: bool) -> ReferenceFrame {
        ReferenceFrame {
            origin_mm: position(frame, [0., 0., z]),
            u: frame.u,
            v: frame.v.map(|v| if up { v } else { -v }),
        }
    }
    /// Thread faces keep their indices (stock caps become hex annuli or the
    /// shoulder); six walls follow, then a head's top.
    pub(super) fn faces(&self, frame: PlaneFrame) -> Result<Vec<Face>, GeometryError> {
        let properties = self.radial().face_properties(frame, None)?;
        let mut faces = self.radial().faces(frame)?;
        let area = self.hex_area();
        for (level, &cap) in self.caps.iter().enumerate() {
            let p = properties
                .iter()
                .find(|f| f.face_index == cap)
                .ok_or_else(|| unsupported("Missing thread cap measure"))?;
            let hole = I::new(p.area_mm2[0], p.area_mm2[1]);
            let center = p.centroid_mm.map(|[a, b]| I::new(a, b));
            let remaining = area - hole;
            let z = self.span_mm[level];
            let axis_point = position(frame, [0., 0., z]).map(I::point);
            let c: [I; 3] =
                std::array::from_fn(|k| (area * axis_point[k] - hole * center[k]) / remaining);
            if !remaining.valid() || remaining.lo <= 0. || c.iter().any(|v| !v.valid()) {
                return Err(unsupported(
                    "A pipe hex stock plane measure is not resolved",
                ));
            }
            let up = self.is_nut() && level == 1;
            let plane = Self::plane_frame(frame, z, up);
            let slot = faces
                .iter_mut()
                .find(|f| f.index == cap as u32)
                .ok_or_else(|| unsupported("Missing thread cap face"))?;
            *slot = Face {
                index: cap as u32,
                plane: Some((plane.origin_mm, plane.normal())),
                center_mm: c.map(|v| v.lo + (v.hi - v.lo) * 0.5),
                area_mm2: remaining.lo + (remaining.hi - remaining.lo) * 0.5,
            };
        }
        let base = self.thread_face_count() as u32;
        let side = self.across_flats_mm / 3_f64.sqrt();
        let height = self.span_mm[1] - self.span_mm[0];
        for k in 0..6 {
            let wall = self.wall_frame(frame, k);
            faces.push(Face {
                index: base + k as u32,
                plane: Some((wall.origin_mm, wall.normal())),
                center_mm: wall.origin_mm,
                area_mm2: side * height,
            });
        }
        if !self.is_nut() {
            let top = Self::plane_frame(frame, self.span_mm[1], true);
            faces.push(Face {
                index: base + 6,
                plane: Some((top.origin_mm, top.normal())),
                center_mm: top.origin_mm,
                area_mm2: area.lo + (area.hi - area.lo) * 0.5,
            });
        }
        Ok(faces)
    }
    pub(super) fn mesh(
        &self,
        frame: PlaneFrame,
        options: multicad_rust_kernel::display::Tessellation,
    ) -> Result<Mesh, GeometryError> {
        let mesh_failure = |what: &str, cause: String| {
            failure(
                if cause.contains("Limit") {
                    GeometryErrorKind::Limit
                } else if cause.contains("Cancelled") {
                    GeometryErrorKind::Cancelled
                } else {
                    GeometryErrorKind::Unsupported
                },
                format!("{what}: {cause}"),
            )
        };
        // The thread is traced at its world occurrence and joined to the
        // stock in that chart, so every pose meshes the same retained source.
        let boundary = self.radial().placed_boundary(frame)?;
        let trace = boundary
            .tessellate_with_cap_trace(options, Default::default(), None)
            .map_err(|e| mesh_failure("Thread trace", format!("{e:?}")))?;
        let kind = match self.kind {
            HexStockKind::Nut => Kind::HexNut {
                across_flats_mm: self.across_flats_mm,
            },
            HexStockKind::Head => Kind::HexHead {
                across_flats_mm: self.across_flats_mm,
                height_mm: self.span_mm[1] - self.span_mm[0],
            },
        };
        let joined = joined::join_in_chart(
            ThreadBoundaryMesh::Retained(&trace),
            kind,
            [frame.origin, frame.u, frame.v, frame.normal()],
            joined::MeshOptions {
                tessellation: options,
                maximum_work: 512_000_000,
                maximum_bytes: 64 * 1024 * 1024,
                additional_resident_bytes: boundary.retained_budget_bytes(),
            },
            None,
        )
        .map_err(|e| mesh_failure("Hex stock join", format!("{e:?}")))?;
        let mut mesh = joined.into_mesh();
        // The join numbers walls from the thread face count, then the lower
        // and upper stock planes; those planes keep the caps' own indices.
        let base = self.thread_face_count() as u32;
        let lower = self.caps[0] as u32;
        let upper = if self.is_nut() {
            self.caps[1] as u32
        } else {
            base + 6
        };
        for face in &mut mesh.face_indices {
            if *face == base + 6 {
                *face = lower;
            } else if *face == base + 7 {
                *face = upper;
            }
        }
        crate::rust_display::placed_mesh(Ok(mesh), PlaneFrame::default(), false)
    }
    pub(crate) fn topology(
        &self,
        frame: PlaneFrame,
        owner: EntityId,
        revision: u64,
    ) -> Result<TopologySnapshot, GeometryError> {
        let mut result = self.radial().topology(frame, owner, revision)?;
        if self.is_nut() {
            let inward = |origin: &Origin| Origin {
                source: origin.source,
                role: if origin.role.starts_with("thread/") || origin.role == "thread.root" {
                    format!("receiving.inward/{}", origin.role)
                } else {
                    origin.role.clone()
                },
            };
            for element in &mut result.elements {
                element.origins = element.origins.iter().map(inward).collect();
                element.boundary = element.boundary.iter().map(inward).collect();
            }
        }
        let prefix = self.prefix();
        let wall = |k: usize| self.origin(format!("{prefix}.wall/{k}"));
        let rim = |level: usize, k: usize| self.origin(format!("{prefix}.rim/{level}/{k}"));
        let corner_edge = |k: usize| self.origin(format!("{prefix}.corner/{k}"));
        let corner_vertex =
            |level: usize, k: usize| self.origin(format!("{prefix}.vertex/{level}/{k}"));
        let ordinal = |kind: TopologyKind, elements: &[TopologyElement]| {
            elements
                .iter()
                .filter(|e| e.key.kind == kind)
                .map(|e| e.key.ordinal + 1)
                .max()
                .unwrap_or(0)
        };
        let first_edge = ordinal(TopologyKind::Edge, &result.elements);
        let first_vertex = ordinal(TopologyKind::Vertex, &result.elements);
        let base = self.thread_face_count() as u32;
        let mut level_face: [Option<Origin>; 2] = [None, None];
        for (level, &cap) in self.caps.iter().enumerate() {
            let element = result
                .elements
                .iter_mut()
                .find(|e| e.key.kind == TopologyKind::Face && e.key.ordinal == cap as u32)
                .ok_or_else(|| unsupported("Missing thread cap topology"))?;
            element.boundary.extend((0..6).map(|k| rim(level, k)));
            let up = self.is_nut() && level == 1;
            let plane = Self::plane_frame(frame, self.span_mm[level], up);
            element.geometry = Some(AnalyticGeometry::Plane { frame: plane });
            for anchor in &mut element.frames {
                anchor.frame = plane;
            }
            level_face[level] = element.origins.iter().next().cloned();
        }
        if !self.is_nut() {
            level_face[1] = Some(self.origin("head.top"));
        }
        let [Some(lower_face), Some(upper_face)] = level_face else {
            return Err(unsupported("The pipe hex stock planes have no topology"));
        };
        let level_face = [lower_face, upper_face];
        let mut faces = Vec::new();
        for k in 0..6 {
            let plane = self.wall_frame(frame, k);
            faces.push(TopologyElement {
                key: TopologyKey {
                    kind: TopologyKind::Face,
                    ordinal: base + k as u32,
                },
                origins: BTreeSet::from([wall(k)]),
                boundary: BTreeSet::from([
                    rim(0, k),
                    rim(1, k),
                    corner_edge(k),
                    corner_edge((k + 1) % 6),
                ]),
                frames: vec![AnchoredFrame {
                    origin: wall(k),
                    frame: plane,
                }],
                geometry: Some(AnalyticGeometry::Plane { frame: plane }),
            });
        }
        if !self.is_nut() {
            let top = Self::plane_frame(frame, self.span_mm[1], true);
            faces.push(TopologyElement {
                key: TopologyKey {
                    kind: TopologyKind::Face,
                    ordinal: base + 6,
                },
                origins: BTreeSet::from([level_face[1].clone()]),
                boundary: (0..6).map(|k| rim(1, k)).collect(),
                frames: vec![AnchoredFrame {
                    origin: level_face[1].clone(),
                    frame: top,
                }],
                geometry: Some(AnalyticGeometry::Plane { frame: top }),
            });
        }
        let mut edges = Vec::new();
        for (level, face) in level_face.iter().enumerate() {
            let z = self.span_mm[level];
            for k in 0..6 {
                edges.push(TopologyElement {
                    key: TopologyKey {
                        kind: TopologyKind::Edge,
                        ordinal: first_edge + (6 * level + k) as u32,
                    },
                    origins: BTreeSet::from([rim(level, k)]),
                    boundary: BTreeSet::from([face.clone(), wall(k)]),
                    frames: vec![],
                    geometry: Some(AnalyticGeometry::Line {
                        a_mm: position(frame, self.corner(k, z)),
                        b_mm: position(frame, self.corner(k + 1, z)),
                    }),
                });
            }
        }
        for k in 0..6 {
            edges.push(TopologyElement {
                key: TopologyKey {
                    kind: TopologyKind::Edge,
                    ordinal: first_edge + 12 + k as u32,
                },
                origins: BTreeSet::from([corner_edge(k)]),
                boundary: BTreeSet::from([wall((k + 5) % 6), wall(k)]),
                frames: vec![],
                geometry: Some(AnalyticGeometry::Line {
                    a_mm: position(frame, self.corner(k, self.span_mm[0])),
                    b_mm: position(frame, self.corner(k, self.span_mm[1])),
                }),
            });
        }
        let mut vertices = Vec::new();
        for (level, face) in level_face.iter().enumerate() {
            for k in 0..6 {
                vertices.push(TopologyElement {
                    key: TopologyKey {
                        kind: TopologyKind::Vertex,
                        ordinal: first_vertex + (6 * level + k) as u32,
                    },
                    origins: BTreeSet::from([corner_vertex(level, k)]),
                    boundary: BTreeSet::from([face.clone(), wall((k + 5) % 6), wall(k)]),
                    frames: vec![],
                    geometry: Some(AnalyticGeometry::Point {
                        position_mm: position(frame, self.corner(k, self.span_mm[level])),
                    }),
                });
            }
        }
        result.elements.extend(faces);
        result.elements.extend(edges);
        result.elements.extend(vertices);
        result
            .validate()
            .map_err(|e| failure(GeometryErrorKind::InvalidResult, e))?;
        Ok(result)
    }
    pub(super) fn display_edges(
        &self,
        frame: PlaneFrame,
        deflection: f64,
    ) -> Result<Vec<DisplayEdge>, GeometryError> {
        let mut hex = self.topology(frame, self.feature, 0)?;
        hex.elements.retain(|e| {
            e.key.kind == TopologyKind::Edge
                && matches!(e.geometry, Some(AnalyticGeometry::Line { .. }))
        });
        crate::rust_display::edges(None, &hex, deflection)
    }
    fn stock(&self, frame: PlaneFrame) -> HexStock {
        HexStock {
            kind: self.kind,
            origin: frame.origin,
            normal: frame.normal(),
            reference: frame.u,
            across_flats_mm: self.across_flats_mm,
            span_mm: self.span_mm,
        }
    }
    pub(super) fn supports_step(&self, frame: PlaneFrame, tolerance: f64) -> bool {
        self.radial().placed_boundary(frame).is_ok_and(|boundary| {
            multicad_rust_kernel::step::write_pipe_hex_hardware(
                &boundary,
                self.stock(frame),
                tolerance,
            )
            .is_ok()
        })
    }
    pub(super) fn append_step(
        &self,
        writer: &mut multicad_rust_kernel::step::Writer,
        frame: PlaneFrame,
        tolerance: f64,
    ) -> Result<(), multicad_rust_kernel::step::StepError> {
        let boundary = self.radial().placed_boundary(frame).map_err(|cause| {
            if cause.kind == GeometryErrorKind::Limit {
                multicad_rust_kernel::step::StepError::Limit
            } else {
                multicad_rust_kernel::step::StepError::Unsupported
            }
        })?;
        writer.append_pipe_hex_hardware(&boundary, self.stock(frame), tolerance)
    }
}
