//! Double-rounded hex nuts and fully threaded hex-head screws. The literal
//! rounded thread keeps its five bands; the hex stock adds six planar walls,
//! annular caps or a seat, and a head top. Neither form claims STEP yet.
use super::rounded_thread_body::RoundedThreadBody;
use super::*;
use crate::rust_prism_topology::position;
use multicad_core::references::{
    AnalyticGeometry, AnchoredFrame, Origin, ReferenceFrame, TopologyElement, TopologyKey,
    TopologyKind,
};
use multicad_rust_kernel::{
    rounded_thread::{HexForm, HexStock, Spec},
    rounded_thread_mesh,
};
use std::{collections::BTreeSet, f64::consts::FRAC_PI_3};

type Point = [f64; 3];
const OPERATION: &str = "Rust rounded hex hardware";

pub(crate) struct RoundedHexBody {
    pub(crate) geometry: HexStock,
    pub(super) feature: EntityId,
    inner: Rc<RoundedThreadBody>,
}

fn failure(kind: GeometryErrorKind, message: impl Into<String>) -> GeometryError {
    error(OPERATION, kind, message)
}

impl Solid {
    /// A hex nut whose through bore is the literal rounded cutter.
    pub fn rust_rounded_hex_nut(
        spec: Spec,
        across_flats_mm: f64,
        minimum_wall_mm: f64,
        feature: EntityId,
        frame: PlaneFrame,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        if !minimum_wall_mm.is_finite()
            || minimum_wall_mm < 0.
            || across_flats_mm / 2. - spec.major_radius_mm
                < minimum_wall_mm.max(policy.minimum_feature_mm)
        {
            return Err(failure(
                GeometryErrorKind::Input,
                "Hex nut flats leave insufficient wall around the rounded thread",
            ));
        }
        Self::rounded_hex(spec, across_flats_mm, HexForm::Nut, feature, frame, policy)
    }

    /// A fully threaded screw: thread `0 <= z <= L` in the thread frame, hex
    /// head `-height <= z <= 0` seated at the thread start.
    pub fn rust_rounded_hex_head_screw(
        spec: Spec,
        across_flats_mm: f64,
        head_height_mm: f64,
        feature: EntityId,
        frame: PlaneFrame,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        policy.length(head_height_mm)?;
        if across_flats_mm / 2. - spec.major_radius_mm < policy.minimum_feature_mm {
            return Err(failure(
                GeometryErrorKind::Input,
                "Hex head seat is below the native minimum feature around the rounded thread",
            ));
        }
        Self::rounded_hex(
            spec,
            across_flats_mm,
            HexForm::Head {
                height_mm: head_height_mm,
            },
            feature,
            frame,
            policy,
        )
    }

    fn rounded_hex(
        spec: Spec,
        across_flats_mm: f64,
        form: HexForm,
        feature: EntityId,
        frame: PlaneFrame,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        let inner = Self::rust_rounded_thread(spec, feature, frame, policy)?;
        policy.length(across_flats_mm)?;
        let SolidBody::RoundedThread(inner, _) = inner.body else {
            unreachable!()
        };
        let geometry = HexStock::new(inner.geometry, across_flats_mm, form).map_err(|e| {
            failure(
                GeometryErrorKind::Input,
                format!("Invalid rounded hex stock: {e:?}"),
            )
        })?;
        let corner = across_flats_mm / 3_f64.sqrt();
        frame.checked_metric_extent(
            [corner, corner, spec.length_mm + geometry.hex_height_mm()],
            policy,
        )?;
        let body = RoundedHexBody {
            geometry,
            feature,
            inner,
        };
        let bounds = body.bounds(frame)?;
        policy.point(bounds.min)?;
        policy.point(bounds.max)?;
        Ok(Self {
            body: SolidBody::RoundedHex(Rc::new(body), frame),
            policy,
            provenance: vec![],
            rust_origins: vec![],
            retained_placement: None,
        })
    }

    pub(crate) fn rounded_hex_body(&self) -> Option<(&RoundedHexBody, PlaneFrame)> {
        match &self.body {
            SolidBody::RoundedHex(body, frame) => Some((body, *frame)),
            _ => None,
        }
    }

    /// Exact source stock for typed consumers.
    pub fn rust_rounded_hex_data(&self) -> Option<(HexStock, PlaneFrame)> {
        self.rounded_hex_body()
            .map(|(body, frame)| (body.geometry, frame))
    }
}

fn cross(a: Point, b: Point) -> Point {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn dot(a: Point, b: Point) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

impl RoundedHexBody {
    fn origin(&self, role: impl Into<String>) -> Origin {
        Origin {
            source: self.feature,
            role: role.into(),
        }
    }
    fn is_nut(&self) -> bool {
        matches!(self.geometry.form(), HexForm::Nut)
    }
    fn prefix(&self) -> &'static str {
        if self.is_nut() { "nut" } else { "head" }
    }
    pub(super) fn transformed_frame(
        &self,
        frame: PlaneFrame,
        policy: GeometryPolicy,
    ) -> Result<PlaneFrame, GeometryError> {
        let corner = self.geometry.across_flats_mm() / 3_f64.sqrt();
        frame.checked_metric_extent(
            [
                corner,
                corner,
                self.geometry.length_mm() + self.geometry.hex_height_mm(),
            ],
            policy,
        )?;
        Ok(frame)
    }
    pub(super) fn volume_mm3(&self) -> f64 {
        self.geometry.volume_mm3()
    }
    fn hex_corner(&self, k: usize, z: f64) -> Point {
        let [x, y] = self.geometry.hex_vertices()[k % 6];
        [x, y, z]
    }
    pub(super) fn bounds(&self, frame: PlaneFrame) -> Result<Bounds, GeometryError> {
        let [low, high] = self.geometry.hex_interval_mm();
        let mut min = [f64::INFINITY; 3];
        let mut max = [f64::NEG_INFINITY; 3];
        for z in [low, high] {
            for k in 0..6 {
                let p = position(frame, self.hex_corner(k, z));
                for i in 0..3 {
                    min[i] = min[i].min(p[i].next_down());
                    max[i] = max[i].max(p[i].next_up());
                }
            }
        }
        if !self.is_nut() {
            let thread = self.inner.bounds(frame)?;
            for i in 0..3 {
                min[i] = min[i].min(thread.min[i]);
                max[i] = max[i].max(thread.max[i]);
            }
        }
        Ok(Bounds { min, max })
    }
    pub(super) fn contains(&self, frame: PlaneFrame, point: Point) -> bool {
        let n = frame.normal();
        let rows = [cross(frame.v, n), cross(n, frame.u), n];
        let det = dot(frame.u, rows[0]);
        let delta = std::array::from_fn(|i| point[i] - frame.origin[i]);
        self.geometry
            .contains(rows.map(|r| dot(r, delta) / det), self.inner.tolerance())
    }
    fn wall_frame(&self, frame: PlaneFrame, side: usize) -> ReferenceFrame {
        let angle = FRAC_PI_3 * (side as f64 + 1.);
        let [low, high] = self.geometry.hex_interval_mm();
        let (s, c) = angle.sin_cos();
        let half = self.geometry.across_flats_mm() / 2.;
        let n = frame.normal();
        // u along the flat (counterclockwise), v along the axis: u × v is the
        // outward wall normal.
        let local_u = [-s, c, 0.];
        let direction = |p: Point| -> Point {
            std::array::from_fn(|i| frame.u[i] * p[0] + frame.v[i] * p[1] + n[i] * p[2])
        };
        ReferenceFrame {
            origin_mm: position(frame, [half * c, half * s, 0.5 * (low + high)]),
            u: direction(local_u),
            v: direction([0., 0., 1.]),
        }
    }
    fn plane_frame(&self, frame: PlaneFrame, z: f64, up: bool) -> ReferenceFrame {
        ReferenceFrame {
            origin_mm: position(frame, [0., 0., z]),
            u: frame.u,
            v: frame.v.map(|v| if up { v } else { -v }),
        }
    }
    pub(super) fn faces(&self, frame: PlaneFrame) -> Result<Vec<Face>, GeometryError> {
        let mut faces = self.inner.faces(frame)?;
        let properties = self.inner.properties()?;
        let area = self.geometry.hex_area_mm2();
        let length = self.geometry.length_mm();
        // Annular caps (nut) or the seat (screw) keep the hex outside the
        // thread section; the screw tip stays the thread's own cap.
        let annular: &[usize] = if self.is_nut() { &[0, 1] } else { &[0] };
        for &i in annular {
            let (center, hole_area) = properties.faces[i];
            let remaining = area - hole_area;
            let up = if self.is_nut() { i == 1 } else { true };
            let plane = self.plane_frame(frame, i as f64 * length, up);
            faces[i] = Face {
                index: i as u32,
                plane: Some((plane.origin_mm, plane.normal())),
                center_mm: position(
                    frame,
                    [
                        -hole_area * center[0] / remaining,
                        -hole_area * center[1] / remaining,
                        i as f64 * length,
                    ],
                ),
                area_mm2: remaining,
            };
        }
        let height = self.geometry.hex_height_mm();
        let side = self.geometry.across_flats_mm() / 3_f64.sqrt();
        for k in 0..6 {
            let wall = self.wall_frame(frame, k);
            faces.push(Face {
                index: HexStock::wall_face(k),
                plane: Some((wall.origin_mm, wall.normal())),
                center_mm: wall.origin_mm,
                area_mm2: side * height,
            });
        }
        if let HexForm::Head { height_mm } = self.geometry.form() {
            let top = self.plane_frame(frame, -height_mm, false);
            faces.push(Face {
                index: 13,
                plane: Some((top.origin_mm, top.normal())),
                center_mm: top.origin_mm,
                area_mm2: area,
            });
        }
        Ok(faces)
    }
    fn tessellation_error(e: rounded_thread_mesh::Error) -> GeometryError {
        match e {
            rounded_thread_mesh::Error::Cancelled => failure(
                GeometryErrorKind::Cancelled,
                "Rounded hex hardware mesh superseded",
            ),
            rounded_thread_mesh::Error::Mesh(e) => failure(
                GeometryErrorKind::InvalidResult,
                format!("Rounded hex hardware mesh: {e:?}"),
            ),
        }
    }
    pub(super) fn mesh(
        &self,
        frame: PlaneFrame,
        options: multicad_rust_kernel::display::Tessellation,
    ) -> Result<Mesh, GeometryError> {
        let mesh = match rounded_thread_mesh::tessellate_hex(&self.geometry, options) {
            Ok(mesh) => Ok(mesh),
            Err(rounded_thread_mesh::Error::Mesh(e)) => Err(e),
            Err(e) => return Err(Self::tessellation_error(e)),
        };
        crate::rust_display::placed_mesh(mesh, frame, false)
    }
    pub(super) fn volume_error_bound(
        &self,
        options: multicad_rust_kernel::display::Tessellation,
    ) -> Result<f64, GeometryError> {
        rounded_thread_mesh::hex_volume_error_bound(&self.geometry, options)
            .map_err(|e| failure(GeometryErrorKind::Limit, format!("{e:?}")))
    }
    pub(crate) fn topology(
        &self,
        frame: PlaneFrame,
        owner: EntityId,
        revision: u64,
    ) -> Result<TopologySnapshot, GeometryError> {
        let mut result = self.inner.topology(frame, owner, revision)?;
        let prefix = self.prefix();
        if self.is_nut() {
            // The same meridian has the opposite material side in a nut.
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
        let [low, high] = self.geometry.hex_interval_mm();
        // Hex levels: index 0 is the lower plane, 1 the upper plane.
        let level_face: [Origin; 2] = if self.is_nut() {
            [self.origin("cap.start"), self.origin("cap.end")]
        } else {
            [self.origin("head.top"), self.origin("cap.start")]
        };
        if !self.is_nut() {
            result.elements[0].origins.insert(self.origin("head.seat"));
            let seat = self.plane_frame(frame, 0., true);
            result.elements[0].geometry = Some(AnalyticGeometry::Plane { frame: seat });
            for anchor in &mut result.elements[0].frames {
                anchor.frame = seat;
            }
        }
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
        let level_index = |level: usize| -> Option<usize> {
            // The upper hex level of a nut is cap 1, its lower cap 0; a head's
            // upper level is the seat (cap 0) and its lower level the new top.
            match (self.is_nut(), level) {
                (true, l) => Some(l),
                (false, 1) => Some(0),
                (false, _) => None,
            }
        };
        for level in 0..2 {
            if let Some(face) = level_index(level) {
                result.elements[face]
                    .boundary
                    .extend((0..6).map(|k| rim(level, k)));
            }
        }
        let mut faces = Vec::new();
        for k in 0..6 {
            let plane = self.wall_frame(frame, k);
            faces.push(TopologyElement {
                key: TopologyKey {
                    kind: TopologyKind::Face,
                    ordinal: HexStock::wall_face(k),
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
            let top = self.plane_frame(frame, low, false);
            faces.push(TopologyElement {
                key: TopologyKey {
                    kind: TopologyKind::Face,
                    ordinal: 13,
                },
                origins: BTreeSet::from([level_face[0].clone()]),
                boundary: (0..6).map(|k| rim(0, k)).collect(),
                frames: vec![AnchoredFrame {
                    origin: level_face[0].clone(),
                    frame: top,
                }],
                geometry: Some(AnalyticGeometry::Plane { frame: top }),
            });
        }
        let levels = [low, high];
        let mut edges = Vec::new();
        for level in 0..2 {
            for k in 0..6 {
                edges.push(TopologyElement {
                    key: TopologyKey {
                        kind: TopologyKind::Edge,
                        ordinal: first_edge + (6 * level + k) as u32,
                    },
                    origins: BTreeSet::from([rim(level, k)]),
                    boundary: BTreeSet::from([level_face[level].clone(), wall(k)]),
                    frames: vec![],
                    geometry: Some(AnalyticGeometry::Line {
                        a_mm: position(frame, self.hex_corner(k, levels[level])),
                        b_mm: position(frame, self.hex_corner(k + 1, levels[level])),
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
                // Corner k joins the walls on either side of it.
                boundary: BTreeSet::from([wall((k + 5) % 6), wall(k)]),
                frames: vec![],
                geometry: Some(AnalyticGeometry::Line {
                    a_mm: position(frame, self.hex_corner(k, low)),
                    b_mm: position(frame, self.hex_corner(k, high)),
                }),
            });
        }
        let mut vertices = Vec::new();
        for level in 0..2 {
            for k in 0..6 {
                vertices.push(TopologyElement {
                    key: TopologyKey {
                        kind: TopologyKind::Vertex,
                        ordinal: first_vertex + (6 * level + k) as u32,
                    },
                    origins: BTreeSet::from([corner_vertex(level, k)]),
                    boundary: BTreeSet::from([
                        level_face[level].clone(),
                        wall((k + 5) % 6),
                        wall(k),
                    ]),
                    frames: vec![],
                    geometry: Some(AnalyticGeometry::Point {
                        position_mm: position(frame, self.hex_corner(k, levels[level])),
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
    pub(in crate::geometry) fn display_edges(
        &self,
        frame: PlaneFrame,
        deflection: f64,
    ) -> Result<Vec<DisplayEdge>, GeometryError> {
        let mut edges = self.inner.display_edges(frame, deflection)?;
        let mut hex = self.topology(frame, self.feature, 0)?;
        hex.elements.retain(|e| {
            e.key.kind == TopologyKind::Edge
                && matches!(e.geometry, Some(AnalyticGeometry::Line { .. }))
        });
        edges.extend(crate::rust_display::edges(None, &hex, deflection)?);
        Ok(edges)
    }
}
#[cfg(test)]
mod tests;
