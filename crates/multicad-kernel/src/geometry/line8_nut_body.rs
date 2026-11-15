//! Actual retained prismatic nut stock and policy-qualified original tap phase.
//! Stock frame and bore offset stay separate source expressions.
use super::*;
use crate::HelicalThread;
use crate::geometry::internal_thread_treatment_body::interval::Interval as I;
use multicad_core::references::{
    AnalyticGeometry, AnchoredFrame, TopologyElement, TopologyKey, TopologyKind,
};
use multicad_rust_kernel::{
    authored_flat_thread::{AuthoredThread, SourceOptions},
    authored_internal_treatment::FacePropertyOptions,
    axial_pipe_thread::authored::StockClipPlacement,
    line8_nut::{
        Boundary, Classification, Edge as SourceEdge, FaceMeasure, Nut, Relief,
        Vertex as SourceVertex,
    },
};
use std::collections::BTreeSet;
type Point = [f64; 3];
mod mating;
mod mesh;
#[cfg(test)]
mod tests;
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Line8NutCut {
    pub profile: HelicalThread,
    pub frame: PlaneFrame,
    pub bore_from_start_mm: f64,
    pub relief: Relief,
    pub minimum_wall_mm: f64,
}
#[derive(Clone, Copy, Debug)]
pub struct Line8NutTransfer {
    pub legacy_position_error_mm: f64,
    pub legacy_normal_error_rad: f64,
    pub phase_position_error_mm: f64,
    pub phase_derivative_error_mm: f64,
    pub phase_normal_error_rad: f64,
    pub axis_alignment_mm: f64,
}
#[derive(Clone)]
pub(crate) struct Line8NutBody {
    pub(super) feature: EntityId,
    source: Nut,
    boundary: Boundary,
    stock: StockClipPlacement,
    outer: RigidTransform,
    transfer: Line8NutTransfer,
    volume: [f64; 2],
    measures: Vec<FaceMeasure>,
    bounds: Bounds,
    requested_wall_mm: f64,
}
fn fail(c: multicad_rust_kernel::line8_nut::Error) -> GeometryError {
    use multicad_rust_kernel::line8_nut::Error as E;
    error(
        "retained Line8 nut",
        match c {
            E::Invalid => GeometryErrorKind::Input,
            E::Limit => GeometryErrorKind::Limit,
            E::Cancelled => GeometryErrorKind::Cancelled,
            _ => GeometryErrorKind::Unsupported,
        },
        format!("Actual retained nut source proof failed: {c:?}"),
    )
}
fn unresolved() -> GeometryError {
    fail(multicad_rust_kernel::line8_nut::Error::Unresolved)
}
fn source_error(c: multicad_rust_kernel::authored_flat_thread::Error) -> GeometryError {
    use multicad_rust_kernel::authored_flat_thread::Error as E;
    error(
        "retained Line8 tap",
        match c {
            E::Invalid => GeometryErrorKind::Input,
            E::Limit => GeometryErrorKind::Limit,
            E::Cancelled => GeometryErrorKind::Cancelled,
            _ => GeometryErrorKind::Unsupported,
        },
        format!("Original authored tap proof failed: {c:?}"),
    )
}
fn placement_error(c: multicad_rust_kernel::axial_pipe_thread::Error) -> GeometryError {
    use multicad_rust_kernel::axial_pipe_thread::Error as E;
    error(
        "retained Line8 placement",
        match c {
            E::Invalid => GeometryErrorKind::Input,
            E::Limit => GeometryErrorKind::Limit,
            E::Cancelled => GeometryErrorKind::Cancelled,
            _ => GeometryErrorKind::Unsupported,
        },
        format!("Literal stock placement proof failed: {c:?}"),
    )
}
fn check(interrupt: Option<&dyn Fn() -> bool>) -> Result<(), GeometryError> {
    if interrupt.is_some_and(|f| f()) {
        Err(fail(multicad_rust_kernel::line8_nut::Error::Cancelled))
    } else {
        Ok(())
    }
}
fn iv(v: [f64; 2]) -> I {
    I::new(v[0], v[1])
}
fn norm(a: [I; 3]) -> I {
    let n = a.into_iter().fold(I::ZERO, |s, v| s + v.square());
    I::new(n.lo.max(0.), n.hi).sqrt()
}
fn unit(a: [I; 3]) -> Result<[I; 3], GeometryError> {
    let n = norm(a);
    if n.lo <= 0. || !n.hi.is_finite() {
        return Err(unresolved());
    }
    Ok(a.map(|v| v / n))
}
fn cofactor(m: [[I; 3]; 3]) -> [[I; 3]; 3] {
    let pairs = [(1, 2), (2, 0), (0, 1)];
    std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            let (a, b) = pairs[i];
            let (c, d) = pairs[j];
            m[a][c] * m[b][d] - m[a][d] * m[b][c]
        })
    })
}
fn placement(
    matrix: [[f64; 3]; 3],
    translation: Point,
) -> Result<StockClipPlacement, GeometryError> {
    StockClipPlacement::new(matrix, translation).map_err(placement_error)
}
fn matrix_norm(p: StockClipPlacement) -> I {
    (I::ONE + I::point(p.metric_error())).sqrt()
}
fn condition(p: StockClipPlacement) -> I {
    (I::ONE + I::point(p.metric_error())) / (I::ONE - I::point(p.metric_error()))
}
fn qualified(t: Line8NutTransfer, p: GeometryPolicy) -> bool {
    (I::point(t.legacy_position_error_mm) + I::point(t.phase_position_error_mm)).hi
        <= p.linear_tolerance_mm
        && (I::point(t.legacy_normal_error_rad) + I::point(t.phase_normal_error_rad)).hi
            <= p.linear_tolerance_mm
        && t.phase_derivative_error_mm <= p.linear_tolerance_mm
}
impl Solid {
    pub fn rust_cut_line8_nut(
        &self,
        cut: Line8NutCut,
        feature: EntityId,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, GeometryError> {
        let Line8NutCut {
            profile,
            frame: requested,
            bore_from_start_mm,
            relief,
            minimum_wall_mm,
        } = cut;
        policy.validate()?;
        requested.checked(policy)?;
        for value in [
            profile.major_radius_mm,
            profile.pitch_mm,
            profile.length_mm,
            profile.depth_mm,
            profile.crest_width_mm,
            minimum_wall_mm,
        ] {
            policy.length(value)?;
        }

        check(interrupt)?;
        let (stock_body, frame) = self
            .linear_body()
            .ok_or_else(|| fail(multicad_rust_kernel::line8_nut::Error::Unsupported))?;
        let sweep = stock_body
            .straight_source()
            .ok_or_else(|| fail(multicad_rust_kernel::line8_nut::Error::Unsupported))?;
        if self.policy != policy
            || policy.boolean_fuzzy_mm != 0.
            || feature.is_nil()
            || feature == stock_body.feature
            || profile.radial_taper != 0.
            || profile.starts != 1
        {
            return Err(fail(multicad_rust_kernel::line8_nut::Error::Unsupported));
        }
        let spec = multicad_rust_kernel::helical_thread::Spec {
            major_radius_mm: profile.major_radius_mm,
            pitch_mm: profile.pitch_mm,
            length_mm: profile.length_mm,
            depth_mm: profile.depth_mm,
            crest_width_mm: profile.crest_width_mm,
            lower_flank_deg: profile.lower_flank_deg,
            upper_flank_deg: profile.upper_flank_deg,
            starts: profile.starts,
            left_hand: profile.left_hand,
            root_radius_mm: None,
        };
        let thread = AuthoredThread::new_with_interrupt(
            spec,
            SourceOptions {
                position_limit_mm: policy.linear_tolerance_mm / 8.,
                normal_angle_limit_rad: policy.linear_tolerance_mm / 8.,
            },
            interrupt,
        )
        .map_err(source_error)?;
        let axes = [frame.normal(), frame.u, frame.v];
        let matrix = std::array::from_fn(|i| axes.map(|v| v[i]));
        let stock = placement(matrix, frame.origin)?;
        let vectors = placement(matrix, [0.; 3])?;
        let u = vectors
            .inverse_point_bounds(requested.u.map(|v| [v, v]))
            .map_err(placement_error)?;
        let v = vectors
            .inverse_point_bounds(requested.v.map(|v| [v, v]))
            .map_err(placement_error)?;
        let phase = [
            [u[0][0].midpoint(u[0][1]), v[0][0].midpoint(v[0][1])],
            [u[1][0].midpoint(u[1][1]), v[1][0].midpoint(v[1][1])],
        ];
        let source = Nut::new(
            sweep.clone(),
            bore_from_start_mm,
            thread,
            phase,
            relief,
            interrupt,
        )
        .map_err(fail)?;
        let wall = (I::point(source.minimum_thread_wall_lower_bound_mm())
            * (I::ONE - I::point(stock.metric_error())).sqrt())
        .lo;
        if wall < minimum_wall_mm {
            return Err(unresolved());
        }
        let axis = stock
            .point_bounds([[bore_from_start_mm, bore_from_start_mm], [0.; 2], [0.; 2]])
            .map_err(placement_error)?;
        let origin_error = std::array::from_fn(|i| iv(axis[i]) - I::point(requested.origin[i]));
        let axis_alignment = norm(origin_error).hi;
        let old = [requested.u, requested.v, requested.normal()];
        let b: [[I; 3]; 3] = std::array::from_fn(|i| old.map(|a| I::point(a[i])));
        let a: [[I; 3]; 3] = std::array::from_fn(|i| {
            [
                I::point(matrix[i][0]) * I::point(phase[0][0])
                    + I::point(matrix[i][1]) * I::point(phase[1][0]),
                I::point(matrix[i][0]) * I::point(phase[0][1])
                    + I::point(matrix[i][1]) * I::point(phase[1][1]),
                I::point(matrix[i][2]),
            ]
        });
        let delta: [[I; 3]; 3] =
            std::array::from_fn(|i| std::array::from_fn(|j| a[i][j] - b[i][j]));
        let domain = [
            I::new(-profile.major_radius_mm, profile.major_radius_mm),
            I::new(-profile.major_radius_mm, profile.major_radius_mm),
            I::new(0., profile.length_mm),
        ];
        let position = norm(std::array::from_fn(|i| {
            (0..3).fold(origin_error[i], |sum, j| sum + delta[i][j] * domain[j])
        }))
        .hi;
        let operator = delta
            .into_iter()
            .flatten()
            .fold(I::ZERO, |sum, v| sum + I::point(v.magnitude()));
        let tau = I::new(
            std::f64::consts::TAU.next_down(),
            std::f64::consts::TAU.next_up(),
        );
        let derivative = I::point(profile.major_radius_mm)
            * (I::point(3.)
                + tau * (I::point(2.) + I::point(profile.length_mm) / I::point(profile.pitch_mm)))
            + I::point(profile.length_mm);
        let requested_placement =
            placement(std::array::from_fn(|i| old.map(|a| a[i])), requested.origin)?;
        let ca = cofactor(a);
        let cb = cofactor(b);
        let normal_difference = (0..3)
            .flat_map(|i| (0..3).map(move |j| (ca[i][j] - cb[i][j]).magnitude()))
            .fold(I::ZERO, |sum, x| sum + I::point(x));
        let phase_normal = (I::point(4.) * normal_difference
            / (I::ONE - I::point(requested_placement.metric_error())))
        .hi;
        let original = thread.transfer_bounds();
        let transfer = Line8NutTransfer {
            legacy_position_error_mm: (matrix_norm(requested_placement)
                * I::point(original.maximum_position_error_mm))
            .hi,
            legacy_normal_error_rad: (condition(requested_placement)
                * I::point(original.maximum_normal_angle_error_rad))
            .hi,
            phase_position_error_mm: position,
            phase_derivative_error_mm: (operator * derivative).hi,
            phase_normal_error_rad: phase_normal,
            axis_alignment_mm: axis_alignment,
        };
        if !qualified(transfer, policy) {
            return Err(unresolved());
        }
        let boundary = source.boundary(131072, interrupt).map_err(fail)?;
        let volume = source
            .volume_bounds(1e-7, 2048, interrupt)
            .map_err(fail)?
            .bounds_mm3();
        let mut measures = Vec::with_capacity(17);
        for i in 0..17 {
            check(interrupt)?;
            measures.push(
                source
                    .face_measure(i, FacePropertyOptions::default(), interrupt)
                    .map_err(fail)?,
            );
        }
        let mut body = Line8NutBody {
            feature,
            source,
            boundary,
            stock,
            outer: RigidTransform::default(),
            transfer,
            volume,
            measures,
            requested_wall_mm: minimum_wall_mm,
            bounds: Bounds {
                min: [0.; 3],
                max: [0.; 3],
            },
        };
        let bounds = body.point_bounds(body.source.bounds_mm())?;
        policy.point(bounds.map(|v| v[0]))?;
        policy.point(bounds.map(|v| v[1]))?;
        body.bounds = Bounds {
            min: bounds.map(|v| v[0]),
            max: bounds.map(|v| v[1]),
        };
        body.validate_precision(policy, interrupt)?;
        check(interrupt)?;
        Ok(Self {
            body: SolidBody::Line8Nut(Rc::new(body)),
            policy,
            provenance: vec![],
            rust_origins: vec![],
            retained_placement: None,
        })
    }
    pub(crate) fn line8_nut_body(&self) -> Option<&Line8NutBody> {
        if let SolidBody::Line8Nut(b) = &self.body {
            Some(b)
        } else {
            None
        }
    }
    pub fn rust_line8_nut_source(&self) -> Option<&Nut> {
        self.line8_nut_body().map(|b| &b.source)
    }
    pub fn rust_line8_nut_requested_wall_mm(&self) -> Option<f64> {
        self.line8_nut_body().map(|b| b.requested_wall_mm)
    }
    pub fn rust_line8_nut_transfer(&self) -> Option<Line8NutTransfer> {
        self.line8_nut_body().map(|b| b.transfer)
    }
    pub fn rust_line8_nut_contains(&self, p: Point) -> Result<bool, GeometryError> {
        self.policy.point(p)?;
        let body = self.line8_nut_body().ok_or_else(unresolved)?;
        match body
            .source
            .classify_box(body.inverse(p)?, None)
            .map_err(fail)?
        {
            Classification::Inside => Ok(true),
            Classification::Outside => Ok(false),
            Classification::Undetermined => Err(unresolved()),
        }
    }
}
impl Line8NutBody {
    fn validate_precision(
        &self,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<(), GeometryError> {
        for &vertex in self.boundary.vertices() {
            check(interrupt)?;
            let point = self
                .source
                .vertex_bounds(vertex, 2048, interrupt)
                .map_err(fail)?
                .position_mm();
            let world = self.point_bounds(point)?;
            let midpoint = world.map(|v| v[0].midpoint(v[1]));
            if norm(std::array::from_fn(|i| {
                iv(world[i]) - I::point(midpoint[i])
            }))
            .hi > policy.linear_tolerance_mm / 8.
            {
                return Err(unresolved());
            }
        }
        // Integral enclosures can be wider than a source vertex enclosure.
        // Check the actual centroid cache before publishing an infallible
        // faces/topology descriptor, including after every retained pose.
        for measure in &self.measures {
            check(interrupt)?;
            self.point_bounds(measure.centroid_mm())?;
        }
        Ok(())
    }
    fn outer(&self) -> StockClipPlacement {
        placement(self.outer.rotation, self.outer.translation_mm)
            .expect("qualified immutable nut outer placement")
    }
    fn point_bounds(&self, mut p: [[f64; 2]; 3]) -> Result<[[f64; 2]; 3], GeometryError> {
        let shifted = iv(p[0]) + I::point(self.source.bore_from_start_mm());
        p[0] = [shifted.lo, shifted.hi];
        let p = self.stock.point_bounds(p).map_err(placement_error)?;
        self.outer().point_bounds(p).map_err(placement_error)
    }
    fn inverse(&self, p: Point) -> Result<[[f64; 2]; 3], GeometryError> {
        let p = self
            .outer()
            .inverse_point_bounds(p.map(|v| [v, v]))
            .map_err(placement_error)?;
        let mut p = self
            .stock
            .inverse_point_bounds(p)
            .map_err(placement_error)?;
        let x = iv(p[0]) - I::point(self.source.bore_from_start_mm());
        p[0] = [x.lo, x.hi];
        Ok(p)
    }
    fn vector(&self, p: Point) -> Point {
        let a = self
            .stock
            .matrix()
            .map(|r| r.into_iter().zip(p).map(|(a, b)| a * b).sum());
        self.outer.vector(a)
    }
    fn point(&self, p: Point) -> Point {
        self.point_bounds(p.map(|v| [v, v]))
            .expect("qualified source point")
            .map(|v| v[0].midpoint(v[1]))
    }
    pub(super) fn volume_mm3(&self) -> f64 {
        let v = iv(self.volume)
            * iv(self.stock.determinant_bounds())
            * iv(self.outer().determinant_bounds());
        v.lo.midpoint(v.hi)
    }
    pub(super) fn bounds(&self) -> Bounds {
        self.bounds
    }
    pub(super) fn faces(&self) -> Vec<Face> {
        let area = iv(self.stock.area_scale_bounds()) * iv(self.outer().area_scale_bounds());
        self.measures
            .iter()
            .map(|m| {
                let a = iv(m.area_mm2()) * area;
                let center = self
                    .point_bounds(m.centroid_mm())
                    .expect("qualified source face centroid")
                    .map(|v| v[0].midpoint(v[1]));
                Face {
                    index: m.face(),
                    plane: self.plane(m.face()).map(|f| (f.origin_mm, f.normal())),
                    center_mm: center,
                    area_mm2: a.lo.midpoint(a.hi),
                }
            })
            .collect()
    }
    fn plane(&self, face: u32) -> Option<ReferenceFrame> {
        let bounds = self.source.bounds_mm();
        let (origin, u, v) = match face {
            0 => ([0., 0., 0.], [1., 0., 0.], [0., -1., 0.]),
            4 | 10 => ([0., 0., bounds[2][1]], [1., 0., 0.], [0., 1., 0.]),
            15 => (
                [0., 0., self.source.relief().start_mm],
                [1., 0., 0.],
                [0., 1., 0.],
            ),
            8 | 9 => (
                [bounds[0][if face == 8 { 0 } else { 1 }], 0., 0.],
                [0., 1., 0.],
                [0., 0., if face == 8 { -1. } else { 1. }],
            ),
            1 | 2 | 3 | 5 | 6 | 7 => {
                let a = self.source.stock().base()[0][face as usize];
                let b = self.source.stock().base()[0][(face as usize + 1) % 8];
                let dy = b[0] - a[0];
                let dz = b[1] - a[1];
                let length = dy.hypot(dz);
                (
                    [0., a[0], a[1]],
                    [1., 0., 0.],
                    [0., -dy / length, -dz / length],
                )
            }
            _ => return None,
        };
        Some(ReferenceFrame {
            origin_mm: self.point(origin),
            u: self.vector(u),
            v: self.vector(v),
        })
    }
    fn origin(&self, role: String) -> Origin {
        Origin {
            source: self.feature,
            role: format!("profile.line8.retained-nut/{role}"),
        }
    }
    fn edge_origin(&self, edge: SourceEdge) -> Origin {
        self.origin(format!("edge/{edge:?}"))
    }
    fn vertex_origin(&self, vertex: SourceVertex) -> Origin {
        self.origin(format!("vertex/{vertex:?}"))
    }
    pub(crate) fn topology(&self, owner: EntityId, revision: u64) -> TopologySnapshot {
        let mut elements = Vec::new();
        for face in self.boundary.faces() {
            let origin = self.origin(format!("face/{}", face.index()));
            let frame = self.plane(face.index());
            let boundary = face
                .loops()
                .iter()
                .flatten()
                .flat_map(|e| {
                    let mut out = vec![self.edge_origin(e.edge)];
                    out.extend(
                        e.edge
                            .endpoints()
                            .expect("owned source edge")
                            .map(|v| self.vertex_origin(v)),
                    );
                    out
                })
                .collect();
            elements.push(TopologyElement {
                key: TopologyKey {
                    kind: TopologyKind::Face,
                    ordinal: face.index(),
                },
                origins: BTreeSet::from([origin.clone()]),
                boundary,
                frames: frame
                    .into_iter()
                    .map(|frame| AnchoredFrame {
                        origin: origin.clone(),
                        frame,
                    })
                    .collect(),
                geometry: frame.map(|frame| AnalyticGeometry::Plane { frame }),
            });
        }
        for (i, &edge) in self.boundary.edges().iter().enumerate() {
            elements.push(TopologyElement {
                key: TopologyKey {
                    kind: TopologyKind::Edge,
                    ordinal: i as u32,
                },
                origins: BTreeSet::from([self.edge_origin(edge)]),
                boundary: edge
                    .endpoints()
                    .expect("owned source edge")
                    .map(|v| self.vertex_origin(v))
                    .into_iter()
                    .collect(),
                frames: vec![],
                geometry: None,
            });
        }
        for (i, &vertex) in self.boundary.vertices().iter().enumerate() {
            let p = self
                .source
                .vertex_bounds(vertex, 2048, None)
                .expect("qualified immutable source vertex")
                .position_mm();
            let position = self
                .point_bounds(p)
                .expect("qualified source vertex world bounds")
                .map(|v| v[0].midpoint(v[1]));
            elements.push(TopologyElement {
                key: TopologyKey {
                    kind: TopologyKind::Vertex,
                    ordinal: i as u32,
                },
                origins: BTreeSet::from([self.vertex_origin(vertex)]),
                boundary: BTreeSet::new(),
                frames: vec![],
                geometry: Some(AnalyticGeometry::Point {
                    position_mm: position,
                }),
            });
        }
        TopologySnapshot {
            owner,
            revision,
            elements,
        }
    }
    pub(crate) fn planar_loops(&self, index: u32) -> Vec<Vec<(u32, bool)>> {
        if self.plane(index).is_none() {
            return vec![];
        }
        let Some(face) = self.boundary.faces().get(index as usize) else {
            return vec![];
        };
        face.loops()
            .iter()
            .map(|ring| {
                ring.iter()
                    .map(|c| {
                        (
                            self.boundary
                                .edges()
                                .iter()
                                .position(|e| *e == c.edge)
                                .expect("owned source edge index")
                                as u32,
                            !c.forward,
                        )
                    })
                    .collect()
            })
            .collect()
    }
}
impl Line8NutBody {
    pub(in crate::geometry) fn transformed(
        &self,
        outer: RigidTransform,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        let mut body = self.clone();
        let prior = self.outer();
        let calibration = placement(outer.rotation, outer.translation_mm)?;
        body.outer = RigidTransform {
            rotation: std::array::from_fn(|i| {
                std::array::from_fn(|j| {
                    (0..3)
                        .map(|k| outer.rotation[i][k] * self.outer.rotation[k][j])
                        .sum()
                })
            }),
            translation_mm: std::array::from_fn(|i| {
                (0..3)
                    .map(|k| outer.rotation[i][k] * self.outer.translation_mm[k])
                    .sum::<f64>()
                    + outer.translation_mm[i]
            }),
        };
        body.outer.validate(policy)?;
        placement(body.outer.rotation, body.outer.translation_mm)?;
        let exact: [[I; 3]; 3] = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                (0..3).fold(I::ZERO, |s, k| {
                    s + I::point(outer.rotation[i][k]) * I::point(self.outer.rotation[k][j])
                })
            })
        });
        let difference: [[I; 3]; 3] = std::array::from_fn(|i| {
            std::array::from_fn(|j| exact[i][j] - I::point(body.outer.rotation[i][j]))
        });
        let translation: [I; 3] = std::array::from_fn(|i| {
            (0..3).fold(I::point(outer.translation_mm[i]), |s, k| {
                s + I::point(outer.rotation[i][k]) * I::point(self.outer.translation_mm[k])
            }) - I::point(body.outer.translation_mm[i])
        });
        let mut domain = self.source.bounds_mm();
        let shifted = iv(domain[0]) + I::point(self.source.bore_from_start_mm());
        domain[0] = [shifted.lo, shifted.hi];
        let base_bounds = self.stock.point_bounds(domain).map_err(placement_error)?;
        let base_axis = self
            .stock
            .point_bounds([[self.source.bore_from_start_mm(); 2], [0.; 2], [0.; 2]])
            .map_err(placement_error)?;
        let residual = |bounds: [[f64; 2]; 3]| {
            norm(std::array::from_fn(|i| {
                (0..3).fold(translation[i], |s, j| s + difference[i][j] * iv(bounds[j]))
            }))
            .hi
        };
        let operator = difference
            .into_iter()
            .flatten()
            .fold(I::ZERO, |s, v| s + I::point(v.magnitude()));
        let spec = self.source.thread().spec();
        let tau = I::new(
            std::f64::consts::TAU.next_down(),
            std::f64::consts::TAU.next_up(),
        );
        let derivative = matrix_norm(self.stock)
            * (I::point(3. * spec.major_radius_mm)
                * (I::point(3.)
                    + tau * (I::point(2.) + I::point(spec.length_mm) / I::point(spec.pitch_mm)))
                + I::point(spec.length_mm));
        let c0 = cofactor(exact);
        let c1 = cofactor(body.outer.rotation.map(|r| r.map(I::point)));
        let normal_error = (0..3)
            .flat_map(|i| (0..3).map(move |j| (c0[i][j] - c1[i][j]).magnitude()))
            .fold(I::ZERO, |s, v| s + I::point(v));
        let lower = (I::ONE - I::point(prior.metric_error()))
            * (I::ONE - I::point(calibration.metric_error()));
        let normal_error = I::point(4.) * normal_error / lower;
        let scale = matrix_norm(calibration);
        let condition = condition(calibration);
        body.transfer.legacy_position_error_mm =
            (I::point(self.transfer.legacy_position_error_mm) * scale).hi;
        body.transfer.legacy_normal_error_rad =
            (I::point(self.transfer.legacy_normal_error_rad) * condition).hi;
        body.transfer.phase_position_error_mm = (I::point(self.transfer.phase_position_error_mm)
            * scale
            + I::point(residual(base_bounds)))
        .hi;
        body.transfer.phase_derivative_error_mm =
            (I::point(self.transfer.phase_derivative_error_mm) * scale + operator * derivative).hi;
        body.transfer.phase_normal_error_rad =
            (I::point(self.transfer.phase_normal_error_rad) * condition + normal_error).hi;
        body.transfer.axis_alignment_mm =
            (I::point(self.transfer.axis_alignment_mm) * scale + I::point(residual(base_axis))).hi;
        if !qualified(body.transfer, policy) {
            return Err(unresolved());
        }
        let bounds = body.point_bounds(body.source.bounds_mm())?;
        policy.point(bounds.map(|v| v[0]))?;
        policy.point(bounds.map(|v| v[1]))?;
        body.bounds = Bounds {
            min: bounds.map(|v| v[0]),
            max: bounds.map(|v| v[1]),
        };
        body.validate_precision(policy, None)?;
        Ok(body)
    }
    pub(in crate::geometry) fn sample_edge(
        &self,
        index: u32,
        q: multicad_rust_kernel::display::Tessellation,
    ) -> Option<Result<Vec<Point>, multicad_rust_kernel::display::MeshError>> {
        use multicad_rust_kernel::display::MeshError;
        let key = *self.boundary.edges().get(index as usize)?;
        Some((|| {
            if !q.deflection.is_finite()
                || q.deflection <= 0.
                || !(0.01..=1.).contains(&q.angle)
                || q.triangle_limit == 0
            {
                return Err(MeshError::InvalidOptions);
            }
            let scale = matrix_norm(self.stock) * matrix_norm(self.outer());
            let m = I::point(
                self.source
                    .curve_second_derivative_bound_mm(key)
                    .map_err(|_| MeshError::UnresolvedBoundary)?,
            ) * scale;
            let mut n = 1usize;
            while (m / I::point(8.) / I::point(n as f64) / I::point(n as f64)).hi
                > q.deflection / 4.
            {
                n = n.checked_mul(2).ok_or(MeshError::BudgetExceeded)?;
                if n >= q.triangle_limit.min(131072) {
                    return Err(MeshError::BudgetExceeded);
                }
            }
            if (n + 1) * 2048 > 16_000_000 {
                return Err(MeshError::BudgetExceeded);
            }
            let mut out = Vec::with_capacity(n + 1);
            let mut rounding = 0_f64;
            for i in 0..=n {
                let t = i as f64 / n as f64;
                let point = self
                    .source
                    .curve_bounds(key, [t, t], 2048, None)
                    .map_err(|_| MeshError::UnresolvedBoundary)?
                    .position_mm();
                let point = self
                    .point_bounds(point)
                    .map_err(|_| MeshError::UnresolvedBoundary)?;
                let mid = point.map(|v| v[0].midpoint(v[1]));
                rounding =
                    rounding.max(norm(std::array::from_fn(|i| iv(point[i]) - I::point(mid[i]))).hi);
                out.push(mid);
            }
            if (m / I::point(8.) / I::point(n as f64) / I::point(n as f64) + I::point(rounding)).hi
                > q.deflection
            {
                return Err(MeshError::UnresolvedBoundary);
            }
            Ok(out)
        })())
    }
}
