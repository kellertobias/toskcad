//! Policy-qualified derived external NPT material retaining the exact original
//! stock. Original nominal thread metadata is not relabelled as exact geometry.
use super::*;
type Point = [f64; 3];
use crate::geometry::internal_thread_treatment_body::interval::Interval as I;
use multicad_core::references::{
    AnalyticGeometry, AnchoredFrame, TopologyElement, TopologyKey, TopologyKind,
};
use multicad_rust_kernel::axial_pipe_thread::{
    Error as Cause,
    authored::{
        AuthoredSpec, AuthoredThread, StockClipClassification, StockClipFace,
        StockClipMeasureOptions, StockClipPlacement,
    },
    external_stock::{
        Boundary, Curve as SourceCurve, ExternalStockCut, Face as SourceFace, FaceMeasure,
        Vertex as SourceVertex,
    },
};
use std::collections::BTreeSet;
mod curves;
mod mesh;
#[cfg(test)]
mod tests;
const OPERATION: &str = "derived external NPT stock cut";
#[derive(Clone, Copy, Debug)]
pub struct ExternalNptStockTransfer {
    pub position_error_mm: f64,
    pub derivative_error_mm: f64,
    pub axis_alignment_mm: f64,
}
#[derive(Clone)]
pub(crate) struct ExternalNptStockBody {
    pub(super) feature: EntityId,
    source: ExternalStockCut,
    boundary: Boundary,
    volume: [f64; 2],
    faces: Vec<Face>,
    measures: Vec<FaceMeasure>,
    vertices: Vec<Point>,
    pose: RigidTransform,
    transfer: ExternalNptStockTransfer,
}
fn failure(cause: Cause) -> GeometryError {
    error(
        OPERATION,
        match cause {
            Cause::Invalid => GeometryErrorKind::Input,
            Cause::Limit => GeometryErrorKind::Limit,
            Cause::Cancelled => GeometryErrorKind::Cancelled,
            _ => GeometryErrorKind::Unsupported,
        },
        format!("Exact-stock derived source proof failed: {cause:?}"),
    )
}
fn checkpoint(interrupt: Option<&dyn Fn() -> bool>) -> Result<(), GeometryError> {
    if interrupt.is_some_and(|f| f()) {
        Err(failure(Cause::Cancelled))
    } else {
        Ok(())
    }
}
fn pose_point(pose: RigidTransform, p: Point) -> Point {
    let v = pose.vector(p);
    std::array::from_fn(|i| v[i] + pose.translation_mm[i])
}
fn norm(v: [I; 3]) -> I {
    {
        let n = v[0].square() + v[1].square() + v[2].square();
        I::new(n.lo.max(0.), n.hi).sqrt()
    }
}
fn iv(v: [f64; 2]) -> I {
    I::new(v[0], v[1])
}
impl Solid {
    /// Actual stock datums own the untouched tail. Thread radius and source
    /// placement changes are explicit, measured and separately policy gated.
    pub fn rust_cut_external_npt_stock(
        &self,
        spec: AuthoredSpec,
        feature: EntityId,
        frame: PlaneFrame,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, GeometryError> {
        policy.validate()?;
        frame.checked(policy)?;
        checkpoint(interrupt)?;
        let SolidBody::AxialCylinder(stock, stock_frame, _) = &self.body else {
            return Err(failure(Cause::Unsupported));
        };
        if self.policy != policy
            || policy.boolean_fuzzy_mm != 0.
            || feature.is_nil()
            || feature.uuid_bits() == stock.feature()
            || frame.u != stock_frame.u
            || frame.v != stock_frame.v
        {
            return Err(failure(Cause::Unsupported));
        }
        let thread =
            AuthoredThread::new(spec, policy.linear_tolerance_mm, policy.minimum_feature_mm)
                .map_err(failure)?
                .stock_clip_plan(Default::default(), interrupt)
                .map_err(failure)?;
        let source = ExternalStockCut::new(thread, **stock, Default::default(), interrupt)
            .map_err(failure)?;
        let axes = [stock_frame.u, stock_frame.v, stock_frame.normal()];
        let pose = if self.axial_cylinder().is_some() {
            if stock_frame.u != [1., 0., 0.] || stock_frame.v != [0., 1., 0.] {
                return Err(failure(Cause::Unsupported));
            }
            RigidTransform::default()
        } else {
            let zero = Rational::new(0, 1).expect("zero");
            if stock.circle().center() != [zero; 2] || stock.bottom() != zero {
                return Err(failure(Cause::Unsupported));
            }
            RigidTransform {
                rotation: std::array::from_fn(|i| axes.map(|axis| axis[i])),
                translation_mm: stock_frame.origin,
            }
        };
        pose.validate(policy)?;
        let calibration =
            StockClipPlacement::new(pose.rotation, pose.translation_mm).map_err(failure)?;
        let axis = calibration
            .point_bounds(source.stock_axis_base_bounds_mm())
            .map_err(failure)?;
        let alignment = norm(std::array::from_fn(|i| {
            iv(axis[i]) - I::point(frame.origin[i])
        }))
        .hi;
        let metric = (I::point(1.) + I::point(calibration.metric_error())).sqrt();
        let transfer = ExternalNptStockTransfer {
            position_error_mm: (I::point(alignment)
                + metric * I::point(source.original_source_position_transfer_mm()))
            .hi,
            derivative_error_mm: (metric
                * I::point(source.original_source_derivative_transfer_mm()))
            .hi,
            axis_alignment_mm: alignment,
        };
        // D1 is an absolute original normalized coedge/patch error in mm. Use
        // the same linear budget conservatively; this does not waive normals.
        if transfer.position_error_mm > policy.linear_tolerance_mm
            || transfer.derivative_error_mm > policy.linear_tolerance_mm
        {
            return Err(error(
                OPERATION,
                GeometryErrorKind::Unsupported,
                format!("Transfer exceeds policy: {transfer:?}"),
            ));
        }
        let bounds = calibration
            .point_bounds(source.bounds_mm())
            .map_err(failure)?;
        policy.point(bounds.map(|v| v[0]))?;
        policy.point(bounds.map(|v| v[1]))?;
        let boundary = source
            .boundary(Default::default(), interrupt)
            .map_err(failure)?;
        let mut vertices = Vec::with_capacity(boundary.vertices().len());
        for &key in boundary.vertices() {
            checkpoint(interrupt)?;
            let point = source
                .vertex_bounds(key, Default::default(), interrupt)
                .map_err(failure)?
                .position_mm();
            let p = point.map(|v| v[0].midpoint(v[1]));
            let e = norm(std::array::from_fn(|i| iv(point[i]) - I::point(p[i]))).hi;
            if e > policy.linear_tolerance_mm / 8. {
                return Err(failure(Cause::Unresolved));
            }
            vertices.push(p);
        }
        let extent = source.stock_radius_mm().max(source.stock_length_mm());
        let measures_options = StockClipMeasureOptions {
            area_width_mm2: (policy.linear_tolerance_mm * extent / 32.).max(1e-7),
            moment_width_mm3: (policy.linear_tolerance_mm * extent * extent / 32.).max(1e-7),
            centroid_width_mm: policy.linear_tolerance_mm,
            ..Default::default()
        };
        let mut faces = Vec::with_capacity(boundary.faces().len());
        let mut measures = Vec::with_capacity(boundary.faces().len());
        let mut area = 0_f64;
        for (i, entry) in boundary.faces().iter().enumerate() {
            checkpoint(interrupt)?;
            let measure = source
                .face_measure(entry.face(), measures_options, interrupt)
                .map_err(|cause| {
                    let mut error = failure(cause);
                    error.message = format!("Face {:?}: {}", entry.face(), error.message);
                    error
                })?;
            area = (area + measure.area_mm2()[1]).next_up();
            faces.push(Face {
                index: i as u32,
                plane: plane(entry.face(), &source).map(|f| (f.origin_mm, f.normal())),
                center_mm: measure.centroid_mm().map(|v| v[0].midpoint(v[1])),
                area_mm2: measure.area_mm2()[0].midpoint(measure.area_mm2()[1]),
            });
            measures.push(measure);
        }
        let volume = source
            .volume_bounds(
                (policy.linear_tolerance_mm * area / 32.).max(1e-6),
                Default::default(),
                interrupt,
            )
            .map_err(failure)?
            .source_volume_mm3();
        checkpoint(interrupt)?;
        Ok(Self {
            body: SolidBody::ExternalNptStock(Rc::new(ExternalNptStockBody {
                feature,
                source,
                boundary,
                volume,
                faces,
                measures,
                vertices,
                pose,
                transfer,
            })),
            policy,
            provenance: vec![],
            rust_origins: vec![],
            retained_placement: None,
        })
    }
    pub(crate) fn external_npt_stock_body(&self) -> Option<&ExternalNptStockBody> {
        if let SolidBody::ExternalNptStock(body) = &self.body {
            Some(body)
        } else {
            None
        }
    }
    /// Distinct derived material; the legacy original-clipped getter stays None.
    pub fn rust_external_npt_stock_source(&self) -> Option<&ExternalStockCut> {
        self.external_npt_stock_body().map(|b| &b.source)
    }
    pub fn rust_external_npt_stock_transfer(&self) -> Option<ExternalNptStockTransfer> {
        self.external_npt_stock_body().map(|b| b.transfer)
    }
    pub fn rust_external_npt_stock_contains(&self, point: Point) -> Result<bool, GeometryError> {
        self.policy.point(point)?;
        let body = self
            .external_npt_stock_body()
            .ok_or_else(|| failure(Cause::Unsupported))?;
        let local = body
            .calibration()
            .inverse_point_bounds(point.map(|x| [x, x]))
            .map_err(failure)?;
        match body
            .source
            .classify_box(local, Default::default(), None)
            .map_err(failure)?
            .classification()
        {
            StockClipClassification::Inside => Ok(true),
            StockClipClassification::Outside => Ok(false),
            StockClipClassification::Undetermined => Err(failure(Cause::Unresolved)),
        }
    }
}
fn plane(face: SourceFace, source: &ExternalStockCut) -> Option<ReferenceFrame> {
    let (z, sign) = match face {
        SourceFace::Thread(StockClipFace::StartCap) => (source.stock_bottom_bounds_mm(), -1.),
        SourceFace::Shoulder => (source.thread_join_bounds_mm(), -1.),
        SourceFace::TailCap => (source.stock_top_bounds_mm(), 1.),
        _ => return None,
    };
    let center = source.stock().circle().center().map(Rational::as_f64);
    Some(ReferenceFrame {
        origin_mm: [center[0], center[1], z[0].midpoint(z[1])],
        u: [1., 0., 0.],
        v: [0., sign, 0.],
    })
}
impl ExternalNptStockBody {
    fn calibration(&self) -> StockClipPlacement {
        StockClipPlacement::new(self.pose.rotation, self.pose.translation_mm)
            .expect("qualified immutable derived source pose")
    }
    pub(in crate::geometry) fn transformed(
        &self,
        outer: RigidTransform,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        let mut body = self.clone();
        let calibration =
            StockClipPlacement::new(outer.rotation, outer.translation_mm).map_err(failure)?;
        body.pose = RigidTransform {
            rotation: std::array::from_fn(|i| {
                std::array::from_fn(|j| {
                    (0..3)
                        .map(|k| outer.rotation[i][k] * self.pose.rotation[k][j])
                        .sum()
                })
            }),
            translation_mm: pose_point(outer, self.pose.translation_mm),
        };
        body.pose.validate(policy)?;
        let placed = body
            .calibration()
            .point_bounds(body.source.bounds_mm())
            .map_err(failure)?;
        policy.point(placed.map(|v| v[0]))?;
        policy.point(placed.map(|v| v[1]))?;
        let norm = (I::point(1.) + I::point(calibration.metric_error())).sqrt();
        body.transfer.position_error_mm = (I::point(body.transfer.position_error_mm) * norm).hi;
        body.transfer.derivative_error_mm = (I::point(body.transfer.derivative_error_mm) * norm).hi;
        body.transfer.axis_alignment_mm = (I::point(body.transfer.axis_alignment_mm) * norm).hi;
        if body.transfer.position_error_mm > policy.linear_tolerance_mm
            || body.transfer.derivative_error_mm > policy.linear_tolerance_mm
        {
            return Err(failure(Cause::Unresolved));
        }
        Ok(body)
    }
    pub(super) fn volume_mm3(&self) -> f64 {
        let d = self.calibration().determinant_bounds();
        let v = iv(self.volume) * iv(d);
        v.lo.midpoint(v.hi)
    }
    pub(super) fn bounds(&self) -> Bounds {
        let v = self
            .calibration()
            .point_bounds(self.source.bounds_mm())
            .expect("qualified derived stock bounds");
        Bounds {
            min: v.map(|b| b[0]),
            max: v.map(|b| b[1]),
        }
    }
    fn frame(&self, f: ReferenceFrame) -> ReferenceFrame {
        ReferenceFrame {
            origin_mm: pose_point(self.pose, f.origin_mm),
            u: self.pose.vector(f.u),
            v: self.pose.vector(f.v),
        }
    }
    pub(super) fn faces(&self) -> Vec<Face> {
        let scale = self.calibration().area_scale_bounds();
        self.faces
            .iter()
            .zip(&self.measures)
            .map(|(f, m)| {
                let mut f = f.clone();
                f.center_mm = pose_point(self.pose, f.center_mm);
                let a = iv(m.area_mm2()) * iv(scale);
                f.area_mm2 = a.lo.midpoint(a.hi);
                f.plane = self
                    .boundary
                    .faces()
                    .get(f.index as usize)
                    .and_then(|v| plane(v.face(), &self.source))
                    .map(|p| self.frame(p))
                    .map(|f| (f.origin_mm, f.normal()));
                f
            })
            .collect()
    }
    fn origin(&self, role: String) -> Origin {
        Origin {
            source: self.feature,
            role: format!("pipe.npt.derived-stock/{role}"),
        }
    }
    fn edge_origin(&self, key: SourceCurve) -> Origin {
        self.origin(format!("edge/{key:?}"))
    }
    fn vertex_origin(&self, key: SourceVertex) -> Origin {
        self.origin(format!("vertex/{key:?}"))
    }
    pub(crate) fn topology(&self, owner: EntityId, revision: u64) -> TopologySnapshot {
        let mut elements = Vec::new();
        for (i, face) in self.boundary.faces().iter().enumerate() {
            let origin = self.origin(format!("face/{:?}", face.face()));
            let plane = plane(face.face(), &self.source).map(|f| self.frame(f));
            let boundary = face
                .loops()
                .iter()
                .flatten()
                .flat_map(|c| {
                    let mut v = vec![self.edge_origin(c.curve)];
                    v.extend(c.curve.endpoints().map(|p| self.vertex_origin(p)));
                    v
                })
                .collect();
            elements.push(TopologyElement {
                key: TopologyKey {
                    kind: TopologyKind::Face,
                    ordinal: i as u32,
                },
                origins: BTreeSet::from([origin.clone()]),
                boundary,
                frames: plane
                    .into_iter()
                    .map(|frame| AnchoredFrame {
                        origin: origin.clone(),
                        frame,
                    })
                    .collect(),
                geometry: plane.map(|frame| AnalyticGeometry::Plane { frame }),
            });
        }
        for (i, &edge) in self.boundary.curves().iter().enumerate() {
            elements.push(TopologyElement {
                key: TopologyKey {
                    kind: TopologyKind::Edge,
                    ordinal: i as u32,
                },
                origins: BTreeSet::from([self.edge_origin(edge)]),
                boundary: edge
                    .endpoints()
                    .map(|v| self.vertex_origin(v))
                    .into_iter()
                    .collect(),
                frames: vec![],
                geometry: None,
            });
        }
        for (i, &vertex) in self.boundary.vertices().iter().enumerate() {
            elements.push(TopologyElement {
                key: TopologyKey {
                    kind: TopologyKind::Vertex,
                    ordinal: i as u32,
                },
                origins: BTreeSet::from([self.vertex_origin(vertex)]),
                boundary: BTreeSet::new(),
                frames: vec![],
                geometry: Some(AnalyticGeometry::Point {
                    position_mm: pose_point(self.pose, self.vertices[i]),
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
        let Some(face) = self.boundary.faces().get(index as usize) else {
            return vec![];
        };
        if plane(face.face(), &self.source).is_none() {
            return vec![];
        }
        face.loops()
            .iter()
            .map(|ring| {
                ring.iter()
                    .map(|c| {
                        (
                            self.boundary
                                .curves()
                                .iter()
                                .position(|v| *v == c.curve)
                                .expect("physical derived edge") as u32,
                            !c.forward,
                        )
                    })
                    .collect()
            })
            .collect()
    }
}
