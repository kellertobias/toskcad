//! Placed occupied material of the unchanged local clipped NPT source.
//! Local STEP export earns independent literal material and source
//! correspondence. Imported geometry and placed STEP admission remain separate.
use super::*;
use multicad_core::references::{
    AnalyticGeometry, AnchoredFrame, TopologyElement, TopologyKey, TopologyKind,
};
use multicad_rust_kernel::{
    axial_pipe_thread::{
        Error as Cause,
        authored::{
            StockClipClassification, StockClipCurve, StockClipFace, StockClipFaceMeasure,
            StockClipMeasureOptions, StockClipPlan, StockClipVertex,
        },
    },
    display::{MeshError, Tessellation},
};
use std::collections::BTreeSet;
mod mesh;
mod placement;
mod step;
#[cfg(test)]
mod tests;
#[derive(Clone)]
pub(crate) struct ClippedNptBody {
    pub(super) feature: EntityId,
    source: StockClipPlan,
    volume: [f64; 2],
    faces: Vec<Face>,
    measures: Vec<StockClipFaceMeasure>,
    vertices: [[f64; 3]; 12],
    pose: RigidTransform,
    step_capability: std::cell::Cell<Option<(u64, bool)>>,
}
fn failure(cause: Cause) -> GeometryError {
    error(
        "clipped authored NPT",
        match cause {
            Cause::Invalid => GeometryErrorKind::Input,
            Cause::Limit => GeometryErrorKind::Limit,
            Cause::Cancelled => GeometryErrorKind::Cancelled,
            Cause::Unsupported | Cause::Unresolved => GeometryErrorKind::Unsupported,
        },
        format!("Original clipped NPT source proof failed: {cause:?}"),
    )
}
fn checkpoint(interrupt: Option<&dyn Fn() -> bool>) -> Result<(), GeometryError> {
    if interrupt.is_some_and(|f| f()) {
        Err(failure(Cause::Cancelled))
    } else {
        Ok(())
    }
}
impl Solid {
    pub(super) fn from_clipped_npt_source(
        source: StockClipPlan,
        feature: EntityId,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, GeometryError> {
        checkpoint(interrupt)?;
        let b = source.bounds_mm();
        policy.point(b.map(|v| v[0]))?;
        policy.point(b.map(|v| v[1]))?;
        // Qualify the finite source vertices before expensive integral queries.
        let mut vertices = [[0.; 3]; 12];
        for (i, key) in source.vertices().into_iter().enumerate() {
            let point = source
                .vertex_bounds(key, Default::default(), interrupt)
                .map_err(failure)?;
            let bounds = point.position_mm();
            let p = bounds.map(|v| v[0].midpoint(v[1]));
            let error = bounds.into_iter().zip(p).fold(0f64, |sum, (b, x)| {
                (sum + (x - b[0]).abs().max((b[1] - x).abs()).next_up()).next_up()
            });
            if error > policy.linear_tolerance_mm / 8. {
                return Err(failure(Cause::Unresolved));
            }
            policy.point(p)?;
            vertices[i] = p;
        }
        // These facade scalar descriptors have no caller-specified absolute
        // volume/moment precision. Scale their query widths dimensionally;
        // the source centroid still must meet the actual linear policy.
        let radius = source.stock_radius_limit_mm();
        let length = source.source().spec().length_mm;
        let extent = radius.max(length);
        // Active clipping earns stockR > rootEnd > L/32. Use the tighter
        // radius/length proxy as well, rather than unbounded slender L².
        let moment_scale = (extent * extent).min(radius * (radius + length));
        let measures = StockClipMeasureOptions {
            area_width_mm2: (policy.linear_tolerance_mm * extent / 32.).max(1e-7),
            moment_width_mm3: (policy.linear_tolerance_mm * moment_scale / 32.).max(1e-7),
            centroid_width_mm: policy.linear_tolerance_mm,
            ..Default::default()
        };
        // Seven distinct bounded source queries: at most7*16M logical work,
        // each64MiB. They retain only seven finite scalar descriptors here;
        // no mesh/file volume or transformed-area assumption is substituted.
        let mut faces = Vec::with_capacity(7);
        let mut bounds = Vec::with_capacity(7);
        let mut area_upper = 0f64;
        for (i, face) in source.physical_faces().into_iter().enumerate() {
            checkpoint(interrupt)?;
            let measure = source
                .face_measure(face, measures, interrupt)
                .map_err(failure)?;
            area_upper = (area_upper + measure.area_mm2()[1]).next_up();
            bounds.push(measure);
            faces.push(Face {
                index: i as u32,
                plane: plane(face, source.source().spec().length_mm)
                    .map(|f| (f.origin_mm, f.normal())),
                center_mm: measure.centroid_mm().map(|v| v[0].midpoint(v[1])),
                area_mm2: measure.area_mm2()[0].midpoint(measure.area_mm2()[1]),
            });
        }
        // The already-earned ACTUAL full source area controls this scalar
        // volume query. A bounding box never substitutes for occupied material.
        let volume_width = (policy.linear_tolerance_mm * area_upper / 32.).max(1e-6);
        let volume = source
            .volume_bounds(volume_width, Default::default(), interrupt)
            .map_err(failure)?
            .source_volume_mm3();
        checkpoint(interrupt)?;
        Ok(Self {
            body: SolidBody::ClippedNpt(Rc::new(ClippedNptBody {
                feature,
                source,
                volume,
                faces,
                measures: bounds,
                vertices,
                pose: RigidTransform::default(),
                step_capability: std::cell::Cell::new(None),
            })),
            policy,
            provenance: vec![],
            rust_origins: vec![],
            retained_placement: None,
        })
    }
    pub(crate) fn clipped_npt_body(&self) -> Option<&ClippedNptBody> {
        match &self.body {
            SolidBody::ClippedNpt(body) => Some(body),
            _ => None,
        }
    }
    /// Original authored material in its unchanged local chart. The separate
    /// placement getter locates it; imported geometry never acquires a recipe.
    pub fn rust_clipped_npt_source(&self) -> Option<&StockClipPlan> {
        self.clipped_npt_body().map(|b| &b.source)
    }
    pub fn rust_clipped_npt_placement(&self) -> Option<RigidTransform> {
        self.clipped_npt_body().map(|b| b.pose)
    }
    /// Placed source-volume enclosure, never independently decoded file mass.
    pub fn rust_clipped_npt_volume_bounds_mm3(&self) -> Option<[f64; 2]> {
        self.clipped_npt_body().map(|b| b.placed_volume())
    }
    /// Retained local source descriptors; use the placement getter separately.
    pub fn rust_clipped_npt_face_measures(&self) -> Option<&[StockClipFaceMeasure]> {
        self.clipped_npt_body().map(|b| b.measures.as_slice())
    }
    pub fn rust_clipped_npt_contains(&self, point_mm: [f64; 3]) -> Result<bool, GeometryError> {
        self.policy.point(point_mm)?;
        let body = self
            .clipped_npt_body()
            .ok_or_else(|| failure(Cause::Unsupported))?;
        let world = body.bounds();
        if (0..3).any(|i| point_mm[i] < world.min[i] || point_mm[i] > world.max[i]) {
            return Ok(false);
        }
        let local = if body.pose == RigidTransform::default() {
            point_mm.map(|x| [x, x])
        } else {
            body.calibration()
                .inverse_point_bounds(point_mm.map(|x| [x, x]))
                .map_err(failure)?
        };
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
fn plane(face: StockClipFace, length: f64) -> Option<ReferenceFrame> {
    match face {
        StockClipFace::StartCap | StockClipFace::EndCap => Some(ReferenceFrame {
            origin_mm: [
                0.,
                0.,
                if face == StockClipFace::StartCap {
                    0.
                } else {
                    length
                },
            ],
            u: [1., 0., 0.],
            v: [
                0.,
                if face == StockClipFace::StartCap {
                    -1.
                } else {
                    1.
                },
                0.,
            ],
        }),
        _ => None,
    }
}
fn face_tag(face: StockClipFace) -> &'static str {
    match face {
        StockClipFace::LowerFlank => "lower",
        StockClipFace::NaturalCrest => "crest.natural",
        StockClipFace::StockCylinder => "crest.stock",
        StockClipFace::UpperFlank => "upper",
        StockClipFace::Root => "root",
        StockClipFace::StartCap => "cap.start",
        StockClipFace::EndCap => "cap.end",
    }
}
impl ClippedNptBody {
    pub(super) fn volume_mm3(&self) -> f64 {
        let volume = self.placed_volume();
        volume[0].midpoint(volume[1])
    }
    pub(super) fn bounds(&self) -> Bounds {
        let b = if self.pose == RigidTransform::default() {
            self.source.bounds_mm()
        } else {
            self.calibration()
                .point_bounds(self.source.bounds_mm())
                .expect("qualified placed source bounds")
        };
        Bounds {
            min: b.map(|v| v[0]),
            max: b.map(|v| v[1]),
        }
    }
    pub(super) fn faces(&self) -> Vec<Face> {
        self.placed_faces()
    }
    fn origin(&self, role: String) -> Origin {
        Origin {
            source: self.feature,
            role: format!("pipe.npt.clipped/{role}"),
        }
    }
    fn edge_origin(&self, key: StockClipCurve) -> Origin {
        self.origin(format!("edge/{key:?}"))
    }
    fn vertex_origin(&self, key: StockClipVertex) -> Origin {
        self.origin(format!("vertex/{key:?}"))
    }
    pub(crate) fn topology(&self, owner: EntityId, revision: u64) -> TopologySnapshot {
        let faces = self.source.physical_faces();
        let edges = self.source.physical_curves();
        let vertices = self.source.vertices();
        let mut elements: Vec<_> = faces
            .into_iter()
            .enumerate()
            .map(|(i, face)| {
                let origin = self.origin(format!("face/{}", face_tag(face)));
                let plane = plane(face, self.source.source().spec().length_mm);
                TopologyElement {
                    key: TopologyKey {
                        kind: TopologyKind::Face,
                        ordinal: i as u32,
                    },
                    origins: BTreeSet::from([origin.clone()]),
                    boundary: BTreeSet::new(),
                    frames: plane
                        .into_iter()
                        .map(|frame| AnchoredFrame {
                            origin: origin.clone(),
                            frame,
                        })
                        .collect(),
                    geometry: plane.map(|frame| AnalyticGeometry::Plane { frame }),
                }
            })
            .collect();
        for (i, key) in edges.into_iter().enumerate() {
            let origin = self.edge_origin(key);
            let ends = key.endpoints().map(|v| self.vertex_origin(v));
            for face in self.source.curve_faces(key) {
                let index = faces
                    .iter()
                    .position(|x| *x == face)
                    .expect("actual source face");
                elements[index].boundary.insert(origin.clone());
                elements[index].boundary.extend(ends.clone());
            }
            elements.push(TopologyElement {
                key: TopologyKey {
                    kind: TopologyKind::Edge,
                    ordinal: i as u32,
                },
                origins: BTreeSet::from([origin]),
                boundary: ends.into_iter().collect(),
                frames: vec![],
                geometry: None,
            });
        }
        for (i, key) in vertices.into_iter().enumerate() {
            elements.push(TopologyElement {
                key: TopologyKey {
                    kind: TopologyKind::Vertex,
                    ordinal: i as u32,
                },
                origins: BTreeSet::from([self.vertex_origin(key)]),
                boundary: BTreeSet::new(),
                frames: vec![],
                geometry: Some(AnalyticGeometry::Point {
                    position_mm: self.vertices[i],
                }),
            });
        }
        let mut snapshot = TopologySnapshot {
            owner,
            revision,
            elements,
        };
        self.place_topology(&mut snapshot);
        snapshot
    }
    pub(crate) fn planar_loops(&self, index: u32) -> Vec<Vec<(u32, bool)>> {
        let faces = self.source.physical_faces();
        let Some(&face) = faces.get(index as usize) else {
            return vec![];
        };
        if !matches!(face, StockClipFace::StartCap | StockClipFace::EndCap) {
            return vec![];
        }
        let edges = self.source.physical_curves();
        vec![
            self.source
                .face_loop(face)
                .coedges()
                .iter()
                .map(|e| {
                    (
                        edges
                            .iter()
                            .position(|x| *x == e.curve)
                            .expect("physical cap edge") as u32,
                        !e.forward,
                    )
                })
                .collect(),
        ]
    }
    pub(super) fn sample_edge(
        &self,
        index: u32,
        options: Tessellation,
    ) -> Option<Result<Vec<[f64; 3]>, MeshError>> {
        let curve = *self.source.physical_curves().get(index as usize)?;
        Some((|| {
            quality(options)?;
            let calibration = self.calibration();
            let result = self
                .source
                .curve_polyline(
                    curve,
                    options.deflection / 4.,
                    options.triangle_limit.min(131072),
                    Default::default(),
                    None,
                )
                .map_err(|e| match e {
                    Cause::Invalid => MeshError::InvalidOptions,
                    Cause::Limit => MeshError::BudgetExceeded,
                    _ => MeshError::UnresolvedBoundary,
                })?;
            self.place_polyline(&result, options, calibration)
        })())
    }
}

fn quality(q: Tessellation) -> Result<(), MeshError> {
    if !q.deflection.is_finite()
        || q.deflection <= 0.
        || !(0.01..=1.).contains(&q.angle)
        || !(1..=5_000_000).contains(&q.triangle_limit)
    {
        Err(MeshError::InvalidOptions)
    } else {
        Ok(())
    }
}
