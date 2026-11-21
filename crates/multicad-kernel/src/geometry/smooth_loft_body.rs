//! Owner-thread smooth loft solids: planar caps, C¹ cubic corner rails and
//! exact ruled walls that keep every station's authored edge origin.
use super::*;
use crate::rust_prism_topology::{Builder, place_frame, polygon_planar_loops, position};
use multicad_core::references::{AnalyticGeometry, AnchoredFrame};
use multicad_rust_kernel::planar_boolean::SmoothLoft;

type Point = [f64; 3];

struct LoftFace {
    plane: Option<ReferenceFrame>,
    anchors: Vec<AnchoredFrame>,
    origins: Vec<Origin>,
    center: Point,
}

pub(crate) struct SmoothLoftBody {
    pub(super) geometry: SmoothLoft,
    pub(super) feature: EntityId,
    faces: Vec<LoftFace>,
}

fn subtract(a: Point, b: Point) -> Point {
    std::array::from_fn(|i| a[i] - b[i])
}

fn norm([x, y, z]: Point) -> f64 {
    x.hypot(y).hypot(z)
}

impl Solid {
    /// A smooth loft through every station. `sections` holds each station's
    /// authored edge origins in corner order.
    pub(super) fn from_smooth_loft(
        geometry: SmoothLoft,
        feature: EntityId,
        placement: PlaneFrame,
        sections: &[Vec<Origin>],
        frames: [ReferenceFrame; 2],
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        let count = geometry.count();
        if sections.len() != geometry.heights().len()
            || sections.iter().any(|edges| edges.len() != count)
        {
            return Err(input("Loft section does not have a source origin per edge"));
        }
        let rings = geometry.rings();
        for ring in &rings {
            for (&a, &b) in ring.iter().zip(ring.iter().cycle().skip(1)) {
                policy.point(position(placement, a))?;
                policy.length(norm(subtract(b, a)))?;
            }
        }
        let bounds = geometry.bounds();
        policy.point(position(placement, bounds[0]))?;
        policy.point(position(placement, bounds[1]))?;

        let mut faces = Vec::with_capacity(count + 2);
        for (cap, mut frame) in frames.into_iter().enumerate() {
            let normal = if cap == 0 {
                [0., 0., -1.]
            } else {
                [0., 0., 1.]
            };
            frame.v = cross(normal, frame.u);
            let role = if cap == 0 { "cap.start" } else { "cap.end" };
            let origin = Origin {
                source: feature,
                role: role.into(),
            };
            let ring = &rings[cap];
            let center = std::array::from_fn(|axis| {
                let anchor = ring[0][axis];
                anchor + ring.iter().map(|p| p[axis] - anchor).sum::<f64>() / ring.len() as f64
            });
            faces.push(LoftFace {
                plane: Some(frame),
                anchors: vec![AnchoredFrame {
                    origin: origin.clone(),
                    frame,
                }],
                origins: vec![origin],
                center,
            });
        }
        let middle = (geometry.heights().len() - 1) / 2;
        for i in 0..count {
            // Every station's edge, and the generated side identity of
            // each, names the one continuous wall.
            let origins = sections
                .iter()
                .flat_map(|edges| {
                    let origin = &edges[i];
                    [
                        Origin {
                            source: feature,
                            role: format!("side/{}/{}", origin.source, origin.role),
                        },
                        origin.clone(),
                    ]
                })
                .collect();
            faces.push(LoftFace {
                plane: None,
                anchors: vec![],
                origins,
                center: geometry.point(i, 0.5, middle, 0.5),
            });
        }
        let body = SmoothLoftBody {
            geometry,
            feature,
            faces,
        };
        body.topology(placement, feature, 0)
            .validate()
            .map_err(|message| error("loft topology", GeometryErrorKind::InvalidResult, message))?;
        Ok(Self {
            body: SolidBody::SmoothLoft(Rc::new(body), placement),
            policy,
            provenance: vec![],
            rust_origins: vec![],
            retained_placement: None,
        })
    }

    pub(super) fn from_imported_smooth_loft(
        geometry: SmoothLoft,
        feature: EntityId,
        placement: PlaneFrame,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        let sections: Vec<Vec<_>> = (0..geometry.heights().len())
            .map(|station| {
                (0..geometry.count())
                    .map(|edge| Origin {
                        source: feature,
                        role: format!("import.loft.section.{station}.edge.{edge}"),
                    })
                    .collect()
            })
            .collect();
        let frames = geometry.rings().map(|ring| ReferenceFrame {
            origin_mm: ring[0],
            u: [1., 0., 0.],
            v: [0., 1., 0.],
        });
        Self::from_smooth_loft(geometry, feature, placement, &sections, frames, policy)
    }

    pub(crate) fn smooth_loft_body(&self) -> Option<(&SmoothLoftBody, PlaneFrame)> {
        match &self.body {
            SolidBody::SmoothLoft(body, frame) => Some((body, *frame)),
            _ => None,
        }
    }
}

impl SmoothLoftBody {
    pub(super) fn bounds(&self, placement: PlaneFrame) -> Bounds {
        // The rail poles enclose every wall by the convex hull property.
        let z = self.geometry.pole_heights();
        let mut bounds = Bounds {
            min: [f64::INFINITY; 3],
            max: [f64::NEG_INFINITY; 3],
        };
        for rail in self.geometry.rails() {
            for (&[x, y], &z) in rail.iter().zip(&z) {
                for (axis, coordinate) in position(placement, [x, y, z]).into_iter().enumerate() {
                    bounds.min[axis] = bounds.min[axis].min(coordinate);
                    bounds.max[axis] = bounds.max[axis].max(coordinate);
                }
            }
        }
        bounds
    }

    pub(super) fn faces(&self, placement: PlaneFrame) -> Result<Vec<Face>, GeometryError> {
        let areas = self
            .geometry
            .cap_areas_mm2()
            .into_iter()
            .map(Ok)
            .chain((0..self.geometry.count()).map(|wall| self.geometry.wall_area_mm2(wall)));
        self.faces
            .iter()
            .zip(areas)
            .enumerate()
            .map(|(index, (face, area))| {
                let area_mm2 = area.map_err(|cause| {
                    error(
                        "loft area",
                        GeometryErrorKind::Limit,
                        format!("Could not bound loft area: {cause:?}"),
                    )
                })?;
                let plane = face.plane.map(|p| place_frame(placement, p));
                Ok(Face {
                    index: index as u32,
                    plane: plane.map(|p| (p.origin_mm, p.normal())),
                    center_mm: position(placement, face.center),
                    area_mm2,
                })
            })
            .collect()
    }

    pub(crate) fn topology(
        &self,
        placement: PlaneFrame,
        owner: EntityId,
        revision: u64,
    ) -> TopologySnapshot {
        let count = self.geometry.count();
        let mut builder =
            Builder::new(self.faces.iter().map(|face| face.origins.clone()).collect());
        for (cap, face) in self.faces[..2].iter().enumerate() {
            builder.face(
                cap,
                face.plane.map(|frame| place_frame(placement, frame)),
                None,
                2..count + 2,
            );
        }
        for i in 0..count {
            let previous = (i + count - 1) % count;
            let next = (i + 1) % count;
            builder.face(i + 2, None, None, [0, 1, previous + 2, next + 2]);
        }
        let rings = self.geometry.rings();
        let point = |point| position(placement, point);
        for (cap, ring) in rings.iter().enumerate() {
            for i in 0..count {
                let [a_mm, b_mm] = [ring[i], ring[(i + 1) % count]].map(point);
                builder.non_face(AnalyticGeometry::Line { a_mm, b_mm }, &[cap, i + 2]);
            }
        }
        // The corner rails are cubic B-splines, not analytic curves.
        for i in 0..count {
            let previous = (i + count - 1) % count;
            builder.non_analytic_edge(&[previous + 2, i + 2]);
        }
        for (cap, ring) in rings.iter().enumerate() {
            for (i, &vertex) in ring.iter().enumerate() {
                let previous = (i + count - 1) % count;
                builder.non_face(
                    AnalyticGeometry::Point {
                        position_mm: point(vertex),
                    },
                    &[cap, previous + 2, i + 2],
                );
            }
        }
        let mut elements = builder.into_elements();
        for (element, face) in elements.iter_mut().zip(&self.faces) {
            element.frames = face
                .anchors
                .iter()
                .map(|anchor| AnchoredFrame {
                    origin: anchor.origin.clone(),
                    frame: place_frame(placement, anchor.frame),
                })
                .collect();
        }
        TopologySnapshot {
            owner,
            revision,
            elements,
        }
    }

    /// Whether imported rails, oriented from the lower cap, carry this
    /// loft's rail poles and station knots, in any corner order.
    pub(super) fn matches_rails(
        &self,
        placement: PlaneFrame,
        readback: &[multicad_rust_kernel::step_read::SplineEdge],
        tolerance: f64,
    ) -> bool {
        let z = self.geometry.pole_heights();
        let heights = self.geometry.heights();
        let mut used = vec![false; readback.len()];
        readback.len() == self.geometry.count()
            && self.geometry.rails().iter().all(|rail| {
                let poles: Vec<_> = rail
                    .iter()
                    .zip(&z)
                    .map(|(&[x, y], &z)| position(placement, [x, y, z]))
                    .collect();
                let found = readback.iter().enumerate().position(|(i, curve)| {
                    let knots = curve.knots();
                    !used[i]
                        && knots.len() == heights.len()
                        && knots.iter().zip(heights).all(|(a, b)| {
                            let span = knots[knots.len() - 1] - knots[0];
                            ((a - knots[0]) / span * heights[heights.len() - 1] - b).abs()
                                <= tolerance
                        })
                        && curve.poles_mm().len() == poles.len()
                        && curve
                            .poles_mm()
                            .iter()
                            .zip(&poles)
                            .all(|(a, b)| norm(subtract(*a, *b)) <= tolerance)
                });
                found.inspect(|&i| used[i] = true).is_some()
            })
    }

    pub(crate) fn planar_loops(&self, face: u32) -> Vec<Vec<(u32, bool)>> {
        if face >= 2 {
            return vec![];
        }
        polygon_planar_loops(std::slice::from_ref(&self.geometry.rings()[0]), 0, face)
    }
}
