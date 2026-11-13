//! Box stock with one oblique cylindrical cut on its retained supports.
//!
//! The stock planes, tool caps and tool cylinder stay analytic. Display and
//! manufacturing triangles are sampled from this boundary and never define it.
use super::*;
use crate::rust_prism_topology::{Builder, place_frame, position};
use multicad_core::references::{
    AnalyticGeometry, AnchoredFrame, ReferenceFrame, TopologyKind, TopologySnapshot,
};
use multicad_rust_kernel::oblique_bore::{self, Cut, ObliqueBore, Support, Tool};
use std::collections::BTreeSet;

fn subtract(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| a[i] - b[i])
}

const OPERATION: &str = "oblique cylindrical cut";
const MISSING_PLANE: &str = "Box stock topology does not publish every support plane";

/// Authored identity of one source face. Each source frame is kept
/// literally as captured and in the stock chart for later placements.
struct Source {
    origins: Vec<Origin>,
    frames: Vec<(AnchoredFrame, ReferenceFrame)>,
}

pub(crate) struct ObliqueBoreBody {
    pub(super) geometry: ObliqueBore,
    pub(super) feature: EntityId,
    /// Stock plane frames in the stock chart, indexed like the kernel.
    planes: [ReferenceFrame; 6],
    stock: [Source; 6],
    /// Tool start cap, end cap and side.
    tool: [Source; 3],
    /// Placement that captured the literal source frames.
    authored: PlaneFrame,
}

fn refusal(kind: GeometryErrorKind, message: &str) -> GeometryError {
    error(OPERATION, kind, message)
}

fn kernel_error(cause: oblique_bore::Error) -> GeometryError {
    match cause {
        oblique_bore::Error::Invalid => refusal(
            GeometryErrorKind::Input,
            "Invalid stock or cylinder for an oblique cut",
        ),
        oblique_bore::Error::Unsupported => refusal(
            GeometryErrorKind::Unsupported,
            "Oblique cylindrical subtraction supports one connected result shell",
        ),
        oblique_bore::Error::Unresolved => refusal(
            GeometryErrorKind::Unsupported,
            "Oblique cylindrical subtraction is unresolved within the certified margin",
        ),
        oblique_bore::Error::Cancelled => refusal(
            GeometryErrorKind::Cancelled,
            "Oblique cylindrical subtraction was cancelled",
        ),
    }
}

fn direction(frame: PlaneFrame, local: [f64; 3]) -> [f64; 3] {
    position(
        PlaneFrame {
            origin: [0.; 3],
            ..frame
        },
        local,
    )
}

fn sources(snapshot: &TopologySnapshot, ordinal: usize, chart: PlaneFrame) -> Source {
    let face = snapshot
        .elements
        .iter()
        .find(|e| e.key.kind == TopologyKind::Face && e.key.ordinal as usize == ordinal);
    let axes = [chart.u, chart.v, chart.normal()];
    let local = |vector: [f64; 3]| axes.map(|axis| dot(vector, axis));
    let frames = face.map_or(vec![], |face| {
        let frames = face.frames.iter().map(|anchor| {
            let charted = ReferenceFrame {
                origin_mm: local(subtract(anchor.frame.origin_mm, chart.origin)),
                u: local(anchor.frame.u),
                v: local(anchor.frame.v),
            };
            (anchor.clone(), charted)
        });
        frames.collect()
    });
    Source {
        origins: face.map_or(vec![], |face| face.origins.iter().cloned().collect()),
        frames,
    }
}

impl Solid {
    /// Subtract a finite cylinder that crosses hole-free box stock obliquely.
    /// The result retains the stock planes and the actual cylinder and caps.
    /// Refusals leave both operands unchanged.
    pub fn subtract_oblique_cylinder(
        &self,
        tool: &Self,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Option<Self>, GeometryError> {
        policy.validate()?;
        if policy.boolean_fuzzy_mm != 0. || policy != self.policy || policy != tool.policy {
            return Err(refusal(
                GeometryErrorKind::Input,
                "An oblique cut requires the exact geometry policy of both operands",
            ));
        }
        self.oblique_cylinder_cut(tool, interrupt)
            .unwrap_or_else(|| {
                Err(refusal(
                    GeometryErrorKind::Unsupported,
                    "An oblique cut requires hole-free box stock and a cylinder oblique to it",
                ))
            })?
            .map(|result| result.retain_authored_boolean_sources(&[self, tool]))
            .transpose()
    }

    /// Hole-free box stock as its exact size, the placement of its minimum
    /// corner chart, that chart in the unplaced source and the source feature.
    fn box_stock(&self) -> Option<([Rational; 3], PlaneFrame, PlaneFrame, EntityId)> {
        if let Some((prism, frame)) = self.rust_prism_frame() {
            return (prism.rounded_end().is_none() && prism.inner_loops().next().is_none()).then(
                || {
                    let feature = EntityId::from_uuid_bits(prism.feature());
                    (prism.dimensions(), frame, PlaneFrame::default(), feature)
                },
            );
        }
        let placed = self.polygon_prism()?;
        let [ring] = placed.prism.rings() else {
            return None;
        };
        let zero = Rational::new(0, 1).expect("zero is an exact length");
        let size: [Rational; 2] =
            std::array::from_fn(|i| ring.iter().map(|p| p[i]).fold(zero, Rational::max));
        let corners = ring;
        // A rectangle from the profile origin: every corner is extreme and
        // consecutive corners share exactly one coordinate.
        let rectangle = corners.len() == 4
            && corners
                .iter()
                .all(|p| (0..2).all(|i| p[i] == zero || p[i] == size[i]))
            && (0..4).all(|i| {
                let [a, b] = [corners[i], corners[(i + 1) % 4]];
                (a[0] == b[0]) != (a[1] == b[1])
            });
        if !rectangle
            || !placed.prism.circular_holes().is_empty()
            || !placed.pockets.is_empty()
            || placed.bore_treatments.is_some()
        {
            return None;
        }
        let height = placed.prism.height();
        let frame = placed.frame;
        Some(if placed.reverse_sweep {
            // Swapping the profile axes keeps a right-handed chart whose
            // travel follows the reversed sweep from the same exact origin.
            let chart = PlaneFrame {
                origin: frame.origin,
                u: frame.v,
                v: frame.u,
            };
            let local = PlaneFrame {
                origin: [0.; 3],
                u: [0., 1., 0.],
                v: [1., 0., 0.],
            };
            ([size[1], size[0], height], chart, local, placed.feature)
        } else {
            let local = PlaneFrame::default();
            ([size[0], size[1], height], frame, local, placed.feature)
        })
    }

    /// Topology of the stock in its unplaced source chart.
    fn box_stock_source(&self, owner: EntityId) -> Result<TopologySnapshot, String> {
        if let Some(prism) = self.rust_prism() {
            return crate::rust_prism_topology::snapshot(prism, PlaneFrame::default(), owner, 0);
        }
        let mut local = self.polygon_prism().expect("box stock").clone();
        local.frame = PlaneFrame::default();
        local.snapshot(owner, 0)
    }

    /// `None` when the operands belong to another provider: the stock is not
    /// a plain box or the tool axis is parallel to a stock edge.
    pub(in crate::geometry) fn oblique_cylinder_cut(
        &self,
        tool: &Self,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        let (cylinder, tool_frame) = tool.axial_cylinder_frame()?;
        let (exact_size, frame, local_chart, feature) = self.box_stock()?;
        let size = exact_size.map(Rational::as_f64);
        let height = cylinder.height().as_f64();
        let stock_axes = [frame.u, frame.v, frame.normal()];
        let local = |vector: [f64; 3]| stock_axes.map(|axis| dot(vector, axis));
        // The orthogonal provider owns every travel within this chart noise
        // of a stock edge direction. The same bound snaps a single component
        // so a tool parallel to one face pair keeps exact generators.
        let numeric = 32. * f64::EPSILON * size.into_iter().fold(1., f64::max);
        let mut travel = local(tool_frame.normal());
        let aligned = travel
            .iter()
            .filter(|c| (*c * height).abs() <= numeric)
            .count();
        if aligned > 1 {
            return None;
        }
        let snapped = travel
            .iter()
            .any(|c| *c != 0. && (*c * height).abs() <= numeric);
        for component in &mut travel {
            if (*component * height).abs() <= numeric {
                *component = 0.;
            }
        }
        let base: [f64; 3] = std::array::from_fn(|i| tool_frame.origin[i] - frame.origin[i]);
        // Exact contact is decidable only when this chart reproduces the
        // source operands: binary lengths, a signed-cardinal stock chart, an
        // exact origin difference and an unsnapped travel.
        let binary = |length: Rational| {
            let (numerator, denominator) = length.ratio();
            denominator.count_ones() == 1 && numerator.unsigned_abs() <= 1 << 53
        };
        let radius = cylinder.circle().radius();
        let literal = !snapped
            && exact_size.into_iter().all(binary)
            && binary(radius)
            && binary(cylinder.height())
            && stock_axes
                .as_flattened()
                .iter()
                .all(|c| *c == 0. || c.abs() == 1.)
            && (0..3).all(|i| {
                base[i] + frame.origin[i] == tool_frame.origin[i]
                    && tool_frame.origin[i] - base[i] == frame.origin[i]
            });
        let cutter = Tool {
            base_mm: local(base),
            axes: [local(tool_frame.u), local(tool_frame.v), travel],
            radius_mm: radius.as_f64(),
            height_mm: height,
            literal,
        };
        Some(self.oblique_cut(tool, (size, frame, local_chart, feature), cutter, interrupt))
    }

    fn oblique_cut(
        &self,
        tool: &Self,
        (size, frame, local_chart, feature): ([f64; 3], PlaneFrame, PlaneFrame, EntityId),
        cutter: Tool,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Option<Self>, GeometryError> {
        let geometry =
            match oblique_bore::subtract(size, cutter, interrupt).map_err(kernel_error)? {
                Cut::Unchanged => return self.share().map(Some),
                Cut::Empty => return Ok(None),
                Cut::Body(body) => *body,
            };
        let owner = EntityId::new();
        let failure = |message: String| refusal(GeometryErrorKind::InvalidResult, &message);
        let raw = self.box_stock_source(owner).map_err(failure)?;
        // Whole-solid origins stay on the result; the faces keep the rest.
        let mut plain = self.share()?;
        plain.rust_origins.clear();
        let stock = plain.topology_snapshot(owner, 0)?;
        let cylinder = tool.topology_snapshot(owner, 0)?;
        // Identify each stock plane among the unplaced source faces by its
        // actual support plane, then express that plane in the cut chart.
        let chart_axes = [local_chart.u, local_chart.v, local_chart.normal()];
        let mut planes = [ReferenceFrame {
            origin_mm: [0.; 3],
            u: [0.; 3],
            v: [0.; 3],
        }; 6];
        let mut ordinals = [0; 6];
        for surface in 0..6 {
            let support = Support::Stock {
                axis: surface / 2,
                high: surface % 2 == 1,
            };
            let (point, normal) = geometry.plane(support).expect("stock planes");
            let [point, normal] = [position(local_chart, point), direction(local_chart, normal)];
            let found = raw.elements.iter().find_map(|element| {
                let Some(AnalyticGeometry::Plane { frame }) = element.geometry else {
                    return None;
                };
                let offset = dot(subtract(frame.origin_mm, point), normal);
                (dot(frame.normal(), normal) > 0.5 && offset == 0.).then_some((element.key, frame))
            });
            let (key, plane) = found.ok_or_else(|| failure(MISSING_PLANE.into()))?;
            let chart = |vector: [f64; 3]| chart_axes.map(|axis| dot(vector, axis));
            planes[surface] = ReferenceFrame {
                origin_mm: chart(plane.origin_mm),
                u: chart(plane.u),
                v: chart(plane.v),
            };
            ordinals[surface] = key.ordinal as usize;
        }
        let body = ObliqueBoreBody {
            geometry,
            feature,
            planes,
            stock: ordinals.map(|ordinal| sources(&stock, ordinal, frame)),
            tool: [0, 1, 2].map(|ordinal| sources(&cylinder, ordinal, frame)),
            authored: frame,
        };
        let bounds = body.bounds(frame);
        self.policy.point(bounds.min)?;
        self.policy.point(bounds.max)?;
        Ok(Some(Self {
            body: SolidBody::ObliqueBore(Rc::new(body), frame),
            policy: self.policy,
            provenance: vec![],
            rust_origins: self.rust_origins.clone(),
            retained_placement: None,
        }))
    }

    pub(crate) fn oblique_bore_body(&self) -> Option<(&ObliqueBoreBody, PlaneFrame)> {
        match &self.body {
            SolidBody::ObliqueBore(body, frame) => Some((body, *frame)),
            _ => None,
        }
    }
}

impl ObliqueBoreBody {
    pub(super) fn bounds(&self, frame: PlaneFrame) -> Bounds {
        let normal = frame.normal();
        let ranges: [[f64; 2]; 3] = std::array::from_fn(|i| {
            self.geometry
                .extent([frame.u[i], frame.v[i], normal[i]])
                .map(|value| frame.origin[i] + value)
        });
        Bounds {
            min: ranges.map(|range| range[0]),
            max: ranges.map(|range| range[1]),
        }
    }

    pub(super) fn faces(&self, frame: PlaneFrame) -> Vec<Face> {
        self.geometry
            .faces()
            .iter()
            .enumerate()
            .map(|(index, face)| Face {
                index: index as u32,
                plane: self
                    .geometry
                    .plane(face.support)
                    .map(|(point, normal)| (position(frame, point), direction(frame, normal))),
                center_mm: position(frame, face.center_mm),
                area_mm2: face.area_mm2,
            })
            .collect()
    }

    fn source(&self, support: Support) -> &Source {
        match support {
            Support::Stock { axis, high } => &self.stock[2 * axis + usize::from(high)],
            Support::ToolStart => &self.tool[0],
            Support::ToolEnd => &self.tool[1],
            Support::ToolSide => &self.tool[2],
        }
    }

    /// Plane frame of a planar support in the stock chart, normal outward.
    fn plane(&self, support: Support) -> Option<ReferenceFrame> {
        let (origin_mm, normal) = self.geometry.plane(support)?;
        Some(match support {
            Support::Stock { axis, high } => self.planes[2 * axis + usize::from(high)],
            _ => {
                let [u, v, _] = self.geometry.tool().axes;
                let v = if dot(cross(u, v), normal) < 0. {
                    v.map(|c| -c)
                } else {
                    v
                };
                ReferenceFrame { origin_mm, u, v }
            }
        })
    }

    pub(crate) fn topology(
        &self,
        frame: PlaneFrame,
        owner: EntityId,
        revision: u64,
    ) -> Result<TopologySnapshot, String> {
        let cut = &self.geometry;
        let mut build = Builder::new(
            cut.faces()
                .iter()
                .map(|face| self.source(face.support).origins.clone())
                .collect(),
        );
        let mut neighbors = vec![BTreeSet::new(); cut.faces().len()];
        let mut incident = vec![BTreeSet::new(); cut.vertices().len()];
        for edge in cut.edges() {
            let [a, b] = edge.faces;
            neighbors[a].insert(b);
            neighbors[b].insert(a);
            for end in edge.ends {
                incident[end].extend(edge.faces);
            }
        }
        for (index, face) in cut.faces().iter().enumerate() {
            let plane = self
                .plane(face.support)
                .map(|plane| place_frame(frame, plane));
            build.face(index, plane, None, neighbors[index].iter().copied());
        }
        for edge in cut.edges() {
            let geometry = match edge.curve {
                None => {
                    let [a_mm, b_mm] = edge.ends.map(|end| position(frame, cut.vertices()[end]));
                    AnalyticGeometry::Line { a_mm, b_mm }
                }
                Some(arc) => {
                    let principal = arc.principal();
                    let placed = place_frame(
                        frame,
                        ReferenceFrame {
                            origin_mm: principal.center_mm,
                            u: principal.u,
                            v: principal.v,
                        },
                    );
                    if matches!(edge.supports[0], Support::Stock { .. }) {
                        AnalyticGeometry::EllipticArc {
                            frame: placed,
                            semiaxes_mm: principal.semiaxes_mm,
                            start_rad: principal.start_rad,
                            sweep_rad: principal.sweep_rad,
                        }
                    } else {
                        AnalyticGeometry::CircularArc {
                            frame: placed,
                            radius_mm: cut.tool().radius_mm,
                            start_rad: principal.start_rad,
                            sweep_rad: principal.sweep_rad,
                        }
                    }
                }
            };
            build.non_face(geometry, &edge.faces);
        }
        for (vertex, faces) in cut.vertices().iter().zip(incident) {
            let faces: Vec<usize> = faces.into_iter().collect();
            build.non_face(
                AnalyticGeometry::Point {
                    position_mm: position(frame, *vertex),
                },
                &faces,
            );
        }
        let mut elements = build.into_elements();
        // Source frames are literal at the placement that captured them and
        // follow the stock chart rigidly afterwards.
        for (element, face) in elements.iter_mut().zip(cut.faces()) {
            let anchors = self.source(face.support).frames.iter();
            element.frames = anchors
                .map(|(literal, charted)| AnchoredFrame {
                    origin: literal.origin.clone(),
                    frame: if frame == self.authored {
                        literal.frame
                    } else {
                        place_frame(frame, *charted)
                    },
                })
                .collect();
        }
        let snapshot = TopologySnapshot {
            owner,
            revision,
            elements,
        };
        snapshot.validate()?;
        Ok(snapshot)
    }
}
