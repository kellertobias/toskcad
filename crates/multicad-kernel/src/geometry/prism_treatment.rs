//! Constant-size fillets and symmetric chamfers on any selection of edges
//! of a plain prism: a line-and-arc profile extruded along its normal, such
//! as a box or a sketch region extrusion.
//!
//! Side edges (along the extrusion) are treated in the profile: their
//! corners are rounded or bevelled and the prism is extruded again. Cap
//! edges become exact blend faces on a retained boundary; see [`blend`].
//! Where both are selected at one corner, the cap blends continue around
//! the rounded corner, and at equal sizes the corner blend is a sphere.
mod blend;
mod profile;

use super::authored_sweep::FaceName;
use super::*;
use crate::{CurvedEdge2d, CurvedProfile, EdgeTreatment, TreatmentEdge};
use blend::{Kind, Plan, refused};
use multicad_core::references::{
    AnalyticGeometry, ReferenceFrame, TopologyElement, TopologyKey, TopologyKind, TopologySnapshot,
};
use multicad_rust_kernel::normal_sweep::Edge as SweepEdge;
use profile::{CornerCut, Seg, Source};

/// A plain prism in its local chart: the ring runs counterclockwise about
/// the frame normal and spans `span` along it.
struct Prism {
    frame: PlaneFrame,
    ring: Vec<Seg>,
    span: [f64; 2],
    cap_names: [FaceName; 2],
    side_names: Vec<FaceName>,
    sweep: Option<SweepSource>,
}

/// What a normal sweep needs to be extruded again with rounded corners.
struct SweepSource {
    reverse: bool,
    feature: EntityId,
    origins: Vec<Origin>,
    authored_reversed: bool,
    path_origin: Option<Origin>,
    cap_anchors: Option<[ReferenceFrame; 2]>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Pick {
    /// Edge `i` of the cap at level 0 (bottom) or 1 (top).
    Cap(usize, usize),
    /// The side edge through vertex `k`, where edge `k` starts.
    Side(usize),
}

fn seg(edge: SweepEdge) -> Option<Seg> {
    match edge {
        SweepEdge::Line { a, b } => Some(Seg::Line { a, b }),
        SweepEdge::Arc {
            center,
            a,
            b,
            clockwise,
        } => Some(Seg::Arc {
            c: center,
            r: profile::norm(profile::sub(a, center)),
            a,
            b,
            cw: clockwise,
        }),
        SweepEdge::Circle { .. } => None,
    }
}

fn curved(seg: Seg) -> CurvedEdge2d {
    match seg {
        Seg::Line { a, b } => CurvedEdge2d::Line { a, b },
        Seg::Arc { c, a, b, cw, .. } => CurvedEdge2d::Arc {
            center: c,
            a,
            b,
            clockwise: cw,
        },
    }
}

fn reversed(seg: Seg) -> Seg {
    match seg {
        Seg::Line { a, b } => Seg::Line { a: b, b: a },
        Seg::Arc { c, r, a, b, cw } => Seg::Arc {
            c,
            r,
            a: b,
            b: a,
            cw: !cw,
        },
    }
}

/// The authored name of a source face, with its attachment frame first.
fn source_name(element: Option<&TopologyElement>) -> FaceName {
    let Some(element) = element else {
        return FaceName {
            origins: vec![],
            support: None,
            anchor: None,
        };
    };
    let mut origins: Vec<_> = element.origins.iter().cloned().collect();
    let anchor = element.frames.first();
    if let Some(anchor) = anchor
        && let Some(position) = origins.iter().position(|o| *o == anchor.origin)
    {
        origins.swap(0, position);
    }
    FaceName {
        anchor: anchor
            .filter(|anchor| origins.first() == Some(&anchor.origin))
            .map(|anchor| anchor.frame),
        origins,
        support: None,
    }
}

fn face(topology: &TopologySnapshot, ordinal: u32) -> Option<&TopologyElement> {
    let key = TopologyKey {
        kind: TopologyKind::Face,
        ordinal,
    };
    topology.elements.iter().find(|element| element.key == key)
}

/// The planar source face through `point` with outward `normal`.
fn plane_face(
    topology: &TopologySnapshot,
    point: [f64; 3],
    normal: [f64; 3],
    tolerance: f64,
) -> Option<&TopologyElement> {
    topology.elements.iter().find(|element| {
        element.key.kind == TopologyKind::Face
            && matches!(element.geometry, Some(AnalyticGeometry::Plane { frame })
                if dot(frame.normal(), normal) > 1. - 1e-12
                    && dot(frame.normal(), std::array::from_fn(|i| point[i] - frame.origin_mm[i])).abs() <= tolerance)
    })
}

impl Prism {
    fn local(&self, p: [f64; 3]) -> [f64; 3] {
        let d: [f64; 3] = std::array::from_fn(|i| p[i] - self.frame.origin[i]);
        [
            dot(d, self.frame.u),
            dot(d, self.frame.v),
            dot(d, self.frame.normal()),
        ]
    }

    fn level(&self, z: f64, tolerance: f64) -> Option<usize> {
        (0..2).find(|&level| (z - self.span[level]).abs() <= tolerance)
    }

    /// Which prism edge the analytic edge geometry is.
    fn pick(&self, geometry: AnalyticGeometry, tolerance: f64) -> Option<Pick> {
        let same = |p: [f64; 2], q: [f64; 2]| profile::norm(profile::sub(p, q)) <= tolerance;
        let ends = |seg: Seg, p: [f64; 2], q: [f64; 2]| {
            (same(seg.start(), p) && same(seg.end(), q))
                || (same(seg.start(), q) && same(seg.end(), p))
        };
        match geometry {
            AnalyticGeometry::Line { a_mm, b_mm } => {
                let (a, b) = (self.local(a_mm), self.local(b_mm));
                let (pa, pb) = ([a[0], a[1]], [b[0], b[1]]);
                if (a[2] - b[2]).abs() <= tolerance {
                    let level = self.level(a[2], tolerance)?;
                    let i = self
                        .ring
                        .iter()
                        .position(|&seg| matches!(seg, Seg::Line { .. }) && ends(seg, pa, pb))?;
                    return Some(Pick::Cap(level, i));
                }
                let k = self.ring.iter().position(|seg| same(seg.start(), pa))?;
                let spans_height = same(pa, pb)
                    && self.level(a[2].min(b[2]), tolerance) == Some(0)
                    && self.level(a[2].max(b[2]), tolerance) == Some(1);
                spans_height.then_some(Pick::Side(k))
            }
            AnalyticGeometry::CircularArc {
                frame,
                radius_mm,
                start_rad,
                sweep_rad,
            } => {
                let at = |angle: f64| {
                    let (s, c) = angle.portable_sin_cos();
                    self.local(std::array::from_fn(|i| {
                        frame.origin_mm[i] + radius_mm * (c * frame.u[i] + s * frame.v[i])
                    }))
                };
                let center = self.local(frame.origin_mm);
                let (p, q) = (at(start_rad), at(start_rad + sweep_rad));
                if dot(frame.normal(), self.frame.normal()).abs() < 1. - 1e-12 {
                    return None;
                }
                let level = self.level(center[2], tolerance)?;
                let i = self.ring.iter().position(|&seg| {
                    matches!(seg, Seg::Arc { c, r, .. }
                        if same(c, [center[0], center[1]])
                            && (r - radius_mm).abs() <= tolerance)
                        && ends(seg, [p[0], p[1]], [q[0], q[1]])
                })?;
                Some(Pick::Cap(level, i))
            }
            _ => None,
        }
    }
}

impl Solid {
    /// The prism views of this body: a normal sweep has one; a box cell
    /// has one along each axis.
    fn prism_views(&self, topology: &TopologySnapshot) -> Result<Vec<Prism>, GeometryError> {
        let tolerance = self.policy.linear_tolerance_mm;
        if let Some((body, frame)) = self.normal_sweep_body() {
            let Some(view) = body.treatable_prism() else {
                return Ok(vec![]);
            };
            let Some(ring) = view
                .ring
                .iter()
                .map(|&edge| seg(edge))
                .collect::<Option<Vec<_>>>()
            else {
                return Ok(vec![]);
            };
            // Cap ordinal 0 lies at height 0, the start of a reversed sweep's
            // top.
            let ordinals = if view.reverse { [1, 0] } else { [0, 1] };
            let side_names = (0..ring.len())
                .map(|i| source_name(face(topology, 2 + i as u32)))
                .collect();
            let sweep =
                (self.provenance.is_empty() && self.retained_placement.is_none()).then(|| {
                    SweepSource {
                        reverse: view.reverse,
                        feature: view.feature,
                        origins: view.origins.iter().map(|&o| o.clone()).collect(),
                        authored_reversed: view.authored_reversed,
                        path_origin: view.path_origin.cloned(),
                        cap_anchors: view.cap_anchors,
                    }
                });
            return Ok(vec![Prism {
                frame,
                ring,
                span: view.span,
                cap_names: ordinals.map(|ordinal| source_name(face(topology, ordinal))),
                side_names,
                sweep,
            }]);
        }
        let Some(cells) = self.orthogonal_cells()? else {
            return Ok(vec![]);
        };
        let [stock] = cells.as_slice() else {
            return Ok(vec![]);
        };
        Ok([2, 0, 1]
            .into_iter()
            .map(|along| {
                let [i, j] = [(along + 1) % 3, (along + 2) % 3];
                let axis = |k: usize| std::array::from_fn(|m| if m == k { 1. } else { 0. });
                let frame = PlaneFrame {
                    origin: stock.min,
                    u: axis(i),
                    v: axis(j),
                };
                let size = std::array::from_fn::<f64, 3, _>(|k| stock.max[k] - stock.min[k]);
                let corners = [[0., 0.], [size[i], 0.], [size[i], size[j]], [0., size[j]]];
                let ring: Vec<_> = (0..4)
                    .map(|k| Seg::Line {
                        a: corners[k],
                        b: corners[(k + 1) % 4],
                    })
                    .collect();
                let point = |p: [f64; 3]| -> [f64; 3] {
                    std::array::from_fn(|m| {
                        frame.origin[m]
                            + p[0] * frame.u[m]
                            + p[1] * frame.v[m]
                            + p[2] * frame.normal()[m]
                    })
                };
                let normal = frame.normal();
                let cap_names = [0, 1].map(|level| {
                    let sign = if level == 1 { 1. } else { -1. };
                    let z = if level == 1 { size[along] } else { 0. };
                    source_name(plane_face(
                        topology,
                        point([0., 0., z]),
                        normal.map(|n| sign * n),
                        tolerance,
                    ))
                });
                let side_names = ring
                    .iter()
                    .map(|&seg| {
                        let out = seg.inward(seg.start()).map(|x| -x);
                        let a = seg.start();
                        let world_out =
                            std::array::from_fn(|m| out[0] * frame.u[m] + out[1] * frame.v[m]);
                        source_name(plane_face(
                            topology,
                            point([a[0], a[1], 0.]),
                            world_out,
                            tolerance,
                        ))
                    })
                    .collect();
                Prism {
                    frame,
                    ring,
                    span: [0., size[along]],
                    cap_names,
                    side_names,
                    sweep: None,
                }
            })
            .collect())
    }

    /// Fillet or chamfer any selection of edges of a plain prism. `None`
    /// when this body is not a plain prism or an edge is not one of its
    /// prism edges, so that the caller reports its own refusal.
    pub(crate) fn rust_treat_prism_edges(
        &self,
        edges: &[TreatmentEdge],
        treatment: EdgeTreatment,
        feature: EntityId,
        policy: GeometryPolicy,
        topology: &TopologySnapshot,
    ) -> Option<Result<Self, GeometryError>> {
        if self.policy != policy || policy.boolean_fuzzy_mm != 0. {
            return None;
        }
        match self.local_box_treatment_source() {
            Err(cause) => return Some(Err(cause)),
            Ok(Some((source, placement))) => {
                let local = match source.topology_snapshot(topology.owner, topology.revision) {
                    Ok(local) => local,
                    Err(cause) => return Some(Err(cause)),
                };
                return source
                    .rust_treat_prism_edges(edges, treatment, feature, policy, &local)
                    .map(|result| result.and_then(|solid| solid.transformed(placement)));
            }
            Ok(None) => {}
        }
        let views = match self.prism_views(topology) {
            Ok(views) => views,
            Err(cause) => return Some(Err(cause)),
        };
        let tolerance = policy.linear_tolerance_mm;
        // A box has a view along each axis. Prefer the view that treats the
        // most edges on its caps: two edges meeting at a box corner are then
        // both cap edges of their common face.
        let classified = views.iter().filter_map(|prism| {
            let picks = edges
                .iter()
                .map(|edge| {
                    let geometry = topology
                        .elements
                        .iter()
                        .find(|element| element.key == edge.key)?
                        .geometry?;
                    prism.pick(geometry, tolerance)
                })
                .collect::<Option<Vec<Pick>>>()?;
            let sides = picks.iter().filter(|p| matches!(p, Pick::Side(_))).count();
            Some((sides, prism, picks))
        });
        let (_, prism, picks) = classified.min_by_key(|(sides, ..)| *sides)?;
        Some(treat_prism(
            prism, edges, &picks, treatment, feature, policy,
        ))
    }
}

fn treat_prism(
    prism: &Prism,
    edges: &[TreatmentEdge],
    picks: &[Pick],
    treatment: EdgeTreatment,
    feature: EntityId,
    policy: GeometryPolicy,
) -> Result<Solid, GeometryError> {
    let (kind, size) = match treatment {
        EdgeTreatment::Fillet { radius_mm } => (Kind::Fillet, radius_mm),
        EdgeTreatment::Chamfer { distance_mm } => (Kind::Chamfer, distance_mm),
    };
    let n = prism.ring.len();
    let mut caps: [Vec<Option<Origin>>; 2] = [vec![None; n], vec![None; n]];
    let mut sides: Vec<Option<Origin>> = vec![None; n];
    for (edge, &pick) in edges.iter().zip(picks) {
        let slot = match pick {
            Pick::Cap(level, i) => &mut caps[level][i],
            Pick::Side(k) => &mut sides[k],
        };
        *slot = Some(edge.origin.clone());
    }
    let selected: Vec<bool> = sides.iter().map(Option::is_some).collect();
    let (ring, sources, names, cap_origins) = if selected.contains(&true) {
        let cut = match kind {
            Kind::Fillet => CornerCut::Round(size),
            Kind::Chamfer => CornerCut::Bevel(size),
        };
        let cut_ring = profile::cut_corners(&prism.ring, &selected, cut, policy.minimum_feature_mm)
            .map_err(|message| refused(capitalized(message)))?;
        let mut ring = Vec::with_capacity(cut_ring.len());
        let mut sources = Vec::with_capacity(cut_ring.len());
        let mut names = Vec::with_capacity(cut_ring.len());
        let mut cap_origins: [Vec<Option<Vec<Origin>>>; 2] = [vec![], vec![]];
        for (seg, source) in cut_ring {
            ring.push(seg);
            sources.push(source);
            match source {
                Source::Edge(i) => {
                    names.push(prism.side_names[i].clone());
                    for level in 0..2 {
                        cap_origins[level].push(caps[level][i].clone().map(|o| vec![o]));
                    }
                }
                Source::Corner(k) => {
                    let side = sides[k].clone().expect("selected corner");
                    names.push(FaceName {
                        origins: vec![side.clone()],
                        support: None,
                        anchor: None,
                    });
                    for level in 0..2 {
                        let around = [&caps[level][(k + n - 1) % n], &caps[level][k]];
                        cap_origins[level].push(match around {
                            [Some(a), Some(b)] => Some(vec![a.clone(), b.clone(), side.clone()]),
                            [None, None] => None,
                            _ if kind == Kind::Fillet => {
                                return Err(refused(
                                    "A rounded edge ends at a corner where only one of the two other edges is selected; select both of them or neither",
                                ));
                            }
                            _ => None,
                        });
                    }
                }
            }
        }
        (ring, sources, names, cap_origins)
    } else {
        (
            prism.ring.clone(),
            (0..n).map(Source::Edge).collect(),
            prism.side_names.clone(),
            caps.map(|level| level.into_iter().map(|o| o.map(|o| vec![o])).collect()),
        )
    };
    let any_cap = cap_origins.iter().flatten().any(Option::is_some);
    if !any_cap && let Some(sweep) = &prism.sweep {
        // Kept edges keep their authored identity; new corners take the
        // treatment origin.
        let origins = sources
            .iter()
            .map(|&source| match source {
                Source::Edge(i) => sweep.origins[i].clone(),
                Source::Corner(k) => sides[k].clone().expect("selected corner"),
            })
            .collect();
        return extrude_again(prism, sweep, &ring, origins, policy);
    }
    blend::build(
        &Plan {
            frame: prism.frame,
            ring,
            span: prism.span,
            kind,
            size,
            caps: cap_origins,
            cap_names: prism.cap_names.clone(),
            side_names: names,
        },
        feature,
        policy,
    )
}

fn extrude_again(
    prism: &Prism,
    sweep: &SweepSource,
    ring: &[Seg],
    mut origins: Vec<Origin>,
    policy: GeometryPolicy,
) -> Result<Solid, GeometryError> {
    let outer: Vec<CurvedEdge2d> = if sweep.authored_reversed {
        origins.reverse();
        ring.iter()
            .rev()
            .map(|&seg| curved(reversed(seg)))
            .collect()
    } else {
        ring.iter().map(|&seg| curved(seg)).collect()
    };
    let height = prism.span[1] - prism.span[0];
    let solid = Solid::rust_extrude_normal_profile_with_path(
        &CurvedProfile {
            frame: prism.frame,
            outer,
            holes: vec![],
        },
        if sweep.reverse { -height } else { height },
        sweep.feature,
        &origins,
        sweep.path_origin.clone(),
        policy,
    )?;
    Ok(solid.with_local_cap_anchors(sweep.cap_anchors))
}

fn capitalized(message: &str) -> String {
    let mut chars = message.chars();
    chars
        .next()
        .map(|first| first.to_uppercase().chain(chars).collect())
        .unwrap_or_default()
}
