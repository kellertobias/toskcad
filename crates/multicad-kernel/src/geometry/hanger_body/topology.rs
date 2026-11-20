//! References come from the sealed finite source graph, never mesh ordinals.
use super::*;
use crate::geometry::internal_thread_treatment_body::interval::Interval as I;
use multicad_core::references::{
    AnalyticGeometry, AnchoredFrame, TopologyElement, TopologyKey, TopologyKind,
};
use source::{Component, CurveKey, End, FaceKey, Knot, VertexKey};
use std::collections::{BTreeMap, BTreeSet};

fn component(c: Component) -> &'static str {
    match c {
        Component::Wood => "wood",
        Component::Metric => "metric",
    }
}
fn end(e: End) -> &'static str {
    match e {
        End::Start => "start",
        End::End => "end",
    }
}
fn knot(k: Knot) -> &'static str {
    match k {
        Knot::CrestEnd => "crest-end",
        Knot::UpperRoot => "upper-root",
        Knot::LowerRoot => "lower-root",
        Knot::CrestStart => "crest-start",
    }
}
fn band(b: multicad_rust_kernel::authored_internal_treatment::ThreadBandRole) -> &'static str {
    use multicad_rust_kernel::authored_internal_treatment::ThreadBandRole::*;
    match b {
        Lower => "lower",
        Crest => "crest",
        Upper => "upper",
        Root => "root",
    }
}
fn face_tag(k: FaceKey) -> String {
    match k {
        FaceKey::MetricTip => "metric-tip".into(),
        FaceKey::MetricShoulder => "metric-shoulder".into(),
        FaceKey::WoodShoulder => "wood-shoulder".into(),
        FaceKey::WoodTip => "wood-tip".into(),
        FaceKey::PlainWall => "plain-wall".into(),
        FaceKey::MetricThread { start, band: b } => format!("metric/{start}/{}", band(b)),
        FaceKey::WoodThread { start, band: b } => format!("wood/{start}/{}", band(b)),
    }
}
fn vertex_tag(k: VertexKey) -> String {
    match k {
        VertexKey::Plain { end: e } => format!("plain/{}", end(e)),
        VertexKey::Thread {
            component: c,
            end: e,
            period,
            knot: k,
        } => format!("{}/{}/{period}/{}", component(c), end(e), knot(k)),
    }
}
fn curve_tag(k: CurveKey) -> String {
    match k {
        CurveKey::PlainCircle { end: e } => format!("plain-circle/{}", end(e)),
        CurveKey::Rail {
            component: c,
            period,
            knot: k,
        } => format!("{}/rail/{period}/{}", component(c), knot(k)),
        CurveKey::Cap {
            component: c,
            end: e,
            period,
            band: b,
        } => format!("{}/cap/{}/{period}/{}", component(c), end(e), band(b)),
    }
}
fn origin(feature: EntityId, role: String) -> Origin {
    Origin {
        source: feature,
        role: format!("hanger/{role}"),
    }
}
/// Every bounded helical band also carries the shared `thread/` role that
/// every other modeled thread publishes, qualified by its component.
fn face_origins(feature: EntityId, key: FaceKey) -> BTreeSet<Origin> {
    let mut origins = BTreeSet::from([origin(feature, face_tag(key))]);
    if matches!(
        key,
        FaceKey::MetricThread { .. } | FaceKey::WoodThread { .. }
    ) {
        origins.insert(Origin {
            source: feature,
            role: format!("thread/{}", face_tag(key)),
        });
    }
    origins
}
pub(super) fn refuse(message: impl Into<String>) -> GeometryError {
    error(
        "Rust hanger references",
        GeometryErrorKind::Unsupported,
        message,
    )
}
pub(super) fn world_point(
    frame: PlaneFrame,
    point: [[f64; 2]; 3],
    maximum_error: f64,
    policy: GeometryPolicy,
) -> Result<[f64; 3], GeometryError> {
    let p = point.map(|[lo, hi]| I::new(lo, hi));
    let axes = [frame.u, frame.v, frame.normal()];
    let world: [I; 3] = std::array::from_fn(|i| {
        (0..3).fold(I::point(frame.origin[i]), |s, k| {
            s + I::point(axes[k][i]) * p[k]
        })
    });
    let center = world.map(|x| x.lo.midpoint(x.hi));
    let error = world
        .into_iter()
        .zip(center)
        .fold(I::ZERO, |sum, (bound, x)| {
            sum + I::point((bound - I::point(x)).magnitude())
        });
    if error.hi > maximum_error {
        return Err(refuse(
            "Original source point exceeds requested placement precision",
        ));
    }
    policy.point(center)?;
    Ok(center)
}
impl HangerBody {
    pub(crate) fn topology(
        &self,
        frame: PlaneFrame,
        owner: EntityId,
        revision: u64,
        policy: GeometryPolicy,
    ) -> Result<TopologySnapshot, GeometryError> {
        self.material(frame, policy)?;
        if owner.is_nil() {
            return Err(input("Hanger topology needs its owner identity"));
        }
        let b = &self.boundary;
        let vertices: BTreeMap<_, _> = b
            .vertices()
            .iter()
            .map(|&k| (k, origin(self.feature, format!("vertex/{}", vertex_tag(k)))))
            .collect();
        let edges: BTreeMap<_, _> = b
            .curves()
            .iter()
            .map(|c| {
                (
                    c.key,
                    origin(self.feature, format!("edge/{}", curve_tag(c.key))),
                )
            })
            .collect();
        let mut elements =
            Vec::with_capacity(vertices.len() + edges.len() + self.source.face_count());
        for index in 0..self.source.face_count() {
            let key = self
                .source
                .face_key(index)
                .expect("actual source face index");
            let id = origin(self.feature, face_tag(key));
            let boundary = b
                .face_loops(key)
                .expect("sealed source face")
                .into_iter()
                .flatten()
                .map(|e| edges[&e.curve].clone())
                .collect();
            let station = match key {
                FaceKey::MetricTip => Some((
                    -(I::point(self.source.metric().full_length_mm())
                        + I::point(self.source.spec().plain_length_mm)),
                    false,
                )),
                FaceKey::MetricShoulder => {
                    Some((-I::point(self.source.spec().plain_length_mm), false))
                }
                FaceKey::WoodShoulder => Some((I::ZERO, false)),
                FaceKey::WoodTip => Some((I::point(self.source.wood().full_length_mm()), true)),
                _ => None,
            };
            let plane = station
                .map(|(z, up)| -> Result<ReferenceFrame, GeometryError> {
                    let frame = ReferenceFrame {
                        origin_mm: world_point(
                            frame,
                            [[0., 0.], [0., 0.], [z.lo, z.hi]],
                            policy.linear_tolerance_mm / 8.,
                            policy,
                        )?,
                        u: frame.u,
                        v: if up { frame.v } else { frame.v.map(|x| -x) },
                    };
                    frame.validate().map_err(refuse)?;
                    Ok(frame)
                })
                .transpose()?;
            elements.push(TopologyElement {
                key: TopologyKey {
                    kind: TopologyKind::Face,
                    ordinal: index as u32,
                },
                origins: face_origins(self.feature, key),
                boundary,
                frames: plane
                    .map(|frame| AnchoredFrame { origin: id, frame })
                    .into_iter()
                    .collect(),
                geometry: plane.map(|frame| AnalyticGeometry::Plane { frame }),
            });
        }
        let mut corner_edges: BTreeMap<VertexKey, BTreeSet<Origin>> = BTreeMap::new();
        for (index, c) in b.curves().iter().enumerate() {
            let id = edges[&c.key].clone();
            for v in c.vertices {
                corner_edges.entry(v).or_default().insert(id.clone());
            }
            elements.push(TopologyElement {
                key: TopologyKey {
                    kind: TopologyKind::Edge,
                    ordinal: index as u32,
                },
                origins: BTreeSet::from([id]),
                boundary: c
                    .faces
                    .map(|f| origin(self.feature, face_tag(f)))
                    .into_iter()
                    .collect(),
                frames: vec![],
                geometry: None,
            });
        }
        for (index, &k) in b.vertices().iter().enumerate() {
            let point = b
                .vertex_bounds(k, None)
                .map_err(super::super::hanger_material::failure)?;
            elements.push(TopologyElement {
                key: TopologyKey {
                    kind: TopologyKind::Vertex,
                    ordinal: index as u32,
                },
                origins: BTreeSet::from([vertices[&k].clone()]),
                boundary: corner_edges
                    .remove(&k)
                    .expect("sealed source vertex incidence"),
                frames: vec![],
                geometry: Some(AnalyticGeometry::Point {
                    position_mm: world_point(
                        frame,
                        point,
                        policy.linear_tolerance_mm / 8.,
                        policy,
                    )?,
                }),
            });
        }
        let snapshot = TopologySnapshot {
            owner,
            revision,
            elements,
        };
        snapshot.validate().map_err(refuse)?;
        Ok(snapshot)
    }
}
