//! Solid integration of the actual paired source and its finite annular boundary.
use super::internal_thread_treatment_body::interval::Interval as B;
use super::*;
use multicad_core::references::{
    AnalyticGeometry, AnchoredFrame, TopologyElement, TopologyKey, TopologyKind,
};
use multicad_rust_kernel::{
    authored_internal_treatment::{BoundaryConstraint, BoundaryPointRole},
    wood_insert::{self as source, BoundaryVertexKey, FaceKey},
};
use std::collections::{BTreeMap, BTreeSet};
use std::sync::OnceLock;
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WoodInsertData {
    pub source: source::Insert,
    pub frame: PlaneFrame,
}
pub(crate) struct WoodInsertBody {
    source: source::Insert,
    pub(super) feature: EntityId,
    boundary: OnceLock<source::Boundary>,
}
fn failure(cause: source::Error) -> GeometryError {
    error(
        "Rust wood insert solid",
        match cause {
            source::Error::Invalid => GeometryErrorKind::Input,
            source::Error::Limit => GeometryErrorKind::Limit,
            source::Error::Cancelled => GeometryErrorKind::Cancelled,
            source::Error::Unsupported | source::Error::Unresolved => {
                GeometryErrorKind::Unsupported
            }
        },
        format!("Actual insert source qualification failed: {cause:?}"),
    )
}
fn refused(message: impl Into<String>) -> GeometryError {
    error(
        "Rust wood insert topology",
        GeometryErrorKind::Unsupported,
        message,
    )
}
impl Solid {
    pub fn rust_wood_insert(
        source: source::Insert,
        feature: EntityId,
        frame: PlaneFrame,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, GeometryError> {
        if feature.is_nil() {
            return Err(input("Insert needs its original feature identity"));
        }
        WoodInsertMaterial::new(source, frame, policy, interrupt)?;
        let boundary = source
            .boundary_with_interrupt(source::BoundaryOptions::default(), interrupt)
            .map_err(failure)?;
        let body = WoodInsertBody {
            source,
            feature,
            boundary: OnceLock::from(boundary),
        };
        if interrupt.is_some_and(|f| f()) {
            return Err(failure(source::Error::Cancelled));
        }
        Ok(Self {
            body: SolidBody::WoodInsert(Rc::new(body), frame),
            policy,
            provenance: vec![],
            rust_origins: vec![],
            retained_placement: None,
        })
    }
    pub(crate) fn wood_insert_body(&self) -> Option<(&WoodInsertBody, PlaneFrame)> {
        match &self.body {
            SolidBody::WoodInsert(body, frame) => Some((body, *frame)),
            _ => None,
        }
    }
    pub fn rust_wood_insert_data(&self) -> Option<WoodInsertData> {
        self.wood_insert_body().map(|(body, frame)| WoodInsertData {
            source: body.source,
            frame,
        })
    }
}
fn side(key: BoundaryVertexKey) -> &'static str {
    match key {
        BoundaryVertexKey::Outer(_) => "outer",
        BoundaryVertexKey::Inner(_) => "inner",
    }
}
fn constraint(c: BoundaryConstraint) -> String {
    match c {
        BoundaryConstraint::Rail(i) => format!("r{i}"),
        BoundaryConstraint::Row(i) => format!("a{i}"),
        BoundaryConstraint::Station(i, o) => format!("z{i}.{o}"),
    }
}
fn point_tag(role: BoundaryPointRole) -> String {
    let key = |k: multicad_rust_kernel::authored_internal_treatment::BoundaryPointKey| {
        let mut c = k.0.map(constraint);
        c.sort();
        format!("{}-{}", c[0], c[1])
    };
    match role {
        BoundaryPointRole::Source(k) => format!("p{}", key(k)),
        BoundaryPointRole::Tail { original, end } => format!("t{end}.{}", key(original)),
        BoundaryPointRole::StockRim { end } => format!("s{end}"),
    }
}
fn face_tag(key: FaceKey) -> String {
    use multicad_rust_kernel::authored_internal_treatment::ThreadBandRole;
    match key {
        FaceKey::StartCap => "cap.start".into(),
        FaceKey::EndCap => "cap.end".into(),
        FaceKey::OuterThread { start, band } | FaceKey::InnerThread { start, band } => format!(
            "{}/start/{start}/{}",
            if matches!(key, FaceKey::OuterThread { .. }) {
                "outer"
            } else {
                "inner"
            },
            match band {
                ThreadBandRole::Lower => "lower-flank",
                ThreadBandRole::Crest => "crest",
                ThreadBandRole::Upper => "upper-flank",
                ThreadBandRole::Root => "root",
            }
        ),
    }
}
fn origin(feature: EntityId, role: String) -> Origin {
    Origin {
        source: feature,
        role: format!("wood-insert/{role}"),
    }
}
fn world_point(
    frame: PlaneFrame,
    point: [B; 3],
    policy: GeometryPolicy,
) -> Result<[f64; 3], GeometryError> {
    let axes = [frame.u, frame.v, frame.normal()];
    let p: [B; 3] = std::array::from_fn(|i| {
        (0..3).fold(B::point(frame.origin[i]), |s, k| {
            s + B::point(axes[k][i]) * point[k]
        })
    });
    if p.into_iter().any(|v| !v.valid()) {
        return Err(refused("Source vertex placement exceeds numerical range"));
    }
    let center = p.map(|v| v.lo.midpoint(v.hi));
    let error = p.into_iter().zip(center).fold(B::ZERO, |s, (b, x)| {
        s + B::point((b - B::point(x)).magnitude())
    });
    if error.hi > policy.linear_tolerance_mm / 8. {
        return Err(refused(
            "Source vertex exceeds requested placement precision",
        ));
    }
    policy.point(center)?;
    Ok(center)
}
impl WoodInsertBody {
    pub(super) fn material(
        &self,
        frame: PlaneFrame,
        policy: GeometryPolicy,
    ) -> Result<WoodInsertMaterial, GeometryError> {
        WoodInsertMaterial::new(self.source, frame, policy, None)
    }
    pub(super) fn boundary(&self) -> &source::Boundary {
        self.boundary
            .get()
            .expect("immutable source graph qualified before Solid construction")
    }
    pub(crate) fn topology(
        &self,
        frame: PlaneFrame,
        owner: EntityId,
        revision: u64,
        policy: GeometryPolicy,
    ) -> Result<TopologySnapshot, GeometryError> {
        self.material(frame, policy)?;
        if owner.is_nil() {
            return Err(input("Insert topology needs its owner identity"));
        }
        let b = self.boundary();
        let (nv, ne) = b.topology_counts();
        let nf = self.source.face_count() as usize;
        if nv + ne + nf > 4096 {
            return Err(refused(
                "Insert reference handoff exceeds its descriptor budget",
            ));
        }
        let vertices: BTreeMap<_, _> = b
            .vertices()
            .map(|(k, p)| {
                (
                    k,
                    origin(
                        self.feature,
                        format!("{}/vertex/{}", side(k), point_tag(p.role())),
                    ),
                )
            })
            .collect();
        if vertices.values().collect::<BTreeSet<_>>().len() != nv {
            return Err(refused(
                "Distinct source vertices have ambiguous reference identities",
            ));
        }
        let mut edges = BTreeMap::new();
        let mut used = BTreeSet::new();
        for (key, _) in b.curves() {
            let ends = b.curve_endpoints(key).expect("qualified curve endpoints");
            let mut tags = ends.map(|k| vertices[&k].role.clone());
            tags.sort();
            let mut faces = b
                .curve_faces(key)
                .expect("qualified source uses")
                .map(face_tag);
            faces.sort();
            let role = format!("edge/{}.{}/{}|{}", faces[0], faces[1], tags[0], tags[1]);
            let id = origin(self.feature, role);
            if id.role.len() > 256 || !used.insert(id.clone()) {
                return Err(refused(
                    "Finite source curves need distinct semantic reference roles",
                ));
            }
            edges.insert(key, id);
        }
        let mut elements = Vec::with_capacity(nv + ne + nf);
        for index in 0..nf {
            let face = self
                .source
                .face_key(index as u32)
                .expect("source face role");
            let id = origin(self.feature, face_tag(face));
            let mut boundary = BTreeSet::new();
            for edge in b
                .face_loops(face)
                .expect("qualified face loops")
                .iter()
                .flatten()
            {
                boundary.insert(edges[&edge.curve].clone());
            }
            let plane = if matches!(face, FaceKey::StartCap | FaceKey::EndCap) {
                let z = if face == FaceKey::StartCap {
                    0.
                } else {
                    self.source.full_length_mm()
                };
                let support = ReferenceFrame {
                    origin_mm: world_point(frame, [B::ZERO, B::ZERO, B::point(z)], policy)?,
                    u: frame.u,
                    v: if face == FaceKey::StartCap {
                        frame.v.map(|x| -x)
                    } else {
                        frame.v
                    },
                };
                support.validate().map_err(refused)?;
                Some(support)
            } else {
                None
            };
            elements.push(TopologyElement {
                key: TopologyKey {
                    kind: TopologyKind::Face,
                    ordinal: index as u32,
                },
                origins: BTreeSet::from([id.clone()]),
                boundary,
                frames: plane
                    .map(|p| AnchoredFrame {
                        origin: id,
                        frame: p,
                    })
                    .into_iter()
                    .collect(),
                geometry: plane.map(|frame| AnalyticGeometry::Plane { frame }),
            });
        }
        let mut corners: BTreeMap<BoundaryVertexKey, BTreeSet<Origin>> = BTreeMap::new();
        for (index, (key, _)) in b.curves().enumerate() {
            let id = edges[&key].clone();
            for p in b.curve_endpoints(key).expect("source endpoints") {
                corners.entry(p).or_default().insert(id.clone());
            }
            let boundary = b
                .curve_faces(key)
                .expect("source face incidence")
                .map(|k| origin(self.feature, face_tag(k)))
                .into_iter()
                .collect();
            elements.push(TopologyElement {
                key: TopologyKey {
                    kind: TopologyKind::Edge,
                    ordinal: index as u32,
                },
                origins: BTreeSet::from([id]),
                boundary,
                frames: vec![],
                geometry: None,
            });
        }
        for (index, (key, p)) in b.vertices().enumerate() {
            elements.push(TopologyElement {
                key: TopologyKey {
                    kind: TopologyKind::Vertex,
                    ordinal: index as u32,
                },
                origins: BTreeSet::from([vertices[&key].clone()]),
                boundary: corners
                    .remove(&key)
                    .expect("qualified physical source vertex"),
                frames: vec![],
                geometry: Some(AnalyticGeometry::Point {
                    position_mm: world_point(frame, p.position_bounds().map(B::from), policy)?,
                }),
            });
        }
        let snapshot = TopologySnapshot {
            owner,
            revision,
            elements,
        };
        snapshot.validate().map_err(refused)?;
        Ok(snapshot)
    }
    pub(super) fn display_edges(
        &self,
        frame: PlaneFrame,
        deflection: f64,
        policy: GeometryPolicy,
    ) -> Result<Vec<DisplayEdge>, GeometryError> {
        use crate::geometry::arithmetic::{dot_error, sum_error};
        self.material(frame, policy)?;
        if !deflection.is_finite() || deflection <= 0. {
            return Err(input("Invalid insert curve display quality"));
        }
        let b = self.boundary();
        let mut remaining = 500_000usize;
        let mut output = Vec::with_capacity(b.topology_counts().1);
        let axes = [frame.u, frame.v, frame.normal()];
        for (index, (key, _)) in b.curves().enumerate() {
            if remaining < 2 {
                return Err(failure(source::Error::Limit));
            }
            // The validated near-unit U,V and their stored cross product have
            // Frobenius norm below3. Half the chord allowance remains for
            // actual emitted-point rounding in this raw affine placement.
            let mut line = b
                .curve_polyline(key, deflection / 6., remaining, None)
                .map_err(failure)?;
            remaining -= line.points_mm.len();
            let mut rounding = 0_f64;
            for point in &mut line.points_mm {
                let mut e = 0_f64;
                *point = std::array::from_fn(|i| {
                    let (v, d) = dot_error(axes.map(|a| a[i]), *point);
                    let (v, s) = sum_error(frame.origin[i], v);
                    e = (e + (d + s).next_up()).next_up();
                    v
                });
                policy.point(*point)?;
                rounding = rounding.max(e);
            }
            if (B::point(3.) * B::point(line.maximum_error_mm) + B::point(rounding)).hi > deflection
            {
                return Err(refused(
                    "Placed source curve exceeds its requested chord bound",
                ));
            }
            output.push(DisplayEdge {
                index: index as u32,
                points_mm: line.points_mm,
            });
        }
        Ok(output)
    }
}

#[cfg(test)]
mod tests;
