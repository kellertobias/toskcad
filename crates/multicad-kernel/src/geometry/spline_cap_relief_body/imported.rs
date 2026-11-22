//! Physical per-face identity for opaque retained spline relief imports.
use super::*;
fn numbers(values: impl IntoIterator<Item = f64>) -> Vec<u64> {
    values
        .into_iter()
        .map(|v| if v == 0. { 0 } else { v.to_bits() })
        .collect()
}
fn curve_key(curve: &Curve, forward: bool, frame: PlaneFrame) -> Vec<u64> {
    let mut key = vec![u64::from(forward)];
    match curve {
        Curve::Line { a_mm, b_mm } => {
            key.push(0);
            key.extend(numbers(a_mm.iter().chain(b_mm).copied()));
        }
        Curve::Profile { edge, z_mm } => {
            key.push(1);
            key.extend(numbers([*z_mm]));
            match edge {
                Edge::Analytic(AnalyticEdge::Line { a, b }) => {
                    key.push(0);
                    key.extend(numbers(a.iter().chain(b).copied()));
                }
                Edge::Analytic(AnalyticEdge::Circle {
                    center,
                    radius_mm,
                    clockwise,
                }) => {
                    key.push(1);
                    key.push(u64::from(*clockwise));
                    key.extend(numbers(center.iter().copied().chain([*radius_mm])));
                }
                Edge::Analytic(AnalyticEdge::Arc {
                    a,
                    b,
                    center,
                    clockwise,
                }) => {
                    key.push(2);
                    key.push(u64::from(*clockwise));
                    key.extend(numbers(a.iter().chain(b).chain(center).copied()));
                }
                Edge::Rational(c) => {
                    key.push(4);
                    key.push(u64::from(c.reversed()));
                    key.push(c.poles().len() as u64);
                    key.extend(numbers(c.poles().iter().flatten().copied()));
                    key.extend(numbers(c.weights().iter().copied()));
                }
                Edge::Cubic(c) => {
                    key.push(3);
                    key.push(u64::from(c.reversed()));
                    key.push(c.poles().len() as u64);
                    key.extend(numbers(c.poles().iter().flatten().copied()));
                    key.push(c.knots().len() as u64);
                    key.extend(numbers(c.knots().iter().copied()));
                }
            }
        }
    }
    // Literal source coefficients plus complete placement, rather than a
    // rounded world-only sample, distinguish changed opaque finite surfaces.
    key.extend(numbers(
        frame.origin.into_iter().chain(frame.u).chain(frame.v),
    ));
    key
}
fn face_stamp(geometry: &Relief, index: usize, frame: PlaneFrame) -> ContentHash {
    let graph = geometry.boundary();
    let face = &graph.faces[index];
    let mut wires = face
        .loops
        .iter()
        .map(|ring| {
            let mut curves = ring
                .iter()
                .map(|e| curve_key(&graph.edges[e.edge].curve, e.forward, frame))
                .collect::<Vec<_>>();
            if let Some((i, _)) = curves.iter().enumerate().min_by_key(|(_, s)| *s) {
                curves.rotate_left(i);
            }
            curves
        })
        .collect::<Vec<_>>();
    let outside = wires.remove(0);
    wires.sort();
    let mut values = match &face.surface {
        Surface::Plane { z_mm, upward } => {
            let mut v = vec![0, u64::from(*upward)];
            v.extend(numbers([*z_mm]));
            v
        }
        Surface::Extrusion {
            edge,
            span_mm,
            inward,
        } => {
            let mut v = vec![1, u64::from(*inward)];
            v.extend(numbers(*span_mm));
            v.extend(curve_key(
                &Curve::Profile {
                    edge: edge.clone(),
                    z_mm: span_mm[0],
                },
                true,
                frame,
            ));
            v
        }
    };
    for ring in std::iter::once(outside).chain(wires) {
        values.push(ring.len() as u64);
        for curve in ring {
            values.push(curve.len() as u64);
            values.extend(curve);
        }
    }
    ContentHash::of(
        &values
            .into_iter()
            .flat_map(u64::to_le_bytes)
            .collect::<Vec<_>>(),
    )
}

impl Solid {
    pub(in crate::geometry) fn from_spline_cap_relief_import(
        geometry: Relief,
        frame: PlaneFrame,
        feature: EntityId,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        let labels = (0..geometry.boundary().faces.len())
            .map(|i| Labels {
                origins: BTreeSet::from([Origin {
                    source: feature,
                    role: format!("imported-spline-cap.v1/{}", face_stamp(&geometry, i, frame)),
                }]),
                anchors: vec![],
            })
            .collect();
        let body = SplineCapReliefBody::new(
            geometry,
            feature,
            labels,
            policy,
            source_precision::Receipt::default(),
        )?;
        frame.checked_metric_bounds(body.bounds(PlaneFrame::default()), policy)?;
        let bounds = body.bounds(frame);
        policy.point(bounds.min)?;
        policy.point(bounds.max)?;
        Ok(Self {
            body: SolidBody::SplineCapRelief(Rc::new(body), frame),
            policy,
            provenance: vec![],
            rust_origins: vec![],
            retained_placement: None,
        })
    }
}

pub(in crate::geometry) fn signature(geometry: &Relief, frame: PlaneFrame) -> ContentHash {
    let mut faces = (0..geometry.boundary().faces.len())
        .map(|i| face_stamp(geometry, i, frame).to_string())
        .collect::<Vec<_>>();
    faces.sort();
    ContentHash::of(faces.join("/").as_bytes())
}
