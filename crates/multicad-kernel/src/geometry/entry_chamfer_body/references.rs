//! One logical unchanged source rim, never a duplicate alias on fragments.
use super::*;
use source::{Face as Role, PhaseConstraint, SourceCurve};
use std::collections::BTreeMap;
#[derive(Clone, Debug)]
pub struct EntryChamferEdge {
    pub(super) pieces: Vec<usize>,
    pub(super) original: Option<(u32, usize)>,
    pub(super) rail: Option<String>,
}
impl EntryChamferEdge {
    pub fn source_edges(&self) -> &[usize] {
        &self.pieces
    }
}
pub(super) struct ReferenceView {
    pub edges: Vec<EntryChamferEdge>,
    pub raw_to_edge: Vec<usize>,
    pub thread_faces: BTreeMap<(u32, u32), String>,
    pub corners: BTreeMap<usize, String>,
}
fn refuse() -> GeometryError {
    error(
        "Entry source references",
        GeometryErrorKind::Unsupported,
        "Original finite source rim ownership is not uniquely proved",
    )
}
impl ReferenceView {
    pub fn new(body: &Chamfer, graph: &BoundaryTopology) -> Result<Self, GeometryError> {
        let tool = body.source_thread();
        let mut groups = BTreeMap::<String, Vec<usize>>::new();
        let mut corners = BTreeMap::new();
        let cap = if body.terminal_relief().is_some() {
            "thread.end-relief/shoulder"
        } else {
            "cap.end"
        };
        for (i, edge) in graph.edges().iter().enumerate() {
            if let SourceCurve::EndProfile { start, band, phase } = edge.source() {
                groups
                    .entry(helical_body::phase_role(tool, start, band))
                    .or_default()
                    .push(i);
                for (vertex, phase) in edge.vertices().into_iter().zip(phase) {
                    if let PhaseConstraint::ProfileKnot(k @ 1..=4) = phase.constraint() {
                        let role =
                            format!("{cap}/corner/{}", helical_body::phase_role(tool, start, k));
                        if corners
                            .insert(vertex, role.clone())
                            .is_some_and(|old| old != role)
                        {
                            return Err(refuse());
                        }
                    }
                }
            }
        }
        let mut thread_faces = BTreeMap::new();
        for face in graph.faces() {
            if let Role::ThreadBand { start, band } = face.role() {
                let roles = face
                    .loops()
                    .iter()
                    .flatten()
                    .filter_map(|u| {
                        if let SourceCurve::EndProfile { start, band, .. } =
                            graph.edges()[u.edge].source()
                        {
                            Some(helical_body::phase_role(tool, start, band))
                        } else {
                            None
                        }
                    })
                    .collect::<BTreeSet<_>>();
                let mut roles = roles.into_iter();
                let role = roles.next().ok_or_else(refuse)?;
                if roles.next().is_some() {
                    return Err(refuse());
                }
                thread_faces.insert((start, band), role);
            }
        }
        let mut edges = Vec::new();
        let mut raw_to_edge = vec![usize::MAX; graph.edges().len()];
        for i in 0..graph.edges().len() {
            if raw_to_edge[i] != usize::MAX {
                continue;
            }
            let (pieces, original) = if let SourceCurve::EndProfile { start, band, .. } =
                graph.edges()[i].source()
            {
                let group = groups
                    .remove(&helical_body::phase_role(tool, start, band))
                    .ok_or_else(refuse)?;
                let starts=group.iter().copied().filter(|&id|{
                    if let SourceCurve::EndProfile {band,phase,..}=graph.edges()[id].source(){matches!(phase[0].constraint(),PhaseConstraint::ProfileKnot(k@1..=4)if k==band)}else{false}
                }).collect::<Vec<_>>();
                let [first] = starts.as_slice() else {
                    return Err(refuse());
                };
                let SourceCurve::EndProfile { start, band, .. } = graph.edges()[*first].source()
                else {
                    unreachable!()
                };
                let mut remaining = group.into_iter().collect::<BTreeSet<_>>();
                let mut ordered = vec![];
                let mut next = *first;
                loop {
                    if !remaining.remove(&next) {
                        return Err(refuse());
                    }
                    ordered.push(next);
                    if remaining.is_empty() {
                        break;
                    }
                    let end = graph.edges()[next].vertices()[1];
                    let candidates = remaining
                        .iter()
                        .copied()
                        .filter(|&id| graph.edges()[id].vertices()[0] == end)
                        .collect::<Vec<_>>();
                    let [id] = candidates.as_slice() else {
                        return Err(refuse());
                    };
                    next = *id;
                }
                let SourceCurve::EndProfile { phase, .. } =
                    graph.edges()[*ordered.last().expect("nonempty source chain")].source()
                else {
                    unreachable!()
                };
                let expected = if band == 4 { 1 } else { band + 1 };
                if phase[1].constraint() != PhaseConstraint::ProfileKnot(expected) {
                    return Err(refuse());
                }
                (ordered, Some((start, band)))
            } else {
                (vec![i], None)
            };
            let id = edges.len();
            for &raw in &pieces {
                raw_to_edge[raw] = id;
            }
            let mut rail = None;
            if let SourceCurve::WallRail { start, phase, .. } = graph.edges()[i].source()
                && let PhaseConstraint::ProfileKnot(k @ 1..=4) = phase.constraint()
            {
                let bounds = graph.curve_bounds(i, [0., 0.]).map_err(failure)?;
                let [x, y] = [0, 1].map(|k| bounds[k].lower.abs().max(bounds[k].upper.abs()));
                let radius = ((x * x).next_up() + (y * y).next_up())
                    .next_up()
                    .sqrt()
                    .next_up();
                // At this fixed original phase, the complete helical rail is
                // unchanged iff its radius fits the entry cone already at0.
                if radius <= body.entry_radius_mm() {
                    rail = Some(format!(
                        "thread/{}/leading_rail",
                        helical_body::phase_role(tool, start, k)
                    ));
                }
            }
            edges.push(EntryChamferEdge {
                pieces,
                original,
                rail,
            });
        }
        Ok(Self {
            edges,
            raw_to_edge,
            thread_faces,
            corners,
        })
    }
}
