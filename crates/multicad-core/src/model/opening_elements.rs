//! The 3D door and window elements shown inside wall openings (CAD-429).
//! They are derived from the same resolved openings as the plan symbols
//! (`ResolvedWall::plan_symbols_at`): a door's frame and its leaf standing
//! open on the swing face at the hinge jamb, exactly where the plan draws
//! the leaf line; a window's frame on the two plan pane lines and its
//! glazing between them. These are display elements, not manufacturing
//! geometry: the wall solid and its opening cut are unchanged.
use super::wall::{Hinge, OpeningKind, ResolvedOpening, ResolvedWall, Swing, WallPlan};
use crate::document::EntityId;

/// The largest frame member width and the door leaf thickness, in mm.
pub const OPENING_FRAME_MM: f64 = 50.;
pub const DOOR_LEAF_MM: f64 = 40.;
/// Glazing thickness in mm.
pub const GLAZING_MM: f64 = 8.;

/// What one element piece represents.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum OpeningPieceRole {
    Frame,
    Leaf,
    Glazing,
}
impl OpeningPieceRole {
    pub fn label(self) -> &'static str {
        match self {
            Self::Frame => "frame",
            Self::Leaf => "leaf",
            Self::Glazing => "glazing",
        }
    }
}

/// One box-shaped piece in world millimetres. Corner `i` lies at the high
/// end of the segment direction when bit 0 is set, toward the wall's left
/// normal when bit 1 is set and at the top when bit 2 is set.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OpeningPiece {
    pub role: OpeningPieceRole,
    pub corners: [[f64; 3]; 8],
}

/// The display element of one door or window.
#[derive(Clone, Debug, PartialEq)]
pub struct OpeningElement {
    pub opening: EntityId,
    pub kind: OpeningKind,
    pub pieces: Vec<OpeningPiece>,
}

impl ResolvedWall {
    /// The display elements of every door and window, in source order.
    pub fn opening_elements(&self, plan: &WallPlan) -> Vec<OpeningElement> {
        self.openings
            .iter()
            .filter_map(|o| self.opening_element(plan, o))
            .collect()
    }
    fn opening_element(&self, plan: &WallPlan, o: &ResolvedOpening) -> Option<OpeningElement> {
        let seg = *plan.segments.get(o.segment)?;
        let base = self.base_z_mm;
        let piece = |role, along: [f64; 2], across: [f64; 2], up: [f64; 2]| OpeningPiece {
            role,
            corners: std::array::from_fn(|i| {
                let (a, c, z) = (along[i & 1], across[(i >> 1) & 1], up[(i >> 2) & 1]);
                [
                    seg.start[0] + seg.direction[0] * a + seg.normal[0] * c,
                    seg.start[1] + seg.direction[1] * a + seg.normal[1] * c,
                    base + z,
                ]
            }),
        };
        let (a, b) = (o.offset_mm, o.offset_mm + o.width_mm);
        let (s, h) = (o.sill_mm, o.sill_mm + o.height_mm);
        let f = OPENING_FRAME_MM.min(o.width_mm / 6.).min(o.height_mm / 6.);
        if !(f > 0.) {
            return None;
        }
        let (right, left) = (plan.right_mm, plan.left_mm);
        let mut pieces = vec![];
        match &o.kind {
            OpeningKind::Door { hinge, swing } => {
                // The frame lines the jambs and head through the wall.
                let across = [right, left];
                pieces.push(piece(OpeningPieceRole::Frame, [a, a + f], across, [s, h]));
                pieces.push(piece(OpeningPieceRole::Frame, [b - f, b], across, [s, h]));
                pieces.push(piece(
                    OpeningPieceRole::Frame,
                    [a + f, b - f],
                    across,
                    [h - f, h],
                ));
                // The leaf stands open a quarter turn on the swing face from
                // the hinge jamb, as long as the plan leaf line.
                let t = DOOR_LEAF_MM.min(o.width_mm / 4.);
                let along = match hinge {
                    Hinge::Start => [a, a + t],
                    Hinge::End => [b - t, b],
                };
                let across = match swing {
                    Swing::Left => [left, left + o.width_mm],
                    Swing::Right => [right - o.width_mm, right],
                };
                pieces.push(piece(OpeningPieceRole::Leaf, along, across, [s, h - f]));
            }
            OpeningKind::Window => {
                // The frame sits between the two plan pane lines.
                let mid = (left + right) / 2.;
                let gap = (left - right).abs() / 6.;
                let across = [mid - gap, mid + gap];
                pieces.push(piece(OpeningPieceRole::Frame, [a, a + f], across, [s, h]));
                pieces.push(piece(OpeningPieceRole::Frame, [b - f, b], across, [s, h]));
                pieces.push(piece(
                    OpeningPieceRole::Frame,
                    [a + f, b - f],
                    across,
                    [h - f, h],
                ));
                pieces.push(piece(
                    OpeningPieceRole::Frame,
                    [a + f, b - f],
                    across,
                    [s, s + f],
                ));
                let g = GLAZING_MM.min(gap) / 2.;
                pieces.push(piece(
                    OpeningPieceRole::Glazing,
                    [a + f, b - f],
                    [mid - g, mid + g],
                    [s + f, h - f],
                ));
            }
        }
        Some(OpeningElement {
            opening: o.id,
            kind: o.kind.clone(),
            pieces,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::super::wall::Justification;
    use super::*;

    fn wall(points: &[[f64; 2]], openings: Vec<ResolvedOpening>) -> ResolvedWall {
        ResolvedWall {
            junctions: vec![],
            layers: None,
            points: points.to_vec(),
            closed: false,
            base_z_mm: 100.,
            thickness_mm: 200.,
            height_mm: 2500.,
            justification: Justification::Center,
            openings,
            cutouts: vec![],
        }
    }
    fn door(hinge: Hinge, swing: Swing) -> ResolvedOpening {
        ResolvedOpening {
            id: EntityId::new(),
            kind: OpeningKind::Door { hinge, swing },
            segment: 0,
            offset_mm: 1000.,
            width_mm: 900.,
            height_mm: 2100.,
            sill_mm: 0.,
        }
    }
    fn bounds(p: &OpeningPiece) -> ([f64; 3], [f64; 3]) {
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for c in p.corners {
            for i in 0..3 {
                lo[i] = lo[i].min(c[i]);
                hi[i] = hi[i].max(c[i]);
            }
        }
        (lo, hi)
    }
    fn close(a: [f64; 2], b: [f64; 2]) -> bool {
        (a[0] - b[0]).abs() < 1e-9 && (a[1] - b[1]).abs() < 1e-9
    }

    #[test]
    fn door_frame_lines_the_opening_and_the_leaf_stands_on_the_plan_leaf_line() {
        for (hinge, swing) in [
            (Hinge::Start, Swing::Left),
            (Hinge::End, Swing::Left),
            (Hinge::Start, Swing::Right),
            (Hinge::End, Swing::Right),
        ] {
            // A wall running along +Y, so along/across are not world X/Y.
            let w = wall(&[[0., 0.], [0., 6000.]], vec![door(hinge, swing)]);
            let plan = w.plan().unwrap();
            let elements = w.opening_elements(&plan);
            assert_eq!(elements.len(), 1);
            let e = &elements[0];
            assert_eq!(e.opening, w.openings[0].id);
            let frames: Vec<_> = e
                .pieces
                .iter()
                .filter(|p| p.role == OpeningPieceRole::Frame)
                .map(bounds)
                .collect();
            assert_eq!(frames.len(), 3);
            // Jambs and head stay inside the cut: y 1000–1900, through the
            // 200 mm wall (x −100…100), z 100…2200 above the wall base.
            for (lo, hi) in &frames {
                assert!(lo[0] >= -100. - 1e-9 && hi[0] <= 100. + 1e-9);
                assert!(lo[1] >= 1000. - 1e-9 && hi[1] <= 1900. + 1e-9);
                assert!(lo[2] >= 100. - 1e-9 && hi[2] <= 2200. + 1e-9);
            }
            assert!(frames.iter().any(|(_, hi)| (hi[2] - 2200.).abs() < 1e-9));
            let leaf = e
                .pieces
                .iter()
                .find(|p| p.role == OpeningPieceRole::Leaf)
                .unwrap();
            let (lo, hi) = bounds(leaf);
            assert!((lo[2] - 100.).abs() < 1e-9 && (hi[2] - 2150.).abs() < 1e-9);
            // The plan leaf line runs from the hinge on the swing face to the
            // open end; the 3D leaf has an edge on exactly that line.
            let line = w.plan_symbols(&plan).lines[0];
            let edges: Vec<[[f64; 2]; 2]> = (0..8)
                .flat_map(|i| (0..8).map(move |j| (i, j)))
                .filter(|(i, j)| i < j && (i ^ j) == 2 && i & 4 == 0)
                .map(|(i, j)| {
                    let p = |k: usize| [leaf.corners[k][0], leaf.corners[k][1]];
                    [p(i), p(j)]
                })
                .collect();
            assert!(
                edges
                    .iter()
                    .any(|[p, q]| (close(*p, line[0]) && close(*q, line[1]))
                        || (close(*p, line[1]) && close(*q, line[0]))),
                "{hinge:?} {swing:?}: leaf {edges:?} vs plan {line:?}"
            );
        }
    }

    #[test]
    fn window_frame_and_glazing_sit_between_the_plan_pane_lines() {
        let window = ResolvedOpening {
            id: EntityId::new(),
            kind: OpeningKind::Window,
            segment: 0,
            offset_mm: 2000.,
            width_mm: 1200.,
            height_mm: 1200.,
            sill_mm: 900.,
        };
        let w = wall(&[[0., 0.], [6000., 0.]], vec![window]);
        let plan = w.plan().unwrap();
        let e = &w.opening_elements(&plan)[0];
        assert_eq!(e.kind, OpeningKind::Window);
        let roles: Vec<_> = e.pieces.iter().map(|p| p.role).collect();
        assert_eq!(
            roles,
            vec![
                OpeningPieceRole::Frame,
                OpeningPieceRole::Frame,
                OpeningPieceRole::Frame,
                OpeningPieceRole::Frame,
                OpeningPieceRole::Glazing,
            ]
        );
        // The plan draws the panes at y = ±33.3 (a sixth of 200 mm).
        let panes: Vec<f64> = w.plan_symbols(&plan).lines[..2]
            .iter()
            .map(|l| l[0][1])
            .collect();
        for p in e.pieces.iter().map(bounds) {
            let (lo, hi) = p;
            assert!(lo[0] >= 2000. - 1e-9 && hi[0] <= 3200. + 1e-9);
            assert!(lo[1] >= panes[0].min(panes[1]) - 1e-9);
            assert!(hi[1] <= panes[0].max(panes[1]) + 1e-9);
            assert!(lo[2] >= 1000. - 1e-9 && hi[2] <= 2200. + 1e-9);
        }
        let (lo, hi) = bounds(&e.pieces[4]);
        assert!((hi[1] - lo[1] - GLAZING_MM).abs() < 1e-9);
        assert!((lo[0] - 2050.).abs() < 1e-9 && (hi[0] - 3150.).abs() < 1e-9);
        assert!((lo[2] - 1050.).abs() < 1e-9 && (hi[2] - 2150.).abs() < 1e-9);
    }
}
