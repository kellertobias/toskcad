//! Rectilinear cuts that stay clear of existing circular cuts.
//!
//! When every tool is separated from the closed bounds of each blind pocket
//! and through bore by more than the linear tolerance, removing the tools and
//! removing the round cuts commute exactly:
//! `(stock - round) - tools = (stock - tools) - round`. The plain stock is cut
//! by the exact grid Boolean and every round cut is identified again on the
//! cut stock with its original source identity. The restored mouths, depths
//! and ends must agree with the originals, so a tool that changes a mouth
//! plane can never move a pocket.
use super::*;

fn separated(a: Box3, b: Bounds, tolerance: f64) -> bool {
    (0..3).any(|i| a.max[i] + tolerance < b.min[i] || b.max[i] + tolerance < a.min[i])
}

/// Conservative world bounds of the closed cylinder from `a` to `b`.
fn cylinder_bounds(a: Point, b: Point, radius: f64) -> Box3 {
    Box3 {
        min: std::array::from_fn(|i| a[i].min(b[i]) - radius),
        max: std::array::from_fn(|i| a[i].max(b[i]) + radius),
    }
}

fn moved(actual: Point, expected: Point, tolerance: f64) -> bool {
    actual
        .iter()
        .zip(expected)
        .any(|(a, b)| (a - b).abs() > tolerance)
}

/// Uncut stock, its round cuts, and the blind and through cut origins.
type RoundCutStock = (Solid, Vec<RoundCut>, Vec<Origin>, Vec<Origin>);

struct RoundCut {
    /// World mouth center and bottom center of a blind pocket, or the
    /// entry and exit centers of a through bore.
    ends: [Point; 2],
    through: bool,
    radius: f64,
}

impl Solid {
    /// The stock without its circular cuts and the world placement of each
    /// cut. Blind cuts come first, then through bores, as origins are stored.
    fn round_cut_stock(&self) -> Option<RoundCutStock> {
        let (body, cuts, pocket_origins, through_origins) = match &self.body {
            SolidBody::Orthogonal(body, frame)
                if body.has_curved_cuts()
                    && body.slots.is_empty()
                    && body.imported_roundoff_mm.is_none() =>
            {
                let world = |point| transform_point(*frame, point);
                let faces = &body.boundary.faces;
                let cuts = body
                    .pockets
                    .iter()
                    .map(|pocket| {
                        let normal = faces[pocket.face].normal;
                        let bottom = std::array::from_fn(|i| {
                            pocket.center_mm[i] - normal[i] * pocket.depth_mm
                        });
                        RoundCut {
                            ends: [world(pocket.center_mm), world(bottom)],
                            through: false,
                            radius: pocket.radius_mm,
                        }
                    })
                    .chain(body.through_bores.iter().map(|bore| RoundCut {
                        ends: bore.centers_mm.map(world),
                        through: true,
                        radius: bore.radius_mm,
                    }))
                    .collect();
                let mut stock = (**body).clone();
                stock.pockets = Rc::from([]);
                stock.pocket_origins = Rc::from([]);
                stock.through_bores = Rc::from([]);
                stock.through_origins = Rc::from([]);
                (
                    SolidBody::Orthogonal(Rc::new(stock), *frame),
                    cuts,
                    body.pocket_origins.to_vec(),
                    body.through_origins.to_vec(),
                )
            }
            SolidBody::PolygonPrism(stock)
                if !stock.pockets.is_empty()
                    && stock.prism.circular_holes().is_empty()
                    && stock.bore_treatments.is_none() =>
            {
                let sign = if stock.reverse_sweep { -1. } else { 1. };
                let height = stock.prism.height().as_f64();
                let cuts = stock
                    .pockets
                    .iter()
                    .map(|pocket| {
                        let [x, y] = pocket.center_mm;
                        let ends = [pocket.mouth_z_mm(height), pocket.floor_z_mm(height)]
                            .map(|z| transform_point(stock.frame, [x, y, sign * z]));
                        RoundCut {
                            ends,
                            through: false,
                            radius: pocket.radius_mm,
                        }
                    })
                    .collect();
                let origins = stock.pocket_origins.to_vec();
                let mut stock = (**stock).clone();
                stock.pockets = Rc::from([]);
                stock.pocket_origins = Rc::from([]);
                (
                    SolidBody::PolygonPrism(Rc::new(stock)),
                    cuts,
                    origins,
                    Vec::new(),
                )
            }
            _ => return None,
        };
        Some((
            Self {
                body,
                policy: self.policy,
                provenance: vec![],
                rust_origins: self.rust_origins.clone(),
                retained_placement: None,
            },
            cuts,
            pocket_origins,
            through_origins,
        ))
    }

    /// The uncut stock when every circular cut is separated from `probe` by
    /// more than `tolerance`. Inside `probe` the stock is then the identical
    /// point set, so any query confined to `probe` may use its planar cells.
    pub(in crate::geometry) fn stock_clear_of_round_cuts(
        &self,
        probe: Bounds,
        tolerance: f64,
    ) -> Option<Self> {
        let (uncut, cuts, ..) = self.round_cut_stock()?;
        cuts.iter()
            .all(|cut| {
                separated(
                    cylinder_bounds(cut.ends[0], cut.ends[1], cut.radius),
                    probe,
                    tolerance,
                )
            })
            .then_some(uncut)
    }

    /// `None` leaves bodies without round cuts, with slots or authored
    /// overrides, and tools that reach a round cut to the other providers.
    pub(in crate::geometry) fn subtract_clear_of_round_cuts(
        &self,
        tools: &[Self],
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        if tools.is_empty()
            || self.has_authored_metadata_override()
            || policy != self.policy
            || policy.boolean_fuzzy_mm != 0.
        {
            return None;
        }
        let (uncut, cuts, pocket_origins, through_origins) = self.round_cut_stock()?;
        let tolerance = policy.linear_tolerance_mm;
        let supports: Vec<_> = cuts
            .iter()
            .map(|cut| cylinder_bounds(cut.ends[0], cut.ends[1], cut.radius))
            .collect();
        for tool in tools {
            let reach = tool.bounds().ok()?;
            if !supports
                .iter()
                .all(|&support| separated(support, reach, tolerance))
            {
                return None;
            }
        }
        Some((|| {
            let Some(result) =
                uncut.orthogonal_boolean_with_tools(tools, BooleanOperation::Subtract, policy)?
            else {
                return Ok(None);
            };
            let (cut, cut_frame) = result
                .orthogonal_body()
                .expect("rectilinear Boolean result");
            debug_assert_eq!(cut_frame, PlaneFrame::default());
            let overshoot = policy.minimum_feature_mm.max(8. * tolerance);
            let mut pockets = Vec::new();
            let mut bores = Vec::new();
            for round in &cuts {
                let travel = subtract(round.ends[1], round.ends[0]);
                let length = dot(travel, travel).sqrt();
                let delta = travel.map(|x| overshoot * x / length);
                let start = subtract(round.ends[0], delta);
                if round.through {
                    let end = std::array::from_fn(|i| round.ends[1][i] + delta[i]);
                    let restored = cut
                        .identify_through(start, end, round.radius, policy)
                        .map_err(circular_pocket_body::pocket_error)?;
                    if restored
                        .centers_mm
                        .iter()
                        .zip(round.ends)
                        .any(|(&actual, expected)| moved(actual, expected, tolerance))
                    {
                        return Err(unsupported(
                            "Rectilinear cut changes the ends of an existing bore",
                        ));
                    }
                    bores.push(restored);
                } else {
                    let restored = cut
                        .identify_blind(start, round.ends[1], round.radius, policy)
                        .map_err(circular_pocket_body::pocket_error)?;
                    if moved(restored.center_mm, round.ends[0], tolerance)
                        || (restored.depth_mm - length).abs() > tolerance
                    {
                        return Err(unsupported(
                            "Rectilinear cut changes the mouth of an existing pocket",
                        ));
                    }
                    pockets.push(restored);
                }
            }
            result
                .with_orthogonal_round_cuts(pockets, pocket_origins, bores, through_origins)
                .map(Some)
        })())
    }
}
