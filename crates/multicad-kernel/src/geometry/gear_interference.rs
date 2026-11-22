//! The interference solid of two placed straight Rust gears.
//!
//! Spur gears and internal rings are normal extrusions of exact sections on
//! parallel, codirected axes, so their common material is the exact common
//! region of the two sections swept over the common axial span. The region
//! comes from `clearance::intersect_sections`; line and arc pieces stay
//! lines and arcs, and every involute piece becomes its certified cubic
//! B-spline, as in STEP exchange. Each connected region is one extruded
//! solid, and several regions form one exact compound.
use super::spur_gear_body::coaxial_placements;
use super::*;
use crate::{BSpline2d, CurvedEdge2d, CurvedProfile};
use multicad_core::references::Origin;
use multicad_rust_kernel::spur_gear::{
    Edge,
    clearance::{self, Piece, Placement, Section},
};

const OPERATION: &str = "Boolean";

/// Spline error of the exported involute pieces, as for STEP flanks.
const FLANK_ERROR_MM: f64 = 1e-8;

struct Member<'a> {
    section: Section<'a>,
    frame: PlaneFrame,
    height_mm: f64,
    support: Box<dyn Fn([f64; 3]) -> f64 + 'a>,
}

fn member(solid: &Solid) -> Option<Member<'_>> {
    if let Some((body, frame)) = solid.spur_gear_body() {
        let profile = &body.profile;
        if profile.hub().is_some() {
            return None;
        }
        return Some(Member {
            section: Section::Spur(profile),
            frame,
            height_mm: profile.height_mm(),
            support: Box::new(|direction| profile.support(direction)),
        });
    }
    let (body, frame) = solid.internal_ring_body()?;
    let ring = &body.ring;
    Some(Member {
        section: Section::Ring(ring),
        frame,
        height_mm: ring.height_mm(),
        support: Box::new(|direction| ring.support(direction)),
    })
}

fn failure(kind: GeometryErrorKind, message: impl Into<String>) -> GeometryError {
    error(OPERATION, kind, message)
}

/// One exact piece in the common chart, as a sketch curve.
fn sketch(piece: Piece, placement: Placement) -> Result<CurvedEdge2d, GeometryError> {
    use multicad_rust_kernel::normal_sweep::Edge as Planar;
    let at = |p: [f64; 2]| placement.point(p);
    Ok(match piece.edge {
        Edge::Planar(Planar::Line { a, b }) => CurvedEdge2d::Line { a: at(a), b: at(b) },
        Edge::Planar(Planar::Arc {
            center,
            a,
            b,
            clockwise,
        }) => CurvedEdge2d::Arc {
            center: at(center),
            a: at(a),
            b: at(b),
            clockwise,
        },
        Edge::Planar(Planar::Circle { .. }) => {
            return Err(failure(
                GeometryErrorKind::InvalidResult,
                "Rust gear interference returned an unsplit circle",
            ));
        }
        Edge::Involute(flank) => {
            let cubic = flank.cubic_spline(FLANK_ERROR_MM).map_err(|cause| {
                failure(
                    GeometryErrorKind::Limit,
                    format!("Rust gear interference flank spline: {cause:?}"),
                )
            })?;
            CurvedEdge2d::BSpline(BSpline2d {
                degree: 3,
                poles: cubic.poles_mm.into_iter().map(at).collect(),
                knots: cubic.knots,
                multiplicities: cubic.multiplicities.into_iter().map(u32::from).collect(),
                weights: None,
            })
        }
    })
}

fn start(edge: &CurvedEdge2d) -> [f64; 2] {
    match edge {
        CurvedEdge2d::Line { a, .. } | CurvedEdge2d::Arc { a, .. } => *a,
        CurvedEdge2d::BSpline(spline) => spline.poles[0],
        CurvedEdge2d::Circle { .. } => unreachable!("pieces are never circles"),
    }
}

/// A closed sketch wire. Both sections solve each crossing on their own
/// exact edge; the wire takes the start of the following piece, which the
/// crossing tolerance keeps within a nanometre of the end it replaces.
fn wire(pieces: &[Piece], placements: [Placement; 2]) -> Result<Vec<CurvedEdge2d>, GeometryError> {
    let mut edges = pieces
        .iter()
        .map(|&piece| sketch(piece, placements[piece.section]))
        .collect::<Result<Vec<_>, _>>()?;
    let starts: Vec<_> = edges.iter().map(start).collect();
    for (i, edge) in edges.iter_mut().enumerate() {
        let end = starts[(i + 1) % starts.len()];
        match edge {
            CurvedEdge2d::Line { b, .. } | CurvedEdge2d::Arc { b, .. } => *b = end,
            CurvedEdge2d::BSpline(spline) => {
                *spline.poles.last_mut().expect("cubic poles") = end;
            }
            CurvedEdge2d::Circle { .. } => unreachable!("pieces are never circles"),
        }
    }
    Ok(edges)
}

impl Solid {
    /// The exact common material of two placed straight Rust gears whose
    /// pair certificate proved a common interior point. `None` leaves the
    /// pair to the typed interference refusal.
    pub(in crate::geometry) fn straight_gear_interference(
        &self,
        other: &Self,
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        let members = [member(self)?, member(other)?];
        let frames = [members[0].frame, members[1].frame];
        let heights = [members[0].height_mm, members[1].height_mm];
        let [first, second] = members;
        let sections = [first.section, second.section];
        let placements = match coaxial_placements(
            (OPERATION, "Rust gear interference"),
            policy,
            frames,
            [
                (first.height_mm, first.support),
                (second.height_mm, second.support),
            ],
        ) {
            Ok(placements) => placements,
            Err(_) => return None,
        };
        // The common axial span in the first frame.
        let normal = frames[0].normal();
        let lift = (0..3)
            .map(|i| (frames[1].origin[i] - frames[0].origin[i]) * normal[i])
            .sum::<f64>();
        let low = lift.max(0.);
        let high = (lift + heights[1]).min(heights[0]);
        if high - low <= policy.linear_tolerance_mm {
            return None;
        }
        let regions = match clearance::intersect_sections(
            sections,
            placements,
            policy.linear_tolerance_mm,
            None,
        ) {
            Ok(Some(regions)) => regions,
            Ok(None) | Err(clearance::Error::Invalid) => return None,
            Err(clearance::Error::Limit) => {
                return Some(Err(failure(
                    GeometryErrorKind::Limit,
                    "Rust gear interference exceeded its crossing budget",
                )));
            }
            Err(clearance::Error::Cancelled) => {
                return Some(Err(failure(
                    GeometryErrorKind::Cancelled,
                    "Rust gear interference superseded",
                )));
            }
        };
        Some(self.extrude_regions(&regions, placements, frames[0], [low, high], policy))
    }

    fn extrude_regions(
        &self,
        regions: &[clearance::Region],
        placements: [Placement; 2],
        frame: PlaneFrame,
        [low, high]: [f64; 2],
        policy: GeometryPolicy,
    ) -> Result<Option<Self>, GeometryError> {
        if regions.is_empty() {
            return Ok(None);
        }
        let feature = match self.rust_topology() {
            Some(Ok(topology)) => topology.owner,
            Some(Err(error)) => return Err(error),
            None => EntityId::new(),
        };
        let normal = frame.normal();
        let base = PlaneFrame {
            origin: std::array::from_fn(|i| frame.origin[i] + normal[i] * low),
            ..frame
        };
        let mut solids = Vec::with_capacity(regions.len());
        for (index, region) in regions.iter().enumerate() {
            let profile = CurvedProfile {
                frame: base,
                outer: wire(&region.outer, placements)?,
                holes: region
                    .holes
                    .iter()
                    .map(|hole| wire(hole, placements))
                    .collect::<Result<_, _>>()?,
            };
            let pieces = region.outer.iter().chain(region.holes.iter().flatten());
            let origins: Vec<_> = pieces
                .enumerate()
                .map(|(i, piece)| Origin {
                    source: feature,
                    role: format!("gear/interference/{index}/{}/{i}", piece.section),
                })
                .collect();
            let region_feature = EntityId::new();
            solids.push(Self::rust_certified_region_extrusion(
                &profile,
                high - low,
                region_feature,
                &origins,
                policy,
            )?);
        }
        if solids.len() == 1 {
            return Ok(solids.pop());
        }
        let members: Vec<&Self> = solids.iter().collect();
        RustCompound::assemble(&members, feature, policy, None).map(Some)
    }
}
