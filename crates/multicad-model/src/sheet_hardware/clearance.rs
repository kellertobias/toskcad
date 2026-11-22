//! Exact hardware clearance on a formed Rust sheet without a Boolean.
//!
//! Every hardware piece is a prism coaxial with its hole, so a cylinder of
//! the piece's circumradius over its axial span covers it. The sheet is the
//! web prism and one swept cross-section per flange. A plane separates a
//! cylinder from a sheet piece when their support functions along its normal
//! leave no overlap. Contact is allowed, as it removes no volume.
use crate::SheetLayout;
use crate::sheet_forming::{PlacedFlange, section};
use multicad_core::model::{
    CutoutFace, HardwareItem, HardwareSide, ResolvedPattern, ResolvedPreset,
};
use multicad_kernel::CurvedEdge2d;
use multicad_precision::PortableMath as _;
use std::f64::consts::TAU;

type Point = [f64; 3];

fn dot(a: Point, b: Point) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn add(a: Point, b: Point, s: f64) -> Point {
    std::array::from_fn(|i| a[i] + b[i] * s)
}
fn neg(a: Point) -> Point {
    a.map(|c| -c)
}

/// A solid cylinder from `start` running `length` along the unit `axis`.
#[derive(Clone, Copy, Debug)]
pub(super) struct Envelope {
    start: Point,
    axis: Point,
    length: f64,
    radius: f64,
    /// The piece passes through its receiving hole.
    through: bool,
}

impl Envelope {
    fn support(&self, d: Point) -> f64 {
        let along = dot(self.axis, d);
        let across = (dot(d, d) - along * along).max(0.).sqrt();
        dot(self.start, d) + (self.length * along).max(0.) + self.radius * across
    }
}

/// The 2D support of an exact boundary edge in direction `d`.
fn edge_support(edge: &CurvedEdge2d, d: [f64; 2]) -> f64 {
    let at = |p: [f64; 2]| p[0] * d[0] + p[1] * d[1];
    match *edge {
        CurvedEdge2d::Line { a, b } => at(a).max(at(b)),
        CurvedEdge2d::Circle { center, radius_mm } => {
            at(center) + radius_mm * d[0].portable_hypot(d[1])
        }
        CurvedEdge2d::Arc {
            center,
            a,
            b,
            clockwise,
        } => {
            let ends = at(a).max(at(b));
            let (from, to) = if clockwise { (b, a) } else { (a, b) };
            let angle = |p: [f64; 2]| (p[1] - center[1]).portable_atan2(p[0] - center[0]);
            let sweep = (angle(to) - angle(from)).rem_euclid(TAU);
            let toward = (d[1].portable_atan2(d[0]) - angle(from)).rem_euclid(TAU);
            let radius = (a[0] - center[0]).portable_hypot(a[1] - center[1]);
            if toward <= sweep {
                ends.max(at(center) + radius * d[0].portable_hypot(d[1]))
            } else {
                ends
            }
        }
        // Splines bound nothing exactly here; never claim clearance.
        CurvedEdge2d::BSpline(_) => f64::INFINITY,
    }
}

/// One convex-hull support source of the folded sheet.
enum Piece<'a> {
    Web,
    Flange(&'a PlacedFlange, Vec<CurvedEdge2d>),
}

struct Sheet<'a> {
    layout: &'a SheetLayout,
    pieces: Vec<Piece<'a>>,
}

impl<'a> Sheet<'a> {
    /// The sheet's pieces and the index of the one carrying the hardware.
    /// A carrying flange splits into its final leg, a slab like the web,
    /// and the rest of its section.
    fn new(layout: &'a SheetLayout, carrier: Option<usize>) -> (Self, usize) {
        let t = layout.thickness_mm;
        let mut pieces = vec![Piece::Web];
        let mut carried = 0;
        for f in &layout.flanges {
            let section = section(&f.flange, t);
            if carrier != Some(f.flange.index) {
                pieces.push(Piece::Flange(f, section.edges));
                continue;
            }
            let leg = |role: &str| {
                [".leg.a", ".leg.b", ".end"]
                    .iter()
                    .any(|suffix| role.ends_with(suffix))
            };
            let rest = section
                .edges
                .iter()
                .zip(&section.roles)
                .filter(|(_, role)| !leg(role))
                .map(|(edge, _)| edge.clone())
                .collect();
            let start = section.leg_start;
            let at = |a: f64, b: f64| -> [f64; 2] {
                std::array::from_fn(|i| start.point[i] + start.tangent[i] * a + start.normal[i] * b)
            };
            let length = f.flange.leg_mm;
            let corners = [at(0., 0.), at(length, 0.), at(length, t), at(0., t)];
            let slab = (0..4)
                .map(|i| CurvedEdge2d::Line {
                    a: corners[i],
                    b: corners[(i + 1) % 4],
                })
                .collect();
            pieces.push(Piece::Flange(f, rest));
            carried = pieces.len();
            pieces.push(Piece::Flange(f, slab));
        }
        (Self { layout, pieces }, carried)
    }

    fn in_plane(&self, d: Point) -> [f64; 2] {
        [dot(d, self.layout.frame.u), dot(d, self.layout.frame.v)]
    }

    fn support(&self, piece: &Piece<'_>, d: Point) -> f64 {
        let frame = self.layout.frame;
        let thickness = self.layout.thickness_direction;
        let origin = dot(frame.origin, d);
        let up = dot(thickness, d);
        match piece {
            Piece::Web => {
                let planar = self.in_plane(d);
                let radial = self
                    .layout
                    .base
                    .outer
                    .iter()
                    .map(|edge| edge_support(edge, planar))
                    .fold(f64::NEG_INFINITY, f64::max);
                origin + radial + (self.layout.thickness_mm * up).max(0.)
            }
            Piece::Flange(f, edges) => {
                let [along, outward] = [f.along, f.outward].map(|axis| {
                    let p = self.in_plane(d);
                    axis[0] * p[0] + axis[1] * p[1]
                });
                let start = self.in_plane(d);
                let start = f.start[0] * start[0] + f.start[1] * start[1];
                let radial = edges
                    .iter()
                    .map(|edge| edge_support(edge, [outward, up]))
                    .fold(f64::NEG_INFINITY, f64::max);
                origin + start + radial + (f.width * along).max(0.)
            }
        }
    }

    fn directions(&self, envelope: &Envelope) -> Vec<Point> {
        let frame = self.layout.frame;
        let world = |[x, y]: [f64; 2]| std::array::from_fn(|i| frame.u[i] * x + frame.v[i] * y);
        let mut axes = vec![
            envelope.axis,
            frame.u,
            frame.v,
            self.layout.thickness_direction,
        ];
        for f in &self.layout.flanges {
            axes.extend([world(f.along), world(f.outward)]);
            let leg = section(&f.flange, self.layout.thickness_mm).leg_start;
            for vector in [leg.tangent, leg.normal] {
                let o = world(f.outward);
                axes.push(std::array::from_fn(|i| {
                    o[i] * vector[0] + self.layout.thickness_direction[i] * vector[1]
                }));
            }
        }
        axes.into_iter()
            .flat_map(|axis| [axis, neg(axis)])
            .collect()
    }

    fn clear(&self, piece: &Piece<'_>, envelope: &Envelope, tolerance: f64) -> bool {
        self.directions(envelope)
            .into_iter()
            .any(|d| self.support(piece, d) + envelope.support(neg(d)) <= tolerance)
    }
}

/// Covering cylinders of one hardware item at one pattern point. `surface`
/// is the face on the hardware side, `rear` the opposite face, and `out`
/// points away from the sheet on the hardware side.
pub(super) fn envelopes(
    preset: &ResolvedPreset,
    item: HardwareItem,
    surface: Point,
    rear: Point,
    out: Point,
) -> Vec<Envelope> {
    let back = neg(out);
    let cylinder = |start, axis, length, radius, through| Envelope {
        start,
        axis,
        length,
        radius,
        through,
    };
    let hexagon = |across_flats: f64| across_flats / 3f64.sqrt();
    match (preset, item) {
        (
            ResolvedPreset::Screw {
                thread,
                length_mm,
                washers,
                ..
            },
            HardwareItem::Screw,
        ) => {
            let d = thread.data();
            let seat = add(surface, out, if *washers { d.washer[2] } else { 0. });
            vec![
                cylinder(seat, out, d.head[1], d.head[0] / 2., false),
                cylinder(seat, back, *length_mm, d.diameter / 2., true),
            ]
        }
        (ResolvedPreset::Screw { thread, nut, .. }, HardwareItem::Washer) => {
            let d = thread.data();
            let mut rings = vec![cylinder(surface, out, d.washer[2], d.washer[1] / 2., false)];
            if *nut {
                rings.push(cylinder(rear, back, d.washer[2], d.washer[1] / 2., false));
            }
            rings
        }
        (
            ResolvedPreset::Screw {
                thread, washers, ..
            },
            HardwareItem::Nut,
        ) => {
            let d = thread.data();
            let lift = if *washers { d.washer[2] } else { 0. };
            vec![cylinder(
                add(rear, back, lift),
                back,
                d.nut[1],
                hexagon(d.nut[0]),
                false,
            )]
        }
        (
            ResolvedPreset::Standoff {
                length_mm,
                across_flats_mm,
                ..
            },
            _,
        ) => vec![cylinder(
            surface,
            out,
            *length_mm,
            hexagon(*across_flats_mm),
            false,
        )],
        (
            ResolvedPreset::Insert {
                body_mm, height_mm, ..
            },
            _,
        ) => vec![cylinder(surface, out, *height_mm, body_mm / 2., false)],
        _ => vec![],
    }
}

/// Whether the hardware of one pattern item clears the formed sheet, or
/// `None` when the support planes cannot prove it either way.
///
/// A piece through its hole stays inside the hole's column in the carrying
/// web or leg slab when it is no wider than the hole. Everything else, including every other
/// sheet piece, must be separated by a plane.
pub(super) fn formed_clearance(
    layout: &SheetLayout,
    pattern: &ResolvedPattern,
    item: HardwareItem,
    tolerance: f64,
) -> Result<Option<bool>, String> {
    let space = layout.space();
    let (frame, n) = space.face_frame(pattern.face, &layout.flanges)?;
    let t = layout.thickness_mm;
    let (out, surface_offset, back_offset) = match pattern.side {
        HardwareSide::Material => (n, t, 0.),
        HardwareSide::Base => (neg(n), 0., t),
    };
    let (sheet, carrier) = Sheet::new(
        layout,
        match pattern.face {
            CutoutFace::Base => None,
            CutoutFace::Flange { index } => Some(index as usize),
        },
    );
    if let CutoutFace::Flange { index } = pattern.face
        && carrier == 0
    {
        return Err(format!("Flange {} does not exist", index + 1));
    }
    let hole = pattern.preset.hole_mm() / 2.;
    for p in &pattern.points {
        let base = frame.point(*p);
        let pieces = envelopes(
            &pattern.preset,
            item,
            add(base, n, surface_offset),
            add(base, n, back_offset),
            out,
        );
        for envelope in pieces {
            for (index, piece) in sheet.pieces.iter().enumerate() {
                if envelope.through && index == carrier {
                    // Inside the slab the piece stays in the hole's column.
                    if envelope.radius > hole + tolerance {
                        return Ok(Some(true));
                    }
                    continue;
                }
                if !sheet.clear(piece, &envelope, tolerance) {
                    return Ok(None);
                }
            }
        }
    }
    Ok(Some(false))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_supports_follow_arc_sweeps_and_cylinder_axes() {
        let quarter = CurvedEdge2d::Arc {
            center: [0., 0.],
            a: [2., 0.],
            b: [0., 2.],
            clockwise: false,
        };
        let diagonal = std::f64::consts::FRAC_1_SQRT_2;
        assert!((edge_support(&quarter, [diagonal, diagonal]) - 2.).abs() < 1e-12);
        // Outside the sweep only the endpoints count.
        assert!((edge_support(&quarter, [-diagonal, -diagonal]) + 2f64.sqrt()).abs() < 1e-12);
        let reversed = CurvedEdge2d::Arc {
            center: [0., 0.],
            a: [0., 2.],
            b: [2., 0.],
            clockwise: true,
        };
        assert!((edge_support(&reversed, [diagonal, diagonal]) - 2.).abs() < 1e-12);
        let washer = Envelope {
            start: [1., 2., 3.],
            axis: [0., 0., 1.],
            length: 0.5,
            radius: 3.5,
            through: false,
        };
        assert!((washer.support([0., 0., 1.]) - 3.5).abs() < 1e-12);
        assert!((washer.support([0., 0., -1.]) + 3.).abs() < 1e-12);
        assert!((washer.support([1., 0., 0.]) - 4.5).abs() < 1e-12);
    }
}
