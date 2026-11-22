//! Inset, adjacent, multi-bend and hemmed flanges, and polygon bases, do not
//! share one cross-section. Join each flange's exact root boundary to one
//! relieved web rather than fusing display solids.
use super::*;
use crate::sheet_forming::{FlatFrame, section, shape_edges};
use multicad_core::references::{AnchoredFrame, ReferenceFrame};
use multicad_kernel::rust_kernel::{
    normal_sweep::{Edge, signed_area},
    warped_union::{SheetFlange, SheetFlangeOpening, SheetWebEdge, formed_sheet_with_overhangs},
};
use std::collections::BTreeMap;

pub(super) fn build(
    layout: &mut SheetLayout,
    base_origins: &[Origin],
    policy: GeometryPolicy,
) -> Result<Solid, String> {
    let thickness = layout.thickness_mm;
    let count = layout.flanges.len();
    if !(1..=4).contains(&count) {
        return Err("Rust formed sheet supports up to four flanges".into());
    }
    let edge_count =
        layout.base.outer.len() + layout.base.holes.iter().map(Vec::len).sum::<usize>();
    if base_origins.len() != edge_count {
        return Err("Rust formed sheet needs one authored origin per base edge".into());
    }
    let mut origins = Vec::new();
    let feature = layout.forming_feature;
    let authored = |role| Origin {
        source: feature,
        role,
    };
    let mut source = |origin| {
        origins.push(origin);
        origins.len() - 1
    };
    let mut outer = Vec::new();
    let mut developed = Vec::new();
    let web = folded_web(layout)?;
    let mut overhangs = vec![[0.; 2]; count];
    for (index, (edge, web)) in layout.base.outer.iter().zip(web).enumerate() {
        let CurvedEdge2d::Line { a, b } = *edge else {
            unreachable!("folded_web refuses curved base edges")
        };
        let base_origin = &base_origins[index];
        let Some((hinge, f)) = layout
            .flanges
            .iter()
            .enumerate()
            .find(|(_, f)| f.flange.edge == index)
        else {
            outer.push(SheetWebEdge {
                ends: [web.points[0], web.points[1]],
                source: source(base_origin.clone()),
                flange: None,
            });
            developed.push(PlanarEdge::Line { a, b });
            continue;
        };
        overhangs[hinge] = web.overhang;
        for (i, pair) in web.points.windows(2).enumerate() {
            outer.push(SheetWebEdge {
                ends: [pair[0], pair[1]],
                source: source(Origin {
                    source: base_origin.source,
                    role: format!("{}.segment.{i}", base_origin.role),
                }),
                flange: (Some(i) == web.hinge).then_some(hinge),
            });
        }
        let extension = f.flange.developed_mm(thickness, layout.k_factor);
        developed.extend(edge_points(f, a, b, extension)?.windows(2).map(|pair| {
            PlanarEdge::Line {
                a: pair[0],
                b: pair[1],
            }
        }));
    }
    let refusal = crossing(&developed).then(|| {
        "The developed flanges overlap each other or the base in the flat pattern".to_string()
    });
    let mut loops = vec![PlanarLoop {
        outer: true,
        edges: developed,
    }];
    let mut holes = Vec::new();
    let mut hole_origins = base_origins[layout.base.outer.len()..].iter();
    for hole in &layout.base.holes {
        let ring = hole
            .iter()
            .map(|edge| {
                let origin = hole_origins.next().expect("counted base origins").clone();
                Ok((profile_edge(edge)?, source(origin)))
            })
            .collect::<Result<Vec<_>, String>>()?;
        holes.push(oriented(ring, false));
        loops.push(channel::developed_hole(hole)?);
    }
    let mut flange_cuts: Vec<Vec<Opening>> = (0..count).map(|_| Vec::new()).collect();
    let flange_index = |index| {
        layout
            .flanges
            .iter()
            .position(|f| f.flange.index == index as usize)
            .ok_or_else(|| format!("Flange {} does not exist", index + 1))
    };
    for cut in &layout.cutouts {
        let name = cut.source.id();
        let flat = flat_cutout(cut, &layout.flanges, thickness, layout.k_factor)
            .ok_or_else(|| format!("{name} names a missing flange"))?;
        loops.push(channel::developed_hole(&flat)?);
        let outline = shape_edges(&cut.shape);
        let ring = outline
            .iter()
            .enumerate()
            .map(|(i, edge)| {
                Ok((
                    profile_edge(edge)?,
                    source(authored(format!("sheet.{name}.edge.{i}"))),
                ))
            })
            .collect::<Result<Vec<_>, String>>()?;
        match cut.face {
            CutoutFace::Base => {
                let web: Vec<_> = outer.iter().map(|edge| edge.ends).collect();
                if !inside(&outline, &web, policy.minimum_feature_mm) {
                    return Err(format!("{name} is not fully inside its face"));
                }
                holes.push(oriented(ring, false))
            }
            CutoutFace::Flange { index } => {
                flange_cuts[flange_index(index)?].push(Opening::new(layout, cut.face, ring)?)
            }
        }
    }
    for cut in &layout.profile_cuts {
        let flange = match cut.face {
            CutoutFace::Base => None,
            CutoutFace::Flange { index } => Some(flange_index(index)?),
        };
        let outline = if let Some(index) = flange {
            let f = &layout.flanges[index];
            let flat = FlatFrame::flange(f, f.flange.leg_offset_mm(thickness, layout.k_factor));
            super::profile_cut::transform(&cut.outline, |p| flat.coordinates(p), flat.reflected())?
        } else {
            cut.outline.clone()
        };
        if outline.len() != cut.origins.len() {
            return Err("Inset sheet cut needs one authored origin per edge".into());
        }
        let ring = outline
            .iter()
            .zip(&cut.origins)
            .map(|(edge, origin)| Ok((profile_edge(edge)?, source(origin.clone()))))
            .collect::<Result<Vec<_>, String>>()?;
        if let Some(index) = flange {
            flange_cuts[index].push(Opening::new(layout, cut.face, ring)?);
        } else {
            holes.push(oriented(ring, false));
        }
        loops.push(channel::developed_hole(&cut.outline)?);
    }
    let mut authored_frames = BTreeMap::new();
    let thickness_side = dot(layout.frame.normal(), layout.thickness_direction);
    let sections = layout
        .flanges
        .iter()
        .zip(&flange_cuts)
        .map(|(f, cuts)| {
            let profile = section(&f.flange, thickness);
            let edges = profile
                .edges
                .iter()
                .zip(&profile.roles)
                .map(|(edge, role)| {
                    let source = source(authored(role.clone()));
                    if let CurvedEdge2d::Line { a, b } = *edge
                        && (role.ends_with(".leg.a") || role.ends_with(".leg.b"))
                    {
                        let position = |[outward, height]: [f64; 2]| {
                            [
                                f.start[0] + f.outward[0] * outward,
                                f.start[1] + f.outward[1] * outward,
                                thickness_side * height,
                            ]
                        };
                        let [a, b] = [position(a), position(b)];
                        authored_frames.insert(
                            source,
                            AnchoredFrame {
                                origin: authored(format!("side/{feature}/{role}")),
                                frame: ReferenceFrame {
                                    origin_mm: a,
                                    u: unit(std::array::from_fn(|i| b[i] - a[i])),
                                    v: [f.along[0], f.along[1], 0.],
                                },
                            },
                        );
                    }
                    Ok((profile_edge(edge)?, source))
                })
                .collect::<Result<Vec<_>, String>>()?;
            let caps = ["start", "end"].map(|end| {
                source(authored(format!(
                    "sheet.flange.{}.cap.{end}",
                    f.flange.index
                )))
            });
            let mut openings = Vec::new();
            if !cuts.is_empty() {
                let find_wall = |side| {
                    let role = format!("sheet.flange.{}.leg.{side}", f.flange.index);
                    profile
                        .roles
                        .iter()
                        .position(|candidate| candidate == &role)
                        .ok_or("An inset flange opening needs two complete leg walls")
                };
                let walls = [find_wall("a")?, find_wall("b")?];
                openings.extend(cuts.iter().map(|cut| cut.as_sheet_opening(walls)));
            }
            Ok((edges, caps, openings))
        })
        .collect::<Result<Vec<_>, String>>()?;
    let flanges = sections
        .iter()
        .map(|(section, ends, _)| SheetFlange {
            section,
            ends: *ends,
        })
        .collect::<Vec<_>>();
    let openings: Vec<_> = sections
        .iter()
        .map(|(_, _, cuts)| cuts.as_slice())
        .collect();
    let web_caps = [channel::WEB_A, channel::WEB_B].map(|role| source(authored(role.into())));
    let shell = formed_sheet_with_overhangs(
        &outer,
        &holes,
        web_caps,
        &flanges,
        &openings,
        &overhangs,
        thickness.copysign(thickness_side),
        policy.linear_tolerance_mm,
        policy.minimum_feature_mm,
    )
    .map_err(|cause| format!("The inset formed sheet boundary cannot be certified: {cause:?}"))?;
    let solid = Solid::rust_place_formed_sheet(
        shell,
        layout.frame,
        layout.forming_feature,
        &origins,
        web_caps,
        &authored_frames,
        policy,
    )
    .map_err(|cause| cause.to_string())?;
    let contours = PlanarContours {
        frame: ReferenceFrame {
            origin_mm: layout.frame.origin,
            u: layout.frame.u,
            v: layout.frame.v,
        },
        loops,
    };
    contours.validate(policy.linear_tolerance_mm)?;
    match refusal {
        Some(refusal) => layout.developed_refusal = Some(refusal),
        None => layout.developed = Some(contours),
    }
    Ok(solid)
}

/// Whether a closed polyline touches or crosses itself away from the shared
/// ends of consecutive edges. The developed outline of flanges spliced into
/// their base edges is simple exactly when no developed region overlaps
/// another or the base.
fn crossing(outline: &[PlanarEdge]) -> bool {
    let segments: Vec<_> = outline
        .iter()
        .map(|edge| match *edge {
            PlanarEdge::Line { a, b } => Some([a, b]),
            PlanarEdge::Arc { .. } => None,
        })
        .collect();
    let Some(segments) = segments.into_iter().collect::<Option<Vec<_>>>() else {
        return false;
    };
    let n = segments.len();
    (0..n).any(|i| {
        (i + 1..n).any(|j| {
            let adjacent = j == i + 1 || (i == 0 && j == n - 1);
            let [p, q] = segments[i];
            let [r, s] = segments[j];
            if adjacent {
                // Consecutive edges meet at one end; they may not fold back.
                let shared = if j == i + 1 { (q, p, s) } else { (p, q, r) };
                let (corner, back, next) = shared;
                let u = [back[0] - corner[0], back[1] - corner[1]];
                let v = [next[0] - corner[0], next[1] - corner[1]];
                return (u[0] * v[1] - u[1] * v[0]).abs()
                    <= 1e-9 * u[0].hypot(u[1]) * v[0].hypot(v[1])
                    && u[0] * v[0] + u[1] * v[1] > 0.;
            }
            segments_meet(p, q, r, s)
        })
    })
}

/// Whether a line or circle outline lies strictly inside a web polygon,
/// clear of its boundary by at least `clearance`.
fn inside(outline: &[CurvedEdge2d], web: &[[[f64; 2]; 2]], clearance: f64) -> bool {
    let contains = |p: [f64; 2]| {
        web.iter()
            .filter(|[a, b]| (a[1] > p[1]) != (b[1] > p[1]))
            .filter(|[a, b]| p[0] < a[0] + (p[1] - a[1]) / (b[1] - a[1]) * (b[0] - a[0]))
            .count()
            % 2
            == 1
    };
    let distance = |p: [f64; 2], [a, b]: [[f64; 2]; 2]| {
        let d = [b[0] - a[0], b[1] - a[1]];
        let t = (((p[0] - a[0]) * d[0] + (p[1] - a[1]) * d[1]) / (d[0] * d[0] + d[1] * d[1]))
            .clamp(0., 1.);
        (p[0] - a[0] - t * d[0]).hypot(p[1] - a[1] - t * d[1])
    };
    outline.iter().all(|edge| match *edge {
        CurvedEdge2d::Circle { center, radius_mm } => {
            contains(center)
                && web
                    .iter()
                    .all(|&side| distance(center, side) > radius_mm + clearance)
        }
        CurvedEdge2d::Line { a, b } => {
            contains(a)
                && web.iter().all(|&[p, q]| {
                    !segments_meet(a, b, p, q)
                        && distance(a, [p, q]) > clearance
                        && distance(p, [a, b]) > clearance
                })
        }
        // Arcs only come from authored sketch profiles, which the kernel
        // checks against the web itself.
        _ => true,
    })
}

fn segments_meet(p: [f64; 2], q: [f64; 2], r: [f64; 2], s: [f64; 2]) -> bool {
    let orient = |a: [f64; 2], b: [f64; 2], c: [f64; 2]| {
        (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    };
    let within = |a: [f64; 2], b: [f64; 2], c: [f64; 2]| {
        c[0] >= a[0].min(b[0]) - 1e-9
            && c[0] <= a[0].max(b[0]) + 1e-9
            && c[1] >= a[1].min(b[1]) - 1e-9
            && c[1] <= a[1].max(b[1]) + 1e-9
    };
    let scale = [p, q, r, s]
        .iter()
        .flatten()
        .fold(1f64, |m, v| m.max(v.abs()));
    let eps = 1e-12 * scale * scale;
    let [d1, d2, d3, d4] = [
        orient(r, s, p),
        orient(r, s, q),
        orient(p, q, r),
        orient(p, q, s),
    ];
    if ((d1 > eps && d2 < -eps) || (d1 < -eps && d2 > eps))
        && ((d3 > eps && d4 < -eps) || (d3 < -eps && d4 > eps))
    {
        return true;
    }
    (d1.abs() <= eps && within(r, s, p))
        || (d2.abs() <= eps && within(r, s, q))
        || (d3.abs() <= eps && within(p, q, r))
        || (d4.abs() <= eps && within(p, q, s))
}

/// One base edge of the folded web: its points, the segment welded to its
/// flange root, and how far that root runs past the weld at each end.
struct WebEdge {
    points: Vec<[f64; 2]>,
    hinge: Option<usize>,
    overhang: [f64; 2],
}

/// The folded web outline. A relief as wide as its inset reaches the corner;
/// at a right-angled corner its inward side then runs along the neighboring
/// edge, which the web loses for the relief depth. A neighboring flange keeps
/// its full root there, as an exposed face beside the relief.
fn folded_web(layout: &SheetLayout) -> Result<Vec<WebEdge>, String> {
    const TOLERANCE: f64 = 1e-6;
    let edges = layout
        .base
        .outer
        .iter()
        .map(|edge| match *edge {
            CurvedEdge2d::Line { a, b } => Ok((a, b)),
            _ => Err("Rust formed sheet needs straight base edges".to_string()),
        })
        .collect::<Result<Vec<_>, _>>()?;
    let n = edges.len();
    let flange = |index: usize| layout.flanges.iter().find(|f| f.flange.edge == index);
    let direction = |index: usize| {
        let (a, b) = edges[index];
        let length = (b[0] - a[0]).hypot(b[1] - a[1]);
        ([(b[0] - a[0]) / length, (b[1] - a[1]) / length], length)
    };
    // Depth of a corner-flush relief on `index` at its start or end corner,
    // when it runs along the neighbor across that corner.
    let flush = |index: usize, end: bool| -> Option<f64> {
        let f = flange(index)?;
        let inset = if end {
            f.flange.end_inset_mm
        } else {
            f.flange.start_inset_mm
        };
        let [w, depth] = f.flange.relief_mm?;
        let neighbor = if end {
            (index + 1) % n
        } else {
            (index + n - 1) % n
        };
        let (along, _) = direction(neighbor);
        let toward = if end { [-along[0], -along[1]] } else { along };
        let square = (f.outward[0] * toward[0] + f.outward[1] * toward[1] - 1.).abs() < 1e-9;
        (inset > TOLERANCE && (w - inset).abs() <= TOLERANCE && square && n > 1).then_some(depth)
    };
    (0..n)
        .map(|index| {
            let (a, _) = edges[index];
            let (along, length) = direction(index);
            let clip = [flush((index + n - 1) % n, true), flush((index + 1) % n, false)]
                .map(|depth| depth.unwrap_or(0.));
            if clip[0] + clip[1] >= length - TOLERANCE {
                return Err(format!(
                    "The inset formed sheet boundary cannot be certified: corner reliefs remove all of base edge {}",
                    index + 1
                ));
            }
            let at = |s: f64, q: f64| -> [f64; 2] {
                let outward = [along[1], -along[0]];
                std::array::from_fn(|i| a[i] + along[i] * s + outward[i] * q)
            };
            let Some(f) = flange(index) else {
                return Ok(WebEdge {
                    points: vec![at(clip[0], 0.), at(length - clip[1], 0.)],
                    hinge: None,
                    overhang: [0.; 2],
                });
            };
            let outward = [along[1], -along[0]];
            let sign = f.outward[0] * outward[0] + f.outward[1] * outward[1];
            let at = |s: f64, q: f64| at(s, q * sign);
            let relief = f.flange.relief_mm;
            let insets = [f.flange.start_inset_mm, length - f.flange.end_inset_mm];
            let own_flush = [flush(index, false).is_some(), flush(index, true).is_some()];
            let mut points = Vec::new();
            let mut overhang = [0.; 2];
            if !own_flush[0] {
                points.push(at(clip[0], 0.));
            }
            // Start corner: own relief notch, or a root that begins inside the
            // neighbor's notch.
            let root_start = match relief {
                Some([w, depth]) if f.flange.start_inset_mm > TOLERANCE => {
                    if clip[0] > 0. && insets[0] - w < clip[0] + TOLERANCE {
                        return Err(overlap(index));
                    }
                    if !own_flush[0] {
                        points.push(at(insets[0] - w, 0.));
                    }
                    points.extend([at(insets[0] - w, -depth), at(insets[0], -depth)]);
                    insets[0]
                }
                _ if insets[0] > clip[0] + TOLERANCE => insets[0],
                _ => {
                    overhang[0] = clip[0] - insets[0];
                    clip[0]
                }
            };
            let limit = length - clip[1];
            let root_end = match relief {
                Some(_) if f.flange.end_inset_mm > TOLERANCE => insets[1],
                _ if insets[1] < limit - TOLERANCE => insets[1],
                _ => {
                    overhang[1] = insets[1] - limit;
                    limit
                }
            };
            if root_end - root_start <= TOLERANCE {
                return Err(overlap(index));
            }
            let root = [at(root_start, 0.), at(root_end, 0.)];
            points.extend(root);
            if let Some([w, depth]) = relief
                && f.flange.end_inset_mm > TOLERANCE
            {
                if clip[1] > 0. && insets[1] + w > limit - TOLERANCE {
                    return Err(overlap(index));
                }
                points.extend([at(insets[1], -depth), at(insets[1] + w, -depth)]);
                if !own_flush[1] {
                    points.push(at(insets[1] + w, 0.));
                }
            }
            if !own_flush[1] {
                points.push(at(limit, 0.));
            }
            points.dedup_by(|p, q| (p[0] - q[0]).hypot(p[1] - q[1]) < 1e-9);
            let hinge = points.windows(2).position(|pair| pair == root);
            Ok(WebEdge {
                points,
                hinge,
                overhang,
            })
        })
        .collect()
}

fn overlap(index: usize) -> String {
    format!(
        "The inset formed sheet boundary cannot be certified: corner reliefs overlap beside base edge {}",
        index + 1
    )
}

type Ring = Vec<(Edge, usize)>;

fn oriented(mut ring: Ring, ccw: bool) -> Ring {
    let edges: Vec<_> = ring.iter().map(|&(edge, _)| edge).collect();
    if (signed_area(&edges) > 0.) != ccw {
        ring.reverse();
        ring.iter_mut()
            .for_each(|(edge, _)| *edge = edge.reversed());
    }
    ring
}

struct Opening {
    frame: PlaneFrame,
    through_mm: [f64; 3],
    outline: Ring,
}

impl Opening {
    fn new(layout: &SheetLayout, face: CutoutFace, outline: Ring) -> Result<Self, String> {
        let (frame, normal) = layout.space().face_frame(face, &layout.flanges)?;
        let axes = [layout.frame.u, layout.frame.v, layout.frame.normal()];
        let local = |vector| axes.map(|axis| dot(vector, axis));
        Ok(Self {
            frame: PlaneFrame {
                origin: local(std::array::from_fn(|i| {
                    frame.origin[i] - layout.frame.origin[i]
                })),
                u: local(frame.u),
                v: local(frame.v),
            },
            through_mm: local(scale(normal, layout.thickness_mm)),
            outline: oriented(outline, true),
        })
    }

    fn as_sheet_opening(&self, walls: [usize; 2]) -> SheetFlangeOpening<'_> {
        SheetFlangeOpening {
            walls,
            origin_mm: self.frame.origin,
            u: self.frame.u,
            v: self.frame.v,
            through_mm: self.through_mm,
            outline: &self.outline,
        }
    }
}

fn profile_edge(edge: &CurvedEdge2d) -> Result<Edge, String> {
    match *edge {
        CurvedEdge2d::Line { a, b } => Ok(Edge::Line { a, b }),
        CurvedEdge2d::Arc {
            center,
            a,
            b,
            clockwise,
        } => Ok(Edge::Arc {
            center,
            a,
            b,
            clockwise,
        }),
        CurvedEdge2d::Circle { center, radius_mm } => Ok(Edge::Circle {
            center,
            radius_mm,
            clockwise: true,
        }),
        CurvedEdge2d::BSpline(_) => {
            Err("Inset formed sheet requires straight or circular boundaries".into())
        }
    }
}

fn point(f: &PlacedFlange, along: f64, out: f64) -> [f64; 2] {
    std::array::from_fn(|i| f.start[i] + f.along[i] * along + f.outward[i] * out)
}

/// One web edge and its reliefs, with optional outward flat extension.
pub(super) fn edge_points(
    f: &PlacedFlange,
    a: [f64; 2],
    b: [f64; 2],
    extension: f64,
) -> Result<Vec<[f64; 2]>, String> {
    let at = |along, out| point(f, along, out);
    let mut points = vec![a];
    let relief = f.flange.relief_mm;
    if f.flange.start_inset_mm > 1e-9 {
        if let Some([w, depth]) = relief {
            if w > f.flange.start_inset_mm + 1e-9 {
                return Err(format!(
                    "Flange {} relief is wider than its start inset",
                    f.flange.index + 1
                ));
            }
            points.extend([at(-w, 0.), at(-w, -depth), at(0., -depth)]);
        }
        points.push(at(0., 0.));
    }
    points.extend([at(0., extension), at(f.width, extension)]);
    if f.flange.end_inset_mm > 1e-9 {
        points.push(at(f.width, 0.));
        if let Some([w, depth]) = relief {
            if w > f.flange.end_inset_mm + 1e-9 {
                return Err(format!(
                    "Flange {} relief is wider than its end inset",
                    f.flange.index + 1
                ));
            }
            points.extend([
                at(f.width, -depth),
                at(f.width + w, -depth),
                at(f.width + w, 0.),
            ]);
        }
    }
    points.push(b);
    points.dedup_by(|p, q| (p[0] - q[0]).hypot(p[1] - q[1]) < 1e-9);
    Ok(points)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn outline(points: &[[f64; 2]]) -> Vec<PlanarEdge> {
        (0..points.len())
            .map(|i| PlanarEdge::Line {
                a: points[i],
                b: points[(i + 1) % points.len()],
            })
            .collect()
    }

    #[test]
    fn developed_outlines_cross_only_where_regions_overlap() {
        // A relief notch beside a full developed flange stays simple.
        let notched = [
            [0., 0.],
            [0., -5.],
            [10., -5.],
            [10., 0.],
            [9., 0.],
            [9., 1.5],
            [10., 1.5],
            [10., 8.],
            [0., 8.],
        ];
        assert!(!crossing(&outline(&notched)));
        // Two developed legs folded into a notch reach past each other.
        let overlapping = [
            [0., 0.],
            [12., 0.],
            [12., 10.],
            [9., 10.],
            [9., 4.],
            [5., 4.],
            [5., 9.],
            [7., 9.],
            [7., 3.],
            [3., 3.],
            [3., 10.],
            [0., 10.],
        ];
        assert!(crossing(&outline(&overlapping)));
        // Consecutive edges folding back onto each other touch.
        assert!(crossing(&outline(&[
            [0., 0.],
            [4., 0.],
            [2., 0.],
            [2., 2.]
        ])));
    }

    #[test]
    fn base_openings_must_clear_the_folded_web_boundary() {
        let web = [[0., 0.], [10., 0.], [10., 8.], [0., 8.]];
        let web: Vec<_> = (0..4).map(|i| [web[i], web[(i + 1) % 4]]).collect();
        let circle = |x: f64, r: f64| {
            vec![CurvedEdge2d::Circle {
                center: [x, 4.],
                radius_mm: r,
            }]
        };
        assert!(inside(&circle(5., 2.), &web, 1e-5));
        assert!(!inside(&circle(1., 2.), &web, 1e-5));
        assert!(!inside(&circle(12., 1.), &web, 1e-5));
        let square = |x: f64| {
            let p = [[x, 3.], [x + 2., 3.], [x + 2., 5.], [x, 5.]];
            (0..4)
                .map(|i| CurvedEdge2d::Line {
                    a: p[i],
                    b: p[(i + 1) % 4],
                })
                .collect::<Vec<_>>()
        };
        assert!(inside(&square(4.), &web, 1e-5));
        assert!(!inside(&square(9.), &web, 1e-5));
        assert!(!inside(&square(8.), &web, 1e-5));
    }
}
