//! Hidden-line and section drawing of a scene of exact members.
use super::arrange::{Limit, Strand, arrange, clip, distance, locate};
use super::paths::{Lift, ProfilePath, Segment, add, dot3, scaled};
use super::publish::Output;
use super::*;
use multicad_rust_kernel::{normal_sweep, spur_gear};
use std::collections::HashMap;

pub(in crate::drawing) struct SceneMember {
    pub source: usize,
    pub member: Box<dyn Member>,
    pub hatch: Option<DrawingHatch>,
}

struct Scene<'a> {
    members: &'a [SceneMember],
    /// Retained depth range.
    depth: [f64; 2],
    eps: f64,
    /// Shortest occluding material length that hides a point.
    cover: f64,
}

impl Scene<'_> {
    fn retained(&self, member: usize, paper: P2) -> Result<Vec<[f64; 2]>, GeometryError> {
        Ok(self.members[member]
            .member
            .sight(paper)?
            .into_iter()
            .map(|[a, b]| [a.max(self.depth[0]), b.min(self.depth[1])])
            .filter(|[a, b]| b > a)
            .collect())
    }

    fn hidden(&self, point: P3) -> Result<bool, GeometryError> {
        for member in 0..self.members.len() {
            if self
                .retained(member, [point[0], point[1]])?
                .into_iter()
                .any(|[a, b]| b - a.max(point[2]) > self.cover)
            {
                return Ok(true);
            }
        }
        Ok(false)
    }

    fn inside(&self, member: usize, point: P3) -> Result<bool, GeometryError> {
        Ok(self.members[member]
            .member
            .sight([point[0], point[1]])?
            .into_iter()
            .any(|[a, b]| a < point[2] && point[2] < b))
    }
}

fn paper(p: P3) -> P2 {
    [p[0], p[1]]
}

fn cross3(a: P3, b: P3) -> P3 {
    std::array::from_fn(|i| a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3])
}

/// A range of a border path.
type Stretch = (Arc<dyn Path>, [f64; 2]);

struct Piece {
    strand: usize,
    range: [f64; 2],
    ends: [P2; 2],
    mid: P2,
    hidden: bool,
}

/// Draw the members as one scene. `tolerance` is the linear modelling
/// tolerance of the sources; nothing is returned unless every member answers.
pub(in crate::drawing) fn project(
    members: &[SceneMember],
    options: DrawingProjection,
    policy: GeometryPolicy,
) -> Result<Vec<ProjectedVector>, GeometryError> {
    let tolerance = policy.linear_tolerance_mm;
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    let bounds: Vec<[P3; 2]> = members.iter().map(|m| m.member.bounds()).collect();
    for [a, b] in &bounds {
        for k in 0..3 {
            lo[k] = lo[k].min(a[k]);
            hi[k] = hi[k].max(b[k]);
        }
    }
    let scale = lo.iter().chain(&hi).fold(1_f64, |s, v| s.max(v.abs()));
    if !scale.is_finite() {
        return Err(err(
            GeometryErrorKind::InvalidResult,
            "Unbounded drawing scene member",
        ));
    }
    let eps = (4. * tolerance).max(1e-10 * scale);
    let top = options.section_offset_mm;
    let far = top
        .zip(options.section_depth_mm)
        .map(|(top, depth)| top - depth);
    let scene = Scene {
        members,
        depth: [
            far.unwrap_or(f64::NEG_INFINITY),
            top.unwrap_or(f64::INFINITY),
        ],
        eps,
        cover: 1e-6 * scale,
    };
    let windows = [options.section_u_mm, options.section_v_mm];

    // Limits of the retained region, each with the paths of its border.
    let mut limits: Vec<Limit> = Vec::new();
    let plane = |normal: P3, level: f64| Limit {
        value: Box::new(move |p| dot3(normal, p) - level),
        inward: Some(scaled(normal, -1.)),
        disc: None,
    };
    if let Some(top) = top {
        limits.push(plane([0., 0., 1.], top));
    }
    if let Some(far) = far {
        limits.push(plane([0., 0., -1.], -far));
    }
    let lateral_start = limits.len();
    for (axis, window) in windows.iter().enumerate() {
        if let Some([a, b]) = *window {
            let normal: P3 = std::array::from_fn(|i| f64::from(u8::from(i == axis)));
            limits.push(plane(scaled(normal, -1.), -a));
            limits.push(plane(normal, b));
        }
    }
    let detail_limit = limits.len();
    if let Some((center, radius)) = options.detail {
        limits.push(Limit {
            value: Box::new(move |p| {
                let r = (p[0] - center[0]).hypot(p[1] - center[1]);
                r - radius
            }),
            inward: None,
            disc: Some((center, radius)),
        });
    }

    let flat = options.chord_tolerance_mm / 4.;
    // Hatch regions follow the established sampling of exact section rims.
    let sampled = (options.chord_tolerance_mm / 2.).min(policy.minimum_feature_mm / 4.);
    let mut output = Output::new(options);
    let cut = |index: usize| {
        // A limit removes material when the member reaches beyond it.
        let [a, b] = bounds[index];
        top.is_some_and(|top| b[2] > top + eps)
            || far.is_some_and(|far| a[2] < far - eps)
            || windows.iter().enumerate().any(|(axis, window)| {
                window.is_some_and(|[lo, hi]| a[axis] < lo - eps || b[axis] > hi + eps)
            })
    };

    // Section loops of every member, clipped to the window and detail.
    let all: Vec<&Limit> = limits.iter().collect();
    let lateral = &all[lateral_start..];
    let mut strands: Vec<Strand> = Vec::new();
    let mut sections: Vec<Vec<Strand>> = Vec::new();
    for (index, entry) in members.iter().enumerate() {
        let mut own = Vec::new();
        if let Some(top) = top {
            for path in entry.member.section(top, flat)?.into_iter().flatten() {
                for (ts, pts, _) in clip(&path, lateral, flat, eps) {
                    own.push(Strand {
                        member: index,
                        path: path.clone(),
                        ts,
                        pts,
                        plane: Some(0),
                    });
                }
            }
        }
        sections.push(own);
    }

    let mut seams: Vec<Vec<Stretch>> = vec![Vec::new(); members.len()];
    if !options.section_only {
        for (index, entry) in members.iter().enumerate() {
            // Nothing of a member wholly in front of the cut is retained.
            if top.is_some_and(|top| bounds[index][0][2] >= top - eps) {
                continue;
            }
            let mut paths: Vec<(Option<Kind>, Arc<dyn Path>)> = entry
                .member
                .candidates(flat)?
                .into_iter()
                .filter(|c| c.kind != Kind::Smooth || options.tangent_edges)
                .map(|c| (Some(c.kind), c.path))
                .collect();
            for level in [top, far].into_iter().flatten() {
                for path in entry.member.section(level, flat)?.into_iter().flatten() {
                    paths.push((None, path));
                }
            }
            for (kind, path) in paths {
                for (ts, pts, plane) in clip(&path, &all, flat, eps) {
                    let strand = Strand {
                        member: index,
                        path: path.clone(),
                        ts,
                        pts,
                        plane,
                    };
                    let length: f64 = strand
                        .pts
                        .windows(2)
                        .map(|pair| distance(paper(pair[0]), paper(pair[1])))
                        .sum();
                    if length <= eps {
                        continue;
                    }
                    let (edge, interior) = on_boundary(&scene, &strand, kind, &limits)?;
                    if interior && plane == Some(0) && top.is_some() {
                        // A source edge inside the cut face also bounds the
                        // two coplanar faces of the section.
                        let range = [strand.ts[0], *strand.ts.last().unwrap()];
                        seams[index].push((strand.path.clone(), range));
                    }
                    if edge {
                        strands.push(strand);
                    }
                }
            }
        }
    }

    // Borders of the lateral window and the detail circle.
    let margin = 1. + 0.01 * scale;
    let mut borders: Vec<(Strand, bool)> = Vec::new();
    let level = top.unwrap_or(0.);
    for (axis, window) in windows.iter().enumerate() {
        let Some(window) = *window else { continue };
        let other = windows[1 - axis].unwrap_or([lo[1 - axis] - margin, hi[1 - axis] + margin]);
        for (side, value) in window.into_iter().enumerate() {
            let at = |s: f64| -> P3 {
                let mut p = [0., 0., level];
                p[axis] = value;
                p[1 - axis] = s;
                p
            };
            let path: Arc<dyn Path> = Arc::new(Segment(at(other[0]), at(other[1])));
            let skip = lateral_start + 2 * usize::from(axis == 1 && windows[0].is_some()) + side;
            let rest: Vec<&Limit> = limits
                .iter()
                .enumerate()
                .filter(|(i, _)| *i >= lateral_start && *i != skip)
                .map(|(_, limit)| limit)
                .collect();
            borders.extend(
                border_strands(&path, &rest, flat, eps)
                    .into_iter()
                    .map(|s| (s, true)),
            );
        }
    }
    if let Some((center, radius)) = options.detail
        && top.is_some()
    {
        let circle = spur_gear::Edge::Planar(normal_sweep::Edge::Circle {
            center,
            radius_mm: radius,
            clockwise: false,
        });
        for range in [[0., 0.5], [0.5, 1.]] {
            let path: Arc<dyn Path> = Arc::new(ProfilePath {
                edge: circle,
                range,
                lift: Lift {
                    origin: [0., 0., level],
                    x: [1., 0., 0.],
                    y: [0., 1., 0.],
                },
                straight: false,
            });
            let rest: Vec<&Limit> = limits[lateral_start..detail_limit].iter().collect();
            borders.extend(
                border_strands(&path, &rest, flat, eps)
                    .into_iter()
                    .map(|s| (s, false)),
            );
        }
    }

    // Stretches of each border between the events of every boundary on it.
    let mut section_extra: Vec<Vec<Stretch>> = vec![Vec::new(); members.len()];
    let mut laterals: Vec<(usize, Arc<dyn Path>, bool)> = Vec::new();
    for (border, straight) in &borders {
        let mut stops = vec![border.ts[0], *border.ts.last().unwrap()];
        for strand in strands.iter().chain(sections.iter().flatten()) {
            for end in [strand.pts[0], *strand.pts.last().unwrap()] {
                for i in 0..border.pts.len() - 1 {
                    let t = locate(border, i, paper(end), eps);
                    if distance(paper(border.path.point(t)), paper(end)) <= 2. * eps {
                        stops.push(t);
                    }
                }
            }
        }
        stops.sort_by(f64::total_cmp);
        stops.dedup_by(|a, b| {
            distance(paper(border.path.point(*a)), paper(border.path.point(*b))) <= eps
        });
        // Consecutive stretches with the same answer publish as one curve.
        let mut runs: HashMap<(usize, u8), [f64; 2]> = HashMap::new();
        let flush =
            |key: (usize, u8), range: [f64; 2], laterals: &mut Vec<_>, extra: &mut Vec<Vec<_>>| {
                let (member, class) = key;
                if class == 2 {
                    extra[member].push((border.path.clone(), range));
                } else {
                    let [a, b] = range.map(|t| border.path.point(t));
                    laterals.push((member, Arc::new(Segment(a, b)) as Arc<dyn Path>, class == 1));
                }
            };
        for pair in stops.windows(2) {
            let mid = border.path.point(0.5 * (pair[0] + pair[1]));
            for (member, entry) in members.iter().enumerate() {
                let mut found: Vec<u8> = Vec::new();
                let spans = scene.retained(member, paper(mid))?;
                if let Some(top) = top
                    && (cut(member) || options.section_only)
                    && entry
                        .member
                        .sight(paper(mid))?
                        .into_iter()
                        .any(|[a, b]| a < top - eps && b >= top - eps)
                {
                    found.push(2);
                }
                if *straight
                    && !options.section_only
                    && let Some(last) = spans.last()
                {
                    let hidden = scene.hidden([mid[0], mid[1], last[1]])?;
                    found.push(u8::from(hidden));
                }
                for class in 0..3 {
                    let key = (member, class);
                    if found.contains(&class) {
                        runs.entry(key)
                            .and_modify(|range| range[1] = pair[1])
                            .or_insert([pair[0], pair[1]]);
                    } else if let Some(range) = runs.remove(&key) {
                        flush(key, range, &mut laterals, &mut section_extra);
                    }
                }
            }
        }
        let mut rest: Vec<_> = runs.into_iter().collect();
        rest.sort_by_key(|(key, _)| *key);
        for (key, range) in rest {
            flush(key, range, &mut laterals, &mut section_extra);
        }
    }

    if !options.section_only {
        let stops = arrange(&strands, eps);
        let mut pieces: Vec<Piece> = Vec::new();
        for (index, (strand, stops)) in strands.iter().zip(&stops).enumerate() {
            for pair in stops.windows(2) {
                let mid = strand.path.point(0.5 * (pair[0].0 + pair[1].0));
                let ends = [paper(pair[0].1), paper(pair[1].1)];
                if distance(ends[0], ends[1]).max(distance(ends[0], paper(mid))) <= eps {
                    continue;
                }
                let hidden = scene.hidden(mid)?;
                if hidden && !options.hidden_lines {
                    continue;
                }
                pieces.push(Piece {
                    strand: index,
                    range: [pair[0].0, pair[1].0],
                    ends,
                    mid: paper(mid),
                    hidden,
                });
            }
        }
        // Coincident pieces of one source publish once; visible wins.
        let mut order: Vec<usize> = (0..pieces.len()).collect();
        order.sort_by_key(|&i| pieces[i].hidden);
        let mut kept = vec![false; pieces.len()];
        let mut index: HashMap<[i64; 2], Vec<usize>> = HashMap::new();
        let key = |p: P2| p.map(|v| (v / (16. * eps)).floor() as i64);
        // Pieces are found by their lexicographically first end.
        let first = |piece: &Piece| {
            let [a, b] = piece.ends;
            if (a[0], a[1]) <= (b[0], b[1]) { a } else { b }
        };
        let conic = |piece: &Piece| match strands[piece.strand].path.paper(piece.range) {
            Some(Curve::Ellipse(ellipse)) => Some(ellipse),
            _ => None,
        };
        for i in order {
            let piece = &pieces[i];
            let source = members[strands[piece.strand].member].source;
            let [x, y] = key(first(piece));
            let chord = distance(piece.ends[0], piece.ends[1]);
            let duplicate = (-1..=1)
                .flat_map(|dx| (-1..=1).map(move |dy| [x + dx, y + dy]))
                .any(|cell| {
                    index.get(&cell).is_some_and(|others| {
                        others.iter().any(|&j| {
                            let other: &Piece = &pieces[j];
                            let ends = (distance(other.ends[0], piece.ends[0]) <= 2. * eps
                                && distance(other.ends[1], piece.ends[1]) <= 2. * eps)
                                || (distance(other.ends[0], piece.ends[1]) <= 2. * eps
                                    && distance(other.ends[1], piece.ends[0]) <= 2. * eps);
                            members[strands[other.strand].member].source == source
                                && ends
                                && (distance(other.mid, piece.mid) <= 2. * eps
                                    // The same arc of one conic under another
                                    // parametrisation: its middle lies on the
                                    // same side of the common chord.
                                    || (chord > 2. * eps
                                        && distance(other.mid, piece.mid) < 0.5 * chord
                                        && conic(piece).zip(conic(other)).is_some_and(
                                            |(a, b)| arrange::same_conic(a, b, eps),
                                        )))
                        })
                    })
                });
            if !duplicate {
                kept[i] = true;
                index.entry([x, y]).or_default().push(i);
            }
        }
        // Adjacent pieces of one path with the same role publish together.
        let mut i = 0;
        while i < pieces.len() {
            if !kept[i] {
                i += 1;
                continue;
            }
            let first = &pieces[i];
            let mut range = first.range;
            let mut j = i + 1;
            while j < pieces.len()
                && kept[j]
                && pieces[j].strand == first.strand
                && pieces[j].hidden == first.hidden
                && pieces[j].range[0] == range[1]
            {
                range[1] = pieces[j].range[1];
                j += 1;
            }
            let strand = &strands[first.strand];
            output.curve(
                members[strand.member].source,
                if first.hidden {
                    VectorRole::Hidden
                } else {
                    VectorRole::Visible
                },
                strand.path.as_ref(),
                range,
            )?;
            i = j;
        }
        for (member, path, hidden) in &laterals {
            if !*hidden || options.hidden_lines {
                output.curve(
                    members[*member].source,
                    if *hidden {
                        VectorRole::Hidden
                    } else {
                        VectorRole::Visible
                    },
                    path.as_ref(),
                    [0., 1.],
                )?;
            }
        }
    }

    // Cut faces: the exact boundary, then the hatch of each source's region.
    let mut regions: Vec<(usize, DrawingHatch, Vec<[P2; 2]>)> = Vec::new();
    for (index, entry) in members.iter().enumerate() {
        if top.is_none()
            || !(cut(index) || options.section_only)
            || (!options.section_only && top.is_some_and(|top| bounds[index][0][2] >= top - eps))
        {
            continue;
        }
        let mut boundary = Vec::new();
        for strand in &sections[index] {
            let range = [strand.ts[0], *strand.ts.last().unwrap()];
            output.curve(
                entry.source,
                VectorRole::Section,
                strand.path.as_ref(),
                range,
            )?;
            boundary.extend(output.chords(strand.path.as_ref(), range, sampled));
        }
        for (path, range) in &seams[index] {
            output.curve(entry.source, VectorRole::Section, path.as_ref(), *range)?;
        }
        for (path, range) in &section_extra[index] {
            output.curve(entry.source, VectorRole::Section, path.as_ref(), *range)?;
            boundary.extend(output.chords(path.as_ref(), *range, sampled));
        }
        if let Some(hatch) = entry.hatch
            && !boundary.is_empty()
        {
            match regions
                .iter_mut()
                .find(|(source, ..)| *source == entry.source)
            {
                Some(region) => region.2.extend(boundary),
                None => regions.push((entry.source, hatch, boundary)),
            }
        }
    }
    for (source, hatch, boundary) in regions {
        output.hatch(source, hatch, &boundary, tolerance, 2. * eps)?;
    }
    Ok(output.finish(eps))
}

fn border_strands(path: &Arc<dyn Path>, limits: &[&Limit], flat: f64, eps: f64) -> Vec<Strand> {
    clip(path, limits, flat, eps)
        .into_iter()
        .map(|(ts, pts, _)| Strand {
            member: 0,
            path: path.clone(),
            ts,
            pts,
            plane: None,
        })
        .collect()
}

/// A path lying in a limiting plane is an edge of the retained solid only
/// where material just inside the plane ends at it.
fn on_boundary(
    scene: &Scene<'_>,
    strand: &Strand,
    kind: Option<Kind>,
    limits: &[Limit],
) -> Result<(bool, bool), GeometryError> {
    let Some(inward) = strand.plane.and_then(|index| limits[index].inward) else {
        return Ok((true, false));
    };
    let [a, b] = [strand.ts[0], *strand.ts.last().unwrap()];
    let mid = 0.5 * (a + b);
    let h = 1e-3 * (b - a);
    let [p, q] = [strand.path.point(mid - h), strand.path.point(mid + h)];
    let tangent = add(q, scaled(p, -1.));
    let side = cross3(inward, tangent);
    let length = dot3(side, side).sqrt();
    if length == 0. {
        return Ok((true, false));
    }
    // Step far less into the solid than sideways, so that a wall leaving
    // the plane obliquely still separates the two probes.
    let center = add(strand.path.point(mid), scaled(inward, 4. * scene.eps));
    let at = |sign: f64| {
        scene.inside(
            strand.member,
            add(center, scaled(side, sign * 256. * scene.eps / length)),
        )
    };
    // A surface tangent to the plane keeps its outline along the contact,
    // and a source edge keeps separating the two coplanar faces it bounds
    // (as the native common of a solid and a half-space does). A section
    // loop is an edge only where the retained face ends at it.
    let [left, right] = [at(1.)?, at(-1.)?];
    Ok(match kind {
        Some(Kind::Outline) => (left || right, false),
        Some(Kind::Edge) => (left || right, left && right),
        _ => (left != right, false),
    })
}
