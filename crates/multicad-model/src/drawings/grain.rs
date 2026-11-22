//! Faint display-only wood grain on panel faces seen face-on (CAD-365).
//!
//! Grain follows each panel's canonical grain axis (`GrainAxis` on the panel
//! definition, resolved through its carcass), carried through the panel
//! Part's later Transform features. Lines are clipped to the exact face
//! contours, so holes and cutouts stay clear, and to nearer face-on faces of
//! other bodies. Panels with free grain ("Any") have no direction and get no
//! lines. Nothing here changes geometry or manufacturing data.
use super::*;
use multicad_core::{
    manufacturing::{PlanarEdge, PlanarLoop},
    model::{self, GrainAxis},
};

/// Light gray, so outlines, dimensions and text stay readable.
pub const GRAIN_RGB: [u8; 3] = [196, 196, 196];
pub const GRAIN_LINE_MM: f64 = 0.08;
/// Paper spacing between grain lines.
pub const GRAIN_SPACING_MM: f64 = 2.5;
/// Per-view cap; spacing widens rather than exceeding it.
const MAX_GRAIN_LINES: usize = 6000;
const PARALLEL: f64 = 1. - 1e-6;

type P2 = [f64; 2];
type P3 = [f64; 3];
fn dot(a: P3, b: P3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn sub(a: P3, b: P3) -> P3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn cross(a: P3, b: P3) -> P3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn rotate(r: &[[f64; 3]; 3], v: P3) -> P3 {
    [dot(r[0], v), dot(r[1], v), dot(r[2], v)]
}

/// One face-on face region in view coordinates, with its depth toward the
/// viewer and the polygonized loops (outer and inner) for even-odd tests.
struct Region {
    body: usize,
    depth: f64,
    loops: Vec<Vec<P2>>,
}

/// The panel's world grain direction, or `None` when the body is not a
/// panel, its grain is free, or its frame is unavailable.
fn panel_grain(
    model: &Model,
    document: &Document,
    parameters: &multicad_core::parameters::EvaluatedParameters,
    part: EntityId,
) -> Result<Option<(P3, P3)>, ModelError> {
    let Some(feature) = model::panel_feature(document, part) else {
        return Ok(None);
    };
    let source = model::read_panel(document, feature).map_err(|e| error(Some(feature), e))?;
    let resolved =
        model::resolve_panel(document, parameters, &source).map_err(|e| error(Some(feature), e))?;
    let Some(frame) = model.snapshot().panel_frames.get(&feature).copied() else {
        return Ok(None);
    };
    let mut grain = match resolved.grain {
        GrainAxis::U => frame.u,
        GrainAxis::V => frame.v,
        GrainAxis::Free => return Ok(None),
    };
    let mut normal = cross(frame.u, frame.v);
    // Later Transform features of the same Part move or rotate the panel;
    // the grain turns with it, in the Part's feature order.
    if let Some(EntityData::Part { features, .. }) = document.get(part).map(|e| &e.data) {
        for id in features.iter().skip_while(|f| **f != feature).skip(1) {
            let is_transform = matches!(
                document.get(*id).map(|e| &e.data),
                Some(EntityData::Feature { definition, .. })
                    if definition.type_id == model::TRANSFORM_TYPE
            );
            if !is_transform {
                continue;
            }
            let transform =
                model::read_transform(document, *id).map_err(|e| error(Some(*id), e))?;
            let pose = transform
                .parameters
                .resolve(parameters)
                .map_err(|e| error(Some(*id), e))?;
            grain = rotate(&pose.rotation, grain);
            normal = rotate(&pose.rotation, normal);
        }
    }
    Ok(Some((grain, normal)))
}

fn polygon(edges: &[PlanarEdge], chord_mm: f64) -> Vec<P2> {
    let mut points = vec![];
    for edge in edges {
        match *edge {
            PlanarEdge::Line { a, .. } => points.push(a),
            PlanarEdge::Arc {
                a,
                center,
                sweep_rad,
                ..
            } => {
                let r = (a[0] - center[0]).portable_hypot(a[1] - center[1]);
                let step = if r > chord_mm {
                    2. * (1. - chord_mm / r).clamp(-1., 1.).portable_acos()
                } else {
                    std::f64::consts::FRAC_PI_4
                };
                let n = ((sweep_rad.abs() / step.max(1e-3)).ceil() as usize).clamp(2, 256);
                let start = (a[1] - center[1]).portable_atan2(a[0] - center[0]);
                for i in 0..n {
                    let t = start + sweep_rad * i as f64 / n as f64;
                    points.push([
                        center[0] + r * t.portable_cos(),
                        center[1] + r * t.portable_sin(),
                    ]);
                }
            }
        }
    }
    points
}

/// Parameters along `origin + t·dir` where the line crosses the loops.
fn crossings(origin: P2, dir: P2, loops: &[Vec<P2>]) -> Vec<f64> {
    let normal = [-dir[1], dir[0]];
    let side = |p: P2| (p[0] - origin[0]) * normal[0] + (p[1] - origin[1]) * normal[1];
    let along = |p: P2| (p[0] - origin[0]) * dir[0] + (p[1] - origin[1]) * dir[1];
    let mut out = vec![];
    for ring in loops {
        for i in 0..ring.len() {
            let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
            let (sa, sb) = (side(a), side(b));
            // Half-open rule counts a vertex on the line once.
            if (sa > 0.) != (sb > 0.) {
                let k = sa / (sa - sb);
                out.push(along(a) + k * (along(b) - along(a)));
            }
        }
    }
    out.sort_by(f64::total_cmp);
    out
}

/// Inside intervals of a line by the even-odd rule.
fn inside(origin: P2, dir: P2, loops: &[Vec<P2>]) -> Vec<[f64; 2]> {
    crossings(origin, dir, loops)
        .as_chunks::<2>()
        .0
        .iter()
        .copied()
        .filter(|c| c[1] - c[0] > 1e-9)
        .collect()
}

fn subtract(spans: Vec<[f64; 2]>, cut: &[[f64; 2]]) -> Vec<[f64; 2]> {
    let mut out = spans;
    for c in cut {
        out = out
            .into_iter()
            .flat_map(|s| {
                let mut parts = vec![];
                if c[0] > s[0] {
                    parts.push([s[0], c[0].min(s[1])]);
                }
                if c[1] < s[1] {
                    parts.push([c[1].max(s[0]), s[1]]);
                }
                parts.into_iter().filter(|p| p[1] - p[0] > 1e-9)
            })
            .collect();
    }
    out
}

impl Model {
    /// Grain vectors for one view; empty when no visible face-on panel has a
    /// fixed grain. Only plain projection views get grain.
    pub(super) fn grain_vectors(
        &self,
        document: &Document,
        view: &DrawingView,
        providers: &DrawingProviders,
    ) -> Result<Vec<DrawingVector>, ModelError> {
        let guided;
        let view = if view.guide.is_some() {
            let parameters = multicad_core::parameters::evaluate(document)
                .map_err(|e| error(Some(view.id), e))?;
            guided = multicad_core::drawings::resolve_guided(document, &parameters, view)
                .map_err(|e| error(Some(view.id), e))?;
            &guided
        } else {
            view
        };
        if view.kind != ViewKind::Projection {
            return Ok(vec![]);
        }
        let parameters =
            multicad_core::parameters::evaluate(document).map_err(|e| error(Some(view.id), e))?;
        let bodies = providers.resolve(document, self, &view.source)?;
        let (u, v) = (view.frame.u, view.frame.v);
        let w = cross(u, v);
        let origin = view.frame.origin_mm;
        let scale = view.scale.numerator as f64 / view.scale.denominator as f64;
        // Finer than the outline tolerance: hole polygons are inscribed, so
        // grain stays out of each hole to within this chord.
        let chord = (0.1 * view.paper_chord_tolerance_mm / scale).max(1e-4);
        let to_view = |p: P3| [dot(sub(p, origin), u), dot(sub(p, origin), v)];
        let mut regions = vec![];
        let mut grains = vec![];
        for (index, body) in bodies.iter().enumerate() {
            let solid = body.shape.solid();
            let grain = panel_grain(self, document, &parameters, body.material_part)?
                .filter(|(_, normal)| dot(*normal, w).abs() > PARALLEL);
            for face in solid.faces().map_err(|e| error(Some(body.identity), e))? {
                let Some((point, normal)) = face.plane else {
                    continue;
                };
                // Faces facing the viewer: on panels they carry grain, on
                // every body they hide what lies behind them.
                if dot(normal, w) < PARALLEL {
                    continue;
                }
                let frame = PlaneFrame {
                    origin: point,
                    u,
                    v,
                };
                let Ok(contours) = solid.planar_face_contours(face.index, frame) else {
                    continue;
                };
                let offset = to_view(point);
                let loops = contours
                    .loops
                    .iter()
                    .map(|l: &PlanarLoop| {
                        polygon(&l.edges, chord)
                            .into_iter()
                            .map(|p| [p[0] + offset[0], p[1] + offset[1]])
                            .collect::<Vec<_>>()
                    })
                    .filter(|l| l.len() >= 3)
                    .collect::<Vec<_>>();
                if loops.is_empty() {
                    continue;
                }
                if let Some((g, _)) = grain {
                    grains.push((regions.len(), [dot(g, u), dot(g, v)]));
                }
                regions.push(Region {
                    body: index,
                    depth: dot(point, w),
                    loops,
                });
            }
        }
        let mut spacing = GRAIN_SPACING_MM / scale;
        let estimate = |spacing: f64| -> f64 {
            grains
                .iter()
                .map(|(r, _)| {
                    let pts = regions[*r].loops.iter().flatten();
                    let (mut lo, mut hi) = ([f64::MAX; 2], [f64::MIN; 2]);
                    for p in pts {
                        lo = [lo[0].min(p[0]), lo[1].min(p[1])];
                        hi = [hi[0].max(p[0]), hi[1].max(p[1])];
                    }
                    (hi[0] - lo[0]).portable_hypot(hi[1] - lo[1]) / spacing
                })
                .sum()
        };
        while estimate(spacing) > MAX_GRAIN_LINES as f64 {
            spacing *= 1.5;
        }
        let mut out = vec![];
        for (r, g2) in grains {
            let len = g2[0].portable_hypot(g2[1]);
            if len < 1e-9 {
                continue;
            }
            let dir = [g2[0] / len, g2[1] / len];
            let across = [-dir[1], dir[0]];
            let region = &regions[r];
            let points = region.loops.iter().flatten();
            let (mut lo, mut hi) = (f64::MAX, f64::MIN);
            for p in points {
                let s = p[0] * across[0] + p[1] * across[1];
                lo = lo.min(s);
                hi = hi.max(s);
            }
            // Lines sit on a view-fixed lattice, so panels sharing an axis line up.
            let mut k = (lo / spacing).floor() + 0.5;
            while k * spacing < hi {
                let s = k * spacing;
                k += 1.;
                let line_origin = [across[0] * s, across[1] * s];
                let mut spans = inside(line_origin, dir, &region.loops);
                for other in regions
                    .iter()
                    .filter(|o| o.depth > region.depth + 1e-6 && o.body != region.body)
                {
                    spans = subtract(spans, &inside(line_origin, dir, &other.loops));
                }
                for [t0, t1] in spans {
                    let a = [line_origin[0] + dir[0] * t0, line_origin[1] + dir[1] * t0];
                    let b = [line_origin[0] + dir[0] * t1, line_origin[1] + dir[1] * t1];
                    out.push(DrawingVector {
                        color_rgb: Some(GRAIN_RGB),
                        source: bodies[region.body].identity,
                        material_part: bodies[region.body].material_part,
                        role: VectorRole::Grain,
                        line_width_mm: GRAIN_LINE_MM,
                        curve: VectorCurve::Line { a, b },
                        preview_mm: vec![a, b],
                    });
                }
            }
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn even_odd_clipping_keeps_holes_clear_and_subtracts_occluders() {
        let square =
            |x0: f64, y0: f64, s: f64| vec![[x0, y0], [x0 + s, y0], [x0 + s, y0 + s], [x0, y0 + s]];
        let loops = vec![square(0., 0., 10.), square(4., 4., 2.)];
        let spans = inside([0., 5.], [1., 0.], &loops);
        assert_eq!(spans, vec![[0., 4.], [6., 10.]]);
        assert_eq!(inside([0., 1.], [1., 0.], &loops), vec![[0., 10.]]);
        assert!(inside([0., 11.], [1., 0.], &loops).is_empty());
        let cut = inside([0., 5.], [1., 0.], &[square(8., 0., 6.)]);
        assert_eq!(subtract(spans, &cut), vec![[0., 4.], [6., 8.]]);
    }
    #[test]
    fn arcs_polygonize_within_the_chord_tolerance() {
        let circle = [PlanarEdge::Arc {
            a: [5., 0.],
            b: [5., 0.],
            center: [0., 0.],
            sweep_rad: std::f64::consts::TAU,
        }];
        let points = polygon(&circle, 0.01);
        assert!(points.len() > 16);
        assert!(
            points
                .iter()
                .all(|p| (p[0].portable_hypot(p[1]) - 5.).abs() < 1e-9)
        );
    }
}
