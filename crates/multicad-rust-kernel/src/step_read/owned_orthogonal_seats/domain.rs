//! Constructive finite planar-domain correspondence. A globally Lipschitz
//! displacement maps original stock wires and circular holes onto the SAME raw
//! projected trims. k<1 proves bijectivity/orientation and bounds inverse D1.
use super::*;
type P = [f64; 2];
type V2 = [I; 2];
#[derive(Clone)]
struct Hole {
    center: P,
    radius: f64,
    delta: [V2; 3],
    bounds: [[f64; 2]; 2],
}
fn add2(a: V2, b: V2) -> V2 {
    std::array::from_fn(|i| a[i] + b[i])
}
fn sub2(a: V2, b: V2) -> V2 {
    std::array::from_fn(|i| a[i] - b[i])
}
fn scale2(a: V2, s: I) -> V2 {
    a.map(|v| v * s)
}
fn norm2(a: V2) -> f64 {
    a.into_iter()
        .fold(I::point(0.), |s, v| s + I::point(v.magnitude()))
        .hi
}
fn abs_lower(a: I) -> f64 {
    if a.lo >= 0. {
        a.lo
    } else if a.hi <= 0. {
        -a.hi
    } else {
        0.
    }
}
fn invert(b: [[I; 2]; 2]) -> std::result::Result<[[I; 2]; 2], Error> {
    let d = b[0][0] * b[1][1] - b[0][1] * b[1][0];
    if !d.valid() || (d.lo <= 0. && d.hi >= 0.) {
        return Err(Error::Unsupported);
    }
    Ok([[b[1][1] / d, -b[0][1] / d], [-b[1][0] / d, b[0][0] / d]])
}
fn local_delta(inverse: [[I; 2]; 2], raw: Axis, delta: V) -> V2 {
    let q = [dot_v(raw.u, delta), dot_v(raw.v, delta)];
    std::array::from_fn(|i| inverse[i][0] * q[0] + inverse[i][1] * q[1])
}
fn segment_gap(center: P, a: P, b: P) -> f64 {
    let mut lower = 0f64;
    for i in 0..2 {
        let lo = a[i].min(b[i]);
        let hi = a[i].max(b[i]);
        let gap = if center[i] < lo {
            I::point(lo) - I::point(center[i])
        } else if center[i] > hi {
            I::point(center[i]) - I::point(hi)
        } else {
            I::point(0.)
        };
        lower = lower.max(gap.lo);
    }
    lower
}
fn grid_delta(
    p: P,
    rings: &[Vec<P>],
    deltas: &BTreeMap<[u64; 2], V2>,
) -> std::result::Result<V2, Error> {
    let key = p.map(|v| if v == 0. { 0 } else { v.to_bits() });
    if let Some(d) = deltas.get(&key) {
        return Ok(*d);
    }
    for ring in rings {
        for (a, b) in ring.iter().zip(ring.iter().cycle().skip(1)) {
            let axis = usize::from(a[0] == b[0]);
            let fixed = 1 - axis;
            if a[fixed] != b[fixed] {
                return Err(Error::Unsupported);
            }
            if p[fixed] == a[fixed]
                && p[axis] >= a[axis].min(b[axis])
                && p[axis] <= a[axis].max(b[axis])
            {
                let t = (I::point(p[axis]) - I::point(a[axis]))
                    / (I::point(b[axis]) - I::point(a[axis]));
                if !t.valid() {
                    return Err(Error::Unsupported);
                }
                let key = |p: P| p.map(|v| if v == 0. { 0 } else { v.to_bits() });
                return Ok(add2(
                    deltas[&key(*a)],
                    scale2(sub2(deltas[&key(*b)], deltas[&key(*a)]), t),
                ));
            }
        }
    }
    Ok([I::point(0.); 2])
}
fn grid_slope(
    nodes: [(P, V2); 2],
    axis: usize,
    span: I,
    rings: &[Vec<P>],
    deltas: &BTreeMap<[u64; 2], V2>,
    budget: &mut Budget<'_>,
) -> std::result::Result<V2, Error> {
    let [(p, dp), (q, dq)] = nodes;
    let fixed = 1 - axis;
    for ring in rings {
        for (a, b) in ring.iter().zip(ring.iter().cycle().skip(1)) {
            budget.charge(32)?;
            let within = |v: f64| v >= a[axis].min(b[axis]) && v <= a[axis].max(b[axis]);
            if a[fixed] == b[fixed]
                && a[axis] != b[axis]
                && p[fixed] == a[fixed]
                && q[fixed] == a[fixed]
                && within(p[axis])
                && within(q[axis])
            {
                let length = I::point(b[axis]) - I::point(a[axis]);
                if length.valid() && (length.lo > 0. || length.hi < 0.) {
                    budget.charge(256)?;
                    let key = |p: P| p.map(|v| if v == 0. { 0 } else { v.to_bits() });
                    // Both exact nodal values use these SAME original wire
                    // endpoints. Their affine interpolation difference divided
                    // by the grid span is (Db-Da)/(b_axis-a_axis), including
                    // reversed wire orientation. This cancellation preserves
                    // shared coefficient dependence without merging grid rows.
                    return Ok(sub2(deltas[&key(*b)], deltas[&key(*a)]).map(|v| v / length));
                }
            }
        }
    }
    // Interior and collar transitions keep the original independent enclosure.
    Ok(sub2(dq, dp).map(|v| v / span))
}
fn gradient(x: V2, y: V2) -> f64 {
    (I::point(norm2(x)) + I::point(norm2(y))).hi
}
pub(super) fn certify(
    view: &View<'_>,
    boundary: &Boundary<'_, '_>,
    source_face: usize,
    raw_face: usize,
    vertices: &BTreeMap<u32, V>,
    chains: &Chains<'_, '_, '_>,
    budget: &mut Budget<'_>,
) -> std::result::Result<(f64, f64), Error> {
    budget.charge(8192)?;
    let face = &view.source.boundary().faces[source_face];
    let id = boundary.face_entities[raw_face];
    let [_, _, surface, forward] = boundary.shell.args(id, "ADVANCED_FACE")? else {
        return Err(Error::Invalid);
    };
    let [_, placement] = boundary.shell.args(reference(surface)?, "PLANE")? else {
        return Err(Error::Invalid);
    };
    let raw = axis(boundary, reference(placement)?, budget)?;
    let major = face
        .normal
        .iter()
        .position(|n| n.abs() == 1.)
        .ok_or(Error::Unsupported)?;
    let u = (major + 1) % 3;
    let v = (major + 2) % 3;
    let matrix = view.placement.matrix();
    let column = |j| std::array::from_fn(|i| I::point(matrix[i][j]));
    let mu = column(u);
    let mv = column(v);
    let cofactor = scale(cross_v(mu, mv), I::point(face.normal[major]));
    let unit_normal = unit_v(cofactor)?;
    let physical_raw = scale(raw.n, I::point(if sense(forward)? { 1. } else { -1. }));
    let cosine = dot_v(unit_normal, physical_raw);
    if cosine.lo <= 0. {
        return Err(Error::Unsupported);
    }
    let inverse = invert([
        [dot_v(raw.u, mu), dot_v(raw.u, mv)],
        [dot_v(raw.v, mu), dot_v(raw.v, mv)],
    ])?;
    let metric = view.placement.metric_error();
    let sigma_min = (I::point(1.) - I::point(metric)).sqrt() * I::point(cosine.lo);
    let sigma_max = (I::point(1.) + I::point(metric)).sqrt();
    if !sigma_min.valid() || sigma_min.lo <= 0. {
        return Err(Error::Unsupported);
    }
    let rings = face
        .rings
        .iter()
        .map(|r| r.iter().map(|p| [p[u], p[v]]).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    let mut deltas = BTreeMap::new();
    let mut xs = Vec::new();
    let mut ys = Vec::new();
    let mut normal_error = 0f64;
    let mut raw_min = [f64::INFINITY; 2];
    let mut raw_max = [f64::NEG_INFINITY; 2];
    let mut displacement = 0f64;
    for (p, xy) in face.rings.iter().flatten().zip(rings.iter().flatten()) {
        budget.charge(2048)?;
        let vertex = view.raw_vertices[&crate::orthogonal_boundary::point_key(*p)];
        let actual = raw_point(boundary, vertex, budget)?;
        let expected = *vertices.get(&vertex).ok_or(Error::Unsupported)?;
        let d = local_delta(inverse, raw, sub_v(actual, expected));
        displacement = displacement.max(norm2(d));
        deltas.insert(xy.map(|x| if x == 0. { 0 } else { x.to_bits() }), d);
        xs.push(xy[0]);
        ys.push(xy[1]);
        normal_error = normal_error.max(dot_v(sub_v(expected, raw.origin), raw.n).magnitude());
        for (i, basis) in [raw.u, raw.v].into_iter().enumerate() {
            let q = dot_v(sub_v(actual, raw.origin), basis);
            raw_min[i] = raw_min[i].min(q.lo);
            raw_max[i] = raw_max[i].max(q.hi);
        }
    }
    let mut holes = Vec::new();
    for (bore_index, bore) in view.source.bores().iter().enumerate() {
        for end in 0..2 {
            if bore.faces[end] != source_face {
                continue;
            }
            budget.charge(8192)?;
            let chain = &chains.bores()[view.raw_chains[bore_index]];
            let raw_end = (end != 0) ^ view.reversed(bore_index);
            let at = if raw_end {
                chain.rim_edge_entities().len() - 1
            } else {
                0
            };
            let rim = chain.rim_edge_entities()[at];
            let (raw_ring, raw_axis) = raw_ring(boundary, rim, budget)?;
            let (a, b, n) = phase(view, chain, bore_index, boundary, budget)?;
            let sign = orientation(raw_axis.n, n)?;
            let radius = view
                .source
                .seats()
                .iter()
                .find(|s| s.bore == bore_index && s.end == (end != 0))
                .map_or(bore.radius_mm, |s| s.radius_mm);
            let source = Ring {
                center: placed(view.placement, bore.centers_mm[end].map(I::point)),
                a: scale(a, I::point(radius)),
                b: scale(b, I::point(radius * sign)),
            };
            let delta = [
                local_delta(inverse, raw, sub_v(raw_ring.center, source.center)),
                local_delta(inverse, raw, sub_v(raw_ring.a, source.a)),
                local_delta(inverse, raw, sub_v(raw_ring.b, source.b)),
            ];
            let center = [bore.centers_mm[end][u], bore.centers_mm[end][v]];
            holes.push(Hole {
                center,
                radius,
                delta,
                bounds: [[0.; 2]; 2],
            });
        }
    }
    // Source winding/strict disk containment is already earned by Seats. Each
    // chosen rectangle is connected, contains its circle and crosses no stock
    // wire; hence its entire collar lies in the same source face component.
    for i in 0..holes.len() {
        let h = &holes[i];
        let mut gap = f64::INFINITY;
        for r in &rings {
            for (a, b) in r.iter().zip(r.iter().cycle().skip(1)) {
                budget.charge(256)?;
                gap = gap.min((I::point(segment_gap(h.center, *a, *b)) - I::point(h.radius)).lo);
            }
        }
        for (j, other) in holes.iter().enumerate() {
            if i == j {
                continue;
            }
            budget.charge(256)?;
            let dx = abs_lower(I::point(h.center[0]) - I::point(other.center[0]));
            let dy = abs_lower(I::point(h.center[1]) - I::point(other.center[1]));
            gap = gap.min((I::point(dx.max(dy)) - I::point(h.radius) - I::point(other.radius)).lo);
        }
        if !gap.is_finite() || gap <= 0. {
            return Err(Error::Unsupported);
        }
        let a = h.radius + gap / 4.;
        if !a.is_finite() || a <= h.radius {
            return Err(Error::Unsupported);
        }
        let bounds = std::array::from_fn(|axis| [h.center[axis] - a, h.center[axis] + a]);
        holes[i].bounds = bounds;
    }
    let mut collar_k = 0f64;
    for h in &holes {
        let mut amin = f64::INFINITY;
        let mut amax = 0f64;
        for axis in 0..2 {
            let left = I::point(h.center[axis]) - I::point(h.bounds[axis][0]);
            let right = I::point(h.bounds[axis][1]) - I::point(h.center[axis]);
            amin = amin.min(left.lo.min(right.lo));
            amax = amax.max(left.hi.max(right.hi));
        }
        xs.extend(h.bounds[0]);
        ys.extend(h.bounds[1]);
        let width = I::point(amin) - I::point(h.radius);
        if width.lo <= 0. {
            return Err(Error::Unsupported);
        }
        let e = (I::point(norm2(h.delta[0]))
            + I::point(norm2(h.delta[1]))
            + I::point(norm2(h.delta[2])))
        .hi;
        let d = (I::point(norm2(h.delta[1])) + I::point(norm2(h.delta[2]))).hi;
        displacement = displacement.max(e);
        // R(theta) is the first rectangle-ray hit. R>=amin, R<=2amax,
        // |R'|<=2amax²/amin. The radial blend is zero at the fixed rectangle
        // and one at the circle; its finite sectors are analytic. Inside the
        // excluded disk the same boundary data extends affinely.
        let radius = I::point(h.radius);
        let k = I::point(e)
            * (I::point(1.) / width
                + I::point(2.) * I::point(amax).square() / (I::point(amin) * width * radius))
            + I::point(d) / radius;
        if !k.valid() {
            return Err(Error::Unsupported);
        }
        collar_k = collar_k.max(k.hi);
    }
    // Prove actual stored rectangles avoid every stock wire and each other,
    // including rounding in their constructed finite coordinates.
    for (i, h) in holes.iter().enumerate() {
        for r in &rings {
            for (a, b) in r.iter().zip(r.iter().cycle().skip(1)) {
                budget.charge(256)?;
                let separated = (0..2).any(|axis| {
                    h.bounds[axis][1] < a[axis].min(b[axis])
                        || h.bounds[axis][0] > a[axis].max(b[axis])
                });
                if !separated {
                    return Err(Error::Unsupported);
                }
            }
        }
        for other in holes.iter().skip(i + 1) {
            budget.charge(256)?;
            if !(0..2).any(|axis| {
                h.bounds[axis][1] < other.bounds[axis][0]
                    || other.bounds[axis][1] < h.bounds[axis][0]
            }) {
                return Err(Error::Unsupported);
            }
        }
    }
    for values in [&mut xs, &mut ys] {
        values.sort_by(f64::total_cmp);
        values.dedup_by(|a, b| *a == *b);
    }
    if xs.len() < 2 || ys.len() < 2 || xs.len() > 48 || ys.len() > 48 {
        return Err(Error::Unsupported);
    }
    let count = xs.len() * ys.len();
    budget.charge(count * 1024)?;
    let mut grid = Vec::new();
    grid.try_reserve_exact(count).map_err(|_| Error::Limit)?;
    for &y in &ys {
        for &x in &xs {
            budget.charge(256)?;
            grid.push(grid_delta([x, y], &rings, &deltas)?);
        }
    }
    let mut pl_k = 0f64;
    for j in 0..ys.len() - 1 {
        for i in 0..xs.len() - 1 {
            budget.charge(2048)?;
            let dx = I::point(xs[i + 1]) - I::point(xs[i]);
            let dy = I::point(ys[j + 1]) - I::point(ys[j]);
            if dx.lo <= 0. || dy.lo <= 0. {
                return Err(Error::Unsupported);
            }
            let a = grid[j * xs.len() + i];
            let b = grid[j * xs.len() + i + 1];
            let c = grid[(j + 1) * xs.len() + i + 1];
            let d = grid[(j + 1) * xs.len() + i];
            let p = [xs[i], ys[j]];
            let q = [xs[i + 1], ys[j]];
            let r = [xs[i + 1], ys[j + 1]];
            let t = [xs[i], ys[j + 1]];
            let mut slopes = [[I::point(0.); 2]; 4];
            for (slope, (nodes, axis, span)) in slopes.iter_mut().zip([
                ([(p, a), (q, b)], 0, dx),
                ([(q, b), (r, c)], 1, dy),
                ([(t, d), (r, c)], 0, dx),
                ([(p, a), (t, d)], 1, dy),
            ]) {
                *slope = grid_slope(nodes, axis, span, &rings, &deltas, budget)?;
            }
            let [bottom, right, top, left] = slopes;
            pl_k = pl_k.max(gradient(bottom, right)).max(gradient(top, left));
        }
    }
    let k = pl_k.max(collar_k);
    if !k.is_finite() || k >= 0.25 {
        return Err(Error::Unsupported);
    }
    // PL displacement is zero throughout each collar box, so these finite
    // supports are disjoint. Clamping outside the grid rectangle is1-Lipschitz.
    // Global ||Du||<=k<1 makes I+u a proper, orientation-preserving bi-Lipschitz
    // homeomorphism. It maps every stock wire and hole exactly; winding then
    // proves complete trimmed-domain coverage, with no folds or gaps.
    let position = (I::point(normal_error) + sigma_max * I::point(displacement)).hi;
    let sine = norm(cross_v(unit_normal, physical_raw));
    let tan = I::point(sine) / I::point(cosine.lo);
    let operator = tan + sigma_max / sigma_min * I::point(k) / (I::point(1.) - I::point(k));
    let span = (0..2)
        .map(|i| (I::point(raw_max[i]) - I::point(raw_min[i])).hi)
        .fold(0., f64::max);
    let derivative = (operator * I::point(span)).hi;
    if !position.is_finite() || !derivative.is_finite() {
        return Err(Error::Unsupported);
    }
    Ok((position, derivative))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        axial_pipe_thread::authored::StockClipPlacement,
        orthogonal::{self, Box3, Operation},
        orthogonal_bore_seats::{Options as SourceOptions, Seat, Seats},
        orthogonal_boundary, orthogonal_pockets,
        step::orthogonal_seats,
    };
    use std::rc::Rc;
    fn source(center: [f64; 2], radius: f64, mouth: f64) -> Seats {
        sources(&[(center, radius, mouth)])
    }
    fn sources(holes: &[([f64; 2], f64, f64)]) -> Seats {
        let stock = Rc::new(
            orthogonal::boolean(
                &[
                    Box3 {
                        min: [0.; 3],
                        max: [40., 40., 5.],
                    },
                    Box3 {
                        min: [0.; 3],
                        max: [40., 5., 40.],
                    },
                ],
                &[],
                Operation::Union,
                1e-7,
            )
            .unwrap()
            .unwrap(),
        );
        let boundary = orthogonal_boundary::build(stock.quads()).unwrap();
        let bores = holes
            .iter()
            .map(|(center, radius, _)| {
                orthogonal_pockets::identify_through_bore(
                    &stock,
                    &boundary,
                    [-0.1, center[0], center[1]],
                    [40.1, center[0], center[1]],
                    *radius,
                    1e-7,
                    1e-5,
                )
                .unwrap()
            })
            .collect::<Vec<_>>();
        let seats = holes
            .iter()
            .enumerate()
            .map(|(bore, (_, _, mouth))| Seat {
                bore,
                end: false,
                radius_mm: *mouth,
                depth_mm: 0.5,
            })
            .collect::<Vec<_>>();
        Seats::new(
            stock,
            Rc::from(bores),
            &seats,
            SourceOptions {
                clearance_mm: 1e-7,
                minimum_feature_mm: 1e-5,
                maximum_work: 16_000_000,
                maximum_bytes: 64 * 1024 * 1024,
            },
            &|| false,
        )
        .unwrap()
    }
    #[test]
    fn multiple_holes_on_each_nonconvex_cap_preserve_disjoint_collars_and_owners() {
        let source = sources(&[
            ([12., 2.5], 0.5, 1.),
            ([22., 2.5], 0.5, 1.),
            ([2.5, 12.], 0.5, 1.),
            ([2.5, 22.], 0.5, 1.),
        ]);
        for (label, pose) in [
            (
                "identity",
                StockClipPlacement::new([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]], [0.; 3])
                    .unwrap(),
            ),
            (
                "arbitrary",
                StockClipPlacement::new(
                    [[0.36, -0.48, 0.8], [0.8, 0.6, 0.], [-0.48, 0.64, 0.6]],
                    [11., -7., 5.],
                )
                .unwrap(),
            ),
        ] {
            let candidate =
                orthogonal_seats::candidate(&source, pose, Default::default(), &|| false).unwrap();
            let document =
                CompoundDocument::parse(candidate.as_str(), Default::default(), None).unwrap();
            let boundary = decode(
                &document,
                document.represented_members()[0].item_entity(),
                Options::default(),
                &|| false,
            )
            .unwrap();
            let material = boundary
                .heal_material(MaterialOptions::default(), &|| false)
                .expect("public four-hole healing must succeed for both retained poses");
            assert_eq!(material.chains().bores().len(), 4);
            let caps = material.chains().bores()[0].cap_face_entities();
            assert!(
                material
                    .chains()
                    .bores()
                    .iter()
                    .all(|b| b.cap_face_entities() == caps)
            );
            let matched = material
                .check_original(
                    &source,
                    pose,
                    candidate.face_entities(),
                    TransferOptions::default(),
                    &|| false,
                )
                .unwrap();
            assert!(matched.source_to_healed_position_error_mm() < 1e-7);
            assert!(matched.source_to_healed_normalized_derivative_error_mm() < 1e-7);
            eprintln!(
                "four-hole {label}: C0={} D1={} material_work={} source_match_work={}",
                matched.source_to_healed_position_error_mm(),
                matched.source_to_healed_normalized_derivative_error_mm(),
                material.work_used(),
                matched.work_used()
            );
            if label == "arbitrary" {
                let bores = material.healed_source().bores();
                assert!(
                    bores
                        .iter()
                        .enumerate()
                        .any(
                            |(i, a)| bores.iter().skip(i + 1).any(|b| (0..3).any(|axis| {
                                let difference =
                                    (a.centers_mm[0][axis] - b.centers_mm[0][axis]).abs();
                                difference > 0. && difference < 1e-12
                            }))
                        ),
                    "the distinct near-aligned retained coordinates must survive healing"
                );
            }
            assert!(material.healed_volume_bounds_mm3()[0] > 0.);
            let qualified =
                orthogonal_seats::qualify(&source, pose, Default::default(), &|| false).unwrap();
            let actual = qualified
                .verify_actual_member(
                    &document,
                    boundary.item_entity(),
                    Default::default(),
                    &|| false,
                )
                .unwrap();
            assert!(actual.literal_binding_work_used() > 0);
            assert!(actual.literal_binding_work_used() <= actual.work_used());

            // These cylinders share the same two cap owners. Actual wall/seat
            // ownership must still distinguish their different finite holes.
            let first_wall = source.boundary().faces.len();
            for offset in [first_wall, first_wall + source.bores().len()] {
                let mut wrong = candidate.face_entities().to_vec();
                wrong.swap(offset, offset + 1);
                assert!(matches!(
                    material.check_original(
                        &source,
                        pose,
                        &wrong,
                        TransferOptions::default(),
                        &|| false
                    ),
                    Err(Error::Unsupported)
                ));
            }
            let exact = MaterialOptions {
                transfer: TransferOptions {
                    maximum_work: material.work_used(),
                    maximum_bytes: material.query_budget_bytes(),
                    ..TransferOptions::default()
                },
                ..MaterialOptions::default()
            };
            assert!(boundary.heal_material(exact, &|| false).is_ok());
            for (work, bytes) in [
                (
                    exact.transfer.maximum_work - 1,
                    exact.transfer.maximum_bytes,
                ),
                (
                    exact.transfer.maximum_work,
                    exact.transfer.maximum_bytes - 1,
                ),
            ] {
                let refused = MaterialOptions {
                    transfer: TransferOptions {
                        maximum_work: work,
                        maximum_bytes: bytes,
                        ..exact.transfer
                    },
                    ..exact
                };
                assert!(matches!(
                    boundary.heal_material(refused, &|| false),
                    Err(Error::Limit)
                ));
            }
            assert!(
                material
                    .check_original(
                        &source,
                        pose,
                        candidate.face_entities(),
                        TransferOptions::default(),
                        &|| false
                    )
                    .is_ok()
            );
        }
    }
    #[test]
    fn multiple_holes_with_clear_disks_but_overlapping_collar_boxes_refuse() {
        // Both circular mouths fit and are strictly disjoint, while their
        // bounding rectangles overlap in both axes. No rectangular collar
        // certificate can be inferred from circle separation alone.
        let source = sources(&[([12., 1.25], 0.5, 1.), ([13.5, 2.75], 0.5, 1.)]);
        let pose =
            StockClipPlacement::new([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]], [0.; 3]).unwrap();
        let candidate =
            orthogonal_seats::candidate(&source, pose, Default::default(), &|| false).unwrap();
        let document =
            CompoundDocument::parse(candidate.as_str(), Default::default(), None).unwrap();
        let boundary = decode(
            &document,
            document.represented_members()[0].item_entity(),
            Options::default(),
            &|| false,
        )
        .unwrap();
        assert_eq!(
            boundary
                .through_chains(Options::default(), &|| false)
                .unwrap()
                .bores()
                .len(),
            2
        );
        assert!(source.volume_bounds_mm3()[0] > 0.);
        assert!(matches!(
            boundary.heal_material(MaterialOptions::default(), &|| false),
            Err(Error::Unsupported)
        ));
        assert!(matches!(
            orthogonal_seats::qualify(&source, pose, Default::default(), &|| false),
            Err(orthogonal_seats::QualificationError::Unsupported)
        ));
    }
    #[test]
    fn pierced_nonconvex_caps_earn_global_grid_and_hole_collar_maps() {
        let source = source([20., 2.5], 0.5, 1.0);
        let pose = StockClipPlacement::new(
            [[0.8, -0.6, 0.], [0.6, 0.8, 0.], [0., 0., 1.]],
            [11., -7., 5.],
        )
        .unwrap();
        let q = orthogonal_seats::qualify(
            &source,
            pose,
            orthogonal_seats::QualificationOptions::default(),
            &|| false,
        )
        .unwrap();
        assert!(q.source_to_healed_position_error_mm() < 1e-7);
        assert!(q.source_to_healed_normalized_derivative_error_mm() < 1e-7);
        let document =
            CompoundDocument::parse(q.as_str(), owned_insert::DocumentOptions::default(), None)
                .unwrap();
        let proof = q
            .verify_actual_member(
                &document,
                document.represented_members()[0].item_entity(),
                orthogonal_seats::VerificationOptions::default(),
                &|| false,
            )
            .unwrap();
        assert!(proof.healed_volume_bounds_mm3()[0] > 0.);
    }
    #[test]
    fn geometrically_valid_hole_without_rectangular_collar_refuses_honestly() {
        // The disk clears the reentrant corner by sqrt(2)-1.2, but no
        // axis-aligned rectangle enclosing it can avoid that corner.
        let source = source([4., 4.], 1.1, 1.2);
        let pose =
            StockClipPlacement::new([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]], [0.; 3]).unwrap();
        let candidate = orthogonal_seats::candidate(
            &source,
            pose,
            orthogonal_seats::Options::default(),
            &|| false,
        )
        .unwrap();
        let document = CompoundDocument::parse(
            candidate.as_str(),
            owned_insert::DocumentOptions::default(),
            None,
        )
        .unwrap();
        let boundary = decode(
            &document,
            document.represented_members()[0].item_entity(),
            Options::default(),
            &|| false,
        )
        .unwrap();
        assert!(source.volume_bounds_mm3()[0] > 0.);
        assert!(matches!(
            boundary.heal_material(MaterialOptions::default(), &|| false),
            Err(Error::Unsupported)
        ));
        assert!(matches!(
            orthogonal_seats::qualify(
                &source,
                pose,
                orthogonal_seats::QualificationOptions::default(),
                &|| false
            ),
            Err(orthogonal_seats::QualificationError::Unsupported)
        ));
    }
    #[test]
    fn ill_conditioned_construction_grid_does_not_publish_a_fold_certificate() {
        let source = source([6.375 + 2f64.powi(-45), 2.5], 0.5, 1.0);
        let pose = StockClipPlacement::new(
            [[0.8, -0.6, 0.], [0.6, 0.8, 0.], [0., 0., 1.]],
            [1_000_011., -7., 5.],
        )
        .unwrap();
        let candidate = orthogonal_seats::candidate(
            &source,
            pose,
            orthogonal_seats::Options::default(),
            &|| false,
        )
        .unwrap();
        let document = CompoundDocument::parse(
            candidate.as_str(),
            owned_insert::DocumentOptions::default(),
            None,
        )
        .unwrap();
        let boundary = decode(
            &document,
            document.represented_members()[0].item_entity(),
            Options::default(),
            &|| false,
        )
        .unwrap();
        // Actual closed source and valid raw schema exist; the finite grid's
        // derivative enclosure cannot earn k<1 for this nearly coincident
        // constructed collar/stock line. No tolerance or gap is widened.
        assert!(source.volume_bounds_mm3()[0] > 0.);
        assert!(matches!(
            orthogonal_seats::qualify(
                &source,
                pose,
                orthogonal_seats::QualificationOptions::default(),
                &|| false
            ),
            Err(orthogonal_seats::QualificationError::Unsupported)
        ));
        assert_eq!(boundary.conical_face_count(), 1);
    }
}

#[test]
fn second_actual_translated_candidate_earns_nonzero_pl_and_collar_correspondence() {
    use crate::{axial_pipe_thread::authored::StockClipPlacement, step::orthogonal_seats};
    let source = super::super::super::tests::source();
    let baseline =
        StockClipPlacement::new([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]], [0.; 3]).unwrap();
    for (distance, accepted) in [(2e-10, true), (2e-9, false), (2e-6, false)] {
        let moved = StockClipPlacement::new(baseline.matrix(), [distance, 0., 0.]).unwrap();
        let candidate = orthogonal_seats::candidate(
            &source,
            moved,
            orthogonal_seats::Options::default(),
            &|| false,
        )
        .unwrap();
        let document = CompoundDocument::parse(
            candidate.as_str(),
            owned_insert::DocumentOptions::default(),
            None,
        )
        .unwrap();
        let boundary = decode(
            &document,
            document.represented_members()[0].item_entity(),
            Options::default(),
            &|| false,
        )
        .unwrap();
        let material = boundary
            .heal_material(MaterialOptions::default(), &|| false)
            .unwrap();
        let matched = material.check_original(
            &source,
            baseline,
            candidate.face_entities(),
            TransferOptions::default(),
            &|| false,
        );
        if accepted {
            let proof = matched.unwrap();
            let c = proof.source_to_raw_correspondence();
            eprintln!(
                "actual second candidate shift={distance} C0={} D1={}",
                c.maximum_position_error_mm(),
                c.maximum_normalized_derivative_error_mm()
            );
            assert!(c.maximum_position_error_mm() >= distance);
            assert!(c.maximum_normalized_derivative_error_mm() > 1e-9);
            assert!(
                proof.source_to_healed_position_error_mm() < 1e-7
                    && proof.source_to_healed_normalized_derivative_error_mm() < 1e-7
            );
        } else {
            assert!(matches!(matched, Err(Error::Unsupported)));
        }
    }
}
