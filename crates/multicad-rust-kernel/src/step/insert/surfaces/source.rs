//! Actual retained-source restriction, sharing every event anchor by ownership.
use super::*;
use crate::authored_internal_treatment::ExactBoundaryRational;
const QUERY_OUTPUT: usize = 64 * 1024;
fn cardinal(point: &mut [f64; 3], r: f64, q: u8, left: bool) {
    let y = if left { -r } else { r };
    let xy = match q {
        0 | 4 => [r, 0.],
        1 => [0., y],
        2 => [-r, 0.],
        _ => [0., -y],
    };
    point[..2].copy_from_slice(&xy);
}
fn quarter(v: ExactBoundaryRational<'_>) -> Option<u8> {
    domains::cardinal_parts(
        v.signum(),
        v.numerator_magnitude_le(),
        v.denominator_magnitude_le(),
    )
}
pub(super) fn original_vertices(
    domain: &DomainPlan<'_>,
    vertices: &mut Vec<SurfaceVertex>,
    w: &mut Work<'_>,
) -> Result<(), Error> {
    let body = domain.boundary();
    for side in [Side::Outer, Side::Inner] {
        let source = body.source_boundary(side);
        for (i, p) in source.points().iter().enumerate() {
            w.charge(1)?;
            vertices.push(SurfaceVertex {
                owner: VertexOwner::Original(vertex_key(side, i as u32)),
                point: p
                    .position_bounds()
                    .map(|b| I::new(b.lower, b.upper).midpoint()),
            });
        }
        for (cell, polygon) in source.cells().iter().enumerate() {
            let exact = domain
                .exact_cell(side, cell)
                .expect("retained original cell query");
            for (corner, id) in polygon.points().iter().enumerate() {
                w.charge(1)?;
                if let Some(q) = quarter(exact.coordinate(corner, 1).unwrap()) {
                    let radius = source.points()[*id as usize].cylindrical_bounds()[0];
                    let r = I::new(radius.lower, radius.upper).midpoint();
                    let vertex = &mut vertices[vertex_index(body, vertex_key(side, *id))];
                    cardinal(&mut vertex.point, r, q, source.source().spec().left_hand);
                }
            }
        }
        for (_, curve) in source.source_curves() {
            w.charge(1)?;
            if let Some(cap) = curve.faces().into_iter().find(|&f| f < 2) {
                let z = if cap == 0 {
                    0.
                } else {
                    body.source().full_length_mm()
                };
                for id in curve.endpoints() {
                    vertices[vertex_index(body, vertex_key(side, id))].point[2] = z;
                }
            }
        }
    }
    Ok(())
}
#[allow(clippy::too_many_arguments)]
pub(super) fn request(
    domain: &DomainPlan<'_>,
    key: RowKey,
    previous: Option<&Request>,
    vertices: &mut Vec<SurfaceVertex>,
    memory: &mut Memory,
    w: &mut Work<'_>,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<Request, Error> {
    let body = domain.boundary();
    let side = key.curve.owner();
    let source = body.source_boundary(side);
    memory.add(QUERY_OUTPUT)?;
    let reserved = memory.used - source.retained_bytes() - SCRATCH;
    let data = source
        .exact_curve_coefficients(key.curve.source_index(), w.remaining(), reserved, interrupt)
        .map_err(domains::failure)?;
    w.charge(data.work_used())?;
    memory.used -= QUERY_OUTPUT;
    memory.add(data.owned_bytes())?;
    let mut base: [N; 3] = std::array::from_fn(|_| N::integer(0));
    let mut delta = base.clone();
    for axis in 0..3 {
        base[axis] = N::view(data.coordinate(0, axis).unwrap(), w)?;
        delta[axis] = N::view(data.coordinate(1, axis).unwrap(), w)?;
    }
    let station = &domain.stations(side)[key.interval as usize..=key.interval as usize + 1];
    let v = [station[0].number(w)?, station[1].number(w)?];
    let t = [
        v[0].sub(&base[2], w)?.div(&delta[2], w)?,
        v[1].sub(&base[2], w)?.div(&delta[2], w)?,
    ];
    let mut start: [N; 3] = std::array::from_fn(|_| N::integer(0));
    let mut end = start.clone();
    for axis in 0..3 {
        start[axis] = base[axis].add(&delta[axis].mul(&t[0], w)?, w)?;
        end[axis] = base[axis].add(&delta[axis].mul(&t[1], w)?, w)?;
    }
    let mut coefficients = [[I::point(0.); 3]; 2];
    for axis in 0..3 {
        coefficients[0][axis] = start[axis].enclosure(w)?;
        coefficients[1][axis] = end[axis].sub(&start[axis], w)?.enclosure(w)?;
    }
    let curve = body.curve(key.curve).expect("actual source row curve");
    let cap = curve.faces().into_iter().find(|&f| f < 2).map(|cap| {
        if cap == 0 {
            0.
        } else {
            body.source().full_length_mm()
        }
    });
    let left = source.source().spec().left_hand;
    let law = Law::from_bounds(coefficients, left);
    let mut ends = [0; 2];
    for i in 0..2 {
        let original = if t[i].cmp(&N::integer(0), w)? == std::cmp::Ordering::Equal {
            Some(curve.endpoints()[0])
        } else if t[i].cmp(&N::integer(1), w)? == std::cmp::Ordering::Equal {
            Some(curve.endpoints()[1])
        } else {
            None
        };
        ends[i] = if let Some(id) = original {
            vertex_index(body, vertex_key(side, id))
        } else if i == 0
            && previous
                .is_some_and(|p| p.key.curve == key.curve && p.key.interval + 1 == key.interval)
        {
            previous.unwrap().ends[1]
        } else {
            let mut point = law.point(i as f64, w)?.map(|v| v.midpoint());
            if let Some(q) = station[i].cardinal() {
                let r = if i == 0 { &start[0] } else { &end[0] }
                    .enclosure(w)?
                    .midpoint();
                cardinal(&mut point, r, q, left);
            }
            if let Some(z) = cap {
                point[2] = z;
            }
            let id = vertices.len();
            w.charge(1)?;
            vertices.push(SurfaceVertex {
                owner: VertexOwner::Station {
                    curve: key.curve,
                    station: key.interval + i as u32,
                },
                point,
            });
            id
        };
    }
    memory.used -= data.owned_bytes();
    Ok(Request {
        key,
        ends,
        law,
        cap,
    })
}
