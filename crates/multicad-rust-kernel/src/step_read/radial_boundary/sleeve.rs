//! Cylindrical stock around one through retained radial cavity. The cavity
//! tables are decoded from the stock's inward faces and annulus holes, then
//! certified by the caller exactly like an independent radial boundary. The
//! stock is one periodic cylinder whose two rims bound the annuli.
use super::*;

/// Finite cylindrical stock: bottom cap centre, unit axis, seam reference.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SleeveStock {
    pub origin: Point,
    pub normal: Point,
    pub reference: Point,
    pub radius_mm: f64,
    pub height_mm: f64,
}
/// Decoded candidate only: the cavity ownership has no material certificate.
#[derive(Clone)]
pub struct ImportedPipeSleeve {
    pub cavity: radial::BoundaryOwnership,
    pub stock: SleeveStock,
}

/// Decode one candidate sleeve. This is not an occupied-solid admission.
pub fn read_pipe_sleeve(input: &str, tolerance: f64) -> Result<ImportedPipeSleeve> {
    let mut members = compound::read_members_by(input, tolerance, true, true, |shell, id, _| {
        recognize(shell, id)?.ok_or(ReadError::Unsupported)
    })?;
    members.pop().ok_or(ReadError::Invalid)
}

fn check(condition: bool) -> Result<()> {
    if condition {
        Ok(())
    } else {
        Err(ReadError::Unsupported)
    }
}
fn dot(a: Point, b: Point) -> f64 {
    (0..3).map(|i| a[i] * b[i]).sum()
}

fn scan(shell: &Shell<'_, '_>, faces: &[u32]) -> Option<u32> {
    if faces.len() > 35 {
        return None;
    }
    let mut cylinders = Vec::new();
    let mut annuli = 0;
    for &face_id in faces {
        let (surface, _, bounds) = face(shell, face_id).ok()?;
        match shell.simple(surface).ok()?.name {
            "CYLINDRICAL_SURFACE" => cylinders.push(face_id),
            "PLANE" if bounds.len() == 2 => annuli += 1,
            _ => {}
        }
    }
    (cylinders.len() == 1 && annuli == 2).then(|| cylinders[0])
}
/// `None` unless the member has exactly one cylindrical face and two planar
/// faces with an outer rim and one hole; other members keep their readers.
pub(in crate::step_read) fn recognize(
    shell: &mut Shell<'_, '_>,
    id: u32,
) -> Result<Option<ImportedPipeSleeve>> {
    if shell.simple(id)?.name != "MANIFOLD_SOLID_BREP" {
        return Ok(None);
    }
    let faces = shell.member_faces(id)?;
    // Discovery only: a member another reader owns is never refused here.
    let Some(wall) = scan(shell, &faces) else {
        return Ok(None);
    };
    let (surface, outward, bounds) = face(shell, wall)?;
    let [Value::String(_), placement, Value::Number(radius)] =
        shell.args(surface, "CYLINDRICAL_SURFACE")?
    else {
        return Err(ReadError::Invalid);
    };
    let radius = *radius;
    let stock_axis = axis(shell, reference(placement)?)?;
    check(outward && radius.is_finite() && radius > 0.)?;
    let [bound] = bounds[..] else {
        return Err(ReadError::Unsupported);
    };
    let wall_edges = ring(shell, bound, true)?;
    check(wall_edges.len() == 4)?;
    let mut rims = Vec::with_capacity(2);
    for &(edge, _) in &wall_edges {
        let [Value::String(_), _, _, curve_id, _] = shell.args(edge, "EDGE_CURVE")? else {
            return Err(ReadError::Invalid);
        };
        let curve_id = reference(curve_id)?;
        if shell.simple(curve_id)?.name == "CIRCLE" {
            let radial::Curve::Circle { axis, radius_mm } = curve(shell, curve_id)? else {
                return Err(ReadError::Invalid);
            };
            check(radius_mm == radius)?;
            rims.push((edge, axis));
        }
    }
    check(rims.len() == 2 && rims[0].0 != rims[1].0)?;
    let tolerance = shell.tolerance;
    let normal = stock_axis.normal;
    check((dot(normal, normal) - 1.).abs() <= 1e-12)?;
    let heights: [f64; 2] = std::array::from_fn(|i| {
        let p = rims[i].1.origin;
        dot(std::array::from_fn(|k| p[k] - stock_axis.origin[k]), normal)
    });
    let (bottom, top) = if heights[0] <= heights[1] {
        (0, 1)
    } else {
        (1, 0)
    };
    let height = heights[top] - heights[bottom];
    check(heights[bottom].abs() <= tolerance && height.is_finite() && height > tolerance)?;
    for (rim, h) in rims.iter().zip(heights) {
        // Each rim centre lies on the stock axis and its plane is normal to it.
        let offset: Point =
            std::array::from_fn(|i| rim.1.origin[i] - stock_axis.origin[i] - normal[i] * h);
        check(dot(offset, offset).sqrt() <= tolerance)?;
        check((dot(rim.1.normal, normal) - 1.).abs() <= 1e-12)?;
    }
    let mut used = [false; 2];
    let mut decoder = Decoder::default();
    for face_id in faces {
        if face_id == wall {
            continue;
        }
        let (surface, forward, bounds) = face(shell, face_id)?;
        match bounds[..] {
            [single] => {
                // Inward stock face of the cavity: restore its own orientation.
                let edges = ring(shell, single, true)?;
                let restored = edges
                    .iter()
                    .rev()
                    .map(|&(e, d)| (e, !d))
                    .collect::<Vec<_>>();
                decoder.push(shell, surface, !forward, &restored)?;
            }
            [a, b] => {
                check(shell.simple(surface)?.name == "PLANE")?;
                let (outer_first, _, _) = shell.bound(a, false)?;
                let (outer, hole) = if outer_first { (a, b) } else { (b, a) };
                let outer = ring(shell, outer, true)?;
                let [(rim, _)] = outer[..] else {
                    return Err(ReadError::Unsupported);
                };
                let index = rims
                    .iter()
                    .position(|r| r.0 == rim)
                    .ok_or(ReadError::Unsupported)?;
                check(!used[index])?;
                used[index] = true;
                let edges = ring(shell, hole, false)?;
                let restored = edges
                    .iter()
                    .rev()
                    .map(|&(e, d)| (e, !d))
                    .collect::<Vec<_>>();
                decoder.push(shell, surface, forward, &restored)?;
            }
            _ => return Err(ReadError::Unsupported),
        }
    }
    check(used == [true; 2])?;
    let cavity = radial::BoundaryOwnership::new(decoder.input, tolerance, None).map_err(failure)?;
    Ok(Some(ImportedPipeSleeve {
        cavity,
        stock: SleeveStock {
            origin: stock_axis.origin,
            normal,
            reference: stock_axis.reference,
            radius_mm: radius,
            height_mm: height,
        },
    }))
}
