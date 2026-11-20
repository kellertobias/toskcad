//! Strict occupied separation, independent of clipping's contact classification.
use super::*;

pub(super) fn empty(
    operands: &[Operand<'_>],
    tolerance: f64,
    work: &mut Work<'_>,
) -> Result<bool, Error> {
    for other in &operands[1..] {
        let mut separated = true;
        'pairs: for a in operands[0].pieces {
            for b in other.pieces {
                if !pair(a, b, tolerance, work)? {
                    separated = false;
                    break 'pairs;
                }
            }
        }
        if separated {
            return Ok(true);
        }
    }
    Ok(false)
}

/// Every cell pair of operand 0 and some other operand is either strictly
/// separated or meets along a supporting plane of one cell: along that unit
/// normal the closed hulls overlap by no more than the rounding scale of
/// their own world coordinates (and never more than `tolerance / 8`). Such a
/// depth is below what the floating input can represent as material, so the
/// regularized intersection is a contact within that declared uncertainty.
/// An exactly representable sliver, however thin, is never a contact.
pub(super) fn contact(
    operands: &[Operand<'_>],
    tolerance: f64,
    work: &mut Work<'_>,
) -> Result<bool, Error> {
    for other in &operands[1..] {
        let mut meets = true;
        'pairs: for a in operands[0].pieces {
            for b in other.pieces {
                if !pair(a, b, tolerance, work)? && !touching(a, b, tolerance, work)? {
                    meets = false;
                    break 'pairs;
                }
            }
        }
        if meets {
            return Ok(true);
        }
    }
    Ok(false)
}

/// Input rounding scale of two cells: the slack `classify_box` grants to
/// plane evaluation at the same world magnitude.
pub(super) fn contact_allowance(a: &ConvexPiece, b: &ConvexPiece, tolerance: f64) -> f64 {
    let magnitude = a
        .bounds
        .iter()
        .chain(&b.bounds)
        .flatten()
        .map(|v| v.abs())
        .fold(1_f64, f64::max);
    (64. * f64::EPSILON * magnitude).min(tolerance / 8.)
}

pub(super) fn touching(
    a: &ConvexPiece,
    b: &ConvexPiece,
    tolerance: f64,
    work: &mut Work<'_>,
) -> Result<bool, Error> {
    let allowance = contact_allowance(a, b, tolerance);
    for plane in a.planes.iter().chain(&b.planes) {
        let axis = plane.normal;
        // The overlap depth scales with |axis|; bound it from below so a
        // rounded unit normal can only make the certificate stricter.
        let length = dot(axis, axis).sqrt();
        if !(length.is_finite() && length > 0.5) {
            continue;
        }
        let left = project(a, axis, work)?;
        let right = project(b, axis, work)?;
        let depth = (left[1] - right[0])
            .next_up()
            .min((right[1] - left[0]).next_up());
        if depth <= allowance * length.next_down() {
            return Ok(true);
        }
    }
    Ok(false)
}

pub(super) fn pair(
    a: &ConvexPiece,
    b: &ConvexPiece,
    tolerance: f64,
    work: &mut Work<'_>,
) -> Result<bool, Error> {
    for plane in a.planes.iter().chain(&b.planes) {
        let axis = plane.normal;
        // Any finite direction can separate two convex vertex hulls. Using
        // published normals only chooses candidate directions; their accuracy
        // never becomes an assumption about a separating support plane.
        let left = project(a, axis, work)?;
        let right = project(b, axis, work)?;
        // Expand BOTH occupied intervals by the existing geometric tolerance.
        // L1(axis) bounds its Euclidean norm without a rounded square root.
        let length = axis.iter().fold(0., |sum, x| (sum + x.abs()).next_up());
        let pad = (tolerance * length).next_up();
        if (left[1] + pad).next_up() < (right[0] - pad).next_down()
            || (right[1] + pad).next_up() < (left[0] - pad).next_down()
        {
            return Ok(true);
        }
    }
    Ok(false)
}

fn project(piece: &ConvexPiece, axis: Point, work: &mut Work<'_>) -> Result<[f64; 2], Error> {
    work.charge(piece.faces.iter().map(Vec::len).sum())?;
    let mut range = [f64::INFINITY, f64::NEG_INFINITY];
    for point in piece.faces.iter().flatten() {
        // Every multiplication and addition is rounded outwards, including
        // subnormal underflow and cancellation. World coordinates are retained;
        // no approximate inverse authors a different local position.
        let mut interval = [0.; 2];
        for i in 0..3 {
            let product = point[i] * axis[i];
            interval[0] = (interval[0] + product.next_down()).next_down();
            interval[1] = (interval[1] + product.next_up()).next_up();
        }
        range[0] = range[0].min(interval[0]);
        range[1] = range[1].max(interval[1]);
    }
    Ok(range)
}

#[cfg(test)]
mod tests {
    use super::*;
    const TOL: f64 = 1e-7;
    fn block(min: Point, max: Point) -> ConvexPiece {
        ConvexPiece::from_box(min, max, TOL).unwrap()
    }
    #[test]
    fn strict_intervals_prove_clearance_but_never_contact_or_tiny_overlap() {
        let a = block([0.; 3], [1.; 3]);
        for (start, separated) in [
            (2., true),
            (-2., true),
            (1., false),
            (1. - 1e-11, false),
            (1. + 1e-11, false),
        ] {
            let b = block([start, 0., 0.], [start + 1., 1., 1.]);
            assert_eq!(pair(&a, &b, TOL, &mut Work::new(None)).unwrap(), separated);
        }
        let before = a.faces.clone();
        let b = block([2.; 3], [3.; 3]);
        let mut limited = Work::new(None);
        limited.remaining = 0;
        assert_eq!(pair(&a, &b, TOL, &mut limited), Err(Error::Limit));
        assert_eq!(
            pair(&a, &b, TOL, &mut Work::new(Some(&|| true))),
            Err(Error::Cancelled)
        );
        assert!(pair(&a, &b, TOL, &mut Work::new(None)).unwrap());
        assert_eq!(a.faces, before);
    }
    #[test]
    fn directed_projection_encloses_cancellation_underflow_and_large_coordinates() {
        let a = block([0.; 3], [1.; 3]);
        for axis in [
            [f64::from_bits(1), -f64::from_bits(1), 0.],
            [1., -1., 1.],
            [0., 0., 0.],
        ] {
            let [low, high] = project(&a, axis, &mut Work::new(None)).unwrap();
            assert!(low <= -axis[1].abs());
            assert!(high >= axis[0].abs());
        }
        let a = block([1e7; 3], [1e7 + 1.; 3]);
        let b = block([1e7 + 3., 1e7, 1e7], [1e7 + 4., 1e7 + 1., 1e7 + 1.]);
        assert!(pair(&a, &b, TOL, &mut Work::new(None)).unwrap());
    }
}
