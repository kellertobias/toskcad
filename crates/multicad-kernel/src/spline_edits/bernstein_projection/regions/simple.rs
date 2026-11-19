use super::*;

fn monotone(piece: &Piece) -> bool {
    (0..2).any(|axis| {
        let mut direction = None;
        for pair in piece.controls[..=piece.degree].windows(2) {
            match (pair[1][axis] / pair[1][2] - pair[0][axis] / pair[0][2]).sign() {
                Some(Equal) => {}
                Some(sign) => {
                    if direction.is_some_and(|d| d != sign) {
                        return false;
                    }
                    direction = Some(sign);
                }
                None => return false,
            }
        }
        direction.is_some()
    })
}
pub(super) fn disjoint(a: &Piece, b: &Piece) -> bool {
    let [a, b] = [a.hull(), b.hull()];
    (0..2).any(|i| a.0[i].hi < b.0[i].lo || b.0[i].hi < a.0[i].lo)
}
pub(super) fn endpoint_separation(a: &Piece, ai: usize, b: &Piece, bi: usize) -> bool {
    let normalized = |p: [Interval; 3]| [p[0] / p[2], p[1] / p[2]];
    let [pa, pb] = [normalized(a.controls[ai]), normalized(b.controls[bi])];
    let endpoint: [Interval; 2] = std::array::from_fn(|i| Interval {
        lo: pa[i].lo.min(pb[i].lo),
        hi: pa[i].hi.max(pb[i].hi),
    });
    let vector = |p: [Interval; 2]| {
        std::array::from_fn::<_, 2, _>(|i| {
            (p[i].lo + p[i].hi) * 0.5 - (endpoint[i].lo + endpoint[i].hi) * 0.5
        })
    };
    let length = |p: [f64; 2]| p[0].hypot(p[1]);
    // A repeated end control is a stationary end derivative; the first
    // distinct control gives the actual leaving direction. The strict
    // control-side test below certifies whichever direction is chosen.
    let adjacent = |piece: &Piece, i: usize| {
        let order: Vec<usize> = if i == 0 {
            (1..=piece.degree).collect()
        } else {
            (0..piece.degree).rev().collect()
        };
        order
            .into_iter()
            .map(|k| vector(normalized(piece.controls[k])))
            .find(|v| length(*v) > 0.)
            .unwrap_or([0.; 2])
    };
    let [x, y] = [adjacent(a, ai), adjacent(b, bi)];
    let [lx, ly] = [length(x), length(y)];
    if lx == 0. || ly == 0. {
        return false;
    }
    let direction: [f64; 2] = std::array::from_fn(|i| x[i] / lx - y[i] / ly);
    let side = |piece: &Piece, endpoint_index, positive| {
        let mut strict = false;
        for (i, p) in piece.controls[..=piece.degree].iter().enumerate() {
            if i == endpoint_index {
                continue;
            }
            let p = normalized(*p);
            let value = (p[0] - endpoint[0]) * Interval::point(direction[0])
                + (p[1] - endpoint[1]) * Interval::point(direction[1]);
            match value.sign() {
                Some(Equal) if i != 0 && i != piece.degree => {}
                Some(Greater) if positive => strict = true,
                Some(Less) if !positive => strict = true,
                _ => return false,
            }
        }
        strict
    };
    side(a, ai, true) && side(b, bi, false)
}

/// Both pieces lie exactly on one supporting line and share a finite run
/// longer than the modeling tolerance. Exact zero cross products certify
/// collinearity; a rounded or merely close pair is never reported.
pub(super) fn collinear_overlap(a: &Piece, b: &Piece, tolerance: f64) -> bool {
    let normalized = |p: [Interval; 3]| [p[0] / p[2], p[1] / p[2]];
    let start = normalized(a.controls[0]);
    let end = normalized(a.controls[a.degree]);
    let d = [end[0] - start[0], end[1] - start[1]];
    let middle = |v: Interval| (v.lo + v.hi) * 0.5;
    let length = middle(d[0]).hypot(middle(d[1]));
    if !length.is_finite() || length <= tolerance {
        return false;
    }
    let mut runs = [[0.; 2]; 2];
    for (run, piece) in runs.iter_mut().zip([a, b]) {
        for control in &piece.controls[..=piece.degree] {
            let q = normalized(*control);
            let r = [q[0] - start[0], q[1] - start[1]];
            if (r[0] * d[1] - r[1] * d[0]).sign() != Some(Equal) {
                return false;
            }
        }
        // Collinear Bernstein controls keep the curve on the line; its end
        // controls are actual curve points, so their span is a covered run.
        let along = |control: [Interval; 3]| {
            let q = normalized(control);
            middle((q[0] - start[0]) * d[0] + (q[1] - start[1]) * d[1]) / length
        };
        let [x, y] = [
            along(piece.controls[0]),
            along(piece.controls[piece.degree]),
        ];
        *run = [x.min(y), x.max(y)];
    }
    runs[0][1].min(runs[1][1]) - runs[0][0].max(runs[1][0]) > tolerance * 2.
}

/// Transverse self contacts and whether an exact collinear retrace exists.
pub(crate) fn self_contacts(
    spline: &BSpline2d,
    policy: GeometryPolicy,
    interrupt: QueryCheckpoint<'_>,
) -> Result<(Vec<[f64; 2]>, bool), GeometryError> {
    let mut budget = Budget {
        remaining: WORK,
        checkpoint: interrupt,
    };
    let mut stack = pieces(spline, spline.poles[0], policy, interrupt)?;
    stack.reverse();
    let mut leaves = Vec::new();
    while let Some(piece) = stack.pop() {
        budget.step()?;
        if monotone(&piece) {
            leaves.push(piece);
        } else {
            let [a, b] = piece.split()?;
            stack.push(b);
            stack.push(a);
        }
        if leaves.len() + stack.len() > 512 {
            return Err(failure(
                GeometryErrorKind::Limit,
                "Spline regularity proof has too many spans",
            ));
        }
    }
    let mut contacts: Vec<[f64; 2]> = Vec::new();
    let domain = [
        leaves[0].range[0],
        leaves.last().expect("source spans").range[1],
    ];
    let closed = spline.poles.first() == spline.poles.last();
    for i in 0..leaves.len() {
        for j in i + 1..leaves.len() {
            let mut pairs = vec![(leaves[i].clone(), leaves[j].clone())];
            while let Some((a, b)) = pairs.pop() {
                budget.step()?;
                if disjoint(&a, &b) {
                    continue;
                }
                let adjacent = a.range[1] == b.range[0] && endpoint_separation(&a, a.degree, &b, 0);
                let seam = closed
                    && a.range[0] == domain[0]
                    && b.range[1] == domain[1]
                    && endpoint_separation(&a, 0, &b, b.degree);
                if adjacent || seam {
                    continue;
                }
                if collinear_overlap(&a, &b, policy.linear_tolerance_mm) {
                    return Ok((contacts, true));
                }
                if let Some(point) = super::contact::transverse(
                    a.clone(),
                    b.clone(),
                    spline.poles[0],
                    policy,
                    &mut budget,
                )? {
                    if contacts.iter().any(|previous| {
                        (point[0] - previous[0]).hypot(point[1] - previous[1])
                            <= policy.linear_tolerance_mm * 2.
                    }) {
                        return Err(uncertain());
                    }
                    contacts.push(point);
                    continue;
                }
                if a.depth <= b.depth {
                    for child in a.split()? {
                        pairs.push((child, b.clone()));
                    }
                } else {
                    for child in b.split()? {
                        pairs.push((a.clone(), child));
                    }
                }
            }
        }
    }
    Ok((contacts, false))
}
