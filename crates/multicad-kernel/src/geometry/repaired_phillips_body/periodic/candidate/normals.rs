//! Fresh whole smooth-patch normal cones against actual emitted facets.
//! No source tangent or phase is inferred by voting at triangle corners.
use super::*;
fn cross(a: [I; 3], b: [I; 3]) -> [I; 3] {
    std::array::from_fn(|i| a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3])
}
fn dot(a: [I; 3], b: [I; 3]) -> I {
    (0..3).fold(I::ZERO, |s, i| s + a[i] * b[i])
}
fn within(actual: [I; 3], source: [I; 3], angle: f64) -> bool {
    let d = dot(actual, source);
    let squared = cross(actual, source)
        .into_iter()
        .fold(I::ZERO, |s, x| s + x.square());
    // Sum of real squares is nonnegative even when directed addition widens
    // its zero lower endpoint below zero. Only the upper norm is needed.
    let c = I::new(0., squared.hi).sqrt();
    let ratio = c / d;
    d.lo > 0. && ratio.valid() && ratio.hi <= angle
}
fn directions(points: [Point; 3], error: f64, radius: f64) -> Result<[I; 2], GeometryError> {
    let units = points.map(|p| {
        let norm = (I::point(p[0]).square() + I::point(p[1]).square()).sqrt();
        [I::point(p[0]) / norm, I::point(p[1]) / norm]
    });
    let mut lower = 1_f64;
    for a in units {
        for b in units {
            lower = lower.min((a[0] * b[0] + a[1] * b[1]).lo);
        }
    }
    // Every convex combination of these unit rays has norm²>=min pair dot.
    // Positive pair dots establish the complete open half-plane. Normalizing
    // that convex cone encloses every intervening real angular ray.
    if lower <= 0. || error >= radius / 4. {
        return Err(unsupported(
            "Periodic facet has no bounded source angular cover",
        ));
    }
    let norm = I::new(I::point(lower).sqrt().lo, 1.);
    let delta = (I::point(2.) * I::point(error) / I::point(radius)).hi;
    Ok(std::array::from_fn(|i| {
        let interval = I::new(
            units.iter().map(|u| u[i].lo).fold(1., f64::min),
            units.iter().map(|u| u[i].hi).fold(-1., f64::max),
        ) / norm;
        let enlarged = interval + I::new(-delta, delta);
        I::new(enlarged.lo.max(-1.), enlarged.hi.min(1.))
    }))
}
pub(super) fn certify(
    mesh: &DisplayMesh,
    source: &OwnedPeriodicMaterial,
    coefficients: [[[f64; 2]; 2]; 4],
    error: f64,
    angle: f64,
    maximum: usize,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<usize, GeometryError> {
    let mut remaining = maximum;
    let mut charge = |n| {
        checkpoint(interrupt)?;
        remaining = remaining.checked_sub(n).ok_or_else(exhausted)?;
        Ok::<_, GeometryError>(())
    };
    for (tri, &face) in mesh
        .indices
        .as_chunks::<3>()
        .0
        .iter()
        .zip(&mesh.face_indices)
    {
        charge(128)?;
        let points = tri.map(|id| std::array::from_fn(|i| mesh.vertices[3 * id as usize + i]));
        let edges = [1, 2]
            .map(|j| std::array::from_fn(|i| I::point(points[j][i]) - I::point(points[0][i])));
        let actual = cross(edges[0], edges[1]);
        let face = source
            .face_source(face as usize)
            .expect("candidate remaps every facet to the original owned atlas");
        match face {
            PhillipsMaterialFaceSource::ReversedCavity { .. } => {} // traced source plus fresh translation proof in producer
            PhillipsMaterialFaceSource::ThreadTip
            | PhillipsMaterialFaceSource::JoinAnnulus
            | PhillipsMaterialFaceSource::TopAnnulus => {
                let sign = if matches!(face, PhillipsMaterialFaceSource::TopAnnulus) {
                    -1.
                } else {
                    1.
                };
                if !within(actual, [I::ZERO, I::ZERO, I::point(sign)], angle) {
                    return Err(unsupported("Periodic planar facet normal is unresolved"));
                }
            }
            PhillipsMaterialFaceSource::ThreadSide { face: band } => {
                let er = directions(points, error, source.thread().radius_nodes_mm()[1])?;
                let [a, b] = coefficients[band].map(|v| I::new(v[0], v[1]));
                let normal = [er[0] - a * er[1], er[1] + a * er[0], b];
                if !within(actual, normal, angle) {
                    #[cfg(test)]
                    eprintln!(
                        "normal refusal band{band} angle{angle} er{er:?} coeff{a:?}/{b:?} actual{actual:?} expected{normal:?} points{points:?}"
                    );
                    return Err(unsupported(
                        "Periodic whole-band facet normal exceeds requested angle",
                    ));
                }
            }
            PhillipsMaterialFaceSource::HeadMeridian { segment } => {
                let head = source.head().view();
                let radius = head.segment_certificate(segment).unwrap().minimum_radius_mm;
                let er = directions(points, error, radius)?;
                let z = I::new(
                    points.iter().map(|p| p[2]).fold(f64::INFINITY, f64::min),
                    points
                        .iter()
                        .map(|p| p[2])
                        .fold(f64::NEG_INFINITY, f64::max),
                ) + I::new(-error, error);
                let mut pending = vec![(0., 1., 0_u8)];
                while let Some((a, b, depth)) = pending.pop() {
                    charge(257)?;
                    // Existing standalone evaluator includes its immutable constructor
                    // prefix. It is already retained by this exact owned source;
                    // credit that prefix, debit only this fresh257-unit query.
                    let bounds = head
                        .evaluate_with_interrupt(
                            segment,
                            [a, b],
                            head.work_used() + 257,
                            64 * 1024 * 1024,
                            interrupt,
                        )
                        .map_err(|e| match e {
                            multicad_rust_kernel::head_profile::RepairedHeadError::Cancelled => {
                                failure(
                                    GeometryErrorKind::Cancelled,
                                    "Periodic head normal query cancelled",
                                )
                            }
                            multicad_rust_kernel::head_profile::RepairedHeadError::Limit => {
                                exhausted()
                            }
                            _ => unsupported("Periodic head normal source query is unresolved"),
                        })?;
                    let source_z = I::new(bounds.position_mm[1][0], bounds.position_mm[1][1])
                        - I::point(source.head_z_shift_mm());
                    if source_z.hi < z.lo || source_z.lo > z.hi {
                        continue;
                    }
                    let [dr, dz] = bounds.derivative_mm.map(|v| I::new(v[0], v[1]));
                    if within(actual, [dz * er[0], dz * er[1], -dr], angle) {
                        continue;
                    }
                    let mid = a + (b - a) * 0.5;
                    if depth >= 32 || !(a < mid && mid < b) {
                        return Err(unsupported(
                            "Repaired whole-meridian facet normal exceeds requested angle",
                        ));
                    }
                    pending.push((a, mid, depth + 1));
                    pending.push((mid, b, depth + 1));
                }
            }
        }
    }
    Ok(maximum - remaining)
}
