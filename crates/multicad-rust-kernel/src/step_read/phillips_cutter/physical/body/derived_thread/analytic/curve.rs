//! Whole-coefficient analytic factor and tensor transfer bounds.
use super::*;
mod rounded;
fn choose(n: usize, k: usize) -> f64 {
    (0..k).fold(1., |a, i| a * (n - i) as f64 / (i + 1) as f64)
}
fn mul(a: [I; 2], b: [I; 2]) -> [I; 2] {
    [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]]
}
fn exp_factor_bounds(
    trig: &rotation::Rotation,
    angle: [I; 2],
    radius: [I; 2],
    work: &mut Work<'_>,
) -> Result<(Vec<[I; 2]>, [f64; 2])> {
    work.step(32768)?;
    let midpoint = angle[0] / I::point(2.) + angle[1] / I::point(2.);
    let delta = angle[1] - angle[0];
    let rm = radius[0] / I::point(2.) + radius[1] / I::point(2.);
    let dr = radius[1] - radius[0];
    let rotation = trig.sin_cos(midpoint, work)?;
    let mut coefficients = vec![[I::point(0.); 2]; 17];
    let mut power = [I::point(1.), I::point(0.)];
    for order in 0..=15 {
        for axis in 0..2 {
            coefficients[order][axis] = coefficients[order][axis] + rm * power[axis];
            coefficients[order + 1][axis] = coefficients[order + 1][axis] + dr * power[axis];
        }
        power = mul(power, [I::point(0.), delta / I::point((order + 1) as f64)]);
    }
    let mut ordinary = vec![[I::point(0.); 2]; 17];
    for (k, c) in coefficients.into_iter().enumerate() {
        for (i, a) in ordinary.iter_mut().enumerate().take(k + 1) {
            let weight = I::point(choose(k, i)) * I::point((-0.5_f64).powi((k - i) as i32));
            for axis in 0..2 {
                a[axis] = a[axis] + c[axis] * weight;
            }
        }
    }
    let poles = (0..=16)
        .map(|i| {
            let mut p = [I::point(0.); 2];
            for (k, c) in ordinary.iter().enumerate().take(i + 1) {
                let weight = I::point(choose(i, k)) / I::point(choose(16, k));
                for axis in 0..2 {
                    p[axis] = p[axis] + c[axis] * weight;
                }
            }
            mul(rotation, p)
        })
        .collect();
    let half = I::point(delta.magnitude()) / I::point(2.);
    let mut power = I::point(1.);
    let mut previous = power;
    let mut factorial = I::point(1.);
    for k in 1..=16 {
        previous = power;
        power = power * half;
        factorial = factorial * I::point(k as f64);
    }
    let remainder = power / factorial;
    let derivative = I::point(dr.magnitude()) * remainder
        + I::point(radius[0].magnitude().max(radius[1].magnitude()))
            * I::point(delta.magnitude())
            * previous
            / (factorial / I::point(16.));
    Ok((
        poles,
        [
            (I::point(radius[0].magnitude().max(radius[1].magnitude())) * remainder).hi,
            derivative.hi,
        ],
    ))
}
fn factor_transfer(actual: &[[f64; 2]], candidate: &[[I; 2]]) -> [f64; 2] {
    let mut difference = Vec::with_capacity(17);
    for (i, &proposed) in candidate.iter().enumerate() {
        let mut elevated = [I::point(0.); 2];
        for (j, &p) in actual.iter().enumerate() {
            if j <= i && i - j <= 5 {
                let weight = I::point(choose(11, j) * choose(5, i - j)) / I::point(choose(16, i));
                for axis in 0..2 {
                    elevated[axis] = elevated[axis] + I::point(p[axis]) * weight;
                }
            }
        }
        difference.push(std::array::from_fn(|a| elevated[a] - proposed[a]));
    }
    fn refined_bound(p: Vec<[I; 2]>, depth: usize) -> f64 {
        if depth == 0 {
            return p
                .into_iter()
                .map(|a| (a[0].square() + a[1].square()).sqrt().hi)
                .fold(0_f64, f64::max);
        }
        let mut q = p;
        let n = q.len();
        let mut a = vec![[I::point(0.); 2]; n];
        let mut b = a.clone();
        for k in 0..n {
            a[k] = q[0];
            b[n - 1 - k] = q[n - 1 - k];
            for i in 0..n - 1 - k {
                q[i] = std::array::from_fn(|j| (q[i][j] + q[i + 1][j]) * I::point(0.5));
            }
        }
        refined_bound(a, depth - 1).max(refined_bound(b, depth - 1))
    }
    let derivative = difference
        .windows(2)
        .map(|p| std::array::from_fn(|a| I::point(16.) * (p[1][a] - p[0][a])))
        .collect();
    [refined_bound(difference, 3), refined_bound(derivative, 3)]
}

pub(super) struct Factors {
    a: Vec<[f64; 2]>,
    au: Vec<f64>,
    b: [f64; 2],
    bv: f64,
    z: [f64; 2],
    major: f64,
}
pub(super) fn factors(
    thread: &crate::helical_thread::Thread,
    reference: &crate::step::helical::ThreadSplines,
    rotation: &rotation::Rotation,
    periodic: Option<&PeriodicThreadSource>,
    work: &mut Work<'_>,
) -> Result<Factors> {
    let rate = if periodic.is_some() {
        I::point(thread.spec().pitch_mm) / (I::point(2.) * rotation.pi(work)?)
    } else {
        I::point(thread.helical_rate().abs())
    };
    let hand = I::point(thread.helical_rate().signum());
    let [_, start, upper, lower, crest, _] = thread.profile_breaks();
    let phase = periodic.map_or_else(
        || [start, upper, lower, crest, start + thread.spec().pitch_mm].map(I::point),
        |p| p.phases(),
    );
    let major = I::point(thread.spec().major_radius_mm);
    let root = I::point(thread.root_radius_mm());
    let mut radius_bound = thread.spec().major_radius_mm;
    let mut all_a = Vec::with_capacity(reference.bands.len());
    let mut all_au = Vec::with_capacity(reference.bands.len());
    let mut b = [0_f64; 2];

    let mut bv = 0_f64;
    let mut z = [0_f64; 2];
    for (band, spline) in reference.bands.iter().enumerate() {
        let mut a = [0_f64; 2];
        let mut au = 0_f64;
        let last = spline.knots.len() - 1;
        for cell in 0..last {
            let endpoints = [cell, cell + 1].map(|i| {
                let ratio = (I::point(spline.knots[i]) - I::point(spline.knots[0]))
                    / (I::point(spline.knots[last]) - I::point(spline.knots[0]));
                phase[band] + ratio * (phase[band + 1] - phase[band])
            });
            let radii = if let Some(source) = periodic {
                [cell, cell + 1].map(|i| {
                    let ratio = (I::point(spline.knots[i]) - I::point(spline.knots[0]))
                        / (I::point(spline.knots[last]) - I::point(spline.knots[0]));
                    source.radii(band, ratio)
                })
            } else {
                endpoints.map(|p| match band {
                    0 => major + I::point(thread.slope_at_phase(0., 1)) * (p - I::point(start)),
                    1 => root,
                    2 => major + I::point(thread.slope_at_phase(0., 3)) * (p - I::point(crest)),
                    _ => major,
                })
            };
            let angle = endpoints.map(|p| -hand * p / rate);
            let (candidate, remainder, cell_radius, cell_du) =
                if let Some(root) = thread.root_circle().filter(|_| band == 1) {
                    let angles = [
                        I::point(spline.knots[cell]),
                        I::point(spline.knots[cell + 1]),
                    ];
                    let factor = rounded::bounds(root, angles, rate, hand, work)?;
                    (factor.poles, factor.remainder, factor.radius, factor.du)
                } else {
                    let cell_radius = radii[0].magnitude().max(radii[1].magnitude());
                    let (candidate, remainder) = exp_factor_bounds(rotation, angle, radii, work)?;
                    let du = (I::point((radii[1] - radii[0]).magnitude())
                        + I::point(cell_radius) * I::point((angle[1] - angle[0]).magnitude()))
                    .hi;
                    (candidate, remainder, cell_radius, du)
                };
            radius_bound = radius_bound.max(cell_radius);
            let transfer = factor_transfer(&spline.poles[11 * cell..=11 * cell + 11], &candidate);
            for k in 0..2 {
                a[k] = a[k].max((I::point(transfer[k]) + I::point(remainder[k])).hi);
            }
            au = au.max(cell_du);
        }
        all_a.push(a);
        all_au.push(au);
    }
    for cell in 0..reference.axial.knots.len() - 1 {
        let heights = [
            reference.heights[11 * cell],
            reference.heights[11 * cell + 11],
        ];
        let angle = heights.map(|h| hand * I::point(h) / rate);
        let (candidate, remainder) = exp_factor_bounds(rotation, angle, [I::point(1.); 2], work)?;
        let transfer = factor_transfer(
            &reference.axial.poles[11 * cell..=11 * cell + 11],
            &candidate,
        );
        for k in 0..2 {
            b[k] = b[k].max((I::point(transfer[k]) + I::point(remainder[k])).hi);
        }
        bv = bv.max((angle[1] - angle[0]).magnitude());
        let difference = std::array::from_fn::<_, 12, _>(|i| {
            I::point(reference.heights[11 * cell + i])
                - (I::point(heights[0])
                    + I::point(i as f64) / I::point(11.)
                        * (I::point(heights[1]) - I::point(heights[0])))
        });
        z[0] = z[0].max(
            difference
                .iter()
                .map(|d| d.magnitude())
                .fold(0_f64, f64::max),
        );
        z[1] = z[1].max(
            difference
                .windows(2)
                .map(|p| (I::point(11.) * (p[1] - p[0])).magnitude())
                .fold(0_f64, f64::max),
        );
    }

    Ok(Factors {
        a: all_a,
        au: all_au,
        b,
        bv,
        z,
        major: radius_bound,
    })
}

pub(super) fn literal_transfer(
    r: &Raw,
    reference: &crate::step::helical::ThreadSplines,
    index: usize,
    frame: &Target,
    shell: &Shell<'_, '_>,
    work: &mut Work<'_>,
) -> Result<[f64; 3]> {
    let [origin, u, v, n] = frame.columns;
    let mut distance = 0_f64;
    let mut derivative = [0_f64; 2];
    // Inferred and supplied paths both retain the actual raw cell map. A
    // mismatched external reference refuses before indexing or allocating any
    // coefficient-difference buffer; no metadata replacement is performed.
    let radial = reference.bands.get(index).ok_or(ReadError::Unsupported)?;
    {
        let (width, height) = (r.ids.len(), r.ids.first().map_or(0, |row| row.len()));
        require_matching_grid(
            width,
            height,
            r.ids.iter().all(|row| row.len() == height),
            radial.poles.len(),
            reference.axial.poles.len(),
        )?;
        let mut differences = Vec::with_capacity(width * height);
        for (i, row) in r.ids.iter().enumerate() {
            for (j, &id) in row.iter().enumerate() {
                work.step(200)?;
                let ar = radial.poles[i].map(I::point);
                let az = reference.axial.poles[j].map(I::point);
                let xy = [ar[0] * az[0] - ar[1] * az[1], ar[0] * az[1] + ar[1] * az[0]];
                let point = std::array::from_fn::<_, 3, _>(|k| {
                    I::point(origin[k])
                        + I::point(u[k]) * xy[0]
                        + I::point(v[k]) * xy[1]
                        + I::point(n[k])
                            * (I::point(reference.heights[j]) - I::point(frame.z_shift))
                });
                let actual = shell.point(id)?.map(|x| I::new(x.next_down(), x.next_up()));
                let d = std::array::from_fn(|k| actual[k] - point[k]);
                distance = distance.max(reduce::l1(d));
                differences.push(d);
            }
        }
        for i in 0..width {
            for j in 0..height {
                let a = differences[i * height + j];
                if i + 1 < width {
                    derivative[0] = derivative[0].max(reduce::l1(std::array::from_fn(|k| {
                        I::point(11.) * (differences[(i + 1) * height + j][k] - a[k])
                    })));
                }
                if j + 1 < height {
                    derivative[1] = derivative[1].max(reduce::l1(std::array::from_fn(|k| {
                        I::point(11.) * (differences[i * height + j + 1][k] - a[k])
                    })));
                }
            }
        }
    }

    Ok([distance, derivative[0], derivative[1]])
}

pub(super) fn world_transfer(
    index: usize,
    frame: &Target,
    literal: [f64; 3],
    g: &Factors,
    work: &mut Work<'_>,
) -> Result<[f64; 3]> {
    work.step(256)?;
    let a = g.a[index];
    let au = g.au[index];
    let b = g.b;
    let bv = g.bv;
    let [_, u, v, n] = frame.columns;
    let dot =
        |a: Point, b: Point| (0..3).fold(I::point(0.), |s, k| s + I::point(a[k]) * I::point(b[k]));
    // Gershgorin bounds the actual XY frame's largest Gram eigenvalue;
    // no approximate rigid-transform assumption is used in this transfer.
    let scale = (I::point(dot(u, u).hi.max(dot(v, v).hi)) + I::point(dot(u, v).magnitude())).sqrt();
    let nz = dot(n, n).sqrt();
    let major = I::point(g.major);
    let error = [
        I::point(a[0]) + (major + I::point(a[0])) * I::point(b[0]),
        I::point(a[1]) + (I::point(au) + I::point(a[1])) * I::point(b[0]),
        I::point(a[0]) * I::point(bv) + (major + I::point(a[0])) * I::point(b[1]),
    ];
    let bounds = std::array::from_fn(|k| {
        (I::point(literal[k]) + scale * error[k] + nz * I::point(g.z[usize::from(k == 2)])).hi
    });
    if bounds.into_iter().all(|x| x.is_finite() && x >= 0.) {
        Ok(bounds)
    } else {
        Err(ReadError::Unsupported)
    }
}

// Whole ORIGINAL raw entry-rim coefficients, with the already-earned actual
// A-angle factors. B is fixed at its ACTUAL retained entry coefficient here;
// interior axial-cell errors are deliberately not substituted for that rim.
pub(super) struct EntryLaw<'l> {
    pub rotation: &'l rotation::Rotation,
    pub factors: &'l Factors,
    pub thread: &'l crate::helical_thread::Thread,
    pub source_z: f64,
}
pub(super) fn entry_transfer(
    raw: &Raw,
    reference: &crate::step::helical::ThreadSplines,
    index: usize,
    frame: &Target,
    shell: &Shell<'_, '_>,
    law: EntryLaw<'_>,
    work: &mut Work<'_>,
) -> Result<[f64; 2]> {
    let EntryLaw {
        rotation,
        factors,
        thread,
        source_z,
    } = law;
    work.step(16384 + raw.ids.len().saturating_mul(512))?;
    let radial = reference.bands.get(index).ok_or(ReadError::Unsupported)?;
    let axial = reference
        .axial
        .poles
        .get(raw.entry)
        .ok_or(ReadError::Unsupported)?;
    let z = *reference
        .heights
        .get(raw.entry)
        .ok_or(ReadError::Unsupported)?;
    require(raw.ids.len() == radial.poles.len())?;
    let [origin, u, v, n] = frame.columns;
    let mut differences: Vec<[I; 3]> = Vec::with_capacity(raw.ids.len());
    let mut literal = 0_f64;
    let mut literal_du = 0_f64;
    for (row, pole) in raw.ids.iter().zip(&radial.poles) {
        let id = *row.get(raw.entry).ok_or(ReadError::Unsupported)?;
        let a = pole.map(I::point);
        let b = axial.map(I::point);
        let xy = mul(a, b);
        let expected = std::array::from_fn::<_, 3, _>(|k| {
            I::point(origin[k])
                + I::point(u[k]) * xy[0]
                + I::point(v[k]) * xy[1]
                + I::point(n[k]) * (I::point(z) - I::point(frame.z_shift))
        });
        let actual = shell.point(id)?.map(|x| I::new(x.next_down(), x.next_up()));
        let difference = std::array::from_fn(|k| actual[k] - expected[k]);
        literal = literal.max(reduce::l1(difference));
        if let Some(previous) = differences.last() {
            literal_du = literal_du.max(reduce::l1(std::array::from_fn(|k| {
                I::point(11.) * (difference[k] - previous[k])
            })));
        }
        differences.push(difference);
    }
    let angle = I::point(thread.helical_rate().signum()) * I::point(source_z)
        / I::point(thread.helical_rate().abs());
    let exact_b = rotation.sin_cos(angle, work)?;
    let b_error = ((I::point(axial[0]) - exact_b[0]).square()
        + (I::point(axial[1]) - exact_b[1]).square())
    .sqrt();
    let a = *factors.a.get(index).ok_or(ReadError::Unsupported)?;
    let au = *factors.au.get(index).ok_or(ReadError::Unsupported)?;
    let dot =
        |p: Point, q: Point| (0..3).fold(I::point(0.), |s, k| s + I::point(p[k]) * I::point(q[k]));
    let scale = (I::point(dot(u, u).hi.max(dot(v, v).hi)) + I::point(dot(u, v).magnitude())).sqrt();
    let normal = dot(n, n).sqrt();
    let error = [
        I::point(a[0]) + (I::point(factors.major) + I::point(a[0])) * b_error,
        I::point(a[1]) + (I::point(au) + I::point(a[1])) * b_error,
    ];
    let bounds = [
        (I::point(literal)
            + scale * error[0]
            + normal * I::point((I::point(z) - I::point(source_z)).magnitude()))
        .hi,
        (I::point(literal_du) + scale * error[1]).hi,
    ];
    require(bounds.into_iter().all(|x| x.is_finite() && x >= 0.))?;
    Ok(bounds)
}

#[cfg(test)]
mod tests;

fn require_matching_grid(
    width: usize,
    height: usize,
    coherent: bool,
    expected_width: usize,
    expected_height: usize,
) -> Result<()> {
    if width > 0 && height > 0 && coherent && width == expected_width && height == expected_height {
        Ok(())
    } else {
        Err(ReadError::Unsupported)
    }
}
