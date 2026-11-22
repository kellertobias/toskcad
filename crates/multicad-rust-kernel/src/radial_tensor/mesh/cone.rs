//! Exact cone sectors parameterized by their actual source rails.
//! Rectification to common height rows has a separately charged physical
//! distance bound; it does not replace the retained rail parameterizations.
use super::*;
use spline::Spline;
pub(super) struct Cone {
    pub face: usize,
    left: Spline,
    right: Spline,
    fractions: [I; 2],
    pub closure_error: f64,
    normal_reserve: f64,
}
pub(super) struct Jet {
    pub position: [I; 3],
    pub du: [I; 3],
    pub dv: [I; 3],
    pub normal: [I; 3],
    pub radial: [I; 2],
    pub angular_rate: [I; 2],
    pub hessian_norms: [I; 3],
}
#[derive(Clone, Copy)]
pub(super) struct Row {
    parameter: I,
    point: [I; 3],
    norm: I,
    gap: I,
    speed: [I; 2],
    acceleration: [I; 2],
}
pub(super) struct Rows {
    range: [f64; 2],
    whole: Row,
    middle: Row,
    ends: [Row; 2],
    jumps: I,
}
impl Cone {
    pub fn new(
        body: &RadialBoundary,
        side: usize,
        rails: &[Spline],
        work: &mut Work<'_>,
    ) -> Result<Self, Error> {
        let caps = &body.ordered.caps;
        let chart = &caps.axial.chart;
        let left = rails[side].clone();
        let right = rails[(side + 1) % rails.len()].clone();
        let mut fractions = [I::point(0.); 2];
        for (end, fraction) in fractions.iter_mut().enumerate() {
            work.step()?;
            let v = end as f64;
            let p = left.bounds([v, v], work)?[0];
            let q = right.bounds([v, v], work)?[0];
            let middle = caps::ends(body.input(), caps.rims[end][side][0]).1;
            let m = chart.project(body.input().vertices[middle]);
            let gap = ordering::angle([p[0], p[1]], [q[0], q[1]])?.ok_or(Error::Unresolved)?;
            let first = ordering::angle([p[0], p[1]], [m[0], m[1]])?.ok_or(Error::Unresolved)?;
            *fraction = (first / gap).positive()?;
            (I::point(1.) - *fraction).positive()?;
        }
        // Every original rail has strictly monotone physical Z. At a
        // shared auxiliary row v, its true height differs from H*v by E.
        // For any angle between the two row rays, continuity and the
        // constructor's same-height strict ordering give a point in the
        // actual sector at some height between those two rail heights.
        // Moving along the cone therefore costs at most sqrt(1+k^2)*E.
        // Normal direction depends on angle, not on this height shift.
        // The remaining finite support/endpoint sewing collar is charged
        // separately; no approximate boundary is declared exact here.
        let axial = I::point(
            left.axial_error(chart.nominal_height, work)?
                .max(right.axial_error(chart.nominal_height, work)?),
        );
        let metric = (I::point(1.) + chart.slope.square()).sqrt();
        let sewing = I::point(4.)
            * I::point(body.error_mm)
            * (I::point(1.) + I::point(chart.slope.magnitude()));
        let closure_error = (metric * axial + sewing).checked()?.hi;
        let normal_reserve =
            angular_reserve(I::point(chart.radius_range.lo), chart.slope, body.error_mm)?;
        Ok(Self {
            face: caps.sides[side],
            left,
            right,
            fractions,
            closure_error,
            normal_reserve,
        })
    }
    pub fn jet(
        &self,
        body: &RadialBoundary,
        uv: [I; 2],
        second: bool,
        work: &mut Work<'_>,
    ) -> Result<Jet, Error> {
        let row = self.row(uv[1], work)?;
        self.jet_at_row(body, uv[0], second, &row, work)
    }
    pub fn row(&self, parameter: I, work: &mut Work<'_>) -> Result<Row, Error> {
        work.charge(128)?;
        let range = [parameter.lo.max(0.), parameter.hi.min(1.)];
        let [p, dp, ddp] = self.left.derivatives(range, work)?;
        let [q, dq, ddq] = self.right.derivatives(range, work)?;
        let gap = ordering::angle([p[0], p[1]], [q[0], q[1]])?.ok_or(Error::Unresolved)?;
        let p2 = (p[0].square() + p[1].square()).positive()?;
        let q2 = (q[0].square() + q[1].square()).positive()?;
        let speed_p = (p[0] * dp[1] - p[1] * dp[0]) / p2;
        let speed_q = (q[0] * dq[1] - q[1] * dq[0]) / q2;
        let accel_p = (p[0] * ddp[1] - p[1] * ddp[0]) / p2
            - I::point(2.) * speed_p * (p[0] * dp[0] + p[1] * dp[1]) / p2;
        let accel_q = (q[0] * ddq[1] - q[1] * ddq[0]) / q2
            - I::point(2.) * speed_q * (q[0] * dq[0] + q[1] * dq[1]) / q2;
        Ok(Row {
            parameter,
            point: p,
            norm: p2.sqrt().positive()?,
            gap,
            speed: [speed_p, speed_q],
            acceleration: [accel_p, accel_q],
        })
    }
    pub fn rows(&self, range: [f64; 2], work: &mut Work<'_>) -> Result<Rows, Error> {
        let middle = range[0] + (range[1] - range[0]) * 0.5;
        Ok(Rows {
            range,
            whole: self.row(I::new(range[0], range[1]), work)?,
            middle: self.row(I::point(middle), work)?,
            ends: [
                self.row(I::point(range[0]), work)?,
                self.row(I::point(range[1]), work)?,
            ],
            jumps: self.left.angular_jumps(range, work)? + self.right.angular_jumps(range, work)?,
        })
    }
    pub fn jet_at_row(
        &self,
        body: &RadialBoundary,
        u: I,
        second: bool,
        row: &Row,
        work: &mut Work<'_>,
    ) -> Result<Jet, Error> {
        work.charge(128)?;
        let chart = &body.ordered.caps.axial.chart;
        let uv = [u, row.parameter];
        let p = row.point;
        let gap = row.gap;
        let [speed_p, speed_q] = row.speed;
        let [accel_p, accel_q] = row.acceleration;
        let delta = self.fractions[1] - self.fractions[0];
        let lambda = self.fractions[0] + delta * uv[1];
        let two = I::point(2.);
        let (w, wu, wv) = if second {
            (
                lambda + (two * uv[0] - I::point(1.)) * (I::point(1.) - lambda),
                two * (I::point(1.) - lambda),
                two * (I::point(1.) - uv[0]) * delta,
            )
        } else {
            (two * uv[0] * lambda, two * lambda, two * uv[0] * delta)
        };
        let phase = w * gap;
        let (sine, cosine) = (phase.sin(), phase.cos());
        let norm = row.norm;
        let radial = [
            (p[0] * cosine - p[1] * sine) / norm,
            (p[0] * sine + p[1] * cosine) / norm,
        ];
        let angular = [-radial[1], radial[0]];
        let theta_u = wu * gap;
        let theta_v = speed_p + wv * gap + w * (speed_q - speed_p);
        let theta_uv =
            (if second { -two * delta } else { two * delta }) * gap + wu * (speed_q - speed_p);
        let theta_vv = accel_p + two * wv * (speed_q - speed_p) + w * (accel_q - accel_p);
        let height = I::point(chart.nominal_height);
        let z = height * uv[1];
        let radius = (chart.radius + chart.slope * z).positive()?;
        let rp = chart.slope * height;
        let hessian_norms = [
            radius * theta_u.square(),
            ((rp * theta_u + radius * theta_uv).square() + (radius * theta_u * theta_v).square())
                .sqrt(),
            ((two * rp * theta_v + radius * theta_vv).square()
                + (radius * theta_v.square()).square())
            .sqrt(),
        ];
        let local = [radius * radial[0], radius * radial[1], z];
        let du = [
            radius * theta_u * angular[0],
            radius * theta_u * angular[1],
            I::point(0.),
        ];
        let dv = [
            chart.slope * height * radial[0] + radius * theta_v * angular[0],
            chart.slope * height * radial[1] + radius * theta_v * angular[1],
            height,
        ];
        let normal = [radial[0], radial[1], -chart.slope];
        let vector = |p: [I; 3]| {
            std::array::from_fn(|k| {
                chart.axes[0][k] * p[0] + chart.axes[1][k] * p[1] + chart.axes[2][k] * p[2]
            })
        };
        let origin_delta = std::array::from_fn(|k| {
            I::point(chart.bottom_point[k]) - I::point(chart.axis_origin[k])
        });
        let station = dot(origin_delta, chart.axes[2]);
        let p = vector([local[0], local[1], local[2] + station]);
        let position = std::array::from_fn(|k| I::point(chart.axis_origin[k]) + p[k]);
        let jet = Jet {
            position,
            du: vector(du),
            dv: vector(dv),
            normal: vector(normal),
            radial,
            angular_rate: [theta_u, theta_v],
            hessian_norms,
        };
        for value in jet
            .position
            .into_iter()
            .chain(jet.du)
            .chain(jet.dv)
            .chain(jet.normal)
        {
            value.checked()?;
        }
        Ok(jet)
    }
    pub fn normal_budget(&self, requested: f64) -> Result<f64, Error> {
        Ok((I::point(requested) - I::point(self.normal_reserve))
            .positive()?
            .lo)
    }
    pub fn split_at(&self, range: [f64; 2]) -> f64 {
        let middle = range[0] + (range[1] - range[0]) * 0.5;
        let a = super::choose_split(&self.left.knots, range);
        let b = super::choose_split(&self.right.knots, range);
        if (a - middle).abs() < (b - middle).abs() {
            a
        } else {
            b
        }
    }
    /// Compare the complete source sector cell to the union of its two
    /// emitted triangles. A frozen-angle frustum quadrilateral is planar
    /// and convex; its circle sag is quadratic in angular width. This
    /// avoids a parameter-correspondence bound's spurious first-order
    /// mixed error on long, straight cone generators.
    pub fn cell_error(
        &self,
        body: &RadialBoundary,
        u_range: [f64; 2],
        rows: &Rows,
        vertices: [Point; 4],
        second: bool,
        work: &mut Work<'_>,
    ) -> Result<(f64, [I; 3]), Error> {
        work.charge(512)?;
        let chart = &body.ordered.caps.axial.chart;
        let range = [u_range, rows.range];
        let uv = range.map(|d| I::new(d[0], d[1]));
        let whole = self.jet_at_row(body, uv[0], second, &rows.whole, work)?;
        let middle = rows.middle.parameter.lo;
        let left = self.jet_at_row(body, I::point(range[0][0]), second, &rows.middle, work)?;
        let right = self.jet_at_row(body, I::point(range[0][1]), second, &rows.middle, work)?;
        let angle =
            (left.angular_rate[0] * (I::point(range[0][1]) - I::point(range[0][0]))).positive()?;
        if angle.hi >= std::f64::consts::PI.next_down() {
            return Err(Error::Unresolved);
        }
        let height = I::point(chart.nominal_height);
        let radius = (chart.radius + chart.slope * height * uv[1]).positive()?;
        let maximum = I::point(radius.hi);
        // At a fixed u, integrating the whole-cell angular-rate bound
        // controls replacement by that ray at the middle height.
        let excursion = (uv[1] - I::point(middle)).checked()?.magnitude();
        let twist = maximum * I::point(whole.angular_rate[1].magnitude()) * I::point(excursion);
        let sag = maximum * (I::point(1.) - (angle / I::point(2.)).cos());
        let delta = std::array::from_fn(|k| {
            I::point(chart.bottom_point[k]) - I::point(chart.axis_origin[k])
        });
        let station = dot(delta, chart.axes[2]);
        let mut corner: f64 = 0.;
        for (i, (ray, v)) in [
            (left.radial, range[1][0]),
            (right.radial, range[1][0]),
            (right.radial, range[1][1]),
            (left.radial, range[1][1]),
        ]
        .into_iter()
        .enumerate()
        {
            let z = height * I::point(v);
            let r = chart.radius + chart.slope * z;
            let p = std::array::from_fn::<_, 3, _>(|k| {
                I::point(chart.axis_origin[k])
                    + chart.axes[0][k] * r * ray[0]
                    + chart.axes[1][k] * r * ray[1]
                    + chart.axes[2][k] * (station + z)
            });
            let error = (0..3)
                .fold(I::point(0.), |s, k| {
                    s + I::point((p[k] - I::point(vertices[i][k])).magnitude())
                })
                .checked()?;
            corner = corner.max(error.hi);
        }
        // Barycentric maps between the two virtual planar triangles and
        // their actual rounded counterparts cost at most this corner
        // error, in either direction. Clamp to the actual angular range
        // at each height to transfer the twist bound both ways.
        let frozen = (twist + sag + I::point(corner) + I::point(self.closure_error)).checked()?;
        let widths = range.map(|d| I::point(d[1]) - I::point(d[0]));
        let regular = (whole.hessian_norms[0] * widths[0].square()
            + I::point(2.) * whole.hessian_norms[1] * widths[0] * widths[1]
            + whole.hessian_norms[2] * widths[1].square())
            / I::point(8.);
        // Across a C0 rail knot only S_v jumps: S_u is continuous.
        // The barycentric interpolation gap of each hinge (v-k)_+ is
        // <=dv/4. Bound the actual jumps instead of assuming C1 sewing.
        let jumps = rows.jumps;
        let discontinuity = maximum * jumps * widths[1] / I::point(4.);
        let mut actual_corner = 0_f64;
        for (i, uv) in [
            [range[0][0], range[1][0]],
            [range[0][1], range[1][0]],
            [range[0][1], range[1][1]],
            [range[0][0], range[1][1]],
        ]
        .into_iter()
        .enumerate()
        {
            let end = usize::from(i >= 2);
            let p = self
                .jet_at_row(body, I::point(uv[0]), second, &rows.ends[end], work)?
                .position;
            let error = (0..3)
                .fold(I::point(0.), |s, k| {
                    s + I::point((p[k] - I::point(vertices[i][k])).magnitude())
                })
                .checked()?;
            actual_corner = actual_corner.max(error.hi);
        }
        let hessian =
            (regular + discontinuity + I::point(actual_corner) + I::point(self.closure_error))
                .checked()?;
        // Recenter the angular field before enclosing its image. Direct
        // rotation of whole-interval left/right rail boxes loses their
        // common longitudinal dependence. The actual angular derivatives
        // still bound every intervening source point, including piecewise
        // smooth rail knots, so a mean-value enclosure is independent of
        // that interval dependency.
        let middle_u = range[0][0] + (range[0][1] - range[0][0]) * 0.5;
        let center = self.jet_at_row(body, I::point(middle_u), second, &rows.middle, work)?;
        let variation = (0..2)
            .fold(I::point(0.), |sum, axis| {
                let origin = if axis == 0 { middle_u } else { middle };
                sum + I::point(whole.angular_rate[axis].magnitude())
                    * I::point((uv[axis] - I::point(origin)).magnitude())
            })
            .checked()?;
        let offset = I::new(-variation.hi, variation.hi);
        let (sine, cosine) = (offset.sin(), offset.cos());
        let local_normal = [
            center.radial[0] * cosine - center.radial[1] * sine,
            center.radial[0] * sine + center.radial[1] * cosine,
            -chart.slope,
        ];
        let centered = std::array::from_fn::<_, 3, _>(|k| {
            chart.axes[0][k] * local_normal[0]
                + chart.axes[1][k] * local_normal[1]
                + chart.axes[2][k] * local_normal[2]
        });
        let mut normal = whole.normal;
        for k in 0..3 {
            normal[k] = I::new(
                normal[k].lo.max(centered[k].lo),
                normal[k].hi.min(centered[k].hi),
            )
            .checked()?;
        }
        Ok((frozen.hi.min(hessian.hi), normal))
    }
}

/// Whole finite-trim angular collar, including steep cone conditioning.
pub(super) fn angular_reserve(radius: I, slope: I, error_mm: f64) -> Result<f64, Error> {
    // Endpoint and support sewing can move the finite angular trim
    // beyond the ideal retained rail. Reserve that normal angle;
    // no extra requested angle is spent implicitly.
    let lower = radius - I::point(error_mm) * (I::point(8.) + I::point(slope.magnitude()));
    let angular_error = if lower.lo > 0. {
        (I::point(16.) * I::point(error_mm) / I::point(lower.lo))
            .checked()?
            .hi
            .min(std::f64::consts::PI.next_up())
    } else {
        std::f64::consts::PI.next_up()
    };
    // A curve end -> shared vertex -> neighboring curve support ->
    // cone support -> recorded closure collar uses at most four
    // displacements, each <=e. Raw end-curve heights are within e of
    // the exact cap planes. Cone radius variation costs |k|e; the
    // padded 8e covers the radial support and displacement chain.
    // The still wider D<=8e chord bound below is therefore safe.
    // For delta<pi, the
    // chord-angle inequality 2 asin(D/(2r)) (D<=8e), followed by
    // asin(x)<=pi*x/2, is below16e/r. Otherwise the complete azimuth
    // is retained. A unit cone normal has angular speed
    // 1/sqrt(1+k^2), including the steep near-axial case.
    let k = if slope.lo > 0. {
        slope.lo
    } else if slope.hi < 0. {
        -slope.hi
    } else {
        0.
    };
    let normal_reserve = (I::point(angular_error) / (I::point(1.) + I::point(k).square()).sqrt())
        .checked()?
        .hi;
    Ok(normal_reserve)
}
