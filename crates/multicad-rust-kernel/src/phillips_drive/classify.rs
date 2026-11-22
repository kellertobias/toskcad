use super::*;
impl Classification {
    fn and(self, b: Self) -> Self {
        match (self, b) {
            (Self::Outside, _) | (_, Self::Outside) => Self::Outside,
            (Self::Inside, Self::Inside) => Self::Inside,
            _ => Self::UnresolvedBoundary,
        }
    }
    fn or(self, b: Self) -> Self {
        match (self, b) {
            (Self::Inside, _) | (_, Self::Inside) => Self::Inside,
            (Self::Outside, Self::Outside) => Self::Outside,
            _ => Self::UnresolvedBoundary,
        }
    }
    fn not(self) -> Self {
        match self {
            Self::Inside => Self::Outside,
            Self::Outside => Self::Inside,
            Self::UnresolvedBoundary => Self::UnresolvedBoundary,
        }
    }
}
fn positive(v: I) -> Classification {
    if v.lo > 0. {
        Classification::Inside
    } else if v.hi < 0. {
        Classification::Outside
    } else {
        Classification::UnresolvedBoundary
    }
}
fn select(value: I, split: I, left: Classification, right: Classification) -> Classification {
    if value.hi < split.lo {
        left
    } else if value.lo > split.hi {
        right
    } else if left == right {
        left
    } else {
        Classification::UnresolvedBoundary
    }
}
fn absolute(v: I) -> I {
    let lo = if v.lo <= 0. && v.hi >= 0. {
        0.
    } else {
        v.lo.abs().min(v.hi.abs())
    };
    I::new(lo, v.magnitude())
}
impl Corner {
    pub(super) fn coordinates(&self, p: [I; 3]) -> [I; 3] {
        let [x, y] = rotate([p[0], p[1]], (4 - self.quarter_turn) % 4);
        // The origin and extrusion expressions cancel analytically before
        // evaluation. This retains the true normal C-C coordinate without a
        // large intermediate origin or projection by cos(beta).
        let s = self.cos_beta * (x - self.b_half) + self.sin_beta * (p[2] - self.b_depth);
        let parameter = p[2] - I::point(self.domain.corner_start_depth_mm) - s * self.sin_beta;
        [s, y, parameter]
    }
    fn classify(&self, p: [I; 3]) -> Classification {
        let [s, y, parameter] = self.coordinates(p);
        let y = absolute(y);
        let inner = match self.inner {
            SourceInner::ShallowV { sine, cosine } => positive(s * sine - y * cosine),
            SourceInner::Rounded { radius_mm } => {
                let r = I::point(radius_mm);
                positive(s - r).or(positive(r.square() - (s - r).square() - y.square()))
            }
        };
        let outer =
            positive((s - self.join_s) * self.sin_theta - (y - self.join_y) * self.cos_theta);
        let profile = select(y, self.join_y, inner, outer);
        let closure = positive(self.far_s - s);
        let start = positive(parameter);
        let end = positive(
            I::point(self.domain.corner_end_depth_mm)
                - I::point(self.domain.corner_start_depth_mm)
                - parameter,
        );
        profile.and(closure).and(start).and(end)
    }
}
pub(super) fn point(cutter: &Cutter, p: [f64; 3]) -> Result<Classification, Error> {
    if !p.into_iter().all(f64::is_finite) {
        return Err(Error::Invalid);
    }
    let e = &cutter.envelope;
    if p[2] < 0. || p[2] > e.tip.hi || p[0].abs() > e.entry.hi || p[1].abs() > e.entry.hi {
        return Ok(Classification::Outside);
    }
    enclosure(cutter, p.map(I::point))
}
/// Classify the entire literal coordinate enclosure, including placement
/// subtraction error. A straddling surface remains unresolved.
pub(super) fn enclosure(cutter: &Cutter, p: [I; 3]) -> Result<Classification, Error> {
    if p.into_iter().any(|v| !v.valid()) {
        return Err(Error::Invalid);
    }
    let e = &cutter.envelope;
    if p[2].hi < 0.
        || p[2].lo > e.tip.hi
        || p[0].lo > e.entry.hi
        || p[0].hi < -e.entry.hi
        || p[1].lo > e.entry.hi
        || p[1].hi < -e.entry.hi
    {
        return Ok(Classification::Outside);
    }
    let rho = (p[0].square() + p[1].square()).sqrt();
    let r = I::point(cutter.spec.mouth_radius_mm);
    let radial = e.entry - rho;
    let roundover =
        positive(radial).and(positive(radial.square() + (p[2] - r).square() - r.square()));
    let flank = positive(
        (I::point(cutter.spec.mouth_diameter_mm) / I::point(2.) - rho) * e.cos_a - p[2] * e.sin_a,
    );
    let floor = positive((e.tip - p[2]) * e.cos_c - rho * e.sin_c);
    let shape = select(p[2], e.join, roundover, select(p[2], e.flank, flank, floor));
    let mut result = shape.and(positive(p[2])).and(positive(e.tip - p[2]));
    for corner in &cutter.corners {
        if result == Classification::Outside {
            break;
        }
        result = result.and(corner.classify(p).not());
    }
    Ok(result)
}
