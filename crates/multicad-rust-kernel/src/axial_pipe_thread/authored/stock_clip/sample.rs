use super::*;
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StockClipPoint {
    position: [[f64; 2]; 3],
    radius: [f64; 2],
    phase: [f64; 2],
    axial: [f64; 2],
    work: usize,
    bytes: usize,
}
impl StockClipPoint {
    pub fn position_mm(&self) -> [[f64; 2]; 3] {
        self.position
    }
    pub fn radius_mm(&self) -> [f64; 2] {
        self.radius
    }
    /// Unwrapped phase of this chart. Incident periodic curves can use
    /// different lifts at the SAME physical endpoint owner.
    pub fn phase_mm(&self) -> [f64; 2] {
        self.phase
    }
    pub fn axial_mm(&self) -> [f64; 2] {
        self.axial
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.bytes
    }
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StockClipSample {
    point: StockClipPoint,
    normal: [[f64; 2]; 3],
}
impl StockClipSample {
    pub fn point(&self) -> &StockClipPoint {
        &self.point
    }
    /// Outward one-sided unit normal on this source patch, including its edge
    /// limits. This is not a literal facet or manufacturing-normal certificate.
    pub fn normal(&self) -> [[f64; 2]; 3] {
        self.normal
    }
}
pub(super) fn unit_parameter(p: [f64; 2]) -> Result<I, Error> {
    let p = interval(p)?;
    if p.lo < 0. || p.hi > 1. {
        return Err(Error::Invalid);
    }
    Ok(p)
}
fn exact(p: I, x: f64) -> bool {
    p.lo == x && p.hi == x
}
fn mix(a: I, b: I, t: I) -> I {
    if exact(t, 0.) {
        a
    } else if exact(t, 1.) {
        b
    } else {
        a + (b - a) * t
    }
}
struct Evaluated {
    point: StockClipPoint,
    angle: I,
    radius: I,
}
impl StockClipPlan {
    fn evaluate<'a>(
        &self,
        patch: StockClipPatch,
        u: [f64; 2],
        v: [f64; 2],
        options: StockClipOptions,
        interrupt: Option<&'a dyn Fn() -> bool>,
    ) -> Result<(Evaluated, Meter<'a>), Error> {
        let mut w = Meter::new(
            options,
            size_of::<Self>() + size_of::<StockClipSample>(),
            interrupt,
        )?;
        w.charge(self.work)?;
        let u = unit_parameter(u)?;
        let v = unit_parameter(v)?;
        // Fixed arithmetic/trig enclosure work; all charts are affine r/z/q
        // compositions and the angular range is bounded by the sealed source.
        w.charge(512)?;
        let (a, b) = match patch.slab {
            StockClipSlab::Before => (StockClipStation::Start, StockClipStation::Contact),
            StockClipSlab::After => (StockClipStation::Contact, StockClipStation::End),
        };
        let station = if exact(v, 0.) {
            Some(a)
        } else if exact(v, 1.) {
            Some(b)
        } else {
            None
        };
        let z = mix(self.height(a), self.height(b), v);
        let after = patch.slab == StockClipSlab::After;
        let d = I::point(self.source.depth());
        let k = I::point(self.source.spec.radial_taper);
        let limit = I::point(self.stock_radius_limit_mm());
        let root = if station == Some(StockClipStation::Contact) {
            limit - d
        } else {
            self.root_radius + k * z
        };
        let crest = if after || station == Some(StockClipStation::Contact) {
            limit
        } else {
            I::point(self.source.spec.major_radius_mm) + k * z
        };
        let contacts = if after {
            self.contacts(z, station)
        } else {
            [-self.half, self.half]
        };
        let i = band_index(patch.band);
        let qends = [
            self.q0,
            contacts[0],
            contacts[1],
            self.q3,
            self.q0 + I::point(self.source.spec.pitch_mm),
        ];
        let q = mix(qends[i], qends[i + 1], u);
        // True authored periodic identity for physical endpoint evaluation;
        // retain the original unwrapped q in the chart/curve receipt below.
        let angle_q = if patch.band == StockClipBand::Root && exact(u, 1.) {
            self.q0
        } else {
            q
        };
        let (ra, rb) = match patch.band {
            StockClipBand::Lower => (root, crest),
            StockClipBand::Crest => (crest, crest),
            StockClipBand::Upper => (crest, root),
            StockClipBand::Root => (root, root),
        };
        let raw = if matches!(patch.band, StockClipBand::Crest | StockClipBand::Root) {
            ra
        } else {
            mix(ra, rb, u)
        };
        // Source-domain theorems, not a tolerance: every retained patch is
        // between the original positive root and the unchanged stock ceiling.
        let radius = I::new(
            raw.lo.max(self.radius_range.lo),
            raw.hi.min(self.radius_range.hi),
        );
        let angle = (z - angle_q) * (I::point(2.) * pi() / I::point(self.source.spec.pitch_mm));
        let position = [radius * angle.cos(), radius * angle.sin(), z];
        // The sealed source bounds r<=1e7, z<=1e7 and |q|<3p,
        // with L/p<128. Trig returns [-1,1], so these finite coordinate
        // operations need no extra numerical admission after the input checks.
        w.charge(1)?;
        let point = StockClipPoint {
            position: position.map(bounds),
            radius: bounds(radius),
            phase: bounds(q),
            axial: bounds(z),
            work: w.used,
            bytes: size_of::<Self>() + size_of::<StockClipSample>() + SCRATCH,
        };
        Ok((
            Evaluated {
                point,
                angle,
                radius,
            },
            w,
        ))
    }
    pub fn point_bounds(
        &self,
        patch: StockClipPatch,
        u: [f64; 2],
        v: [f64; 2],
        options: StockClipOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<StockClipPoint, Error> {
        self.evaluate(patch, u, v, options, interrupt)
            .map(|(p, _)| p.point)
    }
    pub fn sample_bounds(
        &self,
        patch: StockClipPatch,
        u: [f64; 2],
        v: [f64; 2],
        options: StockClipOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<StockClipSample, Error> {
        let (p, mut w) = self.evaluate(patch, u, v, options, interrupt)?;
        w.charge(128)?;
        let m = I::point(3.).sqrt();
        let k = I::point(self.source.spec.radial_taper);
        // Cancel the earned positive angular/axial chart factors analytically.
        // For r=a(q)+k*z, q=z-p*theta/(2pi), the physical normal is
        // [1, a'(q)*p/(2pi*r), -(a'(q)+k)] in the cylindrical basis.
        // This also handles arbitrarily small positive crest phase widths
        // without squaring their tiny parameter-speed coefficients.
        let (a, nz) = match patch.band {
            StockClipBand::Lower => (m - k, -m),
            StockClipBand::Upper => (-m - k, m),
            StockClipBand::Crest if patch.slab == StockClipSlab::After => {
                (I::point(0.), I::point(0.))
            }
            _ => (I::point(0.), -k),
        };
        let nt = a * I::point(self.source.spec.pitch_mm) / (I::point(2.) * pi() * p.radius);
        // AuthoredThread admission earns root>1024eps*(R+p+L+d).
        // Thus |nt|<1e14, |nz|<2 and the positive normal denominator
        // is >=1 (up to outward ulps), far inside the finite numeric range.
        let norm = (I::point(1.) + nt.square() + nz.square()).sqrt();
        let normal = [
            (p.angle.cos() - nt * p.angle.sin()) / norm,
            (p.angle.sin() + nt * p.angle.cos()) / norm,
            nz / norm,
        ];
        // Every exact coordinate belongs to a unit vector. Intersecting this
        // known range avoids interval dependency producing components >1.
        let out =
            std::array::from_fn(|i| bounds(I::new(normal[i].lo.max(-1.), normal[i].hi.min(1.))));
        w.charge(1)?;
        let mut point = p.point;
        point.work = w.used;
        Ok(StockClipSample { point, normal: out })
    }
    pub fn vertex_bounds(
        &self,
        key: StockClipVertex,
        options: StockClipOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<StockClipPoint, Error> {
        let i = rail_index(key.rail);
        let (band, u) = if i == 0 {
            (StockClipBand::Lower, 0.)
        } else {
            (BANDS[i - 1], 1.)
        };
        let (slab, v) = match key.station {
            StockClipStation::Start => (StockClipSlab::Before, 0.),
            StockClipStation::Contact => (StockClipSlab::Before, 1.),
            StockClipStation::End => (StockClipSlab::After, 1.),
        };
        self.point_bounds(
            StockClipPatch { slab, band },
            [u, u],
            [v, v],
            options,
            interrupt,
        )
    }
    pub fn curve_bounds(
        &self,
        key: StockClipCurve,
        t: [f64; 2],
        options: StockClipOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<StockClipPoint, Error> {
        let (patch, u, v) = match key {
            StockClipCurve::Axial { slab, rail } => (
                StockClipPatch {
                    slab,
                    band: BANDS[rail_index(rail)],
                },
                [0., 0.],
                t,
            ),
            StockClipCurve::Ring { station, band } => {
                let (slab, v) = match station {
                    StockClipStation::Start => (StockClipSlab::Before, 0.),
                    StockClipStation::Contact => (StockClipSlab::Before, 1.),
                    StockClipStation::End => (StockClipSlab::After, 1.),
                };
                (StockClipPatch { slab, band }, t, [v, v])
            }
        };
        self.point_bounds(patch, u, v, options, interrupt)
    }
}
