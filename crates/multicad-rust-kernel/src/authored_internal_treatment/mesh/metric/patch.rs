//! Exact finite source cells and whole-circle interval derivative composition.
use super::*;
pub(super) struct Strip {
    lo: Q,
    hi: Q,
    rail: i64,
    u_lo: Q,
    u_hi: Q,
}
pub(super) struct Patch<'a> {
    source: &'a AuthoredMetricThread,
    period: i64,
    band: usize,
    grid: Grid,
    offset: Q,
    base: Q,
    span: Q,
    lo: Q,
    lead: Q,
    cuts: [Line; 2],
    denominator: Q,
    origin: Q,
}
impl<'a> Patch<'a> {
    pub fn new(
        source: &'a AuthoredMetricThread,
        period: i64,
        band: usize,
        grid: Grid,
        w: &mut Work<'_>,
    ) -> Result<Self, Error> {
        w.step()?;
        let spec = source.spec();
        let pitch = Q::float(spec.pitch_mm, w)?;
        let crest = Q::float(spec.crest_width_mm, w)?;
        let phase = |index: usize, w: &mut Work<'_>| {
            let knot = [
                Knot::Begin,
                Knot::CrestEnd,
                Knot::UpperRoot,
                Knot::LowerRoot,
                Knot::CrestStart,
                Knot::End,
            ][index];
            let [(pn, pd), (cn, cd)] = phase_coefficients(knot);
            pitch
                .mul(&Q::integer(pn, w)?, w)?
                .div(&Q::integer(pd, w)?, w)?
                .add(
                    &crest
                        .mul(&Q::integer(cn, w)?, w)?
                        .div(&Q::integer(cd, w)?, w)?,
                    w,
                )
        };
        let base = phase(band, w)?;
        let end = phase(band + 1, w)?;
        let span = end.sub(&base, w)?;
        let offset = pitch.mul(&Q::integer(period, w)?, w)?;
        let lead = pitch.mul(&Q::integer(i64::from(spec.starts), w)?, w)?;
        let zero = Q::integer(0, w)?;
        let one = Q::integer(1, w)?;
        let cuts = [
            Line::new(one.clone(), lead.clone(), zero.clone(), w)?,
            Line::new(
                zero.sub(&one, w)?,
                zero.sub(&lead, w)?,
                Q::float(spec.length_mm, w)?,
                w,
            )?,
        ];
        let denominator = Q::integer(i64::from(grid.effective_rows()), w)?;
        let origin = if grid.shifted {
            Q::float(sampling::TURN_OFFSET, w)?
        } else {
            zero
        };
        Ok(Self {
            source,
            period,
            band,
            grid,
            lo: base.add(&offset, w)?,
            base,
            span,
            offset,
            lead,
            cuts,
            denominator,
            origin,
        })
    }

    pub fn strip(&self, index: usize, w: &mut Work<'_>) -> Result<Strip, Error> {
        w.step()?;
        let n = self.grid.counts[self.band];
        if index >= n {
            return Err(Error::Invalid);
        }
        let fraction = |i: usize, w: &mut Work<'_>| {
            let f = if i == 0 {
                0.
            } else if i == n {
                1.
            } else if self.band == 2 {
                w.charge(8)?;
                0.5 + ((i as f64 / n as f64 - 0.5) * (2. * std::f64::consts::PI / 3.))
                    .portable_sin()
                    / 3_f64.sqrt()
            } else {
                i as f64 / n as f64
            };
            Q::float(f, w)
        };
        let u_lo = fraction(index, w)?;
        let u_hi = fraction(index + 1, w)?;
        let lo = self.lo.add(&self.span.mul(&u_lo, w)?, w)?;
        let hi = self.lo.add(&self.span.mul(&u_hi, w)?, w)?;
        let rail = self.period * self.grid.total() as i64
            + self.grid.counts[..self.band].iter().sum::<usize>() as i64
            + index as i64;
        Ok(Strip {
            lo,
            hi,
            rail,
            u_lo,
            u_hi,
        })
    }
    pub fn cell(&self, strip: &Strip, row: u32, w: &mut Work<'_>) -> Result<Vec<Point>, Error> {
        let rows = self.grid.effective_rows();
        if row >= rows {
            return Err(Error::Invalid);
        }
        let vl = self.row_parameter(row, w)?;
        let vh = self.row_parameter(row + 1, w)?;
        #[cfg(test)]
        if super::tests::FORCE_CLIP.get() {
            return self.full_clip(strip, vl, vh, w);
        }
        w.charge(4)?;
        let low = strip.lo.enclosure(w)?.1 + self.lead.enclosure(w)?.1 * vl.enclosure(w)?.1;
        let high = strip.hi.enclosure(w)?.1 + self.lead.enclosure(w)?.1 * vh.enclosure(w)?.1;
        let length = self.source.spec().length_mm;
        if high.hi < 0. || low.lo > length {
            return Ok(Vec::new());
        }
        if low.lo > 0. && high.hi < length {
            // This is a proof about the ACTUAL affine source chart z=q+lead*v,
            // not a geometric bounding-box approximation. Strict slab bounds
            // make the unchanged exact clip return precisely these corners.
            return Ok(vec![
                Point {
                    q: strip.lo.clone(),
                    v: vl.clone(),
                    active: 0b0101,
                },
                Point {
                    q: strip.hi.clone(),
                    v: vl,
                    active: 0b0110,
                },
                Point {
                    q: strip.hi.clone(),
                    v: vh.clone(),
                    active: 0b1010,
                },
                Point {
                    q: strip.lo.clone(),
                    v: vh,
                    active: 0b1001,
                },
            ]);
        }
        self.full_clip(strip, vl, vh, w)
    }
    fn row_parameter(&self, row: u32, w: &mut Work<'_>) -> Result<Q, Error> {
        let rows = self.grid.effective_rows();
        if rows.is_power_of_two() {
            // row<=rows<=65536 and origin is0 or2^-17: the binary fraction
            // and sum are EXACT binary64 values (<2 and <=18 significant bits).
            let origin = if self.grid.shifted {
                sampling::TURN_OFFSET
            } else {
                0.
            };
            Q::float(f64::from(row) / f64::from(rows) + origin, w)
        } else {
            Q::integer(i64::from(row), w)?
                .div(&self.denominator, w)?
                .add(&self.origin, w)
        }
    }
    fn full_clip(
        &self,
        strip: &Strip,
        vl: Q,
        vh: Q,
        w: &mut Work<'_>,
    ) -> Result<Vec<Point>, Error> {
        Ok(clip::polygon([strip.lo.clone(), strip.hi.clone(), vl, vh], &self.cuts, w)?.0)
    }

    fn u(&self, q: &Q, strip: &Strip, w: &mut Work<'_>) -> Result<Q, Error> {
        if q.cmp(&strip.lo, w)? == Ordering::Equal {
            return Ok(strip.u_lo.clone());
        }
        if q.cmp(&strip.hi, w)? == Ordering::Equal {
            return Ok(strip.u_hi.clone());
        }
        q.sub(&self.offset, w)?
            .sub(&self.base, w)?
            .div(&self.span, w)
    }
    fn unit(u: &Q, w: &mut Work<'_>) -> Result<I, Error> {
        let zero = Q::integer(0, w)?;
        let one = Q::integer(1, w)?;
        if u.cmp(&zero, w)? == Ordering::Less || u.cmp(&one, w)? == Ordering::Greater {
            return Err(Error::Unresolved);
        }
        if u.cmp(&zero, w)? == Ordering::Equal {
            return Ok(I::point(0.));
        }
        if u.cmp(&one, w)? == Ordering::Equal {
            return Ok(I::point(1.));
        }
        let value = u.enclosure(w)?.1;
        Ok(I::new(value.lo.max(0.), value.hi.min(1.)))
    }
    pub fn field(
        &self,
        points: &[Point],
        strip: &Strip,
        w: &mut Work<'_>,
    ) -> Result<source::Field, Error> {
        w.charge(32)?;
        let mut low = &points[0].q;
        let mut high = low;
        let mut vb = I::new(f64::INFINITY, f64::NEG_INFINITY);
        for p in points {
            if p.q.cmp(low, w)? == Ordering::Less {
                low = &p.q;
            }
            if p.q.cmp(high, w)? == Ordering::Greater {
                high = &p.q;
            }
            let v = p.v.enclosure(w)?.1;
            vb.lo = vb.lo.min(v.lo);
            vb.hi = vb.hi.max(v.hi);
        }
        // u(q) is exact affine increasing. Its exact extrema on the convex
        // polygon occur at q extrema; interval conversion need not be repeated
        // for every non-extreme point or artificially union their roundoff.
        let ub = I::new(
            Self::unit(&self.u(low, strip, w)?, w)?.lo,
            Self::unit(&self.u(high, strip, w)?, w)?.hi,
        );
        let axis = |value: I, i: usize| {
            Jet::constant(I::point(value.lo))
                .add(Jet::variable(i).scale(I::point(value.hi) - I::point(value.lo)))
        };
        let u = axis(ub, 0);
        let v = axis(vb, 1);
        let span = self.span.enclosure(w)?.1;
        let bounds = self
            .source
            .profile_bounds(self.band, pair(ub))
            .map_err(source_error)?;
        let r = interval(bounds.radius_mm);
        let m = interval(bounds.slope);
        let dd = interval(bounds.second_derivative_per_mm);
        // Source owns the complete rho/rho_q/rho_qq enclosure on this exact band.
        // Affine q(u) composition preserves all first/second chain-rule terms.
        let derivative = m * span;
        let second = dd * span.square();
        let rho = Jet {
            value: r,
            first: u.first.map(|x| derivative * x),
            second: [
                derivative * u.second[0] + second * u.first[0].square(),
                derivative * u.second[1] + second * u.first[0] * u.first[1],
                derivative * u.second[2] + second * u.first[1].square(),
            ],
        };
        let q = u.scale(span).add(Jet::constant(
            self.base.add(&self.offset, w)?.enclosure(w)?.1,
        ));
        let lead = self.lead.enclosure(w)?.1;
        let z = q.add(v.scale(lead));
        let tau = I::point(2.)
            * I::new(
                std::f64::consts::PI.next_down(),
                std::f64::consts::PI.next_up(),
            );
        let tau = if self.source.spec().left_hand {
            -tau
        } else {
            tau
        };
        let angle = v.scale(tau);
        Ok(source::Field {
            point: [rho.mul(angle.cos()), rho.mul(angle.sin()), z],
            normal: self.normal(angle.value, lead / tau, r, m, &bounds),
        })
    }
    fn normal(
        &self,
        angle: I,
        rate: I,
        r: I,
        m: I,
        bounds: &crate::authored_metric_thread::ProfileBounds,
    ) -> source::Normal {
        #[cfg(test)]
        if super::tests::OLD_ROOT_NORMAL.get() {
            return source::Normal::Radial {
                angle,
                beta: rate * m / r,
                z: -m,
            };
        }
        if self.band == 2 {
            // H>=1/2 is source-owned and strictly positive. Multiplying the
            // actual gradient by H preserves orientation and angle, while
            // removing the avoidable interval quotients X/H and X/(H*r).
            let radial = interval(bounds.normal_radial);
            let z = interval(bounds.normal_axial);
            source::Normal::ScaledRadial {
                angle,
                radial,
                beta: rate * (-z) / r,
                z,
            }
        } else {
            source::Normal::Radial {
                angle,
                beta: rate * m / r,
                z: -m,
            }
        }
    }

    pub fn identity(
        &self,
        p: &Point,
        strip: &Strip,
        row: u32,
        w: &mut Work<'_>,
    ) -> Result<(Key, Box<[u8]>), Error> {
        let zero = Q::integer(0, w)?;
        let one = Q::integer(1, w)?;
        let origin = if self.grid.shifted {
            Q::float(sampling::TURN_OFFSET, w)?
        } else {
            zero.clone()
        };
        let end = one.add(&origin, w)?;
        let wrap = p.v.cmp(&end, w)? == Ordering::Equal;
        let q = if wrap {
            p.q.add(&self.lead, w)?
        } else {
            p.q.clone()
        };
        let v = if wrap { origin } else { p.v.clone() };
        let z = p.q.add(&self.lead.mul(&p.v, w)?, w)?;
        let shift = if wrap {
            i64::from(self.source.spec().starts) * self.grid.total() as i64
        } else {
            0
        };
        let rail = if p.active & 1 != 0 {
            Some(strip.rail + shift)
        } else if p.active & 2 != 0 {
            Some(strip.rail + 1 + shift)
        } else {
            None
        };
        let rows = self.grid.effective_rows();
        let row = if p.active & 4 != 0 {
            Some(row)
        } else if p.active & 8 != 0 {
            Some(row + 1)
        } else {
            None
        }
        .map(|n| if n == rows { 0 } else { n });
        let station = if p.active & 16 != 0 {
            Some(0)
        } else if p.active & 32 != 0 {
            Some(2)
        } else {
            None
        };
        let key = match (rail, row, station) {
            (Some(r), _, Some(s)) => Key([Constraint::Rail(r), Constraint::Station(s, 0)]),
            (Some(r), Some(v), None) => Key([Constraint::Rail(r), Constraint::Row(v)]),
            (None, Some(v), Some(s)) => Key([Constraint::Row(v), Constraint::Station(s, 0)]),
            _ => return Err(Error::Unresolved),
        };
        // Same immutable continuous source + exact canonical(q,v,z) determines
        // radius uniquely, including the shared root/flank tangent identities.
        let encoded = Q::encode(&[&q, &v, &z], w)?;
        Ok((key, encoded))
    }
    pub fn point(&self, p: &Point, strip: &Strip, w: &mut Work<'_>) -> Result<[I; 3], Error> {
        w.charge(16)?;
        let u = Self::unit(&self.u(&p.q, strip, w)?, w)?;
        let r = interval(
            self.source
                .profile_bounds(self.band, pair(u))
                .map_err(source_error)?
                .radius_mm,
        );
        let origin = if self.grid.shifted {
            Q::float(sampling::TURN_OFFSET, w)?
        } else {
            Q::integer(0, w)?
        };
        let end = Q::integer(1, w)?.add(&origin, w)?;
        let v = if p.v.cmp(&end, w)? == Ordering::Equal {
            origin
        } else {
            p.v.clone()
        };
        let tau = I::point(2.)
            * I::new(
                std::f64::consts::PI.next_down(),
                std::f64::consts::PI.next_up(),
            );
        let angle = v.enclosure(w)?.1
            * if self.source.spec().left_hand {
                -tau
            } else {
                tau
            };
        let z = p.q.add(&self.lead.mul(&p.v, w)?, w)?.enclosure(w)?.1;
        Ok([r * angle.cos(), r * angle.sin(), z])
    }
    pub fn face(&self) -> u32 {
        let n = i64::from(self.source.spec().starts);
        let group = (self.period + [0, 0, 0, 1, 1][self.band]).rem_euclid(n);
        let group = if self.source.spec().left_hand {
            group
        } else {
            (n - group) % n
        };
        2 + 4 * group as u32 + [1, 2, 3, 0, 1][self.band]
    }
}
fn interval(v: [f64; 2]) -> I {
    I::new(v[0], v[1])
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn exact_slab_preflight_preserves_tangencies_one_ulp_and_shifted_seams() {
        for length in [1_f64.next_down(), 1., 1_f64.next_up()] {
            let mut spec = super::super::tests::small_source().spec();
            spec.length_mm = length;
            spec.starts = 2;
            let source = AuthoredMetricThread::new(spec).unwrap();
            let mut w = Work {
                remaining: 16_000_000,
                used: 0,
                interrupt: None,
            };
            for (odd, shifted) in [(false, false), (true, false), (false, true)] {
                let grid = Grid {
                    counts: [1, 1, 4, 1, 1],
                    rows: 16,
                    odd,
                    shifted,
                };
                for period in [-3, -2, -1, 0, 1, 2, 3] {
                    for band in 0..5 {
                        let patch = Patch::new(&source, period, band, grid, &mut w).unwrap();
                        for index in 0..grid.counts[band] {
                            let strip = patch.strip(index, &mut w).unwrap();
                            for row in 0..grid.effective_rows() {
                                let fast = patch.cell(&strip, row, &mut w).unwrap();
                                let vl = patch.row_parameter(row, &mut w).unwrap();
                                let vh = patch.row_parameter(row + 1, &mut w).unwrap();
                                let exact = patch.full_clip(&strip, vl, vh, &mut w).unwrap();
                                assert_eq!(fast.len(), exact.len());
                                assert!(fast.is_empty() || (3..=6).contains(&fast.len()));
                                for (a, b) in fast.iter().zip(&exact) {
                                    assert_eq!(a.q.cmp(&b.q, &mut w).unwrap(), Ordering::Equal);
                                    assert_eq!(a.v.cmp(&b.v, &mut w).unwrap(), Ordering::Equal);
                                    assert_eq!(a.active, b.active);
                                    assert_eq!(
                                        patch.identity(a, &strip, row, &mut w).unwrap(),
                                        patch.identity(b, &strip, row, &mut w).unwrap()
                                    );
                                    assert_eq!(
                                        patch.point(a, &strip, &mut w).unwrap().map(pair),
                                        patch.point(b, &strip, &mut w).unwrap().map(pair)
                                    );
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
#[cfg(test)]
mod guards {
    use super::*;
    #[test]
    fn private_chart_domains_refuse_invalid_nodes_without_changing_source() {
        let source = super::super::tests::small_source();
        let grid = Grid {
            counts: [2, 2, 4, 2, 2],
            rows: 16,
            odd: false,
            shifted: false,
        };
        let mut w = Work {
            remaining: 16_000_000,
            used: 0,
            interrupt: None,
        };
        let patch = Patch::new(&source, 0, 1, grid, &mut w).unwrap();
        assert!(patch.strip(2, &mut w).is_err());
        let strip = patch.strip(0, &mut w).unwrap();
        assert!(patch.cell(&strip, 16, &mut w).is_err());
        assert!(Patch::unit(&Q::integer(-1, &mut w).unwrap(), &mut w).is_err());
        let polygon = patch.cell(&strip, 1, &mut w).unwrap();
        let mut invalid = polygon[0].clone();
        invalid.active = 0;
        assert!(patch.identity(&invalid, &strip, 1, &mut w).is_err());
        // Actual repeated native operation budgets reach the second cut-line
        // construction refusal, preserving the source and its exact supports.
        let used = w.used;
        assert!(used > 0);
        for budget in 1..200 {
            let mut limited = Work {
                remaining: budget,
                used: 0,
                interrupt: None,
            };
            let _ = Patch::new(&source, 0, 1, grid, &mut limited);
        }
    }
}
