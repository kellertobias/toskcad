//! Exact q/v chart, v being one angular revolution. True 2pi cancels from
//! source clipping; it remains explicit in all geometric point/normal proofs.
use super::{
    clip::{self, Line, Point},
    exact::Q,
    owners::{Constraint, Key},
    *,
};
use std::cmp::Ordering;
type Footpoint = (Key, Box<[u8]>, [I; 3]);
pub(super) struct Patch {
    pub band: usize,
    pub period: i64,
    pub subdivisions: usize,
    pub lead: Q,
    pub q0: Q,
    pub q1: Q,
    pub slope: Q,
    pub intercept: Q,
    pub entry: Option<(Q, Q)>,
    pub top: Q,
    offset: Option<Q>,
}
pub(super) struct Polygon {
    pub points: Vec<Point>,
    pub stations: [u8; 2],
    pub kind: Kind,
    pub rail0: i64,
    pub row: u32,
    pub rows: u32,
}
impl Patch {
    pub fn new(
        body: &Cavity,
        period: i64,
        band: usize,
        subdivisions: usize,
        w: &mut Work<'_>,
    ) -> Result<Self, Error> {
        Self::with_offset(body, period, band, subdivisions, false, w)
    }
    pub fn with_offset(
        body: &Cavity,
        period: i64,
        band: usize,
        subdivisions: usize,
        shifted: bool,
        w: &mut Work<'_>,
    ) -> Result<Self, Error> {
        Self::from_chart(
            ThreadChart::treated(body),
            period,
            band,
            subdivisions,
            shifted,
            w,
        )
    }
    pub fn from_chart(
        body: ThreadChart<'_>,
        period: i64,
        band: usize,
        subdivisions: usize,
        shifted: bool,
        w: &mut Work<'_>,
    ) -> Result<Self, Error> {
        w.step()?;
        let b = body.thread.bands()[band];
        let pitch = Q::float(body.thread.spec().pitch_mm, w)?;
        let lead = pitch.mul(&Q::integer(i64::from(body.thread.spec().starts), w)?, w)?;
        let offset = pitch.mul(&Q::integer(period, w)?, w)?;
        let a = Q::float(b.phase_mm[0], w)?;
        let z = Q::float(b.phase_mm[1], w)?;
        let r0 = Q::float(b.radius_mm[0], w)?;
        let r1 = Q::float(b.radius_mm[1], w)?;
        let slope = r1.sub(&r0, w)?.div(&z.sub(&a, w)?, w)?;
        let q0 = a.add(&offset, w)?;
        let q1 = z.add(&offset, w)?;
        let intercept = r0.sub(&slope.mul(&q0, w)?, w)?;
        let entry = body
            .entry
            .map(|e| -> Result<_, Error> {
                let radius = Q::float(e.radius, w)?;
                let k = radius
                    .sub(&Q::float(body.thread.minimum_radius_mm(), w)?, w)?
                    .div(&Q::float(e.depth, w)?, w)?;
                Ok((radius, k))
            })
            .transpose()?;
        let top = Q::float(body.top, w)?;
        Ok(Self {
            band,
            period,
            subdivisions,
            lead,
            q0,
            q1,
            slope,
            intercept,
            entry,
            top,
            offset: shifted
                .then(|| Q::float(sampling::TURN_OFFSET, w))
                .transpose()?,
        })
    }
    pub fn cell(
        &self,
        index: usize,
        row: u32,
        rows: u32,
        kind: Kind,
        w: &mut Work<'_>,
    ) -> Result<Polygon, Error> {
        if rows == 0 || row >= rows || index >= self.subdivisions {
            return Err(Error::Invalid);
        }
        let zero = Q::integer(0, w)?;
        let one = Q::integer(1, w)?;
        let span = self.q1.sub(&self.q0, w)?;
        let fraction = |value: usize, w: &mut Work<'_>| {
            Q::integer(value as i64, w)?.div(&Q::integer(self.subdivisions as i64, w)?, w)
        };
        let ql = self.q0.add(&span.mul(&fraction(index, w)?, w)?, w)?;
        let qh = self.q0.add(&span.mul(&fraction(index + 1, w)?, w)?, w)?;
        let denominator = Q::integer(i64::from(rows), w)?;
        let mut vl = Q::integer(i64::from(row), w)?.div(&denominator, w)?;
        let mut vh = Q::integer(i64::from(row + 1), w)?.div(&denominator, w)?;
        if let Some(offset) = &self.offset {
            vl = vl.add(offset, w)?;
            vh = vh.add(offset, w)?;
        }
        // All cuts retain their positive side. No source radius/phase change.
        let mouth = Line::new(one.clone(), self.lead.clone(), zero.clone(), w)?;
        let top = Line::new(
            zero.sub(&one, w)?,
            zero.sub(&self.lead, w)?,
            self.top.clone(),
            w,
        )?;
        let (cuts, stations) = match (kind, &self.entry) {
            (Kind::Cone, Some((a, k))) => {
                let seam = Line::new(
                    zero.sub(&self.slope.add(k, w)?, w)?,
                    zero.sub(&k.mul(&self.lead, w)?, w)?,
                    a.sub(&self.intercept, w)?,
                    w,
                )?;
                ([mouth, seam], [0, 1])
            }
            (Kind::Wall, Some((a, k))) => {
                let seam = Line::new(
                    self.slope.add(k, w)?,
                    k.mul(&self.lead, w)?,
                    self.intercept.sub(a, w)?,
                    w,
                )?;
                ([seam, top], [1, 2])
            }
            (Kind::Wall, None) => ([mouth, top], [0, 2]),
            _ => return Err(Error::Invalid),
        };
        let (points, _) = clip::polygon([ql, qh, vl, vh], &cuts, w)?;
        let rail0 = self
            .period
            .checked_mul(5)
            .and_then(|v| v.checked_add(self.band as i64))
            .and_then(|v| v.checked_mul(self.subdivisions as i64))
            .and_then(|v| v.checked_add(index as i64))
            .ok_or(Error::Limit)?;
        Ok(Polygon {
            points,
            stations,
            kind,
            rail0,
            row,
            rows,
        })
    }
    pub fn field(
        &self,
        thread: &AuthoredThread,
        polygon: &Polygon,
        w: &mut Work<'_>,
    ) -> Result<source::Field, Error> {
        use jet::Jet;
        let mut q = [f64::INFINITY, f64::NEG_INFINITY];
        let mut v = q;
        for p in &polygon.points {
            let a = p.q.enclosure(w)?.1;
            let b = p.v.enclosure(w)?.1;
            q = [q[0].min(a.lo), q[1].max(a.hi)];
            v = [v[0].min(b.lo), v[1].max(b.hi)];
        }
        let axis = |bounds: [f64; 2], i| {
            Jet::constant(I::point(bounds[0]))
                .add(Jet::variable(i).scale(I::point(bounds[1]) - I::point(bounds[0])))
        };
        let q = axis(q, 0);
        let v = axis(v, 1);
        let z = q.add(v.scale(self.lead.enclosure(w)?.1));
        let m = self.slope.enclosure(w)?.1;
        let (r, normal_z, beta) = match polygon.kind {
            Kind::Wall => {
                let r = q
                    .scale(m)
                    .add(Jet::constant(self.intercept.enclosure(w)?.1));
                (r, -m, (I::point(1.) / thread.angular_gain()) * m / r.value)
            }
            Kind::Cone => {
                let (a, k) = self.entry.as_ref().ok_or(Error::Invalid)?;
                let k = k.enclosure(w)?.1;
                (
                    Jet::constant(a.enclosure(w)?.1).sub(z.scale(k)),
                    k,
                    I::point(0.),
                )
            }
        };
        let tau = I::point(2.)
            * I::new(
                std::f64::consts::PI.next_down(),
                std::f64::consts::PI.next_up(),
            );
        let angle = v.scale(if thread.spec().left_hand { -tau } else { tau });
        Ok(source::Field {
            point: [r.mul(angle.cos()), r.mul(angle.sin()), z],
            normal: source::Normal::Radial {
                angle: angle.value,
                beta,
                z: normal_z,
            },
        })
    }
    /// Canonical source coordinates and identity must undergo the same exact
    /// periodic transfer. v=1,j is v=0,j+starts; globalQ advances by exact lead.
    pub fn owned(
        &self,
        p: &Point,
        polygon: &Polygon,
        starts: u32,
        w: &mut Work<'_>,
    ) -> Result<Footpoint, Error> {
        let zero = Q::integer(0, w)?;
        let one = Q::integer(1, w)?;
        let end = match &self.offset {
            Some(offset) => one.add(offset, w)?,
            None => one,
        };
        let wrap = p.v.cmp(&end, w)? == Ordering::Equal;
        let q = if wrap {
            p.q.add(&self.lead, w)?
        } else {
            p.q.clone()
        };
        let v = if wrap {
            self.offset.clone().unwrap_or(zero)
        } else {
            p.v.clone()
        };
        // Physical r,z are evaluated BEFORE periodic transfer; this exact source
        // identity changes Q and v but leaves both values unchanged.
        let z = p.q.add(&self.lead.mul(&p.v, w)?, w)?;
        let r = match polygon.kind {
            Kind::Wall => self.slope.mul(&p.q, w)?.add(&self.intercept, w)?,
            Kind::Cone => {
                let (a, k) = self.entry.as_ref().ok_or(Error::Invalid)?;
                a.sub(&k.mul(&z, w)?, w)?
            }
        };
        let physical = (4..6)
            .find(|&i| p.active & (1 << i) != 0)
            .map(|i| polygon.stations[i - 4]);
        let shift = if wrap {
            i64::from(starts) * 5 * self.subdivisions as i64
        } else {
            0
        };
        let rail = if p.active & 1 != 0 {
            Some(polygon.rail0 + shift)
        } else if p.active & 2 != 0 {
            Some(polygon.rail0 + 1 + shift)
        } else {
            None
        };
        let row = if p.active & 4 != 0 {
            Some(polygon.row)
        } else if p.active & 8 != 0 {
            Some(polygon.row + 1)
        } else {
            None
        }
        .map(|row| if row == polygon.rows { 0 } else { row });
        let key = match (rail, row, physical) {
            (Some(rail), _, Some(station)) => {
                Key([Constraint::Rail(rail), Constraint::Station(station, 0)])
            }
            (Some(rail), Some(row), None) => Key([Constraint::Rail(rail), Constraint::Row(row)]),
            (None, Some(row), Some(station)) => {
                let band = if station == 1 {
                    (self.period + i64::from(starts) * i64::from(wrap)) * 5 + self.band as i64
                } else {
                    0
                };
                Key([Constraint::Row(row), Constraint::Station(station, band)])
            }
            _ => return Err(Error::Unresolved),
        };
        let encoded = Q::encode(&[&q, &v, &r, &z], w)?;
        Ok((
            key,
            encoded,
            [r.enclosure(w)?.1, z.enclosure(w)?.1, v.enclosure(w)?.1],
        ))
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn source_chart_preserves_full_length_phase_and_shared_periodic_ownership() {
        let body = super::super::tests::source(false, 2, Some(0.4), true);
        let w = &mut Work {
            remaining: 16_000_000,
            used: 0,
            interrupt: None,
        };
        let a = Patch::new(&body, -2, 0, 1, w).unwrap();
        let b = Patch::new(&body, 0, 0, 1, w).unwrap();
        // The same actual source q-rail is represented on opposite sides of
        // the angular cut. Full source values, not nominal coordinates, agree.
        let pa = a.cell(0, 7, 8, Kind::Cone, w).unwrap();
        let pb = b.cell(0, 0, 8, Kind::Cone, w).unwrap();
        let aa = pa
            .points
            .iter()
            .find(|p| p.active & 0b1001 == 0b1001)
            .unwrap();
        let bb = pb
            .points
            .iter()
            .find(|p| p.active & 0b0101 == 0b0101)
            .unwrap();
        let (ka, ca, _) = a.owned(aa, &pa, 2, w).unwrap();
        let (kb, cb, _) = b.owned(bb, &pb, 2, w).unwrap();
        assert_eq!(ka, kb);
        assert_eq!(ca, cb);
        // q=p,j and q=0,j+1 are one exact shared phase endpoint too.
        let a = Patch::new(&body, 0, 4, 1, w).unwrap();
        let b = Patch::new(&body, 1, 0, 1, w).unwrap();
        assert_eq!(a.q1.cmp(&b.q0, w).unwrap(), Ordering::Equal);
        let pa = a.cell(0, 0, 8, Kind::Wall, w).unwrap();
        let pb = b.cell(0, 0, 8, Kind::Wall, w).unwrap();
        let aa = pa
            .points
            .iter()
            .find(|p| p.active & 0b0110 == 0b0110)
            .unwrap();
        let bb = pb
            .points
            .iter()
            .find(|p| p.active & 0b0101 == 0b0101)
            .unwrap();
        let (ka, ca, _) = a.owned(aa, &pa, 2, w).unwrap();
        let (kb, cb, _) = b.owned(bb, &pb, 2, w).unwrap();
        assert_eq!(ka, kb);
        assert_eq!(ca, cb);
        assert_eq!(body.full_length_mm(), 3.3);
    }
}
#[cfg(test)]
mod refusal_tests {
    use super::*;
    #[test]
    fn real_small_entry_source_keeps_unqualified_physical_contact_unresolved() {
        let body = super::super::tests::source(false, 1, Some(0.01), true);
        let w = &mut Work {
            remaining: 16_000_000,
            used: 0,
            interrupt: None,
        };
        let mut found = false;
        for band in [1, 3] {
            let patch = Patch::new(&body, -1, band, 1, w).unwrap();
            for row in 0..64 {
                let polygon = patch.cell(0, row, 64, Kind::Cone, w).unwrap();
                for point in &polygon.points {
                    if point.active & 0b110000 == 0b110000 && point.active & 0b1111 == 0 {
                        assert!(matches!(
                            patch.owned(point, &polygon, 1, w),
                            Err(Error::Unresolved)
                        ));
                        found = true;
                    }
                }
            }
        }
        assert!(found);
        let body = super::super::tests::source(false, 1, None, true);
        let patch = Patch::new(&body, 0, 0, 1, w).unwrap();
        assert!(matches!(
            patch.cell(0, 0, 0, Kind::Wall, w),
            Err(Error::Invalid)
        ));
        assert!(matches!(
            patch.cell(0, 0, 64, Kind::Cone, w),
            Err(Error::Invalid)
        ));
    }
}
#[cfg(test)]
mod cancellation_tests {
    use super::*;
    #[test]
    fn every_exact_cell_checkpoint_cancels_without_changing_source_or_reusing_failed_state() {
        for (entry, kind) in [
            (Some(0.4), Kind::Cone),
            (Some(0.4), Kind::Wall),
            (None, Kind::Wall),
        ] {
            let body = super::super::tests::source(false, 1, entry, true);
            let before = body;
            let calls = std::cell::Cell::new(0);
            let count = || {
                calls.set(calls.get() + 1);
                false
            };
            let mut initialization = Work {
                remaining: 16_000_000,
                used: 0,
                interrupt: None,
            };
            let patch = Patch::new(&body, 0, 0, 1, &mut initialization).unwrap();
            patch
                .cell(
                    0,
                    0,
                    64,
                    kind,
                    &mut Work {
                        remaining: 16_000_000,
                        used: 0,
                        interrupt: Some(&count),
                    },
                )
                .unwrap();
            let checkpoints = calls.get();
            for stop in 0..checkpoints {
                let patch = Patch::new(&body, 0, 0, 1, &mut initialization).unwrap();
                let calls = std::cell::Cell::new(0);
                let cancel = || {
                    calls.set(calls.get() + 1);
                    calls.get() > stop
                };
                assert!(matches!(
                    patch.cell(
                        0,
                        0,
                        64,
                        kind,
                        &mut Work {
                            remaining: 16_000_000,
                            used: 0,
                            interrupt: Some(&cancel)
                        }
                    ),
                    Err(Error::Cancelled)
                ));
            }
            assert_eq!(body, before);
        }
    }
}
