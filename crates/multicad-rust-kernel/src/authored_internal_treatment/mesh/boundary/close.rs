//! Actual source cap/tail incidence, with the SAME exact angular footpoints.
use super::*;
pub(super) fn scalar(bytes: &[u8], index: usize, w: &mut Work<'_>) -> Result<Q, Error> {
    let mut offset = 0usize;
    for _ in 0..index {
        w.step()?;
        let header = bytes.get(offset..offset + 5).ok_or(Error::Unresolved)?;
        let n = usize::from(u16::from_le_bytes([header[1], header[2]]));
        let d = usize::from(u16::from_le_bytes([header[3], header[4]]));
        offset += 5 + n + d;
    }
    w.step()?;
    let header = bytes.get(offset..offset + 5).ok_or(Error::Unresolved)?;
    let n = usize::from(u16::from_le_bytes([header[1], header[2]]));
    let d = usize::from(u16::from_le_bytes([header[3], header[4]]));
    Q::decode(
        bytes
            .get(offset..offset + 5 + n + d)
            .ok_or(Error::Unresolved)?,
        w,
    )
}

impl Graph {
    pub(super) fn attach(
        &mut self,
        index: usize,
        face: u32,
        forward: bool,
        w: &mut Work<'_>,
    ) -> Result<(), Error> {
        w.step()?;
        let curve = self.curves.get_mut(index).ok_or(Error::Unresolved)?;
        if curve.uses[1].is_some() || curve.uses[0].unwrap().forward == forward {
            return Err(Error::Unresolved);
        }
        curve.uses[1] = Some(Use { face, forward });
        Ok(())
    }
    fn projected(
        &mut self,
        source: u32,
        end: u8,
        v: &Q,
        r: &Q,
        z: &Q,
        w: &mut Work<'_>,
    ) -> Result<u32, Error> {
        let v = if v.cmp(&Q::integer(1, w)?, w)? == Ordering::Equal {
            Q::integer(0, w)?
        } else {
            v.clone()
        };
        let token = Q::encode(&[&v, r, z], w)?;
        let key = Key([
            Constraint::Station(4 + end, i64::from(source)),
            Constraint::Row(0),
        ]);
        let rzv = [r.enclosure(w)?.1, z.enclosure(w)?.1, v.enclosure(w)?.1];
        self.register_point(key, token, rzv, Origin::Tail { source, end }, w)
    }
    pub(super) fn closed(mut self, w: &mut Work<'_>) -> Result<Self, Error> {
        let mouth = self.rim(0, w)?;
        let mouth_bytes = mouth.capacity() * size_of::<loops::Coedge>();
        self.temporary += mouth_bytes;
        for edge in mouth {
            self.attach(edge.curve, 0, !edge.forward, w)?;
        }
        self.temporary -= mouth_bytes;
        let top = self.rim(1, w)?;
        let top_bytes = top.capacity() * size_of::<loops::Coedge>();
        self.temporary += top_bytes;
        if let Some(tail) = self.source.tail {
            let base =
                2 + 4 * self.source.thread.spec().starts + u32::from(self.source.entry.is_some());
            let r = Q::float(tail.radius, w)?;
            let z0 = Q::float(tail.station, w)?;
            let z1 = Q::float(self.source.full_length_mm(), w)?;
            for edge in top {
                let curve = &self.curves[edge.curve];
                let original = curve.endpoints;
                let v0 = scalar(&curve.exact, 2, w)?;
                let dv = scalar(&curve.exact, 5, w)?;
                let v1 = v0.add(&dv, w)?;
                let low = [
                    self.projected(original[0], 0, &v0, &r, &z0, w)?,
                    self.projected(original[1], 0, &v1, &r, &z0, w)?,
                ];
                let high = [
                    self.projected(original[0], 1, &v0, &r, &z1, w)?,
                    self.projected(original[1], 1, &v1, &r, &z1, w)?,
                ];
                self.attach(edge.curve, base, !edge.forward, w)?;
                let circle0 = self.register(
                    low,
                    [
                        [r.clone(), z0.clone(), v0.clone()],
                        [r.clone(), z0.clone(), v1.clone()],
                    ],
                    Some(4),
                    Use {
                        face: base,
                        forward: edge.forward,
                    },
                    w,
                )? as usize;
                self.attach(
                    circle0,
                    base + 1,
                    !self.curves[circle0].uses[0].unwrap().forward,
                    w,
                )?;
                let circle1 = self.register(
                    high,
                    [[r.clone(), z1.clone(), v0], [r.clone(), z1.clone(), v1]],
                    Some(5),
                    Use {
                        face: base + 1,
                        forward: edge.forward,
                    },
                    w,
                )? as usize;
                self.attach(
                    circle1,
                    1,
                    !self.curves[circle1].uses[0].unwrap().forward,
                    w,
                )?;
            }
        } else {
            for edge in top {
                self.attach(edge.curve, 1, !edge.forward, w)?;
            }
        }
        self.temporary -= top_bytes;
        self.check(0)?;
        w.step()?;
        Ok(self)
    }
}
