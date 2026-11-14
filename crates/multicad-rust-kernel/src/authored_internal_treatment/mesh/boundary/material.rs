//! Original exact stock circles plus reversed void boundary; no mesh surgery.
use super::*;
impl Graph {
    pub(super) fn material(mut self, material: Material, w: &mut Work<'_>) -> Result<Self, Error> {
        w.step()?;
        if self.source != Source::treated(*material.cavity()) || self.material.is_some() {
            return Err(Error::Invalid);
        }
        if self.curves.iter().any(|curve| curve.uses[1].is_none()) {
            return Err(Error::Unresolved);
        }
        // Original tool cap disks become the INNER wires of actual stock
        // annuli. Every receiving surface reverses its occupied-material side.
        for curve in &mut self.curves {
            w.step()?;
            for side in curve.uses.iter_mut().flatten() {
                side.forward = !side.forward;
            }
        }
        let outer = 2
            + 4 * self.source.thread.spec().starts
            + u32::from(self.source.entry.is_some())
            + 2 * u32::from(self.source.tail.is_some());
        let (n, d) = material.stock().circle().radius().ratio();
        let r = Q::integer(n, w)?.div(&Q::integer(d, w)?, w)?;
        let zero = Q::integer(0, w)?;
        let one = Q::integer(1, w)?;
        for end in 0..2u8 {
            let z = Q::float(
                if end == 0 {
                    0.
                } else {
                    self.source.full_length_mm()
                },
                w,
            )?;
            let token = Q::encode(&[&zero, &r, &z], w)?;
            let key = Key([Constraint::Station(6 + end, 0), Constraint::Row(0)]);
            let rzv = [r.enclosure(w)?.1, z.enclosure(w)?.1, I::point(0.)];
            let vertex = self.register_point(key, token, rzv, Origin::Stock { end }, w)?;
            let forward = (end == 1) != self.source.thread.spec().left_hand;
            let circle = self.register(
                [vertex; 2],
                [
                    [r.clone(), z.clone(), zero.clone()],
                    [r.clone(), z, one.clone()],
                ],
                Some(6 + end),
                Use {
                    face: u32::from(end),
                    forward,
                },
                w,
            )? as usize;
            self.attach(
                circle,
                outer,
                !self.curves[circle].uses[0].unwrap().forward,
                w,
            )?;
        }
        self.material = Some(material);
        self.check(0)?;
        w.step()?;
        Ok(self)
    }
}
