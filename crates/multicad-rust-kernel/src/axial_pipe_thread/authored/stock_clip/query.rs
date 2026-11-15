//! Whole-box queries of the original finite authored radial material.
use super::*;
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StockClipClassification {
    Inside,
    Outside,
    Undetermined,
}
#[derive(Clone, Copy, Debug)]
pub struct StockClipQuery {
    classification: StockClipClassification,
    work: usize,
    bytes: usize,
}
impl StockClipQuery {
    pub fn classification(&self) -> StockClipClassification {
        self.classification
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_budget_bytes(&self) -> usize {
        self.bytes
    }
}
fn angle(x: I, y: I) -> I {
    if x.lo > 0. {
        (y / x).atan()
    } else if x.hi < 0. && y.lo >= 0. {
        (y / x).atan() + pi()
    } else if x.hi < 0. && y.hi <= 0. {
        (y / x).atan() - pi()
    } else if y.lo > 0. {
        pi() / I::point(2.) - (x / y).atan()
    } else if y.hi < 0. {
        -pi() / I::point(2.) - (x / y).atan()
    } else {
        I::new(-pi().hi, pi().hi)
    }
}
impl StockClipPlan {
    /// Encloses the complete original material; this bound alone is never
    /// used as proof that its bounding-box corners contain material.
    pub fn bounds_mm(&self) -> [[f64; 2]; 3] {
        let r = self.stock_radius_limit_mm();
        [[-r, r], [-r, r], [0., self.source.spec.length_mm]]
    }
    /// Strict whole-box radial/finite-height classification. Boundary or
    /// unresolved phase boxes remain Undetermined, never silently material.
    pub fn classify_box(
        &self,
        box_mm: [[f64; 2]; 3],
        options: StockClipOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<StockClipQuery, Error> {
        let b = box_mm.map(interval);
        let [x, y, z] = [b[0]?, b[1]?, b[2]?];
        if [x, y, z].iter().any(|v| v.magnitude() > 1e7) {
            return Err(Error::Limit);
        }
        let retained = size_of::<Self>() + size_of::<StockClipQuery>();
        let mut work = Meter::new(options, retained, interrupt)?;
        work.charge(self.work)?;
        work.charge(256)?;
        let length = self.source.spec.length_mm;
        let radius = (x.square() + y.square()).sqrt();
        let classification =
            if z.hi < 0. || z.lo > length || radius.lo > self.stock_radius_limit_mm() {
                StockClipClassification::Outside
            } else {
                let height = I::new(z.lo.max(0.), z.hi.min(length));
                let k = I::point(self.source.spec.radial_taper);
                let root = self.root_radius + k * height;
                let outer = if radius.hi < root.lo {
                    root
                } else {
                    let p = I::point(self.source.spec.pitch_mm);
                    let phase = height - p * angle(x, y) / (I::point(2.) * pi());
                    let bands = [
                        (self.q0, -self.half),
                        (-self.half, self.half),
                        (self.half, self.q3),
                        (self.q3, self.q0 + p),
                    ];
                    let mut hull = I::new(f64::INFINITY, f64::NEG_INFINITY);
                    // The finite authored domain L/p+2<=128, angle in [-pi,pi],
                    // q0>−p and q3<p earn these fixed whole-domain lift bounds.
                    for turn in -2..=130 {
                        work.charge(64)?;
                        let local = phase - I::point(turn as f64) * p;
                        for (band, (lo, hi)) in bands.into_iter().enumerate() {
                            if local.hi < lo.lo || local.lo > hi.hi {
                                continue;
                            }
                            let q = I::new(local.lo.max(lo.lo), local.hi.min(hi.hi));
                            let r = match band {
                                0 => root + (q - self.q0) / self.slopes[0],
                                1 => I::point(self.source.spec.major_radius_mm) + k * height,
                                2 => root + (self.q3 - q) / self.slopes[1],
                                _ => root,
                            };
                            hull = I::new(hull.lo.min(r.lo), hull.hi.max(r.hi));
                        }
                    }
                    // The clipped finite-height path above and sealed
                    // one-period band partition earn at least one lift in
                    // -2..=130. All operands are bounded (<1e10), so this
                    // hull is nonempty/finite; no arbitrary-box assumption.
                    min(hull, I::point(self.stock_radius_limit_mm()))
                };
                if radius.lo > outer.hi {
                    StockClipClassification::Outside
                } else if z.lo > 0. && z.hi < length && radius.hi < outer.lo {
                    StockClipClassification::Inside
                } else {
                    StockClipClassification::Undetermined
                }
            };
        work.charge(1)?;
        Ok(StockClipQuery {
            classification,
            work: work.used,
            bytes: retained + SCRATCH,
        })
    }
}
#[cfg(test)]
mod tests;
