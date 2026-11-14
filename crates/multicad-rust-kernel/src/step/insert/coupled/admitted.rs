//! Build-time capabilities earned from the SAME retained source decomposition.
use super::*;
pub(super) struct Registry<'d, 's> {
    domain: &'d DomainPlan<'s>,
    keys: Vec<RowKey>,
}
impl<'d, 's> Registry<'d, 's> {
    pub(super) fn prepare(
        domain: &'d DomainPlan<'s>,
        memory: &mut Memory,
        w: &mut Work<'_>,
    ) -> Result<Self, Error> {
        let mut keys: Vec<RowKey> = memory.allocate(domain.slabs().len() * 2, w)?;
        for slab in domain.slabs() {
            for curve in slab.rows {
                let key = RowKey {
                    curve,
                    interval: slab.stations[0],
                };
                let mut low = 0;
                let mut high = keys.len();
                let mut found = false;
                while low < high {
                    w.charge(1)?;
                    let mid = low + (high - low) / 2;
                    match keys[mid].cmp(&key) {
                        std::cmp::Ordering::Less => low = mid + 1,
                        std::cmp::Ordering::Greater => high = mid,
                        std::cmp::Ordering::Equal => {
                            found = true;
                            break;
                        }
                    }
                }
                if !found {
                    w.charge(1 + keys.len() - low)?;
                    keys.insert(low, key);
                }
            }
        }
        Ok(Self { domain, keys })
    }
    pub(super) fn len(&self) -> usize {
        self.keys.len()
    }
    pub(super) fn iter(&self) -> impl Iterator<Item = DomainRow<'d, 's>> + '_ {
        self.keys.iter().map(|&key| DomainRow {
            domain: self.domain,
            key,
        })
    }
    pub(super) fn release(&self, memory: &mut Memory) {
        memory.release(&self.keys);
    }
}
pub(super) struct DomainRow<'d, 's> {
    domain: &'d DomainPlan<'s>,
    key: RowKey,
}
impl<'d, 's> DomainRow<'d, 's> {
    pub(super) fn key(&self) -> RowKey {
        self.key
    }
    fn range(&self, w: &mut Work<'_>) -> Result<[N; 2], Error> {
        Ok([
            self.domain.stations(self.key.curve.owner())[self.key.interval as usize].number(w)?,
            self.domain.stations(self.key.curve.owner())[self.key.interval as usize + 1]
                .number(w)?,
        ])
    }
    pub(super) fn original<'c>(
        &'c self,
        curve: &'c Curve<'s>,
        w: &mut Work<'_>,
    ) -> Result<OriginalRow<'c, 's>, Error> {
        w.charge(1)?;
        if !std::ptr::eq(self.domain.boundary(), curve.source) || self.key.curve != curve.key {
            return Err(Error::Invalid);
        }
        let v = self.range(w)?;
        let t = [
            v[0].sub(&curve.affine[0][2], w)?
                .div(&curve.affine[1][2], w)?,
            v[1].sub(&curve.affine[0][2], w)?
                .div(&curve.affine[1][2], w)?,
        ];
        // Registry is built only from WHOLE source-cell slab incidence. Curve
        // retains the exact same producer law/owner/key. Thus the closed slab
        // endpoints map into that original curve's [0,1], without a float range
        // check, clamp, fitted law or independent reconstruction.
        Ok(OriginalRow {
            curve,
            t,
            scope: Scope {
                domain: self.domain,
                side: self.key.curve.owner(),
                interval: self.key.interval,
            },
        })
    }
    pub(super) fn angular<'u>(
        &self,
        angular: &super::super::angular::AngularCarrier,
        reduced: &'u reduced::Reduced,
        w: &mut Work<'_>,
    ) -> Result<AngularSpan<'u, 'd, 's>, Error> {
        let v = self.range(w)?;
        let piece = angular
            .knots
            .windows(2)
            .position(|k| {
                v[0].0 >= BigRational::from_float(k[0]).unwrap()
                    && v[1].0 <= BigRational::from_float(k[1]).unwrap()
            })
            .ok_or(Error::Unresolved)?;
        let a = N(BigRational::from_float(angular.knots[piece]).unwrap());
        let h = N(BigRational::from_float(angular.knots[piece + 1]).unwrap()).sub(&a, w)?;
        let begin = v[0].sub(&a, w)?.div(&h, w)?;
        let end = v[1].sub(&a, w)?.div(&h, w)?;
        // Sorted Domain stations earn v0<v1; the actual containing-piece lookup
        // earns k0<=v0<v1<=k1. Exact affine division therefore earns0<=a<b<=1.
        Ok(AngularSpan {
            scope: Scope {
                domain: self.domain,
                side: self.key.curve.owner(),
                interval: self.key.interval,
            },
            parent: &reduced.poles[piece],
            errors: reduced.errors,
            begin,
            end,
            width: v[1].sub(&v[0], w)?,
        })
    }
}
pub(super) struct Curve<'s> {
    source: &'s Boundary,
    key: BoundaryCurveKey,
    affine: [[N; 3]; 2],
}
impl<'s> Curve<'s> {
    pub(super) fn key(&self) -> BoundaryCurveKey {
        self.key
    }
    pub(super) fn load(
        row: &DomainRow<'_, 's>,
        memory: &mut Memory,
        w: &mut Work<'_>,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, Error> {
        let source = row.domain.boundary();
        let component = source.source_boundary(row.key.curve.owner());
        memory.add(64 * 1024)?;
        let data = component
            .exact_curve_coefficients(
                row.key.curve.source_index(),
                w.remaining(),
                memory.used - component.retained_bytes() - SCRATCH,
                interrupt,
            )
            .map_err(domains::failure)?;
        w.charge(data.work_used())?;
        memory.used -= 64 * 1024;
        memory.add(data.owned_bytes())?;
        let mut affine: [[N; 3]; 2] =
            std::array::from_fn(|_| std::array::from_fn(|_| N::integer(0)));
        for (i, line) in affine.iter_mut().enumerate() {
            for (axis, c) in line.iter_mut().enumerate() {
                *c = N::view(data.coordinate(i, axis).unwrap(), w)?;
            }
        }
        memory.used -= data.owned_bytes();
        Ok(Self {
            source,
            key: row.key.curve,
            affine,
        })
    }
}
#[derive(Clone, Copy)]
pub(super) struct Scope<'d, 's> {
    domain: &'d DomainPlan<'s>,
    side: Side,
    interval: u32,
}
impl Scope<'_, '_> {
    pub(super) fn same(&self, other: &Scope<'_, '_>) -> bool {
        std::ptr::eq(self.domain, other.domain)
            && self.side == other.side
            && self.interval == other.interval
    }
}
pub(super) struct OriginalRow<'c, 's> {
    curve: &'c Curve<'s>,
    t: [N; 2],
    pub(super) scope: Scope<'c, 's>,
}
impl OriginalRow<'_, '_> {
    pub(super) fn affine(&self) -> &[[N; 3]; 2] {
        &self.curve.affine
    }
    pub(super) fn parameters(&self) -> &[N; 2] {
        &self.t
    }
}
pub(super) struct AngularSpan<'u, 'd, 's> {
    scope: Scope<'d, 's>,
    parent: &'u [[number::Grid; 2]; 11],
    errors: [f64; 2],
    begin: N,
    end: N,
    width: N,
}
impl<'u, 'd, 's> AngularSpan<'u, 'd, 's> {
    pub(super) fn scope(&self) -> Scope<'d, 's> {
        self.scope
    }
    pub(super) fn parent(&self) -> &'u [[number::Grid; 2]; 11] {
        self.parent
    }
    pub(super) fn errors(&self) -> [f64; 2] {
        self.errors
    }
    pub(super) fn begin(&self) -> &N {
        &self.begin
    }
    pub(super) fn end(&self) -> &N {
        &self.end
    }
    pub(super) fn width(&self) -> &N {
        &self.width
    }
}
