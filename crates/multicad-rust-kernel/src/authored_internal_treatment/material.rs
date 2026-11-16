//! Actual stock minus the authored MAX cavity. Stock extents are retained
//! rational source data; equality with the binary64 tool height is exact.
use super::*;
use crate::{Rational, cylinder::AxialCylinder};
use num_rational::BigRational;
/// Stock minus the cavity in a stock-relative chart. Local[0,0,0] is the
/// original cylinder's exact center/bottom, returned by stock_chart_origin.
/// Cavity coordinates and tessellation candidates use this local chart; they
/// do not contain an approximated WCS translation. stock() retains all original
/// rational coordinates. A caller must separately qualify any world placement.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Material {
    stock: AxialCylinder,
    cavity: Cavity,
    minimum_wall_mm: f64,
}
fn exact(value: Rational) -> BigRational {
    let (n, d) = value.ratio();
    BigRational::new(n.into(), d.into())
}
fn binary(value: f64) -> BigRational {
    // These callers have constructor-proved finite values, never arbitrary
    // imported floating operands. from_float is exact for every finite f64.
    BigRational::from_float(value).expect("validated finite source dimension")
}
fn integer(n: i64) -> I {
    let f = n as f64;
    if n.unsigned_abs() <= 1_u64 << 53 {
        I::point(f)
    } else {
        I::new(f.next_down(), f.next_up())
    }
}
fn ratio(value: Rational) -> I {
    let (n, d) = value.ratio();
    integer(n) / integer(d)
}
fn checkpoint(interrupt: Option<&dyn Fn() -> bool>) -> Result<(), Error> {
    if interrupt.is_some_and(|f| f()) {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
impl Material {
    pub fn new(stock: AxialCylinder, cavity: Cavity, minimum_wall_mm: f64) -> Result<Self, Error> {
        Self::new_with_interrupt(stock, cavity, minimum_wall_mm, None)
    }
    pub fn new_with_interrupt(
        stock: AxialCylinder,
        cavity: Cavity,
        minimum_wall_mm: f64,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, Error> {
        checkpoint(interrupt)?;
        if !minimum_wall_mm.is_finite() || minimum_wall_mm <= 0. || minimum_wall_mm > 1e7 {
            return Err(Error::Invalid);
        }
        // At most four bounded rational operands: stored integer ratios have
        // <=64bit operands, each binary64 <=1075bits. Fixed comparisons and one
        // subtraction stay below4400bits, with no unbounded expression graph.
        let end = cavity
            .exact_receiving_extent()
            .map_or_else(|| binary(cavity.full_length_mm()), exact);
        if exact(stock.height()) != end {
            return Err(Error::Unresolved);
        }
        checkpoint(interrupt)?;
        let radius = exact(stock.circle().radius());
        if radius > BigRational::from_integer(10_000_000.into()) {
            return Err(Error::Limit);
        }
        if radius - binary(cavity.maximum_radius_mm()) < binary(minimum_wall_mm) {
            return Err(Error::Invalid);
        }
        checkpoint(interrupt)?;
        Ok(Self {
            stock,
            cavity,
            minimum_wall_mm,
        })
    }
    pub fn stock(&self) -> &AxialCylinder {
        &self.stock
    }
    /// Exact origin of the local stock chart in the retained, unplaced
    /// AxialCylinder coordinate system. This is not an active framed placement.
    pub fn stock_chart_origin(&self) -> [Rational; 3] {
        let [x, y] = self.stock.circle().center();
        [x, y, self.stock.bottom()]
    }
    pub fn cavity(&self) -> &Cavity {
        &self.cavity
    }
    pub fn minimum_wall_mm(&self) -> f64 {
        self.minimum_wall_mm
    }
    pub(super) fn radius_interval(&self) -> I {
        ratio(self.stock.circle().radius())
    }
    pub fn volume_bounds(&self, options: QueryOptions) -> Result<Enclosure, Error> {
        self.volume_bounds_with_interrupt(options, None)
    }
    pub fn volume_bounds_with_interrupt(
        &self,
        options: QueryOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Enclosure, Error> {
        checkpoint(interrupt)?;
        if !options.tolerance_mm3.is_finite()
            || options.tolerance_mm3 <= 0.
            || options.maximum_cells == 0
            || options.maximum_cells > 1_000_000
        {
            return Err(Error::Invalid);
        }
        let pi = I::new(
            std::f64::consts::PI.next_down(),
            std::f64::consts::PI.next_up(),
        );
        let stock = pi * ratio(self.stock.circle().radius()).square() * ratio(self.stock.height());
        let remaining = (I::point(options.tolerance_mm3) - I::point(stock.width())).lo;
        if remaining <= 0. {
            return Err(Error::Unresolved);
        }
        let hole = self.cavity.removed_volume_with_interrupt(
            QueryOptions {
                tolerance_mm3: remaining * 0.5,
                ..options
            },
            interrupt,
        )?;
        let material = stock - I::new(hole.lower, hole.upper);
        if !material.valid() || material.lo <= 0. || material.width() > options.tolerance_mm3 {
            return Err(Error::Unresolved);
        }
        checkpoint(interrupt)?;
        Ok(Enclosure {
            lower: material.lo,
            upper: material.hi,
        })
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        Circle,
        authored_flat_thread::{AuthoredThread, SourceOptions},
        helical_thread::Spec,
    };
    fn cavity() -> Cavity {
        let thread = AuthoredThread::new(
            Spec {
                major_radius_mm: 1.5,
                pitch_mm: 0.5,
                length_mm: 4.,
                depth_mm: 0.25,
                crest_width_mm: 0.125,
                lower_flank_deg: 30.,
                upper_flank_deg: 30.,
                starts: 1,
                left_hand: false,
                root_radius_mm: None,
            },
            SourceOptions {
                position_limit_mm: 1e-7,
                normal_angle_limit_rad: 1e-7,
            },
        )
        .unwrap();
        Cavity::new(
            thread,
            Intent {
                entry_depth_mm: Some(0.5),
                terminal: Some(Terminal {
                    width_mm: 0.5,
                    radial_depth_mm: 0.125,
                }),
            },
        )
        .unwrap()
    }
    fn stock(height: Rational, radius: Rational) -> AxialCylinder {
        AxialCylinder::new(
            7,
            Circle::new([Rational::new(0, 1).unwrap(); 2], radius).unwrap(),
            Rational::new(0, 1).unwrap(),
            height,
        )
        .unwrap()
    }
    #[test]
    fn actual_rational_stock_is_retained_with_independent_net_volume_and_minimum_wall() {
        let stock = stock(Rational::new(4, 1).unwrap(), Rational::new(5, 2).unwrap());
        let cavity = cavity();
        let material = Material::new(stock, cavity, 0.5).unwrap();
        assert_eq!(*material.stock(), stock);
        assert_eq!(*material.cavity(), cavity);
        assert_eq!(material.minimum_wall_mm(), 0.5);
        let hole = cavity.removed_volume_mm3(QueryOptions::default()).unwrap();
        let actual = material.volume_bounds(QueryOptions::default()).unwrap();
        let analytic = std::f64::consts::PI * 6.25 * 4. - hole.lower.midpoint(hole.upper);
        assert!(actual.lower <= analytic && analytic <= actual.upper);
        assert!(actual.lower > 0. && actual.upper - actual.lower <= 1e-7);
        assert!(Material::new(stock, cavity, 1.).is_err());
    }
    #[test]
    fn non_dyadic_actual_height_is_not_silently_identified_with_rounded_tool_extent() {
        let mut spec = cavity().thread.spec();
        spec.length_mm = 3.3;
        let thread = AuthoredThread::new(
            spec,
            SourceOptions {
                position_limit_mm: 1e-7,
                normal_angle_limit_rad: 1e-7,
            },
        )
        .unwrap();
        let source = Cavity::new(
            thread,
            Intent {
                entry_depth_mm: Some(0.5),
                terminal: Some(Terminal {
                    width_mm: 0.5,
                    radial_depth_mm: 0.125,
                }),
            },
        )
        .unwrap();
        let real = stock(Rational::new(33, 10).unwrap(), Rational::new(3, 1).unwrap());
        assert_eq!(real.height().as_f64(), source.full_length_mm());
        assert_eq!(Material::new(real, source, 0.5), Err(Error::Unresolved));
    }
    #[test]
    fn exact_stock_wall_limits_and_real_volume_precision_floor_refuse_atomically() {
        let source = cavity();
        let shape = stock(Rational::new(4, 1).unwrap(), Rational::new(5, 2).unwrap());
        for wall in [0., -1., f64::NAN, f64::INFINITY, 10_000_001.] {
            assert_eq!(Material::new(shape, source, wall), Err(Error::Invalid));
        }
        let huge = stock(
            Rational::new(4, 1).unwrap(),
            Rational::new(10_000_001, 1).unwrap(),
        );
        assert_eq!(Material::new(huge, source, 0.5), Err(Error::Limit));
        let large = Material::new(
            stock(
                Rational::new(4, 1).unwrap(),
                Rational::new(1_000_000, 1).unwrap(),
            ),
            source,
            0.5,
        )
        .unwrap();
        assert_eq!(
            large.volume_bounds(QueryOptions::default()),
            Err(Error::Unresolved)
        );
        let pi = I::new(
            std::f64::consts::PI.next_down(),
            std::f64::consts::PI.next_up(),
        );
        let exact_stock =
            pi * ratio(large.stock.circle().radius()).square() * ratio(large.stock.height());
        // The stock's own enclosure fits, and the hole integral resolves, but
        // final subtraction at this magnitude cannot honor the additional
        // sub-ULP width. This is a genuine final certificate precision refusal.
        assert_eq!(
            large.volume_bounds(QueryOptions {
                tolerance_mm3: exact_stock.width() + 1e-6,
                ..QueryOptions::default()
            }),
            Err(Error::Unresolved)
        );
        let relaxed = large
            .volume_bounds(QueryOptions {
                tolerance_mm3: 1.,
                ..QueryOptions::default()
            })
            .unwrap();
        assert!(relaxed.lower > 0. && relaxed.upper - relaxed.lower <= 1.);
        let unusual = stock(
            Rational::new(4, 1).unwrap(),
            Rational::new(900_000_000_000_000_001, 300_000_000_000_000_000).unwrap(),
        );
        let bounds = Material::new(unusual, source, 0.5)
            .unwrap()
            .volume_bounds(QueryOptions::default())
            .unwrap();
        assert!(bounds.lower > 0. && bounds.upper - bounds.lower <= 1e-7);
        assert_eq!(large.cavity, source);
    }
    #[test]
    fn cancelled_and_limited_queries_never_change_the_retained_source() {
        let material = Material::new(
            stock(Rational::new(4, 1).unwrap(), Rational::new(5, 2).unwrap()),
            cavity(),
            0.5,
        )
        .unwrap();
        let original = material;
        for stop in 0..3 {
            let calls = std::cell::Cell::new(0);
            let cancel = || {
                calls.set(calls.get() + 1);
                calls.get() > stop
            };
            assert_eq!(
                Material::new_with_interrupt(material.stock, material.cavity, 0.5, Some(&cancel)),
                Err(Error::Cancelled)
            );
        }
        for options in [
            QueryOptions {
                tolerance_mm3: 0.,
                ..QueryOptions::default()
            },
            QueryOptions {
                maximum_cells: 0,
                ..QueryOptions::default()
            },
            QueryOptions {
                maximum_cells: 1_000_001,
                ..QueryOptions::default()
            },
        ] {
            assert_eq!(material.volume_bounds(options), Err(Error::Invalid));
        }
        assert_eq!(
            material.volume_bounds(QueryOptions {
                maximum_cells: 1,
                ..QueryOptions::default()
            }),
            Err(Error::Limit)
        );
        let calls = std::cell::Cell::new(0);
        let count = || {
            calls.set(calls.get() + 1);
            false
        };
        let expected = material
            .volume_bounds_with_interrupt(QueryOptions::default(), Some(&count))
            .unwrap();
        let checkpoints = calls.get();
        for stop in [0, 1, checkpoints - 1] {
            let calls = std::cell::Cell::new(0);
            let cancel = || {
                calls.set(calls.get() + 1);
                calls.get() > stop
            };
            assert_eq!(
                material.volume_bounds_with_interrupt(QueryOptions::default(), Some(&cancel)),
                Err(Error::Cancelled)
            );
        }
        assert_eq!(material, original);
        assert_eq!(
            material.volume_bounds(QueryOptions::default()).unwrap(),
            expected
        );
    }
}
