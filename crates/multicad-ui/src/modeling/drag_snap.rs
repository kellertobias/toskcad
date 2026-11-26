//! One drag-snap policy for viewport dimension handles: primitive and panel
//! corners, the primitive move gizmo and the Extrude distance arrow
//! (CAD-402, CAD-403).
//!
//! - The snap step is the floor-grid spacing the viewport currently displays
//!   (`viewport::grid_step`), so the visible grid is the snap grid at every
//!   zoom level.
//! - Holding Shift during a drag releases the grid, like Move and Rotate
//!   (CAD-390). A free drag still lands on [`FREE_RESOLUTION_MM`] so a dragged
//!   value is a short decimal the exact kernel can admit (CAD-401).
//! - Snapped values are bit-exact: `n × step` is formed from the integer `n`
//!   and the decimal digits of the step, never by repeated floating-point
//!   addition, so a 10 mm grid gives exactly `90.0`, never `89.99999999999999`.
//! - Typed values bypass this policy entirely.
use crate::egui;
use multicad_precision::Rational;

/// Resolution of a free (Shift) drag in millimetres.
pub const FREE_RESOLUTION_MM: f64 = 0.001;
/// Largest length or coordinate magnitude a drag can produce. It keeps the
/// far primitive corner well inside the geometry policy's model bounds.
pub const MAX_DRAG_MM: f64 = 1e6;
/// Decimal places accepted for a grid step.
const MAX_PLACES: u32 = 9;

/// A positive decimal step `units / 10^places`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Step {
    units: i64,
    places: u32,
}
impl Step {
    pub const FREE: Self = Self {
        units: 1,
        places: 3,
    };
    /// The exact decimal a finite positive step denotes, if it has at most
    /// nine decimal places (grid steps are `base × 10^k`).
    pub fn new(step_mm: f64) -> Option<Self> {
        if !step_mm.is_finite() || step_mm <= 0. || step_mm > MAX_DRAG_MM {
            return None;
        }
        (0..=MAX_PLACES).find_map(|places| {
            let scaled = step_mm * 10_f64.powi(places as i32);
            let units = scaled.round();
            (units >= 1. && (scaled - units).abs() <= scaled * 1e-12).then_some(Self {
                units: units as i64,
                places,
            })
        })
    }
    fn divisor(self) -> f64 {
        10_f64.powi(self.places as i32)
    }
    pub fn mm(self) -> f64 {
        self.value(1)
    }
    /// `index × step`, rounded once from the exact decimal.
    pub fn value(self, index: i64) -> f64 {
        // Both operands are exact binary64 integers below 2^53, so the
        // division rounds the exact decimal once, as parsing it would.
        (index as f64 * self.units as f64) / self.divisor()
    }
    /// Nearest grid index of a finite value inside the drag range.
    pub fn index(self, value: f64) -> Option<i64> {
        if !value.is_finite() {
            return None;
        }
        let clamped = value.clamp(-MAX_DRAG_MM, MAX_DRAG_MM);
        Some((clamped * self.divisor() / self.units as f64).round() as i64)
    }
    /// The grid index of a value that lies exactly on this grid.
    pub fn exact_index(self, value: f64) -> Option<i64> {
        self.index(value)
            .filter(|index| self.value(*index) == value)
    }
    /// Coarser decimal steps for free drags that would exceed exact admission.
    fn coarser(self) -> Option<Self> {
        if self.places > 0 {
            Some(Self {
                places: self.places - 1,
                ..self
            })
        } else {
            (self.units <= MAX_DRAG_MM as i64).then_some(Self {
                units: self.units * 10,
                ..self
            })
        }
    }
}

/// Snap state of one drag frame.
#[derive(Clone, Copy, Debug, PartialEq, Default)]
pub struct DragSnap {
    /// Spacing of the displayed grid, when it is a usable decimal.
    grid: Option<Step>,
    /// The documented modifier (Shift) is held.
    free: bool,
}
impl DragSnap {
    /// Free drag at the fine resolution, used where no grid is known.
    pub const FREE: Self = Self {
        grid: None,
        free: false,
    };
    pub fn new(grid_mm: Option<f64>, free: bool) -> Self {
        Self {
            grid: grid_mm.and_then(Step::new),
            free,
        }
    }
    /// The same grid with the modifier state of the current frame.
    pub fn with_free(self, free: bool) -> Self {
        Self { free, ..self }
    }
    /// The modifier that releases grid snapping during a drag.
    pub fn free_modifier(modifiers: egui::Modifiers) -> bool {
        modifiers.shift
    }
    /// The grid step applied now, `None` while free.
    pub fn grid_mm(self) -> Option<f64> {
        self.active().map(Step::mm)
    }
    fn active(self) -> Option<Step> {
        self.grid.filter(|_| !self.free)
    }
    /// The quantum every dragged value lands on.
    pub fn step(self) -> Step {
        self.active().unwrap_or(Step::FREE)
    }
    /// Nearest snapped value of a dragged coordinate or signed distance.
    #[cfg(test)]
    pub fn coordinate(self, value: f64) -> Option<f64> {
        let step = self.step();
        step.index(value).map(|index| step.value(index))
    }
    /// Nearest snapped positive length, never below one step.
    pub fn length(self, value: f64) -> Option<f64> {
        let step = self.step();
        step.index(value).map(|index| step.value(index.max(1)))
    }
    /// Nearest snapped nonzero signed distance. A value that would snap to
    /// zero takes one step in the direction of `value` (or of `previous` when
    /// `value` is exactly zero), so a drag never produces an empty extrusion.
    pub fn signed(self, value: f64, previous: f64) -> Option<f64> {
        let step = self.step();
        let index = step.index(value)?;
        let index = if index != 0 {
            index
        } else if value > 0. || value == 0. && previous >= 0. {
            1
        } else {
            -1
        };
        Some(step.value(index))
    }
    /// Feedback naming the active spacing and the modifier.
    pub fn label(self) -> String {
        match (self.grid, self.free) {
            (Some(step), false) => format!(
                "grid {} · Shift for free",
                crate::preferences::format_mm(step.mm())
            ),
            (Some(_), true) => format!(
                "free (Shift) · {}",
                crate::preferences::format_mm(FREE_RESOLUTION_MM)
            ),
            (None, _) => format!(
                "free · {}",
                crate::preferences::format_mm(FREE_RESOLUTION_MM)
            ),
        }
    }
    /// The step to use for a box so that its exact volume stays inside the
    /// kernel's admission range. A grid step is kept as is (a box beyond it
    /// is refused by the caller); a free step coarsens decimal by decimal.
    pub(crate) fn box_step(self, admissible: impl Fn(Step) -> bool) -> Option<Step> {
        let mut step = self.step();
        if self.active().is_some() {
            return admissible(step).then_some(step);
        }
        loop {
            if admissible(step) {
                return Some(step);
            }
            step = step.coarser()?;
        }
    }
}

/// Exact admission of a rectangular box as the Rust kernel evaluates it:
/// every extent and the base area and volume are bounded rationals.
pub(crate) fn box_admissible(size_mm: [f64; 3]) -> bool {
    let exact = |value: f64| Rational::parse_decimal(&value.to_string()).ok();
    let Some([l, w, h]) = size_mm
        .iter()
        .map(|value| exact(*value))
        .collect::<Option<Vec<_>>>()
        .and_then(|values| <[Rational; 3]>::try_from(values).ok())
    else {
        return false;
    };
    size_mm.iter().all(|v| *v > 0. && *v <= MAX_DRAG_MM)
        && l.checked_mul(w)
            .and_then(|area| area.checked_mul(h))
            .is_ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn grid_values_are_bit_exact_decimal_multiples() {
        let snap = DragSnap::new(Some(10.), false);
        assert_eq!(snap.coordinate(91.393_323), Some(90.));
        assert_eq!(snap.length(41.678_413), Some(40.));
        assert_eq!(snap.coordinate(-14.99), Some(-10.));
        let fine = DragSnap::new(Some(0.1), false);
        // Never accumulated: 3 × 0.1 is the decimal 0.3, not 0.30000000000000004.
        assert_eq!(fine.coordinate(0.31), Some(0.3));
        assert_eq!(fine.coordinate(0.31).unwrap().to_string(), "0.3");
        let coarse = DragSnap::new(Some(100.0000000000001), false);
        assert_eq!(coarse.grid_mm(), Some(100.));
        assert_eq!(coarse.coordinate(1249.), Some(1200.));
        assert_eq!(DragSnap::new(Some(2.5), false).coordinate(6.4), Some(7.5));
    }

    #[test]
    fn lengths_and_signed_distances_never_snap_to_zero() {
        let snap = DragSnap::new(Some(10.), false);
        assert_eq!(snap.length(0.), Some(10.));
        assert_eq!(snap.length(-37.), Some(10.));
        assert_eq!(snap.signed(3., 10.), Some(10.));
        assert_eq!(snap.signed(-3., 10.), Some(-10.));
        assert_eq!(snap.signed(0., -20.), Some(-10.));
        assert_eq!(snap.signed(-26., 10.), Some(-30.));
        assert_eq!(DragSnap::FREE.length(0.), Some(FREE_RESOLUTION_MM));
        assert_eq!(snap.coordinate(f64::NAN), None);
        assert_eq!(snap.length(f64::INFINITY), None);
        assert_eq!(snap.coordinate(1e300), Some(MAX_DRAG_MM));
    }

    #[test]
    fn the_modifier_releases_the_grid_to_the_free_resolution() {
        let grid = DragSnap::new(Some(10.), false);
        let free = DragSnap::new(Some(10.), true);
        assert_eq!(grid.grid_mm(), Some(10.));
        assert_eq!(free.grid_mm(), None);
        assert_eq!(free.coordinate(91.393_323_7), Some(91.393));
        assert_eq!(DragSnap::FREE.coordinate(91.393_323_7), Some(91.393));
        assert!(DragSnap::free_modifier(egui::Modifiers::SHIFT));
        assert!(!DragSnap::free_modifier(egui::Modifiers::NONE));
        assert!(!DragSnap::free_modifier(egui::Modifiers::ALT));
        assert!(grid.label().contains("Shift for free"));
        assert!(free.label().starts_with("free (Shift)"));
    }

    #[test]
    fn unusable_grids_fall_back_to_the_free_resolution() {
        for step in [0., -10., f64::NAN, f64::INFINITY, 1e-12, 1e9] {
            assert_eq!(DragSnap::new(Some(step), false).grid_mm(), None, "{step}");
        }
        assert_eq!(Step::new(1. / 3.), None);
    }

    #[test]
    fn box_admission_matches_the_exact_kernel_volume_bound() {
        assert!(box_admissible([90., 40., 20.]));
        assert!(box_admissible([91.393, 41.678, 20.]));
        // The reported drag: three full binary64 decimals overflow the volume.
        assert!(!box_admissible([
            91.393_323_123_456_78,
            41.678_413_987_654_32,
            20.123_456_789_012_34
        ]));
        assert!(!box_admissible([0., 10., 10.]));
        assert!(!box_admissible([2e6, 10., 10.]));
        // A free drag coarsens until the exact volume is admitted.
        let size = [2000.123, 3000.457, 1500.789];
        let step = DragSnap::FREE
            .box_step(|step| box_admissible(size.map(|v| step.value(step.index(v).unwrap()))))
            .unwrap();
        assert!(step.mm() > FREE_RESOLUTION_MM, "{step:?}");
        // A grid step is never coarsened silently.
        let grid = DragSnap::new(Some(10.), false);
        assert_eq!(grid.box_step(|_| false), None);
    }
}
