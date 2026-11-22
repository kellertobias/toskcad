//! Read-only descriptor of an internal thread in cylindrical stock.
use super::*;

/// Exact source data of a cylindrical stock less a helical thread void, for
/// typed analytic consumers such as drawing providers.
///
/// In the thread frame the stock occupies `r <= stock_radius_mm` for
/// `0 <= z <= stock_length_mm`. The void is `r < R(theta, z)` of `thread`
/// for `0 <= z <= thread.spec().length_mm`, which never exceeds the stock
/// length; beyond it, a partial-depth cut keeps the full stock disc as its
/// floor. An overrunning tool is reported as `tool`, with its exact prefix
/// inside stock as `thread`; otherwise both are equal.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InternalThreadData {
    pub thread: Thread,
    pub tool: Thread,
    pub stock_radius_mm: f64,
    pub stock_length_mm: f64,
    /// Thread frame: origin on the entry cap, normal into the stock.
    pub frame: PlaneFrame,
    /// Retained placement of the stock cylinder.
    pub stock_frame: PlaneFrame,
}

impl Solid {
    /// Exact retained geometry of a cylindrical stock with an internal
    /// helical thread. Read-only: other bodies return `None`.
    pub fn rust_internal_thread_data(&self) -> Option<InternalThreadData> {
        let (body, frame) = self.helical_body()?;
        let HelicalForm::Bore(stock) = &body.form else {
            return None;
        };
        Some(InternalThreadData {
            thread: body.geometry,
            tool: stock
                .prefix
                .as_ref()
                .map_or(body.geometry, |prefix| *prefix.tool()),
            stock_radius_mm: stock.radius(),
            stock_length_mm: stock.length(),
            frame,
            stock_frame: stock.frame_in_thread.placed(frame),
        })
    }
}
