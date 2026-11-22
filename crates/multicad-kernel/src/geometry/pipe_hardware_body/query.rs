//! Whole inverse-query boxes, never a midpoint substituted for a posed point.
use super::*;
fn abs(p: I) -> I {
    I::new(
        if p.lo <= 0. && p.hi >= 0. {
            0.
        } else {
            p.lo.abs().min(p.hi.abs())
        },
        p.magnitude(),
    )
}
impl Solid {
    pub fn rust_pipe_hardware_contains(&self, point: [f64; 3]) -> Result<bool, GeometryError> {
        self.policy.point(point)?;
        self.pipe_hardware_body()
            .ok_or_else(|| failure(Cause::Unsupported))?
            .contains(point)
    }
}
impl PipeHardwareBody {
    pub(in crate::geometry) fn contains(&self, point: [f64; 3]) -> Result<bool, GeometryError> {
        let p = self.placement.inverse_point(point.map(I::point))?;
        let (af, head) = match self.source.intent.kind() {
            PipeHardwareKind::HexHead {
                across_flats_mm,
                height_mm,
            } => (across_flats_mm, Some(height_mm)),
            PipeHardwareKind::HexNut { across_flats_mm } => (across_flats_mm, None),
        };
        let length = self.source.intent.thread().spec().length_mm;
        let top = if let Some(h) = head {
            I::point(length) + I::point(h)
        } else {
            I::point(length)
        };
        if p[2].hi < 0. || p[2].lo > top.hi {
            return Ok(false);
        }
        let x = abs(p[0]);
        let y = abs(p[1]);
        let diagonal = I::point(3.).sqrt() * x + y;
        let half = I::point(af) / I::point(2.);
        if y.lo > half.hi || diagonal.lo > af {
            return Ok(false);
        }
        let hex_inside = y.hi < half.lo && diagonal.hi < af;
        if !hex_inside {
            return Err(failure(Cause::Unresolved));
        }
        if head.is_some() && p[2].lo > length && p[2].hi < top.lo {
            return Ok(true);
        }
        if p[2].lo <= 0. || p[2].hi >= length {
            return Err(failure(Cause::Unresolved));
        }
        let bounds = p.map(|x| [x.lo, x.hi]);
        let inside = if let Some(source) = self.source.thread.rust_clipped_npt_source() {
            let query = source
                .classify_box(bounds, Default::default(), None)
                .map_err(failure)?;
            match query.classification() {
                multicad_rust_kernel::axial_pipe_thread::authored::StockClipClassification::Inside=>true,
                multicad_rust_kernel::axial_pipe_thread::authored::StockClipClassification::Outside=>false,
                multicad_rust_kernel::axial_pipe_thread::authored::StockClipClassification::Undetermined=>return Err(failure(Cause::Unresolved)),
            }
        } else {
            let (source, _) = self
                .source
                .thread
                .rust_radial_boundary()
                .expect("sealed actual retained pipe source");
            let query = source
                .classify_box(bounds, Default::default(), None)
                .map_err(radial_failure)?;
            match query.location() {
                multicad_rust_kernel::radial_tensor::Location::Inside => true,
                multicad_rust_kernel::radial_tensor::Location::Outside => false,
                _ => return Err(failure(Cause::Unresolved)),
            }
        };
        Ok(if head.is_some() { inside } else { !inside })
    }
}
