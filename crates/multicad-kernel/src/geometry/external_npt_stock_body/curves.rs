use super::*;
use multicad_rust_kernel::display::{MeshError, Tessellation};
impl ExternalNptStockBody {
    pub(in crate::geometry) fn sample_edge(
        &self,
        index: u32,
        q: Tessellation,
    ) -> Option<Result<Vec<Point>, MeshError>> {
        let curve = *self.boundary.curves().get(index as usize)?;
        Some((|| {
            if !q.deflection.is_finite()
                || q.deflection <= 0.
                || !(0.01..=1.).contains(&q.angle)
                || q.triangle_limit == 0
            {
                return Err(MeshError::InvalidOptions);
            }
            let placement = self.calibration();
            if let SourceCurve::Thread(key) = curve {
                let scale = iv(self.source.radial_scale_bounds()).hi.max(1.);
                let metric = (I::point(1.) + I::point(placement.metric_error()))
                    .sqrt()
                    .hi;
                let poly = self
                    .source
                    .thread()
                    .curve_polyline(
                        key,
                        q.deflection / (4. * scale * metric),
                        q.triangle_limit.min(131072),
                        Default::default(),
                        None,
                    )
                    .map_err(|e| match e {
                        Cause::Invalid => MeshError::InvalidOptions,
                        Cause::Limit => MeshError::BudgetExceeded,
                        _ => MeshError::UnresolvedBoundary,
                    })?;
                let inherited = self
                    .source
                    .work_used()
                    .checked_sub(self.source.thread().work_used())
                    .ok_or(MeshError::BudgetExceeded)?;
                let total = poly
                    .work_used()
                    .checked_add(inherited)
                    .and_then(|v| v.checked_add(poly.points_mm().len().checked_mul(128)?))
                    .ok_or(MeshError::BudgetExceeded)?;
                if total > 16_000_000 {
                    return Err(MeshError::BudgetExceeded);
                }
                let base = self.source.stock_axis_base_bounds_mm().map(iv);
                let radial = iv(self.source.radial_scale_bounds());
                let mut out = Vec::with_capacity(poly.points_mm().len());
                let mut rounding = 0_f64;
                for p in poly.points_mm() {
                    let local = [
                        I::point(p[0]) * radial + base[0],
                        I::point(p[1]) * radial + base[1],
                        I::point(p[2]) + base[2],
                    ];
                    let exact = placement
                        .point_bounds(local.map(|v| [v.lo, v.hi]))
                        .map_err(|_| MeshError::UnresolvedBoundary)?;
                    let literal = exact.map(|v| v[0].midpoint(v[1]));
                    rounding = rounding
                        .max(norm(std::array::from_fn(|i| iv(exact[i]) - I::point(literal[i]))).hi);
                    out.push(literal);
                }
                if (I::point(scale) * I::point(metric) * I::point(poly.source_chord_error_mm())
                    + I::point(rounding))
                .hi > q.deflection
                {
                    return Err(MeshError::UnresolvedBoundary);
                }
                return Ok(out);
            }
            let mut work = self.source.work_used();
            let mut bound = |a: f64, b: f64| -> Result<[[f64; 2]; 3], MeshError> {
                let p = self
                    .source
                    .curve_bounds(curve, [a, b], Default::default(), None)
                    .map_err(|_| MeshError::UnresolvedBoundary)?;
                work = work
                    .checked_add(
                        p.work_used()
                            .checked_sub(self.source.work_used())
                            .ok_or(MeshError::BudgetExceeded)?
                            + 128,
                    )
                    .ok_or(MeshError::BudgetExceeded)?;
                if work > 16_000_000 {
                    return Err(MeshError::BudgetExceeded);
                }
                placement
                    .point_bounds(p.position_mm())
                    .map_err(|_| MeshError::UnresolvedBoundary)
            };
            let first = bound(0., 0.)?.map(|v| v[0].midpoint(v[1]));
            let mut out = vec![first];
            let mut todo = vec![(0., 1., 0usize)];
            while let Some((a, b, depth)) = todo.pop() {
                let enclosure = bound(a, b)?;
                let error = norm(std::array::from_fn(|i| {
                    iv(enclosure[i]) - I::point(out.last().unwrap()[i])
                }))
                .hi;
                // Every point of the curve and its eventual segment lies in
                // this ball; diameter is a conservative chord certificate.
                if error <= q.deflection {
                    let endpoint = bound(b, b)?.map(|v| v[0].midpoint(v[1]));
                    out.push(endpoint);
                    if out.len() > q.triangle_limit.min(131072) {
                        return Err(MeshError::BudgetExceeded);
                    }
                } else {
                    if depth >= 24 {
                        return Err(MeshError::UnresolvedBoundary);
                    }
                    let m = a.midpoint(b);
                    todo.push((m, b, depth + 1));
                    todo.push((a, m, depth + 1));
                }
            }
            Ok(out)
        })())
    }
}
