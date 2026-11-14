//! Numeric authority for one finite span of an immutable parsed 3D curve.
//! This earns no face, closed wire, material, source correspondence or mass.
use super::*;

pub(super) struct ParsedCapSpan<'p, 'i> {
    source: &'p Parsed<'i>,
    curve: &'p SplineCurve<'i>,
    span: usize,
    points: [[I; 3]; 12],
}
impl<'p, 'i> ParsedCapSpan<'p, 'i> {
    pub(super) fn load(
        source: &'p Parsed<'i>,
        curve_entity: u32,
        span: usize,
        cache: &[Option<[I; 3]>],
        w: &mut Work<'_>,
    ) -> Result<Self> {
        // Sorted immutable curve lookup needs at most floor(log2(n))+2
        // comparisons under the pinned slice binary-search implementation.
        w.charge(source.curves().len().max(1).ilog2() as usize + 2)
            .map_err(Error::from)?;
        let curve = source.curve(curve_entity).ok_or(Error::Invalid)?;
        if curve.degree() != 11 || span >= curve.knots().len() - 1 {
            return Err(Error::Unsupported);
        }
        // Parsed earns pole_count=(knot_count-1)*degree+1. Each original
        // coordinate enclosure is read from that SAME arena/cache; no fitting.
        let mut points = [[I::point(0.); 3]; 12];
        for (i, point) in points.iter_mut().enumerate() {
            *point = cached(source, cache, curve.pole_entities()[span * 11 + i], w)?;
        }
        Ok(Self {
            source,
            curve,
            span,
            points,
        })
    }
    pub(super) fn points(&self) -> &[[I; 3]; 12] {
        &self.points
    }
    pub(super) fn angular_lower(
        &self,
        sign: f64,
        stack: &mut Vec<projection::Cell>,
        memory: &mut Memory,
        w: &mut Work<'_>,
    ) -> Result<f64> {
        if sign != 1. && sign != -1. {
            return Err(Error::Invalid);
        }
        let x = P::controls(&self.points.map(|v| v[0]));
        let y = P::controls(&self.points.map(|v| v[1]));
        let h = [x.hull(), y.hull()];
        let scale = h[0].magnitude().max(h[1].magnitude()).max(1.);
        w.charge(24).map_err(Error::from)?;
        let tiny = h.iter().all(|a| a.hi - a.lo <= 64. * f64::EPSILON * scale);
        let proof = if tiny {
            None
        } else {
            let k = x
                .product(y.derivative(w).map_err(Error::from)?, w)
                .map_err(Error::from)?
                .sub(
                    y.product(x.derivative(w).map_err(Error::from)?, w)
                        .map_err(Error::from)?,
                    w,
                )
                .map_err(Error::from)?
                .scale(sign, w)
                .map_err(Error::from)?;
            Some(projection::positive([k; 3], stack, w))
        };
        if matches!(
            proof,
            None | Some(Err(crate::step::insert::Error::Unresolved))
        ) {
            stack.clear();
            let mut tokens = [[""; 2]; 12];
            for (i, xy) in tokens.iter_mut().enumerate() {
                let index = point_index(
                    self.source,
                    self.curve.pole_entities()[self.span * 11 + i],
                    w,
                )?;
                for (axis, token) in xy.iter_mut().enumerate() {
                    *token = self.source.points()[index].coordinates()[axis].spelling();
                }
            }
            crate::step::insert::exact_cap::angular_lower(tokens, sign < 0., memory, w)
                .map_err(Error::from)
        } else {
            Ok(proof
                .expect("ordinary cap proof selected")
                .map_err(Error::from)?
                .0
                .lo)
        }
    }
}

#[cfg(test)]
mod tests;
