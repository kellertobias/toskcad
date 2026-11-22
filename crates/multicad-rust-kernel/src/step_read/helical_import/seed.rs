//! Literal first-period tensor domains, without discarding their raw rows.
use super::*;
pub(in super::super) struct TensorHeader<'s, 'input> {
    pub grids: [&'s [Value<'input>]; 4],
    pub radial: [[f64; 2]; 4],
    pub axial: [f64; 2],
}
impl<'s, 'input> TensorHeader<'s, 'input> {
    pub(in super::super) fn retained_origin(
        &self,
        shell: &Shell<'_, '_>,
        first_face: u32,
        entry: Point,
        normal: Point,
    ) -> Result<Point> {
        let first = face(shell, first_face)?;
        let edges = helical::orient(
            &helical::ring(shell, *first.2.first().ok_or(ReadError::Invalid)?, true)?,
            first.1,
        );
        let [_, _, _, (rail, _)] = edges[..] else {
            return Err(ReadError::Invalid);
        };
        let [Value::String(_), a, _, curve, _] = shell.args(rail, "EDGE_CURVE")? else {
            return Err(ReadError::Invalid);
        };
        let args = shell.args(reference(curve)?, "B_SPLINE_CURVE_WITH_KNOTS")?;
        let controls = refs(args.get(2).ok_or(ReadError::Invalid)?)?;
        let point = shell.point(*controls.first().ok_or(ReadError::Invalid)?)?;
        let [Value::String(_), endpoint] = shell.args(reference(a)?, "VERTEX_POINT")? else {
            return Err(ReadError::Invalid);
        };
        if point == shell.point(reference(endpoint)?)? {
            return Ok(entry);
        }
        // The literal full tensor declares its original datum even when the
        // finite face begins at an interior trim. No unseen tool is inferred.
        let ids = refs(&self.grids[0][0])?;
        let first = shell.point(*ids.first().ok_or(ReadError::Invalid)?)?;
        let offset = dot(sub(first, entry), normal);
        Ok(std::array::from_fn(|i| entry[i] + normal[i] * offset))
    }
    pub(in super::super) fn profiled_major(
        shell: &Shell<'_, '_>,
        fourth_face: u32,
        first_forward: bool,
    ) -> Result<f64> {
        let (_, _, bounds) = face(shell, fourth_face)?;
        let edges = helical::orient(
            &helical::ring(shell, *bounds.first().ok_or(ReadError::Invalid)?, true)?,
            first_forward,
        );
        let (edge, _) = *edges.first().ok_or(ReadError::Invalid)?;
        let [Value::String(_), _, _, geometry, _] = shell.args(edge, "EDGE_CURVE")? else {
            return Err(ReadError::Invalid);
        };
        let [Value::String(_), _, Value::Number(radius)] =
            shell.args(reference(geometry)?, "CIRCLE")?
        else {
            return Err(ReadError::Invalid);
        };
        Ok(*radius)
    }
    pub(in super::super) fn read(shell: &'s Shell<'_, 'input>, faces: &[u32]) -> Result<Self> {
        let mut grids: [&[Value<'input>]; 4] = [&[]; 4];
        let mut radial = [[0.; 2]; 4];
        let mut axial = None;
        for i in 0..4 {
            let id = *faces.get(i).ok_or(ReadError::Unsupported)?;
            let (id, _, _) = face(shell, id)?;
            let [
                Value::String(_),
                Value::Number(a),
                Value::Number(b),
                Value::List(rows),
                _,
                _,
                _,
                _,
                _,
                _,
                r,
                z,
                _,
            ] = shell.args(id, "B_SPLINE_SURFACE_WITH_KNOTS")?
            else {
                return Err(ReadError::Invalid);
            };
            require(*a == DEGREE as f64 && *b == DEGREE as f64 && rows.len() > DEGREE)?;
            radial[i] = domain(r)?;
            let domain = domain(z)?;
            if let Some(previous) = axial {
                require(previous == domain)?;
            } else {
                axial = Some(domain);
            }
            grids[i] = rows;
        }
        Ok(Self {
            grids,
            radial,
            axial: axial.ok_or(ReadError::Invalid)?,
        })
    }
    pub(in super::super) fn radius(
        &self,
        shell: &Shell<'_, '_>,
        origin: Point,
        band: usize,
        last: bool,
    ) -> Result<f64> {
        let rows = self.grids[band];
        let row = &rows[if last { rows.len() - 1 } else { 0 }];
        let ids = refs(row)?;
        let point = shell.point(*ids.first().ok_or(ReadError::Invalid)?)?;
        Ok(length(sub(point, origin)))
    }
}
