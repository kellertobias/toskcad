//! Complete periodic tensor support inventory for the retained Hanger source.
//! Cells still require finite trims, shared literal curves and a closed-shell
//! material proof before they may be serialized as an occupied body.
use super::*;
use crate::authored_hanger_bolt::{Boundary, Component, CurveKey, End, FaceKey};
use crate::authored_internal_treatment::ThreadBandRole;
use crate::step::insert::reduced_unit::fit_reduced_unit;
use num_rational::BigRational;
use num_traits::ToPrimitive;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FieldControl {
    /// Homogeneous coefficient enclosure for the shared polynomial angular
    /// approximation and the ORIGINAL symbolic meridian. Not a point sample.
    pub coefficient: [[f64; 2]; 4],
    pub literal: [f64; 4],
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FieldBounds {
    pub position_error_mm: f64,
    /// Original normalized meridian-u and whole unit physical-turn-v.
    pub derivative_error_mm: [f64; 2],
}
pub struct CellSurface {
    key: MeridianRowKey,
    degree: usize,
    columns: usize,
    controls: Vec<FieldControl>,
    bounds: FieldBounds,
}
impl CellSurface {
    pub fn key(&self) -> MeridianRowKey {
        self.key
    }
    pub fn meridian_degree(&self) -> usize {
        self.degree
    }
    pub fn angular_columns(&self) -> usize {
        self.columns
    }
    /// Meridian-major; angular degree10, shared endpoints between pieces.
    pub fn controls(&self) -> &[FieldControl] {
        &self.controls
    }
    pub fn bounds(&self) -> FieldBounds {
        self.bounds
    }
}
pub struct BodyField<'a> {
    source: SourceField<'a>,
    cells: Vec<CellSurface>,
    knots: [Vec<f64>; 2],
    work: usize,
    bytes: usize,
}
fn enclosure(q: &BigRational, w: &mut Work<'_>) -> Result<I, Error> {
    // Sealed ReducedUnitAngular stores96-bit dyadics, magnitude<2. Conversion
    // and comparisons have<1200-bit operands even for a subnormal proposal.
    w.charge(3 * 20 * 20)?;
    checked_enclosure(q).ok_or(Error::Unresolved)
}
fn checked_enclosure(q: &BigRational) -> Option<I> {
    let f = q.to_f64().filter(|x| x.is_finite())?;
    let lo = f.next_down();
    let hi = f.next_up();
    let a = BigRational::from_float(lo)?;
    let b = BigRational::from_float(hi)?;
    (q >= &a && q <= &b).then_some(I::new(lo, hi))
}

fn pair(x: I) -> [f64; 2] {
    [x.lo, x.hi]
}
fn nested(options: CurveOptions, w: &Work<'_>) -> Result<DomainOptions, Error> {
    if w.remaining() == 0 {
        return Err(Error::Limit);
    }
    Ok(DomainOptions {
        maximum_work: w.remaining(),
        maximum_bytes: options.bytes,
    })
}
impl<'a> BodyField<'a> {
    pub fn prepare(
        boundary: &'a Boundary,
        options: CurveOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, Error> {
        options.validate()?;
        let mut w = Work::new(options.work, interrupt)?;
        w.charge(1)?;
        let mut memory = Memory::new(boundary.retained_bytes() + size_of::<Self>(), options.bytes)?;
        // Fixed interval/rational conversion temporaries, including nested
        // meridian queries; output Vec capacities are accounted separately.
        memory.add(256 * 1024)?;
        let source = SourceField::prepare(boundary, nested(options, &w)?, interrupt)?;
        w.charge(source.work_used())?;
        let mut cells = memory.allocate(0, &mut w)?;
        // Sealed source permits at most128 period cells/component and5bands.
        // Preallocate the full count, including partially trimmed end cells.
        let count = [Component::Metric, Component::Wood]
            .into_iter()
            .map(|c| {
                let r = source.cell_periods(c);
                (r[1] - r[0]) as usize * 5
            })
            .sum();
        memory.release(&cells);
        cells = memory.allocate(count, &mut w)?;
        let mut knots = [Vec::new(), Vec::new()];
        for (component_index, component) in
            [Component::Metric, Component::Wood].into_iter().enumerate()
        {
            let s = boundary.source();
            let (pitch, starts, left, radius) = match component {
                Component::Metric => {
                    let v = s.metric().spec();
                    (
                        v.pitch_mm,
                        v.starts,
                        v.left_hand,
                        s.metric().maximum_radius_mm(),
                    )
                }
                Component::Wood => {
                    let v = s.wood().spec();
                    (
                        v.pitch_mm,
                        v.starts,
                        v.left_hand,
                        s.wood().maximum_radius_mm(),
                    )
                }
            };
            let scale = (p(radius) + p(pitch) * p(f64::from(starts)) + p(1.)).hi;
            let angular_options = CurveOptions {
                position_error_mm: options.position_error_mm / (64. * scale),
                derivative_error_mm: options.derivative_error_mm / (64. * scale),
                ..options
            };
            if angular_options.position_error_mm == 0. || angular_options.derivative_error_mm == 0.
            {
                return Err(Error::Unresolved);
            }
            let angular = fit_reduced_unit(left, angular_options, &mut memory, &mut w)?;
            let mut converted = memory.allocate(angular.coefficients().len() * 11, &mut w)?;
            for (i, piece) in angular.coefficients().iter().enumerate() {
                for (j, xy) in piece.iter().enumerate() {
                    w.charge(64)?;
                    let v = if j == 0 {
                        p(angular.knots()[i])
                    } else if j == 10 {
                        p(angular.knots()[i + 1])
                    } else {
                        p(angular.knots()[i])
                            + (p(angular.knots()[i + 1]) - p(angular.knots()[i])) * p(j as f64)
                                / p(10.)
                    };
                    converted.push(([enclosure(&xy[0], &mut w)?, enclosure(&xy[1], &mut w)?], v));
                }
            }
            // At most two endpoints per five bands plus the root middle.
            // Actual enclosure/literal keys are still compared before reuse.
            let mut xy_cache = memory.allocate(11, &mut w)?;
            knots[component_index] = memory.allocate(angular.knots().len(), &mut w)?;
            knots[component_index].extend_from_slice(angular.knots());
            let range = source.cell_periods(component);
            for curve in boundary.curves() {
                w.charge(1)?;
                let CurveKey::Cap {
                    component: owner,
                    end: End::Start,
                    period,
                    band,
                } = curve.key
                else {
                    continue;
                };
                if owner != component {
                    continue;
                }
                let face = *curve
                    .faces
                    .iter()
                    .find(|f| {
                        matches!(f, FaceKey::MetricThread { .. } | FaceKey::WoodThread { .. })
                    })
                    .ok_or(Error::Invalid)?;
                let pieces: &[BandPiece] = if band == ThreadBandRole::Crest {
                    &[BandPiece::CrestLowPhase, BandPiece::CrestHighPhase]
                } else {
                    &[BandPiece::Whole]
                };
                for &piece in pieces {
                    for j in range[0]..range[1] {
                        w.charge(1)?;
                        let reference = j - i32::from(piece == BandPiece::CrestLowPhase);
                        if reference.rem_euclid(starts as i32) != period as i32 {
                            continue;
                        }
                        let lower = source.meridian_row(
                            face,
                            piece,
                            CellSeam {
                                period: j,
                                upper_turn: false,
                            },
                            nested(options, &w)?,
                            interrupt,
                        )?;
                        w.charge(lower.work_used())?;
                        let upper = source.meridian_row(
                            face,
                            piece,
                            CellSeam {
                                period: j,
                                upper_turn: true,
                            },
                            nested(options, &w)?,
                            interrupt,
                        )?;
                        w.charge(upper.work_used())?;
                        cells.push(net(
                            [lower, upper],
                            &AngularNetData {
                                controls: &converted,
                                knots: angular.knots(),
                            },
                            &mut xy_cache,
                            [angular.position_error(), angular.derivative_error()],
                            options,
                            &mut memory,
                            &mut w,
                        )?);
                    }
                }
            }
            for row in &xy_cache {
                memory.release(&row.controls);
            }
            memory.release(&xy_cache);
            memory.release(&converted);
            angular.release_memory(&mut memory);
        }
        w.charge(1)?;
        memory.used -= 256 * 1024;
        Ok(Self {
            source,
            cells,
            knots,
            work: w.used,
            bytes: memory.used,
        })
    }
    pub fn source_field(&self) -> &SourceField<'a> {
        &self.source
    }
    pub fn cells(&self) -> &[CellSurface] {
        &self.cells
    }
    pub fn angular_knots(&self, component: Component) -> &[f64] {
        &self.knots[usize::from(component == Component::Wood)]
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn retained_bytes(&self) -> usize {
        self.bytes
    }
}
struct AngularNetData<'a> {
    controls: &'a [([I; 2], I)],
    knots: &'a [f64],
}
#[derive(Clone, Copy, PartialEq)]
struct XYKey {
    source: [[f64; 2]; 4],
    literal: [[f64; 2]; 2],
}
struct XYControl {
    coefficient: [[f64; 2]; 2],
    literal: [f64; 2],
}
struct XYRow {
    key: XYKey,
    controls: Vec<XYControl>,
}
fn xy_row<'a>(
    key: XYKey,
    angular: &AngularNetData<'_>,
    cache: &'a mut Vec<XYRow>,
    memory: &mut Memory,
    w: &mut Work<'_>,
) -> Result<&'a XYRow, Error> {
    // At most11 actual source rows. Charge the complete bounded scan even
    // when a match occurs early; no callback is hidden inside an iterator.
    w.charge(16 * cache.len())?;
    if let Some(i) = cache.iter().position(|row| row.key == key) {
        return Ok(&cache[i]);
    }
    let spans = angular.knots.len() - 1;
    let columns = 10 * spans + 1;
    let mut controls = memory.allocate(columns, w)?;
    for column in 0..columns {
        w.charge(96)?;
        let piece = (column / 10).min(spans - 1);
        let k = column - piece * 10;
        let coefficient = if column == 0 {
            [key.source[0], key.source[1]]
        } else if column + 1 == columns {
            [key.source[2], key.source[3]]
        } else {
            angular.controls[piece * 11 + k]
                .0
                .map(|x| pair(hull(key.source[0]) * x))
        };
        let literal = if column == 0 {
            key.literal[0]
        } else if column + 1 == columns {
            key.literal[1]
        } else {
            coefficient.map(|x| hull(x).midpoint())
        };
        controls.push(XYControl {
            coefficient,
            literal,
        });
    }
    // Source constructors retain at most11 actual homogeneous radial rows;
    // a shared radius can only reduce this count. No implicit Vec growth.
    cache.push(XYRow { key, controls });
    Ok(&cache[cache.len() - 1])
}
fn net(
    rows: [MeridianRow; 2],
    angular: &AngularNetData<'_>,
    cache: &mut Vec<XYRow>,
    angular_error: [f64; 2],
    options: CurveOptions,
    memory: &mut Memory,
    w: &mut Work<'_>,
) -> Result<CellSurface, Error> {
    let [lower, upper] = rows;
    let knots = angular.knots;
    let degree = lower.degree();
    let spans = knots.len() - 1;
    let columns = 10 * spans + 1;
    let mut controls = memory.allocate((degree + 1) * columns, w)?;
    let mut radial = 0_f64;
    let mut radial_du = 0_f64;
    let mut weight_du = 0_f64;
    let minimum_weight = lower
        .controls()
        .iter()
        .map(|c| c.source_bounds[3][0])
        .fold(1_f64, f64::min);
    for (row, (a, b)) in lower.controls().iter().zip(upper.controls()).enumerate() {
        let ah = a.source_bounds.map(hull);
        let bh = b.source_bounds.map(hull);
        radial = radial.max((ah[0] / ah[3]).magnitude());
        if row < degree {
            let next = lower.controls()[row + 1].source_bounds.map(hull);
            radial_du = radial_du.max((p(degree as f64) * (next[0] - ah[0])).magnitude());
            weight_du = weight_du.max((p(degree as f64) * (next[3] - ah[3])).magnitude());
        }
        let key = XYKey {
            source: [
                a.source_bounds[0],
                a.source_bounds[1],
                b.source_bounds[0],
                b.source_bounds[1],
            ],
            literal: [[a.literal[0], a.literal[1]], [b.literal[0], b.literal[1]]],
        };
        let xy = xy_row(key, angular, cache, memory, w)?;
        let dz = bh[2] - ah[2];
        w.charge(16)?;
        for column in 0..columns {
            // After immutable XY/v reuse, only one interval multiply/add,
            // one midpoint and coefficient copying remain (<48 scalar steps).
            w.charge(48)?;
            let piece = (column / 10).min(spans - 1);
            let k = column - piece * 10;
            let v = angular.controls[piece * 11 + k].1;
            let xy = &xy.controls[column];
            let expected_z = if column == 0 {
                ah[2]
            } else if column + 1 == columns {
                bh[2]
            } else {
                ah[2] + dz * v
            };
            let literal_z = if column == 0 {
                a.literal[2]
            } else if column + 1 == columns {
                b.literal[2]
            } else {
                expected_z.midpoint()
            };
            controls.push(FieldControl {
                coefficient: [
                    xy.coefficient[0],
                    xy.coefficient[1],
                    pair(expected_z),
                    a.source_bounds[3],
                ],
                literal: [xy.literal[0], xy.literal[1], literal_z, a.literal[3]],
            });
        }
    }
    let wm = p(minimum_weight);
    let source_du = (p(radial_du) / wm + p(radial) * p(weight_du) / wm).hi;
    let bounds = certify(
        &controls,
        NetDomain {
            degree,
            columns,
            knots,
            angular_error,
            radial: [radial, source_du],
        },
        options,
        w,
    )?;
    Ok(CellSurface {
        key: lower.key(),
        degree,
        columns,
        controls,
        bounds,
    })
}
struct NetDomain<'a> {
    degree: usize,
    columns: usize,
    knots: &'a [f64],
    angular_error: [f64; 2],
    radial: [f64; 2],
}
fn certify(
    controls: &[FieldControl],
    domain: NetDomain<'_>,
    options: CurveOptions,
    w: &mut Work<'_>,
) -> Result<FieldBounds, Error> {
    let NetDomain {
        degree,
        columns,
        knots,
        angular_error,
        radial,
    } = domain;
    let mut e = [0_f64; 2]; // XYZ L1 and scalar weight coefficient residuals.
    let mut de = [[0_f64; 2]; 2];
    let mut a = 0_f64;
    let mut da = [0_f64; 2];
    let mut dw = [0_f64; 2];
    let mut literal_dw = [0_f64; 2];
    let mut source_min = f64::INFINITY;
    let mut literal_min = f64::INFINITY;
    let source = |row: usize, col: usize| controls[row * columns + col].coefficient.map(hull);
    let literal = |row: usize, col: usize| controls[row * columns + col].literal.map(p);
    for row in 0..=degree {
        for col in 0..columns {
            w.charge(128)?;
            let h = source(row, col);
            let g = literal(row, col);
            let delta: [I; 4] = std::array::from_fn(|k| h[k] - g[k]);
            source_min = source_min.min(h[3].lo);
            literal_min = literal_min.min(g[3].lo);
            a = a.max(l1([h[0], h[1], h[2]]).hi);
            e[0] = e[0].max(l1([delta[0], delta[1], delta[2]]).hi);
            e[1] = e[1].max(delta[3].magnitude());
            for axis in 0..2 {
                let next = if axis == 0 {
                    (row < degree).then(|| {
                        (
                            source(row + 1, col),
                            literal(row + 1, col),
                            p(degree as f64),
                        )
                    })
                } else if col + 1 < columns {
                    let piece = col / 10;
                    Some((
                        source(row, col + 1),
                        literal(row, col + 1),
                        p(10.) / (p(knots[piece + 1]) - p(knots[piece])),
                    ))
                } else {
                    None
                };
                if let Some((n, m, scale)) = next {
                    let dh: [I; 4] = std::array::from_fn(|k| scale * (n[k] - h[k]));
                    let dg: [I; 4] = std::array::from_fn(|k| scale * (m[k] - g[k]));
                    let d: [I; 4] = std::array::from_fn(|k| dh[k] - dg[k]);
                    da[axis] = da[axis].max(l1([dh[0], dh[1], dh[2]]).hi);
                    dw[axis] = dw[axis].max(dh[3].magnitude());
                    literal_dw[axis] = literal_dw[axis].max(dg[3].magnitude());
                    de[axis][0] = de[axis][0].max(l1([d[0], d[1], d[2]]).hi);
                    de[axis][1] = de[axis][1].max(d[3].magnitude());
                }
            }
        }
    }
    if source_min <= 0. || literal_min <= 0. {
        return Err(Error::Unresolved);
    }
    w.charge(256)?;
    let f = p(a) / p(source_min);
    let numerator = p(e[0]) + f * p(e[1]);
    let position = numerator / p(literal_min) + p(radial[0]) * p(angular_error[0]);
    let derivative = std::array::from_fn(|axis| {
        let df = (p(da[axis]) + f * p(dw[axis])) / p(source_min);
        let quotient = (p(de[axis][0]) + df * p(e[1]) + f * p(de[axis][1])) / p(literal_min)
            + numerator * p(literal_dw[axis]) / p(literal_min).square();
        (quotient
            + if axis == 0 {
                p(radial[1]) * p(angular_error[0])
            } else {
                p(radial[0]) * p(angular_error[1])
            })
        .hi
    });
    let bounds = FieldBounds {
        position_error_mm: position.hi,
        derivative_error_mm: derivative,
    };
    if !bounds.position_error_mm.is_finite()
        || bounds.position_error_mm > options.position_error_mm
        || bounds
            .derivative_error_mm
            .iter()
            .any(|x| !x.is_finite() || *x > options.derivative_error_mm)
    {
        return Err(Error::Unresolved);
    }
    Ok(bounds)
}
#[cfg(test)]
mod tests;
