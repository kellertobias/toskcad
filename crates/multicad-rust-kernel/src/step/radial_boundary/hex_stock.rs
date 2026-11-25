//! Hex prism stock joined to one already qualified radial boundary, written
//! with every boundary value unchanged and read back before it is accepted.
use super::*;
use crate::step_read::{HexStock, HexStockKind};

impl Writer {
    /// Append a pipe nut (stock around a through cavity) or a pipe hex head
    /// (stock on the thread's large-end cap). Rejected appends are transactional.
    pub fn append_pipe_hex_hardware(
        &mut self,
        boundary: &RadialBoundary,
        stock: HexStock,
        tolerance: f64,
    ) -> Result<(), StepError> {
        let input = canonical(boundary.input());
        let mut candidate = StepFile::new();
        let solid = encode_hex(&mut candidate, &input, stock, tolerance)?;
        candidate.finish(vec![solid]);
        verify_hex(&candidate.data, boundary, stock, tolerance)?;
        self.append(|file| encode_hex(file, &input, stock, tolerance))
    }
}
pub fn write_pipe_hex_hardware(
    boundary: &RadialBoundary,
    stock: HexStock,
    tolerance: f64,
) -> Result<String, StepError> {
    super::super::writer::write_one(|writer| {
        writer.append_pipe_hex_hardware(boundary, stock, tolerance)
    })
}
/// Independent readback of actual bytes: the decoded radial tables must equal
/// the expected boundary exactly and the decoded corners the written ones.
pub fn verify_pipe_hex_hardware(
    text: &str,
    expected: &RadialBoundary,
    stock: HexStock,
    tolerance: f64,
) -> Result<(), StepError> {
    verify_hex(text, expected, stock, tolerance)
}
fn verify_hex(
    text: &str,
    expected: &RadialBoundary,
    stock: HexStock,
    tolerance: f64,
) -> Result<(), StepError> {
    let imported = crate::step_read::read_pipe_hex_hardware(text, tolerance).map_err(|e| {
        if e == crate::step_read::ReadError::Limit {
            StepError::Limit
        } else {
            StepError::Unsupported
        }
    })?;
    let corners = stock.corners().ok_or(StepError::InvalidGeometry)?;
    let bits = |p: &Point| p.map(|v| if v == 0. { 0f64.to_bits() } else { v.to_bits() });
    let mut written = corners.iter().flatten().map(bits).collect::<Vec<_>>();
    let mut decoded = imported.corners.iter().map(bits).collect::<Vec<_>>();
    written.sort_unstable();
    decoded.sort_unstable();
    if imported.kind != stock.kind
        || written != decoded
        || canonical(imported.boundary.input()) != canonical(expected.input())
    {
        return Err(StepError::Unsupported);
    }
    Ok(())
}
fn encode_hex(
    file: &mut StepFile,
    input: &BoundaryInput,
    stock: HexStock,
    tolerance: f64,
) -> Result<Entity, StepError> {
    if !tolerance.is_finite() || tolerance <= 0. {
        return Err(StepError::InvalidGeometry);
    }
    let [lower, upper] = stock.corners().ok_or(StepError::InvalidGeometry)?;
    let [low, high] = stock.span_mm;
    if high - low <= tolerance {
        return Err(StepError::InvalidGeometry);
    }
    let normal = stock.normal;
    let height = |p: Point| {
        (0..3)
            .map(|i| (p[i] - stock.origin[i]) * normal[i])
            .sum::<f64>()
    };
    // Nut: the cavity caps at both span ends become annuli. Head: the cap at
    // the lower span end becomes the shoulder; the other cap stays below it.
    let mut annuli = [None; 2];
    for (index, face) in input.faces.iter().enumerate() {
        if let Surface::Plane(a) = input.surfaces[face.surface] {
            let h = height(a.origin);
            let end = if (h - low).abs() <= tolerance {
                0
            } else if (h - high).abs() <= tolerance {
                1
            } else {
                match stock.kind {
                    HexStockKind::Head if h < low - tolerance => continue,
                    _ => return Err(StepError::InvalidGeometry),
                }
            };
            if stock.kind == HexStockKind::Head && end == 1 {
                return Err(StepError::InvalidGeometry);
            }
            if annuli[end].replace(index).is_some() {
                return Err(StepError::InvalidGeometry);
            }
        }
    }
    match (stock.kind, annuli) {
        (HexStockKind::Nut, [Some(_), Some(_)]) | (HexStockKind::Head, [Some(_), None]) => {}
        _ => return Err(StepError::InvalidGeometry),
    }
    let Tables { surfaces, edges } = tables(file, input);
    let nut = stock.kind == HexStockKind::Nut;
    let mut faces = Vec::with_capacity(input.faces.len() + 7);
    let mut holes = Vec::new();
    for (index, f) in input.faces.iter().enumerate() {
        let end = annuli.iter().position(|a| *a == Some(index));
        // Nut material is outside the cavity, so each cavity wire reverses.
        let wire = if nut {
            f.outer
                .iter()
                .rev()
                .map(|&(e, d)| file.oriented_edge(edges[e], !d))
                .collect::<Vec<_>>()
        } else {
            f.outer
                .iter()
                .map(|&(e, d)| file.oriented_edge(edges[e], d))
                .collect::<Vec<_>>()
        };
        if let Some(end) = end {
            // The annulus keeps the cap's own plane. A nut annulus faces the
            // same way as the cavity cap; a head shoulder faces back along
            // the thread, so its sense flips and the rim stays a hole as is.
            let hole = file.bound(&wire, false);
            holes.push((
                faces.len(),
                end,
                hole,
                if nut { f.forward } else { !f.forward },
                f.surface,
            ));
            faces.push(None);
        } else {
            let bound = file.bound(&wire, true);
            faces.push(Some(file.face_on(
                &[bound],
                surfaces[f.surface],
                f.forward != nut,
            )));
        }
    }
    let mut points = Vec::with_capacity(12);
    points.extend(lower);
    points.extend(upper);
    let mut brep = BRep::new(file, points);
    let lower_out: Vec<usize> = (0..6).rev().collect();
    let upper_out: Vec<usize> = (6..12).collect();
    for (slot, end, hole, forward, surface) in holes {
        let ring = if end == 0 { &lower_out } else { &upper_out };
        let outer = brep.bound(ring, true)?;
        faces[slot] = Some(
            brep.file
                .face_on(&[outer, hole], surfaces[surface], forward),
        );
    }
    let mut faces = faces.into_iter().flatten().collect::<Vec<_>>();
    for i in 0..6 {
        let j = (i + 1) % 6;
        let edge = direction(lower[i], lower[j])?;
        let outward = [
            edge[1] * normal[2] - edge[2] * normal[1],
            edge[2] * normal[0] - edge[0] * normal[2],
            edge[0] * normal[1] - edge[1] * normal[0],
        ];
        let bound = brep.bound(&[i, j, 6 + j, 6 + i], true)?;
        let point = brep.vertices[i].cartesian;
        faces.push(brep.file.face(&[bound], point, outward, edge));
    }
    if !nut {
        let bound = brep.bound(&upper_out, true)?;
        let point = brep.vertices[6].cartesian;
        faces.push(brep.file.face(&[bound], point, normal, stock.reference));
    }
    Ok(brep.file.solid(&faces))
}
