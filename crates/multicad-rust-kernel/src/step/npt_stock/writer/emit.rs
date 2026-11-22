use super::*;
use crate::axial_pipe_thread::authored::{StockClipCoedge, StockClipStation};
struct Multiplicities(usize);
impl fmt::Display for Multiplicities {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for i in 0..=self.0 {
            if i > 0 {
                f.write_str(",")?;
            }
            write!(f, "{}", if i == 0 || i == self.0 { 12 } else { 11 })?;
        }
        Ok(())
    }
}
struct Knots(usize);
impl fmt::Display for Knots {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for i in 0..=self.0 {
            if i > 0 {
                f.write_str(",")?;
            }
            write!(f, "{:.16E}", i as f64 / self.0 as f64)?;
        }
        Ok(())
    }
}
struct Rows<'a> {
    points: &'a [Entity],
    height: usize,
}
impl fmt::Display for Rows<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (i, row) in self.points.chunks_exact(self.height).enumerate() {
            if i > 0 {
                f.write_str(",")?;
            }
            write!(f, "{}", Entities(row))?;
        }
        Ok(())
    }
}
fn spline(encoder: &mut Encoder<'_, '_>, points: &[Entity], spans: usize) -> Result<Entity, Error> {
    encoder.entity(format_args!(
        "B_SPLINE_CURVE_WITH_KNOTS('',11,{},.UNSPECIFIED.,.F.,.F.,({}),({}),.UNSPECIFIED.)",
        Entities(points),
        Multiplicities(spans),
        Knots(spans)
    ))
}
fn line_pcurve(
    encoder: &mut Encoder<'_, '_>,
    support: Entity,
    context: Entity,
    side: usize,
) -> Result<Entity, Error> {
    let endpoints = match side {
        0 => [[0., 0.], [0., 1.]],
        1 => [[0., 1.], [1., 1.]],
        2 => [[1., 0.], [1., 1.]],
        _ => [[0., 0.], [1., 0.]],
    };
    let mut points = encoder.allocate(2)?;
    for [u, v] in endpoints {
        points.push(encoder.entity(format_args!("CARTESIAN_POINT('',({u:.16E},{v:.16E}))"))?);
    }
    let curve = encoder.entity(format_args!(
        "B_SPLINE_CURVE_WITH_KNOTS('',1,{},.UNSPECIFIED.,.F.,.F.,(2,2),(0.,1.),.UNSPECIFIED.)",
        Entities(&points)
    ))?;
    let pcurve = pcurve(encoder, support, context, curve)?;
    encoder.memory.release(&points);
    drop(points);
    Ok(pcurve)
}
fn pcurve(
    encoder: &mut Encoder<'_, '_>,
    support: Entity,
    context: Entity,
    curve: Entity,
) -> Result<Entity, Error> {
    let representation = encoder.entity(format_args!(
        "DEFINITIONAL_REPRESENTATION('',({curve}),{context})"
    ))?;
    encoder.entity(format_args!("PCURVE('',{support},{representation})"))
}
fn bound(
    encoder: &mut Encoder<'_, '_>,
    coedges: &[StockClipCoedge],
    plan: &SurfacePlan<'_>,
    edges: &[Entity],
) -> Result<Entity, Error> {
    let mut oriented = encoder.allocate(coedges.len())?;
    for coedge in coedges {
        let index = plan
            .curves()
            .curves()
            .iter()
            .position(|c| c.source_curve() == coedge.curve)
            .expect("original NPT edge owner");
        let edge = edges[index];
        let sense = if coedge.forward { ".T." } else { ".F." };
        oriented.push(encoder.entity(format_args!("ORIENTED_EDGE('',*,*,{edge},{sense})"))?);
    }
    let wire = encoder.entity(format_args!("EDGE_LOOP('',{})", Entities(&oriented)))?;
    encoder.memory.release(&oriented);
    drop(oriented);
    encoder.entity(format_args!("FACE_OUTER_BOUND('',{wire},.T.)"))
}
pub(super) fn body(encoder: &mut Encoder<'_, '_>, plan: &SurfacePlan<'_>) -> Result<Entity, Error> {
    let source = plan.curves().source();
    let owners = plan.curves().vertices();
    let mut xyz = encoder.allocate(12)?;
    let mut vertices = encoder.allocate(12)?;
    for vertex in owners {
        let point = encoder.point(vertex.position_mm())?;
        xyz.push(point);
        vertices.push(encoder.entity(format_args!("VERTEX_POINT('',{point})"))?);
    }
    let mut curve_poles = encoder.allocate(20)?;
    let mut curve_geometry = encoder.allocate(20)?;
    for curve in plan.curves().curves() {
        let count = 11 * curve.spans().len() + 1;
        let mut poles = encoder.allocate(count)?;
        for i in 0..count {
            let span = (i / 11).min(curve.spans().len() - 1);
            let pole = curve.spans()[span].poles_mm()[i - 11 * span];
            let point = if i == 0 || i + 1 == count {
                let owner = curve.source_curve().endpoints()[usize::from(i != 0)];
                let index = owners
                    .iter()
                    .position(|v| v.owner() == owner)
                    .expect("original vertex owner");
                xyz[index]
            } else {
                encoder.point(pole)?
            };
            poles.push(point);
        }
        curve_geometry.push(spline(encoder, &poles, curve.spans().len())?);
        curve_poles.push(poles);
    }
    let mut supports = encoder.allocate(10)?;
    for patch in plan.patches() {
        let dimensions = patch.spans().map(|n| 11 * n + 1);
        let mut poles = encoder.allocate(patch.poles_mm().len())?;
        let coedges = patch.source_patch().coedges();
        for i in 0..dimensions[0] {
            for j in 0..dimensions[1] {
                let boundary = if i == 0 {
                    Some((0, j))
                } else if i + 1 == dimensions[0] {
                    Some((2, j))
                } else if j == 0 {
                    Some((3, i))
                } else if j + 1 == dimensions[1] {
                    Some((1, i))
                } else {
                    None
                };
                let point = if let Some((side, index)) = boundary {
                    let curve = plan
                        .curves()
                        .curves()
                        .iter()
                        .position(|c| c.source_curve() == coedges[side].curve)
                        .expect("source support edge");
                    curve_poles[curve][index]
                } else {
                    encoder.point(patch.poles_mm()[i * dimensions[1] + j])?
                };
                poles.push(point);
            }
        }
        supports.push(encoder.entity(format_args!("B_SPLINE_SURFACE_WITH_KNOTS('',11,11,({}),.UNSPECIFIED.,.F.,.F.,.F.,({}),({}),({}),({}),.UNSPECIFIED.)",Rows {points:&poles,height:dimensions[1]},Multiplicities(patch.spans()[0]),Multiplicities(patch.spans()[1]),Knots(patch.spans()[0]),Knots(patch.spans()[1])))?);
        encoder.memory.release(&poles);
        drop(poles);
    }
    for station in [StockClipStation::Start, StockClipStation::End] {
        let owner = owners
            .iter()
            .find(|v| v.owner().station == station)
            .expect("true cap source station");
        let origin = encoder.point([0., 0., owner.position_mm()[2]])?;
        let normal = encoder.entity(format_args!("DIRECTION('',(0.,0.,1.))"))?;
        let reference = encoder.entity(format_args!("DIRECTION('',(1.,0.,0.))"))?;
        let axis = encoder.entity(format_args!(
            "AXIS2_PLACEMENT_3D('',{origin},{normal},{reference})"
        ))?;
        supports.push(encoder.entity(format_args!("PLANE('',{axis})"))?);
    }
    let context=encoder.entity(format_args!("(GEOMETRIC_REPRESENTATION_CONTEXT(2)PARAMETRIC_REPRESENTATION_CONTEXT()REPRESENTATION_CONTEXT('2D SPACE',''))"))?;
    let mut pcurves = encoder.allocate(20)?;
    pcurves.resize(20, [Entity(0); 2]);
    let mut uses = encoder.allocate(20)?;
    uses.resize(20, 0usize);
    for (face, patch) in plan.patches().iter().enumerate() {
        for (side, coedge) in patch.source_patch().coedges().into_iter().enumerate() {
            let index = plan
                .curves()
                .curves()
                .iter()
                .position(|c| c.source_curve() == coedge.curve)
                .expect("patch edge owner");
            pcurves[index][uses[index]] = line_pcurve(encoder, supports[face], context, side)?;
            uses[index] += 1;
        }
    }
    for (cap, station) in [StockClipStation::Start, StockClipStation::End]
        .into_iter()
        .enumerate()
    {
        for coedge in source.cap_coedges(station).expect("true cap station") {
            let index = plan
                .curves()
                .curves()
                .iter()
                .position(|c| c.source_curve() == coedge.curve)
                .expect("cap edge owner");
            let curve = &plan.curves().curves()[index];
            let mut points = encoder.allocate(curve_poles[index].len())?;
            for i in 0..curve_poles[index].len() {
                let span = (i / 11).min(curve.spans().len() - 1);
                let [x, y, _] = curve.spans()[span].poles_mm()[i - 11 * span];
                points
                    .push(encoder.entity(format_args!("CARTESIAN_POINT('',({x:.16E},{y:.16E}))"))?);
            }
            let geometry = spline(encoder, &points, curve.spans().len())?;
            encoder.memory.release(&points);
            drop(points);
            pcurves[index][uses[index]] = pcurve(encoder, supports[8 + cap], context, geometry)?;
            uses[index] += 1;
        }
    }
    let mut edges = encoder.allocate(20)?;
    for (i, curve) in plan.curves().curves().iter().enumerate() {
        let geometry = curve_geometry[i];
        let curve_surface = encoder.entity(format_args!(
            "SURFACE_CURVE('',{geometry},{},.PCURVE_S1.)",
            Entities(&pcurves[i])
        ))?;
        let endpoints = curve.source_curve().endpoints().map(|owner| {
            vertices[owners
                .iter()
                .position(|v| v.owner() == owner)
                .expect("edge vertex owner")]
        });
        edges.push(encoder.entity(format_args!(
            "EDGE_CURVE('',{},{},{curve_surface},.T.)",
            endpoints[0], endpoints[1]
        ))?);
    }
    let mut faces = encoder.allocate(10)?;
    for (i, patch) in plan.patches().iter().enumerate() {
        let wire = bound(encoder, &patch.source_patch().coedges(), plan, &edges)?;
        let support = supports[i];
        faces.push(encoder.entity(format_args!("ADVANCED_FACE('',({wire}),{support},.F.)"))?);
    }
    for (i, station) in [StockClipStation::Start, StockClipStation::End]
        .into_iter()
        .enumerate()
    {
        let wire = bound(
            encoder,
            &source.cap_coedges(station).expect("true cap"),
            plan,
            &edges,
        )?;
        let support = supports[8 + i];
        let sense = if i == 0 { ".F." } else { ".T." };
        faces.push(encoder.entity(format_args!("ADVANCED_FACE('',({wire}),{support},{sense})"))?);
    }
    let shell = encoder.entity(format_args!("CLOSED_SHELL('',{})", Entities(&faces)))?;
    encoder.entity(format_args!("MANIFOLD_SOLID_BREP('',{shell})"))
}
