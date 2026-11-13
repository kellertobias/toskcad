//! Transfer actual source residual receipts without gratuitous local refinement.
use super::*;
/// The caller supplies bounds earned by the owned source for these exact facets.
/// These values are never inferred from sampled mesh/source comparisons.
pub(in crate::geometry) struct SourceMesh {
    pub mesh: Mesh,
    pub chord_mm: f64,
    pub normal_rad: f64,
}
pub(in crate::geometry) fn qualify(
    frame: PlaneFrame,
    options: Tessellation,
    build: &mut dyn FnMut(Tessellation) -> Result<SourceMesh, GeometryError>,
) -> Result<Mesh, GeometryError> {
    let axes = [frame.u, frame.v, frame.normal()];
    chart(axes)?;
    let source = build(options)?;
    if !source.chord_mm.is_finite()
        || source.chord_mm < 0.
        || source.chord_mm > options.deflection
        || !source.normal_rad.is_finite()
        || source.normal_rad < 0.
        || source.normal_rad > options.angle
    {
        return Err(refused());
    }
    let mut mesh = source.mesh;
    if frame == PlaneFrame::default() {
        super::super::embedding::certify(&mesh)?;
        return Ok(mesh);
    }
    let points = mesh.vertices.as_chunks::<3>().0;
    let world = points
        .iter()
        .map(|p| placed(frame, axes, *p))
        .collect::<Vec<_>>();
    if world
        .iter()
        .any(|(p, e)| !e.is_finite() || p.iter().any(|x| !x.is_finite()))
    {
        return Err(refused());
    }
    let rounding = world.iter().map(|(_, e)| *e).fold(0., f64::max);
    // chart proves both ||A|| and condition(A) <1.00001 for the unchanged
    // actual affine frame. Point rounding is bounded in L1, hence also in L2.
    if up(product(source.chord_mm, 1.00001), rounding) > options.deflection {
        return Err(refused());
    }
    let transfer_angle = (options.angle - product(source.normal_rad, 1.00001)).next_down();
    if !transfer_angle.is_finite() || transfer_angle <= 0. {
        return Err(refused());
    }
    for ids in mesh.indices.as_chunks::<3>().0 {
        let [i, j, k] = ids.map(|i| i as usize);
        let (a, ae) = difference_with_error(points[j], points[i]);
        let (b, be) = difference_with_error(points[k], points[i]);
        let (n, ne) = cross_error(a, ae, b, be);
        let lower = (0..3)
            .map(|i| (n[i].abs() - ne[i]).next_down().max(0.))
            .fold(0., f64::max);
        let length = |v: Point, e: Point| (0..3).fold(0., |s, i| up(s, up(v[i].abs(), e[i])));
        let e = world[i].1.max(world[j].1).max(world[k].1);
        let perturbation = up(
            product(product(4., e), up(length(a, ae), length(b, be))),
            product(4., product(e, e)),
        );
        if !rounding_angle_fits(perturbation, lower, transfer_angle)
            && !normal_transfer(
                axes,
                a,
                ae,
                b,
                be,
                [world[i].0, world[j].0, world[k].0],
                transfer_angle,
            )
        {
            return Err(refused());
        }
    }
    mesh.vertices = world.into_iter().flat_map(|(p, _)| p).collect();
    super::super::embedding::certify(&mesh)?;
    Ok(mesh)
}

#[cfg(test)]
mod tests;
