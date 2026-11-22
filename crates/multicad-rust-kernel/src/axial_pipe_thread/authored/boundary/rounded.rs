//! A complete five-band source graph: no analytic root cone is fabricated.
use super::*;
#[cfg(test)]
mod tests;
pub(super) fn build(
    thread: &AuthoredThread,
    nets: Vec<Net>,
    tensor_error: f64,
    tolerance: f64,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<RawAuthoredBoundary, Error> {
    let count = nets.len();
    let mut graph = Graph(BoundaryInput {
        vertices: Vec::new(),
        curves: Vec::new(),
        edges: Vec::new(),
        surfaces: Vec::new(),
        faces: Vec::new(),
    });
    let mut corners = Vec::<[usize; 2]>::new();
    let mut rails = Vec::new();
    let mut rims = Vec::<[usize; 2]>::new();
    let mut surfaces = Vec::new();
    let mut first = Vec::new();
    let mut previous = Vec::new();
    for (band, net) in nets.into_iter().enumerate() {
        if interrupt.is_some_and(|f| f()) {
            return Err(Error::Cancelled);
        }
        let [w, h] = net.tensor.spans().map(|s| DEGREE * s + 1);
        let points = net.tensor.poles();
        if band == 0 {
            first = points[..h].to_vec();
            let vertices = [graph.vertex(points[0]), graph.vertex(points[h - 1])];
            corners.push(vertices);
            rails.push(graph.edge(
                vertices,
                Curve::Spline {
                    spans: net.tensor.spans()[1],
                    poles: first.clone(),
                },
            ));
        } else {
            identical(&previous[..], &points[..h])?;
        }
        let right_rail = &points[(w - 1) * h..];
        let right = if band + 1 == count {
            identical(&first[..], right_rail)?;
            corners[0]
        } else {
            let vertices = [graph.vertex(right_rail[0]), graph.vertex(right_rail[h - 1])];
            corners.push(vertices);
            rails.push(graph.edge(
                vertices,
                Curve::Spline {
                    spans: net.tensor.spans()[1],
                    poles: right_rail.to_vec(),
                },
            ));
            vertices
        };
        let left = corners[band];
        rims.push(std::array::from_fn(|end| {
            graph.edge(
                [left[end], right[end]],
                Curve::Spline {
                    spans: net.tensor.spans()[0],
                    poles: (0..w).map(|i| points[i * h + end * (h - 1)]).collect(),
                },
            )
        }));
        previous = right_rail.to_vec();
        surfaces.push(graph.surface(Surface::Tensor(net.tensor)));
    }
    for end in 0..2 {
        let z = if end == 0 { 0. } else { thread.spec.length_mm };
        let surface = graph.surface(Surface::Plane(Axis {
            origin: [0., 0., z],
            normal: [0., 0., if end == 0 { -1. } else { 1. }],
            reference: [1., 0., 0.],
        }));
        graph.face(
            surface,
            rims.iter().map(|r| (r[end], true)).collect(),
            end == 1,
            true,
        );
    }
    for band in 0..count {
        graph.face(
            surfaces[band],
            vec![
                (rims[band][0], true),
                (rails[(band + 1) % count], true),
                (rims[band][1], false),
                (rails[band], false),
            ],
            true,
            false,
        );
    }
    let error = (I::point(2.) * I::point(tensor_error)).checked()?.hi;
    error_budget(error, tolerance)?;
    Ok(RawAuthoredBoundary {
        ownership: BoundaryOwnership::new(graph.0, tolerance, interrupt)?,
        source_error_mm: error,
        tensor_error_mm: tensor_error,
        analytic_error_mm: 0.,
        cap_faces: [0, 1],
        band_faces: (2..count + 2).collect(),
        root_faces: vec![4],
    })
}
