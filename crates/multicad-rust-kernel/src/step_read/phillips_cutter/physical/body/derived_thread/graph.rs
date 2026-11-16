//! Exact shared candidate rails and end rows; no coordinate proximity merging.
use super::*;
fn axis(z: f64, sign: f64) -> radial::Axis {
    radial::Axis {
        origin: [0., 0., z],
        normal: [0., 0., sign],
        reference: [1., 0., 0.],
    }
}
pub(super) fn build(
    nets: Vec<Net>,
    stations: [f64; 2],
    work: &mut Work<'_>,
) -> Result<radial::BoundaryInput> {
    let count = nets.len();
    require((3..=32).contains(&count))?;
    let forward = nets[0].forward;
    let mut input = radial::BoundaryInput {
        vertices: Vec::with_capacity(2 * count),
        curves: Vec::with_capacity(3 * count),
        edges: Vec::with_capacity(3 * count),
        surfaces: Vec::with_capacity(count + 2),
        faces: Vec::with_capacity(count + 2),
    };
    for (i, net) in nets.iter().enumerate() {
        work.step(net.poles.len())?;
        let h = 9 * net.spans[1] + 1;
        let next = &nets[(i + 1) % count];
        let nh = 9 * next.spans[1] + 1;
        // The earned common prefix and private Net construction preserve
        // common hand/axial dimensions. Keep the actual NEW literal coefficient
        // equality check: no source closure is inferred from raw tolerances.
        require(net.poles[net.poles.len() - h..] == next.poles[..nh])?;
        input.vertices.extend([net.poles[0], net.poles[h - 1]]);
    }
    // Each original shared rail has exactly one literal derived curve owner.
    // Every surface boundary uses that owner's identical control sequence.
    for (i, net) in nets.iter().enumerate() {
        work.step(1)?;
        let h = 9 * net.spans[1] + 1;
        let w = 9 * net.spans[0] + 1;
        for end in 0..2 {
            let poles = (0..w).map(|u| net.poles[u * h + end * (h - 1)]).collect();
            let curve = input.curves.len();
            input.curves.push(radial::Curve::Spline {
                spans: net.spans[0],
                poles,
            });
            input.edges.push(radial::Edge {
                vertices: [2 * i + end, 2 * ((i + 1) % count) + end],
                curve,
                forward: true,
            });
        }
        let curve = input.curves.len();
        input.curves.push(radial::Curve::Spline {
            spans: net.spans[1],
            poles: net.poles[..h].to_vec(),
        });
        input.edges.push(radial::Edge {
            vertices: [2 * i, 2 * i + 1],
            curve,
            forward: true,
        });
    }
    for (end, &station) in stations.iter().enumerate() {
        let surface = input.surfaces.len();
        input.surfaces.push(radial::Surface::Plane(axis(
            station,
            if end == 0 { -1. } else { 1. },
        )));
        let mut outer = (0..count).map(|i| (3 * i + end, true)).collect::<Vec<_>>();
        if forward != (end == 1) {
            outer = outer.into_iter().rev().map(|(e, s)| (e, !s)).collect();
        }
        input.faces.push(radial::Face {
            surface,
            outer,
            forward: true,
        });
    }
    for (i, net) in nets.into_iter().enumerate() {
        work.step(1)?;
        let next = (i + 1) % count;
        let mut outer = vec![
            (3 * i, true),
            (3 * next + 2, true),
            (3 * i + 1, false),
            (3 * i + 2, false),
        ];
        if !forward {
            outer = outer.into_iter().rev().map(|(e, s)| (e, !s)).collect();
        }
        let surface = input.surfaces.len();
        // Private Net has nonzero spans, exact (9u+1)(9v+1) length,
        // all coordinates finite and abs<1e7, with aggregate262144 controls
        // below RetainedTensor's524288 limit. This repeats an earned invariant,
        // rather than introducing an impossible geometry-error mapping.
        let tensor = RetainedTensor::new(net.spans, net.poles).expect("qualified private Net");
        input.surfaces.push(radial::Surface::Tensor(tensor));
        input.faces.push(radial::Face {
            surface,
            outer,
            forward,
        });
    }
    Ok(input)
}
pub(super) fn bytes(input: &radial::BoundaryInput) -> usize {
    let mut n = std::mem::size_of_val(input)
        + input.vertices.capacity() * std::mem::size_of::<Point>()
        + input.curves.capacity() * std::mem::size_of::<radial::Curve>()
        + input.edges.capacity() * std::mem::size_of::<radial::Edge>()
        + input.surfaces.capacity() * std::mem::size_of::<radial::Surface>()
        + input.faces.capacity() * std::mem::size_of::<radial::Face>();
    for c in &input.curves {
        if let radial::Curve::Spline { poles, .. } = c {
            n = n.saturating_add(poles.capacity() * std::mem::size_of::<Point>());
        }
    }
    for s in &input.surfaces {
        if let radial::Surface::Tensor(t) = s {
            n = n.saturating_add(std::mem::size_of_val(t.poles()) + 64);
        }
    }
    for f in &input.faces {
        n = n.saturating_add(f.outer.capacity() * std::mem::size_of::<(usize, bool)>());
    }
    n
}
