//! Exact smooth lofts: straight cap edges, C¹ cubic B-spline rails and
//! degree (1, 3) B-spline walls that share the rail poles.
use super::*;
use crate::planar_boolean::SmoothLoft;

pub fn write_smooth_loft(
    loft: &SmoothLoft,
    origin: Point,
    u: Point,
    v: Point,
) -> Result<String, StepError> {
    write_one(|writer| writer.append_smooth_loft(loft, origin, u, v))
}

fn list(values: impl IntoIterator<Item = String>) -> String {
    values.into_iter().collect::<Vec<_>>().join(",")
}

pub(super) fn encode(
    file: &mut StepFile,
    placement: &Placement,
    loft: &SmoothLoft,
) -> Result<Entity, StepError> {
    let count = loft.count();
    let points = loft
        .rings()
        .iter()
        .flatten()
        .map(|&[x, y, z]| placement.point_f64([x, y], z))
        .collect();
    let mut brep = BRep::new(file, points);
    let stations = loft.heights().len();
    let knots = list(loft.heights().iter().map(|knot| format!("{knot:.16E}")));
    let multiplicities = list((0..stations).map(|i| {
        if i == 0 || i + 1 == stations {
            "4"
        } else {
            "2"
        }
        .to_owned()
    }));
    let heights = loft.pole_heights();
    // Rails share their end poles with the cap vertices.
    let rails: Vec<Vec<Entity>> = loft
        .rails()
        .iter()
        .enumerate()
        .map(|(corner, rail)| {
            let last = rail.len() - 1;
            rail.iter()
                .zip(&heights)
                .enumerate()
                .map(|(n, (&point, &z))| match n {
                    0 => brep.vertices[corner].cartesian,
                    n if n == last => brep.vertices[corner + count].cartesian,
                    _ => brep.file.point(placement.point_f64(point, z)),
                })
                .collect()
        })
        .collect();
    for (corner, poles) in rails.iter().enumerate() {
        let curve = brep.file.entity(format_args!(
            "B_SPLINE_CURVE_WITH_KNOTS('',3,{},.UNSPECIFIED.,.F.,.F.,({multiplicities}),({knots}),.UNSPECIFIED.)",
            Entities(poles)
        ));
        let [start, end] = [corner, corner + count].map(|i| brep.vertices[i].vertex);
        let edge = brep
            .file
            .entity(format_args!("EDGE_CURVE('',{start},{end},{curve},.T.)"));
        brep.edges.insert([corner, corner + count], edge);
    }
    let mut faces = Vec::with_capacity(count + 2);
    for (cap, ring) in loft.rings().iter().enumerate() {
        let offset = cap * count;
        let mut indices: Vec<_> = (offset..offset + count).collect();
        if (cap == 0) != loft.clockwise() {
            indices.reverse();
        }
        let bound = brep.bound(&indices, true)?;
        let reference = placement.vector(direction(ring[0], ring[1])?);
        let normal = placement.vector([0., 0., if cap == 0 { -1. } else { 1. }]);
        faces.push(
            brep.file
                .face(&[bound], brep.vertices[offset].cartesian, normal, reference),
        );
    }
    for i in 0..count {
        let j = (i + 1) % count;
        let mut indices = [i, j, j + count, i + count];
        if loft.clockwise() {
            indices.reverse();
        }
        let bound = brep.bound(&indices, true)?;
        // Rows run along the section edge. Columns are the rail poles.
        let surface = brep.file.entity(format_args!(
            "B_SPLINE_SURFACE_WITH_KNOTS('',1,3,({},{}),.UNSPECIFIED.,.F.,.F.,.F.,(2,2),({multiplicities}),(0.,1.),({knots}),.UNSPECIFIED.)",
            Entities(&rails[i]),
            Entities(&rails[j])
        ));
        faces.push(brep.file.face_on(&[bound], surface, !loft.clockwise()));
    }
    Ok(brep.file.solid(&faces))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encodes_cubic_rails_and_walls_with_shared_poles_in_both_windings() {
        let rect =
            |w: f64, h: f64, [x, y]: [f64; 2]| vec![[x, y], [x + w, y], [x + w, y + h], [x, y + h]];
        for clockwise in [false, true] {
            let mut sections = vec![
                rect(10., 6., [0., 0.]),
                rect(14., 8., [-2., -1.]),
                rect(10., 6., [0., 0.]),
            ];
            if clockwise {
                sections.iter_mut().for_each(|ring| ring.reverse());
            }
            let loft = SmoothLoft::interpolate(&sections, &[0., 15., 30.], 1e-7).unwrap();
            let step =
                write_smooth_loft(&loft, [1., -2., 5.], [1., 0., 0.], [0., 0.6, 0.8]).unwrap();
            assert_eq!(
                step.matches("=B_SPLINE_SURFACE_WITH_KNOTS('',1,3,").count(),
                4
            );
            assert_eq!(step.matches("=B_SPLINE_CURVE_WITH_KNOTS('',3,").count(), 4);
            assert_eq!(step.matches("=EDGE_CURVE(").count(), 12);
            assert_eq!(step.matches("=ORIENTED_EDGE(").count(), 24);
            assert_eq!(step.matches("=VERTEX_POINT(").count(), 8);
            assert_eq!(step.matches("=PLANE(").count(), 2);
            // Every wall row is exactly the pole list of one rail.
            let poles = |line: &str, prefix: &str| -> Vec<String> {
                let rest = &line[line.find(prefix).unwrap() + prefix.len()..];
                let close = rest.find("),").unwrap();
                let list = rest[..close].trim_start_matches('(');
                vec![list.to_owned()]
            };
            let rails: Vec<_> = step
                .lines()
                .filter(|line| line.contains("=B_SPLINE_CURVE_WITH_KNOTS('',3,"))
                .flat_map(|line| poles(line, "=B_SPLINE_CURVE_WITH_KNOTS('',3,"))
                .collect();
            assert_eq!(rails.len(), 4);
            assert!(rails.iter().all(|rail| rail.split(',').count() == 6));
            for line in step
                .lines()
                .filter(|line| line.contains("=B_SPLINE_SURFACE_WITH_KNOTS('',1,3,"))
            {
                let rest = &line[line.find("1,3,((").unwrap() + 6..];
                let net = &rest[..rest.find("))").unwrap()];
                let rows: Vec<_> = net.split("),(").collect();
                assert_eq!(rows.len(), 2);
                assert!(rows.iter().all(|row| rails.iter().any(|rail| rail == row)));
                assert_ne!(rows[0], rows[1]);
            }
        }
    }
}
