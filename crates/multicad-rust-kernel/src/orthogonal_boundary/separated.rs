//! Several outward orthogonal shells that are separate solids, not cavities.
//!
//! A boundary with more than one shell is either material with voids or a
//! group of independent solids. Each shell here is closed and coherently
//! oriented (the boundary graph pairs every edge once in each direction), has
//! a certified positive enclosed volume, and lies in an enclosing box that is
//! strictly apart from every other shell's box. Strictly apart boxes cannot
//! nest or overlap, so every shell bounds its own solid.
use super::*;
use crate::orthogonal::BooleanError;

type ProofResult<T> = std::result::Result<T, BooleanError>;

/// The face indices of each separate solid, in face order.
pub(crate) fn separated_solids(boundary: &Boundary) -> ProofResult<Vec<Vec<usize>>> {
    box_void::budget(
        boundary.faces.len(),
        boundary
            .faces
            .iter()
            .flat_map(|face| face.rings.iter().map(Vec::len)),
    )?;
    if boundary.shell_count < 2 {
        return Err(BooleanError::Invalid);
    }
    let groups = box_void::components(
        boundary
            .faces
            .iter()
            .map(|face| (face.component, face.rings.as_slice())),
    )
    .ok_or(BooleanError::Invalid)?;
    if groups.len() != boundary.shell_count {
        return Err(BooleanError::Invalid);
    }
    let boxes = groups
        .iter()
        .map(|group| enclosing_box(boundary, group))
        .collect::<ProofResult<Vec<_>>>()?;
    for (i, a) in boxes.iter().enumerate() {
        for b in &boxes[i + 1..] {
            if !apart(*a, *b) {
                return Err(BooleanError::Invalid);
            }
        }
    }
    for (group, [origin, _]) in groups.iter().zip(&boxes) {
        if !positive_volume(boundary, group, *origin) {
            return Err(BooleanError::Invalid);
        }
    }
    Ok(groups)
}

impl Boundary {
    /// One separate shell of this boundary as a single-shell boundary.
    pub(crate) fn shell(&self, faces: &[usize]) -> Self {
        Self {
            faces: faces
                .iter()
                .map(|&index| {
                    let face = &self.faces[index];
                    BoundaryFace {
                        component: face.component,
                        rings: face.rings.clone(),
                        quad_indices: face.quad_indices.clone(),
                        normal: face.normal,
                    }
                })
                .collect(),
            shell_count: 1,
        }
    }
}

fn enclosing_box(boundary: &Boundary, group: &[usize]) -> ProofResult<[Point; 2]> {
    let mut points = group
        .iter()
        .flat_map(|&face| &boundary.faces[face].rings)
        .flatten();
    let first = *points.next().ok_or(BooleanError::Invalid)?;
    Ok(points.fold([first, first], |[min, max], point| {
        [
            std::array::from_fn(|i| min[i].min(point[i])),
            std::array::from_fn(|i| max[i].max(point[i])),
        ]
    }))
}

/// Exact comparison: some axis separates the closed boxes by a gap.
fn apart([a_min, a_max]: [Point; 2], [b_min, b_max]: [Point; 2]) -> bool {
    (0..3).any(|i| b_min[i] > a_max[i] || a_min[i] > b_max[i])
}

/// The divergence-theorem volume of the outward shell, relative to a corner
/// of its box, exceeds a rigorous bound on its floating-point rounding.
fn positive_volume(boundary: &Boundary, group: &[usize], origin: Point) -> bool {
    let mut sum = 0.;
    let mut magnitude = 0.;
    let mut terms = 0usize;
    for &face in group {
        for ring in &boundary.faces[face].rings {
            let local: Vec<Point> = ring
                .iter()
                .map(|point| std::array::from_fn(|i| point[i] - origin[i]))
                .collect();
            for pair in local[1..].windows(2) {
                let [a, b, c] = [local[0], pair[0], pair[1]];
                let products = [
                    a[0] * b[1] * c[2],
                    -a[0] * b[2] * c[1],
                    a[1] * b[2] * c[0],
                    -a[1] * b[0] * c[2],
                    a[2] * b[0] * c[1],
                    -a[2] * b[1] * c[0],
                ];
                for product in products {
                    sum += product;
                    magnitude += product.abs();
                    terms += 1;
                }
            }
        }
    }
    // Every subtraction from the origin, each product and each addition
    // rounds once; (terms + 8) ulps of the absolute sum bounds them all.
    let rounding = (terms as f64 + 8.) * f64::EPSILON * magnitude;
    sum.is_finite() && magnitude.is_finite() && sum > 2. * rounding
}

#[cfg(test)]
mod tests {
    use super::*;

    /// An outward axis-aligned box as six single-ring faces of one shell.
    fn cube(component: u32, min: Point, max: Point) -> Vec<BoundaryFace> {
        let corner = |x: usize, y: usize, z: usize| -> Point {
            [
                [min[0], max[0]][x],
                [min[1], max[1]][y],
                [min[2], max[2]][z],
            ]
        };
        let faces = [
            (
                [
                    corner(0, 0, 0),
                    corner(0, 1, 0),
                    corner(1, 1, 0),
                    corner(1, 0, 0),
                ],
                [0., 0., -1.],
            ),
            (
                [
                    corner(0, 0, 1),
                    corner(1, 0, 1),
                    corner(1, 1, 1),
                    corner(0, 1, 1),
                ],
                [0., 0., 1.],
            ),
            (
                [
                    corner(0, 0, 0),
                    corner(1, 0, 0),
                    corner(1, 0, 1),
                    corner(0, 0, 1),
                ],
                [0., -1., 0.],
            ),
            (
                [
                    corner(0, 1, 0),
                    corner(0, 1, 1),
                    corner(1, 1, 1),
                    corner(1, 1, 0),
                ],
                [0., 1., 0.],
            ),
            (
                [
                    corner(0, 0, 0),
                    corner(0, 0, 1),
                    corner(0, 1, 1),
                    corner(0, 1, 0),
                ],
                [-1., 0., 0.],
            ),
            (
                [
                    corner(1, 0, 0),
                    corner(1, 1, 0),
                    corner(1, 1, 1),
                    corner(1, 0, 1),
                ],
                [1., 0., 0.],
            ),
        ];
        faces
            .into_iter()
            .map(|(ring, normal)| BoundaryFace {
                component,
                rings: vec![ring.to_vec()],
                quad_indices: Vec::new(),
                normal,
            })
            .collect()
    }

    fn reversed(mut faces: Vec<BoundaryFace>) -> Vec<BoundaryFace> {
        for face in &mut faces {
            face.rings[0].reverse();
            face.normal = face.normal.map(|x| -x);
        }
        faces
    }

    fn boundary(shells: Vec<Vec<BoundaryFace>>) -> Boundary {
        Boundary {
            shell_count: shells.len(),
            faces: shells.into_iter().flatten().collect(),
        }
    }

    #[test]
    fn strictly_apart_outward_shells_are_separate_solids() {
        let two = boundary(vec![
            cube(0, [0.; 3], [10., 20., 30.]),
            cube(1, [100., 0., 0.], [105., 5., 5.]),
        ]);
        let groups = separated_solids(&two).unwrap();
        assert_eq!(groups, vec![(0..6).collect::<Vec<_>>(), (6..12).collect()]);
        assert!(two.supports_step());
        let shell = two.shell(&groups[1]);
        assert_eq!(shell.shell_count, 1);
        assert_eq!(shell.faces.len(), 6);
        assert!(shell.supports_step());
    }

    #[test]
    fn cavities_touching_or_overlapping_boxes_and_single_shells_are_refused() {
        let invalid = |shells| {
            assert_eq!(
                separated_solids(&boundary(shells)).err(),
                Some(BooleanError::Invalid)
            )
        };
        // An inward shell inside the material is a void, not a solid.
        invalid(vec![
            cube(0, [0.; 3], [10.; 3]),
            reversed(cube(1, [2.; 3], [3.; 3])),
        ]);
        // An inward shell far away has negative volume.
        invalid(vec![
            cube(0, [0.; 3], [1.; 3]),
            reversed(cube(1, [5.; 3], [6.; 3])),
        ]);
        // Boxes that share a face plane are not strictly apart.
        invalid(vec![
            cube(0, [0.; 3], [1.; 3]),
            cube(1, [1., 0., 0.], [2., 1., 1.]),
        ]);
        // Nested outward shells overlap.
        invalid(vec![cube(0, [0.; 3], [10.; 3]), cube(1, [2.; 3], [3.; 3])]);
        invalid(vec![cube(0, [0.; 3], [1.; 3])]);
        // A miscounted shell total is refused.
        let mut miscounted = boundary(vec![cube(0, [0.; 3], [1.; 3]), cube(1, [5.; 3], [6.; 3])]);
        miscounted.shell_count = 3;
        assert_eq!(
            separated_solids(&miscounted).err(),
            Some(BooleanError::Invalid)
        );
    }
}
