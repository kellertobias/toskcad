//! Equal-distance box bevels with planar, trihedral corner closure.
//! Intersecting selections conservatively retain a positive policy-sized
//! residual after 2d on every stock dimension.
use super::*;
use crate::TreatmentEdge;

fn unsupported() -> GeometryError {
    error(
        "chamfer",
        GeometryErrorKind::Unsupported,
        "Rust intersecting chamfers require complete edges of one rectangular solid with positive 2d residuals",
    )
}
fn checkpoint(interrupt: Option<&dyn Fn() -> bool>) -> Result<(), GeometryError> {
    if interrupt.is_some_and(|f| f()) {
        Err(planar_error(Error::Cancelled))
    } else {
        Ok(())
    }
}
fn coordinate(origin: f64, offset: f64, tolerance: f64) -> Result<f64, GeometryError> {
    let value = origin + offset;
    let part = value - origin;
    let residual = (origin - (value - part)) + (offset - part);
    if !value.is_finite() || residual.abs() > tolerance / 64. {
        return Err(unsupported());
    }
    Ok(value)
}
fn sum(a: f64, b: f64) -> (f64, f64) {
    let value = a + b;
    let part = value - a;
    (value, (a - (value - part)) + (b - part))
}
fn residual_lower(max: f64, min: f64, twice_distance: f64) -> f64 {
    let (length, first) = sum(max, -min);
    let (remaining, second) = sum(length, -twice_distance);
    let error = first.abs() + second.abs();
    if error == 0. {
        remaining
    } else {
        (remaining - error.next_up()).next_down()
    }
}
fn corner_role(mut roles: Vec<&str>) -> String {
    roles.sort_unstable();
    let mut bytes = b"multicad.box-chamfer.corner.v1\0".to_vec();
    for role in roles {
        bytes.extend_from_slice(&(role.len() as u64).to_le_bytes());
        bytes.extend_from_slice(role.as_bytes());
    }
    format!(
        "chamfer.corner/{}",
        multicad_core::project_file::ContentHash::of(&bytes)
    )
}
fn convex(mut faces: Vec<Vec<Point>>, tolerance: f64) -> Result<ConvexPiece, GeometryError> {
    let count = faces.iter().map(Vec::len).sum::<usize>() as f64;
    let middle: Point = std::array::from_fn(|i| faces.iter().flatten().map(|p| p[i] / count).sum());
    for face in &mut faces {
        if dot(
            cross(subtract(face[1], face[0]), subtract(face[2], face[0])),
            subtract(middle, face[0]),
        ) > 0.
        {
            face.reverse();
        }
    }
    ConvexPiece::new(faces, tolerance).map_err(planar_error)
}
fn certificate(
    piece: ConvexPiece,
    tolerance: f64,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<CutMaterial, GeometryError> {
    let pieces = [piece];
    let boundary =
        Union::new_with_interrupt(&pieces, tolerance, interrupt).map_err(planar_error)?;
    CutMaterial::certify_existing_with_interrupt(&boundary, &pieces, tolerance, interrupt)
        .map_err(planar_error)
}

impl Solid {
    pub(crate) fn rust_chamfer_box_edges(
        &self,
        edges: &[TreatmentEdge],
        distance: f64,
        feature: EntityId,
        policy: GeometryPolicy,
        topology: TopologySnapshot,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self, GeometryError> {
        checkpoint(interrupt)?;
        if self.policy != policy
            || policy.boolean_fuzzy_mm != 0.
            || !(2..=12).contains(&edges.len())
        {
            return Err(unsupported());
        }
        if let Some((source, placement)) = self.local_box_treatment_source()? {
            let local = source.topology_snapshot(topology.owner, topology.revision)?;
            return source
                .rust_chamfer_box_edges(edges, distance, feature, policy, local, interrupt)?
                .transformed(placement);
        }
        let cells = self.orthogonal_cells()?.ok_or_else(unsupported)?;
        let [stock] = cells.as_slice() else {
            return Err(unsupported());
        };
        let mut selections = BTreeMap::new();
        for edge in edges {
            let Some(AnalyticGeometry::Line { a_mm: a, b_mm: b }) = topology
                .elements
                .iter()
                .find(|e| e.key == edge.key)
                .and_then(|e| e.geometry)
            else {
                return Err(unsupported());
            };
            let along: Vec<_> = (0..3).filter(|&i| a[i] != b[i]).collect();
            let [along] = along.as_slice() else {
                return Err(unsupported());
            };
            if a[*along].min(b[*along]) != stock.min[*along]
                || a[*along].max(b[*along]) != stock.max[*along]
            {
                return Err(unsupported());
            }
            let mut corner = 0usize;
            for i in 0..3 {
                if i == *along {
                    continue;
                }
                if a[i] != b[i] {
                    return Err(unsupported());
                }
                if a[i] == stock.max[i] {
                    corner |= 1 << i;
                } else if a[i] != stock.min[i] {
                    return Err(unsupported());
                }
            }
            if selections.insert((*along, corner), edge).is_some() {
                return Err(unsupported());
            }
        }
        if edges.len() <= 4
            && selections
                .keys()
                .all(|k| k.0 == selections.keys().next().unwrap().0)
        {
            // The existing constant-size section construction admits a wider
            // parallel-edge envelope; cancellation must not narrow it.
            checkpoint(interrupt)?;
            let result =
                self.rust_chamfer_parallel_box_edges(edges, distance, feature, policy, topology)?;
            checkpoint(interrupt)?;
            return Ok(result);
        }
        policy.length(distance)?;
        for i in 0..3 {
            policy.length(residual_lower(stock.max[i], stock.min[i], 2. * distance))?;
        }
        let tolerance = policy.linear_tolerance_mm;
        let mut certificates = vec![certificate(
            ConvexPiece::from_box(stock.min, stock.max, tolerance).map_err(planar_error)?,
            tolerance,
            interrupt,
        )?];
        let mut origins = vec![None];
        let mut slanted = vec![usize::MAX];
        let mut directions = vec![None];
        let mut edge_lengths = 0.;
        for (&(along, corner), edge) in &selections {
            checkpoint(interrupt)?;
            let [a, b] = [(along + 1) % 3, (along + 2) % 3];
            let mut lower = [stock.min; 3];
            for p in &mut lower {
                for axis in [a, b] {
                    p[axis] = if corner & (1 << axis) != 0 {
                        stock.max[axis]
                    } else {
                        stock.min[axis]
                    };
                }
            }
            lower[1][a] = coordinate(
                lower[1][a],
                if corner & (1 << a) == 0 {
                    distance
                } else {
                    -distance
                },
                tolerance,
            )?;
            lower[2][b] = coordinate(
                lower[2][b],
                if corner & (1 << b) == 0 {
                    distance
                } else {
                    -distance
                },
                tolerance,
            )?;
            let upper = lower.map(|mut p| {
                p[along] = stock.max[along];
                p
            });
            let mut faces = vec![lower.to_vec(), upper.to_vec()];
            faces.extend((0..3).map(|i| {
                let j = (i + 1) % 3;
                vec![lower[i], lower[j], upper[j], upper[i]]
            }));
            certificates.push(certificate(
                convex(faces, tolerance)?,
                tolerance,
                interrupt,
            )?);
            origins.push(Some(edge.origin.clone()));
            slanted.push(3);
            directions.push(Some(along));
            edge_lengths += stock.max[along] - stock.min[along];
        }
        let mut correction = 0.;
        for corner in 0usize..8 {
            let incident: Vec<_> = (0..3)
                .filter_map(|along| selections.get(&(along, corner & !(1 << along))).copied())
                .collect();
            if incident.len() == 2 {
                correction += distance.powi(3) / 3.;
            }
            if incident.len() != 3 {
                continue;
            }
            correction += 2. * distance.powi(3) / 3.;
            let point: Point = std::array::from_fn(|i| {
                if corner & (1 << i) == 0 {
                    stock.min[i]
                } else {
                    stock.max[i]
                }
            });
            let mut v = [point; 4];
            for i in 0..3 {
                v[i + 1][i] = coordinate(
                    point[i],
                    if corner & (1 << i) == 0 {
                        2. * distance
                    } else {
                        -2. * distance
                    },
                    tolerance,
                )?;
            }
            let faces = vec![
                vec![v[1], v[2], v[3]],
                vec![v[0], v[2], v[1]],
                vec![v[0], v[1], v[3]],
                vec![v[0], v[3], v[2]],
            ];
            certificates.push(certificate(
                convex(faces, tolerance)?,
                tolerance,
                interrupt,
            )?);
            origins.push(Some(Origin {
                source: feature,
                role: corner_role(incident.iter().map(|e| e.origin.role.as_str()).collect()),
            }));
            slanted.push(0);
            directions.push(None);
        }
        let operands: Vec<_> = certificates.iter().collect();
        let material =
            CutMaterial::subtract_certified_with_interrupt(&operands, tolerance, interrupt)
                .map_err(planar_error)?
                .ok_or_else(unsupported)?;
        let geometry = material.boundary().clone();
        let expected = (0..3).map(|i| stock.max[i] - stock.min[i]).product::<f64>()
            - 0.5 * distance * distance * edge_lengths
            + correction;
        if geometry.solid_count() != 1
            || geometry.shell_count() != 1
            || (geometry.volume_mm3() - expected).abs() > expected * 1e-8
        {
            return Err(error(
                "chamfer",
                GeometryErrorKind::InvalidResult,
                "Box chamfer material proof disagrees with the independent corner volume",
            ));
        }
        // This provider proves a bijection of its finite convex faces to
        // the six original stock planes and each generated bevel/corner plane.
        // These are authored generated surfaces, not anonymous Boolean fragments.
        let mut counts = vec![0usize; certificates.len()];
        let mut stock_faces = BTreeSet::new();
        let mut labels = Vec::with_capacity(geometry.faces().len());
        let fallback =
            imported_identity::labels_with_prefix(&geometry, feature, "box-chamfer.finite.v1.")?;
        for (ordinal, face) in geometry.faces().iter().enumerate() {
            checkpoint(interrupt)?;
            let source = certificates.get(face.source).ok_or_else(unsupported)?;
            let original = source
                .boundary()
                .faces()
                .get(face.source_face)
                .ok_or_else(unsupported)?;
            counts[face.source] += 1;
            if face.source == 0 {
                if !stock_faces.insert(face.source_face) {
                    return Err(unsupported());
                }
                let point = geometry.vertices()[face.vertices[0]];
                let mut matches = topology.elements.iter().filter(|e| {
                    matches!(e.geometry,
                    Some(AnalyticGeometry::Plane { frame }) if frame.normal()==face.normal
                    && dot(face.normal,subtract(point,frame.origin_mm))==0.)
                });
                let inherited = matches.next().ok_or_else(unsupported)?;
                if matches.next().is_some() {
                    return Err(unsupported());
                }
                labels.push(if inherited.origins.is_empty() {
                    fallback[ordinal].clone()
                } else {
                    (
                        inherited.origins.iter().cloned().collect(),
                        inherited.frames.clone(),
                    )
                });
            } else {
                if original.source_face != slanted[face.source] {
                    return Err(unsupported());
                }
                let origin = origins[face.source].clone().ok_or_else(unsupported)?;
                let point = source.boundary().vertices()[original.vertices[0]];
                let u = if let Some(along) = directions[face.source] {
                    let mut u = [0.; 3];
                    u[along] = 1.;
                    u
                } else {
                    let raw = cross(face.normal, [1., 0., 0.]);
                    let length = dot(raw, raw).sqrt();
                    raw.map(|x| x / length)
                };
                let frame = ReferenceFrame {
                    origin_mm: point,
                    u,
                    v: cross(face.normal, u),
                };
                frame.validate().map_err(|_| unsupported())?;
                labels.push((vec![origin.clone()], vec![AnchoredFrame { origin, frame }]));
            }
        }
        if counts[0] != 6 || counts[1..].iter().any(|&n| n != 1) {
            return Err(unsupported());
        }
        checkpoint(interrupt)?;
        PlanarBooleanBody::build_labeled(
            geometry,
            feature,
            policy,
            Some(Material::Certified(Box::new(material))),
            PlaneFrame::default(),
            |i, _, _| labels[i].clone(),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn rounded_width_cannot_admit_a_subminimum_residual() {
        let max = 1e7;
        let min = -1.0001e-6;
        let required = 1.0002e-6;
        assert!((max - min) - 1e7 > required);
        assert!(residual_lower(max, min, 1e7) < required);
        assert_eq!(residual_lower(10., 0., 2.), 8.);
        assert_eq!(residual_lower(10., 0., 10.), 0.);
    }
    #[test]
    fn authored_role_bytes_have_unambiguous_corner_ownership() {
        let first = vec!["a\0b", "c", "d"];
        let other = vec!["a", "b\0c", "d"];
        assert_eq!(first.join("\0"), other.join("\0"));
        assert_ne!(corner_role(first.clone()), corner_role(other));
        assert_eq!(corner_role(first), corner_role(vec!["d", "c", "a\0b"]));
    }
}
