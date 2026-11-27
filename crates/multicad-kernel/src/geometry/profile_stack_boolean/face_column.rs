//! Join and Cut of a profile extrusion into flat faces of a retained solid.
//!
//! A profile sketched on a flat face of a constructed or imported boundary,
//! such as a box whose edges are all filleted, extrudes along the face
//! normal. Where the tool's column, the profile swept over its span and
//! widened by the minimum feature size, meets the target only on flat faces
//! parallel to the profile plane, the target is material over whole slabs
//! of that column:
//! - every other face's exact bounds are proved apart from the column, so
//!   no other face and none of the parallel faces' edges enter it;
//! - a parallel face whose bounds reach the column therefore holds the whole
//!   profile or none of it, and the target's material, classified by the
//!   boundary's certified containment once per slab, changes exactly at the
//!   faces that hold it.
//!
//! The result is the profile stack of the tool against that slab material,
//! where the material outside the profile is the target's own: walls stand
//! where the result differs from the target, new caps close the tool ends,
//! and each face that holds the profile keeps it as a hole wherever the
//! result no longer changes there. Every other target face is kept as it
//! is. The retained solid validates the boundary, and its integrated volume
//! must equal the target's plus the exact slab sections the tool adds or
//! removes. A profile reaching a curved, slanted or stepped face is refused.
use super::*;

const COLUMN: &str = "flat face extrusion Boolean";

fn column_refused(message: impl Into<String>) -> GeometryError {
    error(COLUMN, GeometryErrorKind::Unsupported, message)
}

fn column_invalid(message: impl Into<String>) -> GeometryError {
    error(COLUMN, GeometryErrorKind::InvalidResult, message)
}

/// Whether the closed region inside `rings` stays more than `margin` away
/// from the axis-aligned rectangle `[lo, hi]`.
fn clear_of(rings: &[Vec<Shape>], lo: P2, hi: P2, margin: f64) -> bool {
    let lo = [lo[0] - margin, lo[1] - margin];
    let hi = [hi[0] + margin, hi[1] + margin];
    let inside = |p: P2| (0..2).all(|i| p[i] >= lo[i] && p[i] <= hi[i]);
    let meets = |shape: &Shape| match *shape {
        Shape::Line { a, b } => {
            // Liang-Barsky clipping of the segment against the rectangle.
            let d = arrangement::sub(b, a);
            let (mut t0, mut t1) = (0_f64, 1_f64);
            for i in 0..2 {
                for (p, q) in [(-d[i], a[i] - lo[i]), (d[i], hi[i] - a[i])] {
                    if p == 0. {
                        if q < 0. {
                            return false;
                        }
                    } else {
                        let r = q / p;
                        if p < 0. {
                            t0 = t0.max(r);
                        } else {
                            t1 = t1.min(r);
                        }
                    }
                }
            }
            t0 <= t1
        }
        Shape::Arc {
            center,
            radius,
            start,
            sweep,
        } => {
            if shape.ends().into_iter().any(inside) {
                return true;
            }
            // Otherwise the arc enters the rectangle through one of its
            // sides.
            let on_arc = |angle: f64| {
                let offset = if sweep >= 0. {
                    (angle - start).rem_euclid(std::f64::consts::TAU)
                } else {
                    (start - angle).rem_euclid(std::f64::consts::TAU)
                };
                offset <= sweep.abs()
            };
            for i in 0..2 {
                let j = 1 - i;
                for side in [lo[i], hi[i]] {
                    let offset = side - center[i];
                    if offset.abs() > radius {
                        continue;
                    }
                    let half = (radius * radius - offset * offset).sqrt();
                    for other in [center[j] - half, center[j] + half] {
                        if other < lo[j] || other > hi[j] {
                            continue;
                        }
                        let mut p = [0.; 2];
                        p[i] = side;
                        p[j] = other;
                        let angle = (p[1] - center[1]).portable_atan2(p[0] - center[0]);
                        if on_arc(angle) {
                            return true;
                        }
                    }
                }
            }
            false
        }
    };
    let center = [0.5 * (lo[0] + hi[0]), 0.5 * (lo[1] + hi[1])];
    !rings.iter().flatten().any(meets)
        && !arrangement::contains(rings.iter().flatten().copied(), center)
}

impl Solid {
    /// Join a profile extrusion onto, or cut it into, flat faces of a
    /// retained solid with authored face names.
    pub(in crate::geometry) fn face_column_boolean(
        &self,
        other: &Self,
        operation: BooleanOperation,
        policy: GeometryPolicy,
    ) -> Option<Result<Option<Self>, GeometryError>> {
        if policy.boolean_fuzzy_mm != 0. || policy != self.policy || policy != other.policy {
            return None;
        }
        let (target, tool) = match operation {
            BooleanOperation::Subtract => (self, other),
            BooleanOperation::Union if self.imported_brep_body().is_some() => (self, other),
            BooleanOperation::Union => (other, self),
            BooleanOperation::Intersect => return None,
        };
        let body = target.imported_brep_body()?;
        body.authored_naming()?;
        let frame = tool
            .normal_sweep_body()
            .map(|(_, frame)| frame)
            .or_else(|| tool.axial_cylinder_frame().map(|(_, frame)| frame))?;
        let chart = Chart::new(frame);
        let extrusion = tool
            .sweep_extrusion(&chart, policy)
            .or_else(|| tool.cylinder_extrusion(&chart, policy))?;
        let topology = target.rust_topology()?.ok()?;
        let feature = self.rust_feature()?;
        Some(
            Stack::column(body, topology, extrusion, operation, chart, policy)
                .and_then(|stack| stack.column_solid(body, feature)),
        )
    }
}

impl Stack {
    fn column(
        body: &ImportedBrepBody,
        topology: TopologySnapshot,
        tool: Extrusion,
        operation: BooleanOperation,
        chart: Chart,
        policy: GeometryPolicy,
    ) -> Result<Self, GeometryError> {
        let tolerance = policy.linear_tolerance_mm;
        let margin = policy.minimum_feature_mm.max(16. * tolerance);
        let [a, b] = tool.span;
        let mut segments = vec![];
        let mut owners = vec![];
        for (r, ring) in tool.rings.iter().enumerate() {
            for (c, &shape) in ring.iter().enumerate() {
                segments.push(Segment { shape, operand: 1 });
                owners.push([1, r, c]);
            }
        }
        let arrangement = Arrangement::new(&segments, tolerance).map_err(|cause| match cause {
            Refusal::Contact => column_refused(
                "The sketch region touches itself tangentially or within the modelling tolerance",
            ),
            Refusal::Limit => error(
                COLUMN,
                GeometryErrorKind::Limit,
                "The sketch region has too many curves for an exact extrusion Boolean",
            ),
        })?;
        // Exact face bounds in the chart, where the column runs along z.
        let solid = &body.solid;
        let moved = solid
            .moved(&|p| chart.local(p), &|d| {
                [
                    dot(d, chart.frame.u),
                    dot(d, chart.frame.v),
                    dot(d, chart.n),
                ]
            })
            .map_err(super::super::imported_brep_body::failure)?;
        let bounds = moved
            .face_bounds_mm()
            .map_err(super::super::imported_brep_body::failure)?;
        let magnitude = [a, b]
            .into_iter()
            .chain(bounds.iter().flat_map(|[lo, hi]| [lo[2], hi[2]]))
            .map(f64::abs)
            .fold(1., f64::max);
        let numeric = 64. * f64::EPSILON * magnitude;
        let reach = "The sketch region reaches a curved, slanted or stepped face of the body, or \
                     comes within the minimum feature size of one; keep it inside one flat face \
                     and clear of its fillets and edges";
        // Parallel flat faces whose bounds reach the column, by level.
        let mut faces_at: Vec<(f64, usize)> = vec![];
        for (index, (face, [lo, hi])) in solid.boundary().faces.iter().zip(&bounds).enumerate() {
            if hi[2] < a - margin || lo[2] > b + margin {
                continue;
            }
            if clear_of(&tool.rings, [lo[0], lo[1]], [hi[0], hi[1]], margin) {
                continue;
            }
            let Surface::Plane(frame) = &face.surface else {
                return Err(column_refused(reach));
            };
            if (dot(frame.z, chart.n).abs() - 1.).abs() > 1e-12 {
                return Err(column_refused(reach));
            }
            let z = chart.local(frame.origin)[2];
            if faces_at.iter().any(|&(l, _)| (l - z).abs() <= numeric) {
                return Err(column_refused(
                    "The sketch region crosses an edge between two flat faces of the body; keep \
                     it inside one face",
                ));
            }
            faces_at.push((z, index));
        }
        let mut levels: Vec<f64> = vec![a, b];
        for &(z, _) in &faces_at {
            if !levels.iter().any(|&l| (l - z).abs() <= numeric) {
                levels.push(z);
            }
        }
        levels.sort_by(f64::total_cmp);
        if levels
            .windows(2)
            .any(|w| w[1] - w[0] < policy.minimum_feature_mm)
            || levels[0] < a - numeric
            || levels[levels.len() - 1] > b + numeric
        {
            return Err(column_refused(
                "The extrusion ends closer than the minimum feature size to a face of the body; \
                 end it exactly on that face or clearly away from it",
            ));
        }
        let level = |z: f64| {
            levels
                .iter()
                .position(|&l| (l - z).abs() <= numeric)
                .expect("span end is a level")
        };
        let slabs = levels.len() - 1;
        // The target's material in the column, slab by slab. It is the same
        // throughout a slab's column, so the first probe whose side is
        // certified decides; a probe on an unrelated support surface may not
        // be.
        if arrangement.pieces.is_empty() {
            return Err(column_invalid("The sketch region has no boundary"));
        }
        let classify = |low: f64, high: f64| -> Result<bool, GeometryError> {
            for fraction in [0.5, 0.375, 0.625, 0.25, 0.75] {
                let z = low + fraction * (high - low);
                for piece in &arrangement.pieces {
                    for t in [0.5, 0.25, 0.75] {
                        let point = chart.world(piece.shape.point(t), z);
                        match body.classify_box([point, point], 0.)? {
                            brep::BoxSide::Inside => return Ok(true),
                            brep::BoxSide::Outside => return Ok(false),
                            brep::BoxSide::Unresolved => {}
                        }
                    }
                }
            }
            Err(column_refused(
                "The body's material around the sketch region could not be classified exactly",
            ))
        };
        let inside = (0..slabs)
            .map(|k| classify(levels[k], levels[k + 1]))
            .collect::<Result<Vec<_>, _>>()?;
        let outside = [
            classify(levels[0] - margin, levels[0])?,
            classify(levels[slabs], levels[slabs] + margin)?,
        ];
        // The material changes exactly at faces holding the whole profile,
        // and on their outward side it ends.
        let mut hosts = BTreeMap::new();
        for (l, _) in levels.iter().enumerate() {
            let below = l.checked_sub(1).map_or(outside[0], |k| inside[k]);
            let above = if l < slabs { inside[l] } else { outside[1] };
            let face = faces_at
                .iter()
                .find(|&&(z, _)| level(z) == l)
                .map(|&(_, index)| index);
            if below == above {
                continue;
            }
            let face = face.ok_or_else(|| {
                column_invalid("The body's material changes in the column without a face")
            })?;
            let brep::Face {
                surface: Surface::Plane(frame),
                same_sense,
                ..
            } = &solid.boundary().faces[face]
            else {
                unreachable!("level faces are planar")
            };
            let outward = dot(frame.z, chart.n) * if *same_sense { 1. } else { -1. };
            if (outward > 0.) != below {
                return Err(column_invalid(
                    "A flat face of the body does not bound its material in the column",
                ));
            }
            hosts.insert(l, face);
        }
        let column = Extrusion {
            rings: vec![],
            sides: vec![],
            span: [levels[0], levels[slabs]],
            topology,
        };
        let spans = [[0, slabs], tool.span.map(level)];
        Ok(Self {
            chart,
            arrangement,
            levels,
            spans,
            operation,
            operands: [column, tool],
            owners,
            policy,
            column: Some(Column {
                inside,
                outside,
                hosts,
            }),
        })
    }

    /// The target's boundary with the column's walls, new caps and holes.
    fn column_solid(
        &self,
        body: &ImportedBrepBody,
        feature: EntityId,
    ) -> Result<Option<Solid>, GeometryError> {
        if (0..self.slabs()).all(|k| {
            self.arrangement
                .pieces
                .iter()
                .all(|p| self.wall(p, k).is_none())
        }) {
            return Err(column_refused(
                "The extrusion lies wholly inside or outside the body and would not change it",
            ));
        }
        // The section the result gains (positive) or loses in each slab.
        let delta = self.volume();
        let patches = self.patches()?;
        let hosts: Vec<Option<usize>> = patches.iter().map(|patch| patch.host).collect();
        let (added, names) = self.assemble(patches, feature)?;
        let naming = body
            .authored_naming()
            .ok_or_else(|| column_invalid("The body has no authored face names"))?;
        let target = body.solid.boundary();
        let (nv, ne) = (target.vertices.len(), target.edges.len());
        let mut vertices = target.vertices.clone();
        vertices.extend(added.vertices);
        let mut edges = target.edges.clone();
        edges.extend(added.edges.into_iter().map(|mut edge| {
            edge.vertices = edge.vertices.map(|v| v + nv);
            edge
        }));
        let mut faces = target.faces.clone();
        let mut face_names = naming.faces.clone();
        for ((mut face, name), host) in added.faces.into_iter().zip(names.faces).zip(hosts) {
            for l in &mut face.loops {
                for (e, _) in &mut l.uses {
                    *e += ne;
                }
                if let Some(v) = &mut l.vertex {
                    *v += nv;
                }
            }
            match host {
                Some(host) => faces[host].loops.extend(face.loops),
                None => {
                    faces.push(face);
                    face_names.push(name);
                }
            }
        }
        let boundary = Boundary {
            vertices,
            edges,
            shells: vec![Shell {
                faces: (0..faces.len()).collect(),
                reversed: false,
            }],
            faces,
        };
        if target.shells.len() != 1 {
            return Err(column_refused(
                "The body has inner voids; this extrusion Boolean supports one closed shell",
            ));
        }
        let tolerance = self.policy.linear_tolerance_mm;
        let solid = brep::Solid::new(boundary, tolerance, None).map_err(|cause| match cause {
            brep::Error::Limit => error(
                COLUMN,
                GeometryErrorKind::Limit,
                "The extrusion Boolean result exceeds the retained topology limit",
            ),
            _ => column_refused(
                "The extrusion Boolean result is not one closed solid; it would split into \
                 several pieces or have faces that only touch",
            ),
        })?;
        let expected = body.solid.volume_mm3() + delta;
        let volume = solid.volume_mm3();
        let bound = solid.volume_error_mm3() + body.solid.volume_error_mm3();
        if (volume - expected).abs() > bound.max(1e-9 * expected.abs()) {
            return Err(column_invalid(format!(
                "The extrusion Boolean boundary integrates to {volume} mm3, not the exact {expected} mm3"
            )));
        }
        Solid::from_authored_brep(
            solid,
            feature,
            AuthoredNaming { faces: face_names },
            self.policy,
        )
        .map(Some)
    }
}
