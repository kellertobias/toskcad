//! Bounds for the actual stored affine columns, including determinant/cofactor.
use super::*;
use crate::geometry::internal_thread_treatment_body::interval::Interval as I;
use multicad_rust_kernel::axial_pipe_thread::authored::StockClipPlacement;
#[derive(Clone, Debug)]
pub struct OrthogonalCountersinkMeasures {
    pub volume_mm3: [f64; 2],
    pub face_area_mm2: Vec<[f64; 2]>,
    pub face_centroid_mm: Vec<[[f64; 2]; 3]>,
    pub maximum_centroid_error_mm: f64,
    pub determinant: [f64; 2],
    pub area_scale: [f64; 2],
}
fn interval(a: [f64; 2]) -> I {
    I::new(a[0], a[1])
}
fn bounds(a: I) -> [f64; 2] {
    [a.lo, a.hi]
}
fn placement_error(c: multicad_rust_kernel::axial_pipe_thread::Error) -> GeometryError {
    use multicad_rust_kernel::axial_pipe_thread::Error as E;
    error(
        "Orthogonal countersink placement",
        match c {
            E::Invalid => GeometryErrorKind::Input,
            E::Limit => GeometryErrorKind::Limit,
            E::Cancelled => GeometryErrorKind::Cancelled,
            _ => GeometryErrorKind::Unsupported,
        },
        format!("Literal stored map proof failed: {c:?}"),
    )
}
pub(super) fn placement(frame: PlaneFrame) -> Result<StockClipPlacement, GeometryError> {
    let axes = [frame.u, frame.v, frame.normal()];
    StockClipPlacement::new(std::array::from_fn(|i| axes.map(|a| a[i])), frame.origin)
        .map_err(placement_error)
}
fn cross(a: [I; 3], b: [I; 3]) -> [I; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn normal(frame: PlaneFrame, n: Point) -> Result<Point, GeometryError> {
    let [a, b, c] = [frame.u, frame.v, frame.normal()].map(|v| v.map(I::point));
    let cofactors = [cross(b, c), cross(c, a), cross(a, b)];
    let mapped: [I; 3] =
        std::array::from_fn(|i| (0..3).fold(I::ZERO, |s, j| s + cofactors[j][i] * I::point(n[j])));
    let length = mapped.iter().fold(I::ZERO, |s, v| s + v.square()).sqrt();
    if !length.valid() || length.lo <= 0. {
        return Err(unsupported("Degenerate actual cofactor normal"));
    }
    Ok(mapped.map(|v| {
        let q = v / length;
        q.lo.midpoint(q.hi)
    }))
}
fn face_scale(frame: PlaneFrame, n: Point) -> I {
    let [a, b, c] = [frame.u, frame.v, frame.normal()].map(|v| v.map(I::point));
    let cofactor = [cross(b, c), cross(c, a), cross(a, b)];
    std::array::from_fn::<_, 3, _>(|i| {
        (0..3).fold(I::ZERO, |s, j| s + cofactor[j][i] * I::point(n[j]))
    })
    .iter()
    .fold(I::ZERO, |s, v| s + v.square())
    .sqrt()
}
impl OrthogonalSeatBody {
    pub(in crate::geometry) fn validate_frame(
        &self,
        frame: PlaneFrame,
        policy: GeometryPolicy,
    ) -> Result<(), GeometryError> {
        frame.checked_metric_bounds(self.original.bounds(PlaneFrame::default()), policy)?;
        let map = placement(frame)?;
        let local = self.source.stock().bounds();
        let world = map
            .point_bounds(std::array::from_fn(|i| [local.min[i], local.max[i]]))
            .map_err(placement_error)?;
        policy.point(world.map(|v| v[0]))?;
        policy.point(world.map(|v| v[1]))?;
        if self.measures(frame)?.maximum_centroid_error_mm > policy.linear_tolerance_mm / 4. {
            return Err(unsupported(
                "Actual transported surface centroid uncertainty exceeds the precision budget",
            ));
        }
        Ok(())
    }
    pub(super) fn world_bounds(&self, frame: PlaneFrame) -> Result<Bounds, GeometryError> {
        let local = self.source.stock().bounds();
        let b = placement(frame)?
            .point_bounds(std::array::from_fn(|i| [local.min[i], local.max[i]]))
            .map_err(placement_error)?;
        Ok(Bounds {
            min: b.map(|v| v[0]),
            max: b.map(|v| v[1]),
        })
    }
    pub(super) fn world_volume_bounds(&self, frame: PlaneFrame) -> Result<[f64; 2], GeometryError> {
        Ok(bounds(
            interval(self.source.volume_bounds_mm3())
                * interval(placement(frame)?.determinant_bounds()),
        ))
    }
    pub(super) fn measures(
        &self,
        frame: PlaneFrame,
    ) -> Result<OrthogonalCountersinkMeasures, GeometryError> {
        let map = placement(frame)?;
        let determinant = map.determinant_bounds();
        let area_scale = map.area_scale_bounds();
        let pi = I::new(PI.next_down(), PI.next_up());
        let mut areas = Vec::new();
        for face in &self.source.boundary().faces {
            let axis = (0..3)
                .find(|&i| face.normal[i] != 0.)
                .expect("qualified original face");
            let [u, v] = [(axis + 1) % 3, (axis + 2) % 3];
            let mut area = I::ZERO;
            for ring in &face.rings {
                for (a, b) in ring.iter().zip(ring.iter().cycle().skip(1)) {
                    area = area + I::point(a[u]) * I::point(b[v]) - I::point(a[v]) * I::point(b[u]);
                }
            }
            areas.push(area * I::point(face.normal[axis]) / I::point(2.));
        }
        for (i, bore) in self.source.bores().iter().enumerate() {
            let radius = I::point(bore.radius_mm);
            for face in bore.faces {
                areas[face] = areas[face] - pi * radius * radius;
            }
            let axis = (0..3)
                .find(|&j| bore.centers_mm[0][j] != bore.centers_mm[1][j])
                .expect("qualified original bore");
            let [a, b] = bore.centers_mm.map(|p| p[axis]);
            let mut length = if a < b {
                I::point(b) - I::point(a)
            } else {
                I::point(a) - I::point(b)
            };
            for s in self.source.seats().iter().filter(|s| s.bore == i) {
                length = length - I::point(s.depth_mm);
            }
            areas.push(I::point(2.) * pi * radius * length);
        }
        for s in self.source.seats() {
            let r = I::point(self.source.bores()[s.bore].radius_mm);
            let outer = I::point(s.radius_mm);
            let d = I::point(s.depth_mm);
            let cap = self.source.bores()[s.bore].faces[usize::from(s.end)];
            areas[cap] = areas[cap] - pi * (outer - r) * (outer + r);
            areas.push(pi * (outer + r) * (d.square() + (outer - r).square()).sqrt());
        }
        let face_area_mm2 = areas
            .into_iter()
            .enumerate()
            .map(|(i, a)| {
                let scale = if let Some(face) = self.source.boundary().faces.get(i) {
                    face_scale(frame, face.normal)
                } else {
                    interval(area_scale)
                };
                let area = a * scale;
                if !area.valid() || area.lo <= 0. {
                    Err(unsupported("Actual transported face area is not resolved"))
                } else {
                    Ok(bounds(area))
                }
            })
            .collect::<Result<_, _>>()?;
        let volume = interval(self.source.volume_bounds_mm3()) * interval(determinant);
        let face_centroid_mm = self.centroid_bounds(frame, map.area_scale_bounds())?;
        let maximum_centroid_error_mm = face_centroid_mm
            .iter()
            .map(|b| {
                b.iter()
                    .fold(I::ZERO, |sum, v| {
                        let mid = v[0].midpoint(v[1]);
                        sum + I::point((interval(*v) - I::point(mid)).magnitude())
                    })
                    .hi
            })
            .fold(0_f64, f64::max);
        Ok(OrthogonalCountersinkMeasures {
            volume_mm3: bounds(volume),
            face_area_mm2,
            face_centroid_mm,
            maximum_centroid_error_mm,
            determinant,
            area_scale,
        })
    }
    fn local_centroids(&self) -> Result<Vec<[I; 3]>, GeometryError> {
        let pi = I::new(PI.next_down(), PI.next_up());
        let mut centers = Vec::new();
        for (index, face) in self.source.boundary().faces.iter().enumerate() {
            let axis = (0..3)
                .find(|&i| face.normal[i] != 0.)
                .expect("qualified original face");
            let [u, v] = [(axis + 1) % 3, (axis + 2) % 3];
            let sign = I::point(face.normal[axis]);
            let mut area = I::ZERO;
            let mut moment = [I::ZERO; 3];
            for ring in &face.rings {
                for (a, b) in ring.iter().zip(ring.iter().cycle().skip(1)) {
                    let c =
                        (I::point(a[u]) * I::point(b[v]) - I::point(a[v]) * I::point(b[u])) * sign;
                    area = area + c / I::point(2.);
                    moment[u] = moment[u] + (I::point(a[u]) + I::point(b[u])) * c / I::point(6.);
                    moment[v] = moment[v] + (I::point(a[v]) + I::point(b[v])) * c / I::point(6.);
                }
            }
            moment[axis] = I::point(face.rings[0][0][axis]) * area;
            for (i, bore) in self.source.bores().iter().enumerate() {
                for end in 0..2 {
                    if bore.faces[end] != index {
                        continue;
                    }
                    let radius = self
                        .source
                        .seats()
                        .iter()
                        .find(|s| s.bore == i && usize::from(s.end) == end)
                        .map_or(bore.radius_mm, |s| s.radius_mm);
                    let lost = pi * I::point(radius) * I::point(radius);
                    area = area - lost;
                    for (j, m) in moment.iter_mut().enumerate() {
                        *m = *m - lost * I::point(bore.centers_mm[end][j]);
                    }
                }
            }
            if !area.valid() || area.lo <= 0. {
                return Err(unsupported("Original finite cap moments are unresolved"));
            }
            centers.push(moment.map(|m| m / area));
        }
        for (i, bore) in self.source.bores().iter().enumerate() {
            let mut points = bore.centers_mm.map(|p| p.map(I::point));
            for s in self.source.seats().iter().filter(|s| s.bore == i) {
                let end = usize::from(s.end);
                let n = self.source.boundary().faces[bore.faces[end]].normal;
                for j in 0..3 {
                    points[end][j] = points[end][j] - I::point(n[j]) * I::point(s.depth_mm);
                }
            }
            centers.push(std::array::from_fn(|j| {
                (points[0][j] + points[1][j]) / I::point(2.)
            }));
        }
        for s in self.source.seats() {
            let bore = self.source.bores()[s.bore];
            let end = usize::from(s.end);
            let r = I::point(bore.radius_mm);
            let outer = I::point(s.radius_mm);
            let depth =
                I::point(s.depth_mm) * (outer + I::point(2.) * r) / (I::point(3.) * (outer + r));
            let normal = self.source.boundary().faces[bore.faces[end]].normal;
            centers.push(std::array::from_fn(|j| {
                I::point(bore.centers_mm[end][j]) - I::point(normal[j]) * depth
            }));
        }
        let stock = self.source.stock().bounds();
        for center in &mut centers {
            for (j, c) in center.iter_mut().enumerate() {
                *c = c
                    .intersect(I::new(stock.min[j], stock.max[j]))
                    .ok_or_else(|| {
                        unsupported("Original surface centroid escaped its retained support hull")
                    })?;
            }
        }
        Ok(centers)
    }
    fn centroid_bounds(
        &self,
        frame: PlaneFrame,
        area: [f64; 2],
    ) -> Result<Vec<[[f64; 2]; 3]>, GeometryError> {
        let map = placement(frame)?;
        let world = self.world_bounds(frame)?;
        // For w in[a,b], mu_w-mu_0=E[(x-mu_0)*(w-(a+b)/2)]/E[w].
        // Each coordinate is bounded by its WHOLE support span*(b-a)/(2a).
        // This bounds azimuth-dependent cone/cylinder area weighting without
        // pretending that mapping the unweighted scalar centroid is exact.
        let variation =
            (I::point(area[1]) - I::point(area[0])) / (I::point(2.) * I::point(area[0]));
        let covariance = std::array::from_fn::<_, 3, _>(|j| {
            ((I::point(world.max[j]) - I::point(world.min[j])) * variation).hi
        });
        self.local_centroids()?
            .into_iter()
            .enumerate()
            .map(|(i, c)| {
                let mut b = map.point_bounds(c.map(bounds)).map_err(placement_error)?;
                if i >= self.source.boundary().faces.len() {
                    for j in 0..3 {
                        b[j] = bounds(interval(b[j]) + I::new(-covariance[j], covariance[j]));
                    }
                }
                Ok(b)
            })
            .collect()
    }
    pub(super) fn world_faces(
        &self,
        local: Vec<Face>,
        frame: PlaneFrame,
    ) -> Result<Vec<Face>, GeometryError> {
        let measures = self.measures(frame)?;
        let map = placement(frame)?;
        local
            .into_iter()
            .map(|f| {
                let point = |p: Point| {
                    map.point_bounds(p.map(|v| [v, v]))
                        .map(|b| b.map(|v| v[0].midpoint(v[1])))
                        .map_err(placement_error)
                };
                let plane = f
                    .plane
                    .map(|(p, n)| Ok((point(p)?, normal(frame, n)?)))
                    .transpose()?;
                let a = measures.face_area_mm2[f.index as usize];
                Ok(Face {
                    plane,
                    center_mm: measures.face_centroid_mm[f.index as usize]
                        .map(|v| v[0].midpoint(v[1])),
                    area_mm2: a[0].midpoint(a[1]),
                    ..f
                })
            })
            .collect()
    }
    pub(super) fn world_mesh(
        &self,
        frame: PlaneFrame,
        options: multicad_rust_kernel::display::Tessellation,
    ) -> Result<Mesh, GeometryError> {
        let map = placement(frame)?;
        let norm = (I::ONE + I::point(map.metric_error())).sqrt().hi;
        let local = multicad_rust_kernel::display::Tessellation {
            deflection: options.deflection / (norm * 4.),
            ..options
        };
        let mut mesh = self.source.mesh(local, &|| false).map_err(failure)?;
        let mut maximum_rounding = 0_f64;
        for p in mesh.vertices.as_chunks_mut::<3>().0 {
            let actual = transform_point(frame, *p);
            let exact = map
                .point_bounds(p.map(|v| [v, v]))
                .map_err(placement_error)?;
            let squared = (0..3).fold(I::ZERO, |s, i| {
                s + (interval(exact[i]) - I::point(actual[i])).square()
            });
            let error = I::new(squared.lo.max(0.), squared.hi).sqrt().hi;
            maximum_rounding = maximum_rounding.max(error);
            *p = actual;
        }
        let error = I::point(local.deflection) * I::point(norm) + I::point(maximum_rounding);
        if !error.valid() || error.hi > options.deflection {
            return Err(unsupported(
                "Actual placed facets exceed the requested whole-source position bound",
            ));
        }
        Ok(Mesh {
            vertices: mesh.vertices,
            indices: mesh.indices,
            face_indices: mesh.face_indices,
        })
    }
}
