//! Occupied contact proof from actual shell values, with a declared planar
//! sewing collar. No authored recipe or added-volume identity admits a shell.
use super::*;
use std::collections::BTreeMap;
mod resources;
use resources::{DELEGATED_WORK, reserve};

#[derive(Clone, Copy, Debug)]
pub struct Options {
    /// Local record/control visits plus separately reserved delegated ceilings.
    /// This is not a measurement of delegated work actually performed.
    pub maximum_work: usize,
    pub maximum_bytes: usize,
    /// Other live caller state, beyond the moved shell counted by this phase.
    pub additional_resident_bytes: usize,
}
impl Default for Options {
    fn default() -> Self {
        Self {
            maximum_work: 96_000_000,
            maximum_bytes: 64 * 1024 * 1024,
            additional_resident_bytes: 0,
        }
    }
}
#[derive(Clone, Copy, Debug)]
pub struct Receipt {
    pub local_work_used: usize,
    pub delegated_work_ceiling_reserved: usize,
    pub aggregate_work_reserved: usize,
    pub retained_source_bytes: usize,
    pub retained_result_bytes_upper_bound: usize,
    /// Conservative live allocation reservation, not measured allocator RSS.
    pub phase_peak_bytes_reserved: usize,
}

/// Qualified sewn/rectified contact representation, not exact occupied
/// material of the unmodified raw supports inside their sewing collars.
/// Original supports remain retained; `join_sewing_mm` records their bound.
pub struct OccupiedHead {
    shell: Ownership,
    thread: RadialBoundary,
    polygon: [[f64; 2]; 6],
    height: [f64; 2],
    planar_sewing_mm: f64,
    join_sewing_mm: f64,
    receipt: Receipt,
}
impl OccupiedHead {
    pub fn shell(&self) -> &Ownership {
        &self.shell
    }
    /// The auxiliary cap closes only the decoded ORIGINAL thread trim. It is
    /// a verification construction, not an extra face in the exported shell.
    pub fn capped_thread(&self) -> &RadialBoundary {
        &self.thread
    }
    pub fn head_polygon_mm(&self) -> &[[f64; 2]; 6] {
        &self.polygon
    }
    pub fn head_height_mm(&self) -> [f64; 2] {
        self.height
    }
    pub fn planar_sewing_mm(&self) -> f64 {
        self.planar_sewing_mm
    }
    /// Outward bound including planar residual, inherited radial sewing and
    /// the complete inherited top axial collar about the shared contact plane.
    pub fn join_sewing_mm(&self) -> f64 {
        self.join_sewing_mm
    }
    pub fn receipt(&self) -> Receipt {
        self.receipt
    }
}
fn invalid<T>() -> Result<T, Error> {
    Err(Error::Invalid)
}
fn plane(input: &Input, face: &Face) -> Result<Axis, Error> {
    match input.surfaces.get(face.surface) {
        Some(Surface::Plane(a)) => {
            if a.origin
                .iter()
                .chain(&a.normal)
                .chain(&a.reference)
                .any(|x| !x.is_finite())
            {
                return Err(Error::Invalid);
            }
            let norm = a.normal[0]
                .portable_hypot(a.normal[1])
                .portable_hypot(a.normal[2]);
            let cross = [
                a.normal[1] * a.reference[2] - a.normal[2] * a.reference[1],
                a.normal[2] * a.reference[0] - a.normal[0] * a.reference[2],
                a.normal[0] * a.reference[1] - a.normal[1] * a.reference[0],
            ];
            let transverse = cross[0].portable_hypot(cross[1]).portable_hypot(cross[2]);
            if !norm.is_finite() || norm == 0. || !transverse.is_finite() || transverse == 0. {
                return Err(Error::Invalid);
            }
            Ok(*a)
        }
        _ => Err(Error::Unsupported),
    }
}
fn ring_vertices(input: &Input, ring: &[(usize, bool)]) -> Result<Vec<usize>, Error> {
    ring.iter()
        .map(|&(e, d)| {
            input
                .edges
                .get(e)
                .map(|e| e.vertices[usize::from(!d)])
                .ok_or(Error::Invalid)
        })
        .collect()
}
/// Prove a single fixed-frame outward head. The actual finite thread is
/// independently requalified before proving the hex halfspaces contain its
/// complete occupied cross-section. The declared collar permits rectifying
/// the inherited uncertain cap to the actual contact plane. Only this sewn
/// reference has disjoint open axial interiors; raw collars are not erased.
/// Volume is not an input to admission.
pub fn certify(
    shell: Ownership,
    tolerance: f64,
    minimum: f64,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<OccupiedHead, Error> {
    certify_with_options(shell, tolerance, minimum, Options::default(), interrupt)
}
pub fn certify_with_options(
    shell: Ownership,
    tolerance: f64,
    minimum: f64,
    options: Options,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<OccupiedHead, Error> {
    if !tolerance.is_finite()
        || tolerance <= 0.
        || tolerance > 1e-4
        || !minimum.is_finite()
        || minimum <= 0.
    {
        return invalid();
    }
    let reservation = reserve(shell.input(), options)?;
    let local_allowance = options
        .maximum_work
        .checked_sub(DELEGATED_WORK)
        .ok_or(Error::Limit)?;
    let mut work = Work {
        remaining: local_allowance,
        interrupt,
    };
    work.step()?;
    let s = shell.input();
    let shoulders = s
        .faces
        .iter()
        .enumerate()
        .filter(|(_, f)| !f.inner.is_empty())
        .collect::<Vec<_>>();
    let [(shoulder_index, shoulder)] = shoulders.as_slice() else {
        return invalid();
    };
    let shoulder_axis = plane(s, shoulder)?;
    if shoulder.forward
        || shoulder_axis.normal != [0., 0., 1.]
        || shoulder_axis.reference != [1., 0., 0.]
        || shoulder_axis.origin[0] != 0.
        || shoulder_axis.origin[1] != 0.
        || shoulder_axis.origin[2] <= minimum
        || shoulder.inner.len() != 1
        || shoulder.outer.len() != 6
    {
        return Err(Error::Unsupported);
    }
    let lower = ring_vertices(s, &shoulder.outer)?;
    let low = shoulder_axis.origin[2];
    let original = |ring: &[(usize, bool)]| {
        ring.iter()
            .all(|&(e, _)| matches!(s.curves[s.edges[e].curve], Curve::Original(_)))
    };
    let straight = |ring: &[(usize, bool)]| {
        ring.iter()
            .all(|&(e, _)| matches!(s.curves[s.edges[e].curve], Curve::Line(_)))
    };
    if !original(&shoulder.inner[0]) {
        return invalid();
    }
    if !straight(&shoulder.outer) {
        return Err(Error::Unsupported);
    }
    let mut thread_faces = Vec::new();
    let mut head_faces = Vec::new();
    for (i, f) in s.faces.iter().enumerate() {
        work.step()?;
        if i == *shoulder_index {
            continue;
        }
        if f.inner.is_empty() && original(&f.outer) {
            thread_faces.push(f);
        } else {
            head_faces.push(f);
        }
    }
    if thread_faces.len() != 5 || head_faces.len() != 7 {
        return Err(Error::Unsupported);
    }
    if head_faces.iter().any(|f| !straight(&f.outer)) {
        return Err(Error::Unsupported);
    }
    let top_candidates=head_faces.iter().filter(|f|matches!(s.surfaces[f.surface],Surface::Plane(a) if a.normal==[0.,0.,1.] && a.origin[2]>low)).copied().collect::<Vec<_>>();
    let [top] = top_candidates.as_slice() else {
        return invalid();
    };
    if !top.forward || !top.inner.is_empty() || top.outer.len() != 6 {
        return invalid();
    }
    let top_axis = plane(s, top)?;
    let high = top_axis.origin[2];
    let axial_extent = I::point(high) - I::point(low);
    if top_axis.origin != [0., 0., high] || top_axis.reference != [1., 0., 0.] {
        return invalid();
    }
    if axial_extent.lo <= minimum.max(32. * tolerance) {
        return Err(Error::Unresolved);
    }
    let top_vertices = ring_vertices(s, &top.outer)?;
    let polygon: [[_; 2]; 6] = std::array::from_fn(|i| {
        let p = s.vertices[top_vertices[i]];
        [p[0], p[1]]
    });
    let mut lower_for_top = [0usize; 6];
    for i in 0..6 {
        work.step()?;
        if s.vertices[top_vertices[i]][2] != high {
            return invalid();
        }
        let matching = lower
            .iter()
            .copied()
            .filter(|&j| s.vertices[j] == [polygon[i][0], polygon[i][1], low])
            .collect::<Vec<_>>();
        let [v] = matching.as_slice() else {
            return invalid();
        };
        lower_for_top[i] = *v;
    }
    let mut planar_sewing = 0f64;
    let mut claimed = BTreeSet::new();
    for i in 0..6 {
        let j = (i + 1) % 6;
        let k = (i + 2) % 6;
        let a = polygon[i];
        let b = polygon[j];
        let c = polygon[k];
        let dx = I::point(b[0]) - I::point(a[0]);
        let dy = I::point(b[1]) - I::point(a[1]);
        let cross = dx * (I::point(c[1]) - I::point(b[1])) - dy * (I::point(c[0]) - I::point(b[0]));
        if cross.lo <= 0. {
            return Err(Error::Unresolved);
        }
        // Local left turns permit a multiply wound star. Every remaining
        // vertex must lie strictly inside EVERY edge halfspace, which proves
        // the complete ordered polygon is the simple convex hull boundary.
        for (k, p) in polygon.iter().enumerate() {
            if k == i || k == j {
                continue;
            }
            let side =
                dx * (I::point(p[1]) - I::point(a[1])) - dy * (I::point(p[0]) - I::point(a[0]));
            if side.lo <= 0. {
                return Err(Error::Unresolved);
            }
        }
        let expected = BTreeSet::from([
            lower_for_top[i],
            lower_for_top[j],
            top_vertices[i],
            top_vertices[j],
        ]);
        let found = head_faces
            .iter()
            .enumerate()
            .filter(|(_, f)| !std::ptr::eq(**f, *top))
            .filter_map(|(fi, f)| {
                let vertices = ring_vertices(s, &f.outer).ok()?;
                (vertices.len() == 4 && vertices.into_iter().collect::<BTreeSet<_>>() == expected)
                    .then_some((fi, *f))
            })
            .collect::<Vec<_>>();
        let [(fi, face)] = found.as_slice() else {
            return invalid();
        };
        if !claimed.insert(*fi) || !face.forward || !face.inner.is_empty() {
            return invalid();
        }
        let axis = plane(s, face)?;
        let norm = dx.square() + dy.square();
        let length = norm.sqrt().positive()?;
        let n = [dy / length, -dx / length, I::point(0.)];
        let raw = axis.normal.map(I::point);
        let actual_norm = (raw[0].square() + raw[1].square() + raw[2].square())
            .sqrt()
            .positive()?;
        let ns = raw.map(|x| x / actual_norm);
        let alignment = n[0] * ns[0] + n[1] * ns[1] + n[2] * ns[2];
        alignment.checked()?;
        if alignment.lo < 1. - 64. * f64::EPSILON || axis.normal[2] != 0. {
            return Err(Error::Unsupported);
        }
        for v in &expected {
            work.step()?;
            let p = s.vertices[*v];
            let residual = (0..3).fold(I::point(0.), |sum, k| {
                sum + (I::point(p[k]) - I::point(axis.origin[k])) * ns[k]
            });
            planar_sewing = planar_sewing.max(residual.magnitude());
        }
    }
    if planar_sewing > tolerance / 16. {
        return Err(Error::Unresolved);
    }
    // Reindex only the actual decoded thread tables. No pole, knot, trim or
    // analytic support is changed; the shoulder's actual inner cycle is capped.
    let mut tfaces = thread_faces
        .into_iter()
        .map(|f| radial::Face {
            surface: f.surface,
            outer: f.outer.clone(),
            forward: f.forward,
        })
        .collect::<Vec<_>>();
    tfaces.push(radial::Face {
        surface: shoulder.surface,
        outer: shoulder.inner[0].clone(),
        forward: true,
    });
    let mut input = radial::BoundaryInput {
        vertices: vec![],
        curves: vec![],
        edges: vec![],
        surfaces: vec![],
        faces: vec![],
    };
    let (mut vertices, mut edges, mut curves, mut surfaces) = (
        BTreeMap::new(),
        BTreeMap::new(),
        BTreeMap::new(),
        BTreeMap::new(),
    );
    for f in tfaces {
        work.step()?;
        let surface = *surfaces.entry(f.surface).or_insert_with(|| {
            let i = input.surfaces.len();
            input.surfaces.push(s.surfaces[f.surface].clone());
            i
        });
        let mut outer = Vec::new();
        for (e, d) in f.outer {
            work.step()?;
            let index = if let Some(&i) = edges.get(&e) {
                i
            } else {
                let source = s.edges[e];
                let mut ends = [0; 2];
                for (k, v) in source.vertices.into_iter().enumerate() {
                    ends[k] = *vertices.entry(v).or_insert_with(|| {
                        let i = input.vertices.len();
                        input.vertices.push(s.vertices[v]);
                        i
                    });
                }
                let curve = if let Some(&i) = curves.get(&source.curve) {
                    i
                } else {
                    let Curve::Original(c) = &s.curves[source.curve] else {
                        return invalid();
                    };
                    if let radial::Curve::Spline { poles, .. } = c {
                        work.charge(poles.len())?;
                    }
                    let i = input.curves.len();
                    input.curves.push(c.clone());
                    curves.insert(source.curve, i);
                    i
                };
                let i = input.edges.len();
                input.edges.push(radial::Edge {
                    vertices: ends,
                    curve,
                    forward: source.forward,
                });
                edges.insert(e, i);
                i
            };
            outer.push((index, d));
        }
        input.faces.push(radial::Face {
            surface,
            outer,
            forward: f.forward,
        });
    }
    let thread = radial::BoundaryOwnership::new(input, tolerance, interrupt)?
        .bind_surfaces(tolerance, interrupt)?
        .certify_axial(tolerance, minimum, interrupt)?
        .certify_caps(interrupt)?
        .certify_ordering(interrupt)?
        .certify_material(tolerance, interrupt)?;
    let caps = thread.cap_axes();
    if caps[0].origin != [0., 0., 0.] || caps[1].origin != [0., 0., low] {
        return Err(Error::Unsupported);
    }
    let top_collar = thread.collars().top_axial_mm;
    let axial_uncertainty = (I::new(top_collar[0], top_collar[1]) - I::point(low))
        .checked()?
        .magnitude();
    let join_sewing = (I::point(planar_sewing)
        + I::point(thread.boundary_error_mm())
        + I::point(axial_uncertainty))
    .checked()?
    .hi;
    let required_height =
        I::point(minimum.max(32. * tolerance)) + I::point(2.) * I::point(join_sewing);
    if axial_extent.lo <= required_height.hi {
        return Err(Error::Unresolved);
    }
    for i in 0..6 {
        work.step()?;
        let a = polygon[i];
        let b = polygon[(i + 1) % 6];
        let normal = [b[1] - a[1], a[0] - b[0], 0.];
        let limit = I::point(normal[0]) * I::point(a[0]) + I::point(normal[1]) * I::point(a[1]);
        // Complete finite-control hull, not an adaptive support heap or samples.
        let support = full_control_support(&thread, normal, &mut work)?;
        let margin = (I::point(minimum) + I::point(join_sewing))
            * (I::point(normal[0]).square() + I::point(normal[1]).square()).sqrt();
        if (limit - I::point(support) - margin).lo <= 0. {
            return Err(Error::Unresolved);
        }
    }
    Ok(OccupiedHead {
        shell,
        thread,
        polygon,
        height: [low, high],
        planar_sewing_mm: planar_sewing,
        join_sewing_mm: join_sewing,
        receipt: Receipt {
            local_work_used: local_allowance - work.remaining,
            delegated_work_ceiling_reserved: DELEGATED_WORK,
            aggregate_work_reserved: DELEGATED_WORK + local_allowance - work.remaining,
            retained_source_bytes: reservation.source,
            retained_result_bytes_upper_bound: reservation.retained,
            phase_peak_bytes_reserved: reservation.peak,
        },
    })
}

fn full_control_support(
    boundary: &RadialBoundary,
    n: [f64; 3],
    work: &mut Work<'_>,
) -> Result<f64, Error> {
    if n[2] != 0. {
        return Err(Error::Unsupported);
    }
    let ni = n.map(I::point);
    let norm = (ni[0].square() + ni[1].square()).sqrt().positive()?;
    let dot = |p: [f64; 3]| {
        (ni[0] * I::point(p[0]) + ni[1] * I::point(p[1]) + ni[2] * I::point(p[2])).checked()
    };
    let mut upper = f64::NEG_INFINITY;
    for p in &boundary.input().vertices {
        work.step()?;
        upper = upper.max(dot(*p)?.hi);
    }
    for curve in &boundary.input().curves {
        work.step()?;
        match curve {
            radial::Curve::Spline { poles, .. } => {
                for p in poles {
                    work.step()?;
                    upper = upper.max(dot(*p)?.hi);
                }
            }
            radial::Curve::Circle { axis, radius_mm } => {
                // Full sphere support conservatively encloses the actual arc,
                // independent of its angular interval and axis parameterization.
                upper = upper.max(
                    (dot(axis.origin)? + I::point(*radius_mm) * norm)
                        .checked()?
                        .hi,
                );
            }
        }
    }
    for surface in &boundary.input().surfaces {
        work.step()?;
        match surface {
            Surface::Tensor(t) => {
                for p in t.poles() {
                    work.step()?;
                    upper = upper.max(dot(*p)?.hi);
                }
            }
            Surface::Cone {
                axis,
                radius_mm,
                semi_angle_rad,
            } => {
                if axis.normal != [0., 0., 1.] && axis.normal != [0., 0., -1.] {
                    return Err(Error::Unsupported);
                }
                let angle = I::point(*semi_angle_rad);
                let slope = (angle.sin() / angle.cos()).checked()?;
                // The actual affine radius attains its extrema at cap planes.
                // Include the entire circular section/core, not only its rims.
                for cap in boundary.cap_axes() {
                    let z = (I::point(cap.origin[2]) - I::point(axis.origin[2]))
                        * I::point(axis.normal[2]);
                    let radius = (I::point(*radius_mm) + slope * z).positive()?;
                    upper = upper.max((dot(axis.origin)? + radius * norm).checked()?.hi);
                }
            }
            Surface::Plane(axis) => {
                upper = upper.max(dot(axis.origin)?.hi);
            }
        }
    }
    (I::point(upper) + I::point(boundary.boundary_error_mm()) * norm)
        .checked()
        .map(|v| v.hi)
}
