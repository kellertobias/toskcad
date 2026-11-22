//! Finite hex-head shell values. This stage preserves the original thread
//! supports and earns multi-loop manifold ownership, separately from material.
use super::*;
use crate::radial_tensor::{self as radial, Axis, RadialBoundary, Surface};
use std::collections::BTreeSet;
pub mod occupied;

#[derive(Clone, Debug, PartialEq)]
pub enum Curve {
    Original(radial::Curve),
    /// Exact finite straight segment, with its actual endpoint values.
    Line([usize; 2]),
}
#[derive(Clone, Debug, PartialEq)]
pub struct Face {
    pub surface: usize,
    pub outer: Vec<(usize, bool)>,
    pub inner: Vec<Vec<(usize, bool)>>,
    pub forward: bool,
}
#[derive(Clone, Debug, PartialEq)]
pub struct Input {
    pub vertices: Vec<[f64; 3]>,
    pub curves: Vec<Curve>,
    pub edges: Vec<radial::Edge>,
    pub surfaces: Vec<Surface>,
    pub faces: Vec<Face>,
}
/// This proof certifies only oriented closed-shell ownership. It deliberately
/// does not grant occupied material, embedding, source matching or STEP admission.
pub struct Ownership {
    input: Input,
}
impl Ownership {
    pub fn input(&self) -> &Input {
        &self.input
    }
    pub fn new(input: Input, interrupt: Option<&dyn Fn() -> bool>) -> Result<Self, Error> {
        let v = input.vertices.len();
        let e = input.edges.len();
        let f = input.faces.len();
        if v > 96 || e > 128 || f > 48 || input.curves.len() > 128 || input.surfaces.len() > 48 {
            return Err(Error::Limit);
        }
        if v == 0 || e == 0 || f == 0 {
            return Err(Error::Invalid);
        }
        // Count retained control storage before any later verification clones.
        // This is an input-family ceiling, not a claim about allocator RSS.
        let mut controls = 0usize;
        for c in &input.curves {
            if let Curve::Original(radial::Curve::Spline { poles, .. }) = c {
                controls = controls.checked_add(poles.len()).ok_or(Error::Limit)?;
            }
        }
        for s in &input.surfaces {
            if let Surface::Tensor(t) = s {
                controls = controls.checked_add(t.poles().len()).ok_or(Error::Limit)?;
            }
        }
        if controls > 524_288 {
            return Err(Error::Limit);
        }
        let mut work = Work {
            remaining: 32_768,
            interrupt,
        };
        for p in &input.vertices {
            work.step()?;
            if p.iter().any(|x| !x.is_finite() || x.abs() > 2e7) {
                return Err(Error::Invalid);
            }
        }
        let mut uses = vec![Vec::new(); e];
        let mut incident = vec![BTreeSet::new(); v];
        let mut links = vec![Vec::new(); v];
        for (i, edge) in input.edges.iter().enumerate() {
            work.step()?;
            let [a, b] = edge.vertices;
            if a >= v || b >= v || a == b || edge.curve >= input.curves.len() {
                return Err(Error::Invalid);
            }
            if let Curve::Line(ends) = &input.curves[edge.curve] {
                let ends = if edge.forward {
                    *ends
                } else {
                    [ends[1], ends[0]]
                };
                if ends != [a, b] {
                    return Err(Error::Invalid);
                }
            }
            incident[a].insert(i);
            incident[b].insert(i);
        }
        let mut chi = 0isize;
        let mut annuli = 0;
        let mut used_surfaces = BTreeSet::new();
        let mut used_curves = BTreeSet::new();
        for (fi, face) in input.faces.iter().enumerate() {
            work.step()?;
            if face.surface >= input.surfaces.len() || face.inner.len() > 1 {
                return Err(Error::Invalid);
            }
            used_surfaces.insert(face.surface);
            chi += 1 - face.inner.len() as isize;
            annuli += face.inner.len();
            let mut face_vertices = BTreeSet::new();
            for ring in std::iter::once(&face.outer).chain(&face.inner) {
                if ring.len() < 3 || ring.len() > 80 {
                    return Err(Error::Invalid);
                }
                for (ri, &(edge, direction)) in ring.iter().enumerate() {
                    work.step()?;
                    let edge_value = input.edges.get(edge).ok_or(Error::Invalid)?;
                    used_curves.insert(edge_value.curve);
                    let [a, b] = edge_value.vertices;
                    let [from, to] = if direction { [a, b] } else { [b, a] };
                    if !face_vertices.insert(from) {
                        return Err(Error::Invalid);
                    }
                    let (next, sense) = ring[(ri + 1) % ring.len()];
                    let next_edge = input.edges.get(next).ok_or(Error::Invalid)?;
                    if to != next_edge.vertices[usize::from(!sense)] {
                        return Err(Error::Invalid);
                    }
                    uses[edge].push((fi, direction));
                    links[to].push((edge, next));
                }
            }
        }
        if annuli != 1
            || v as isize - e as isize + chi != 2
            || used_surfaces.len() != input.surfaces.len()
            || used_curves.len() != input.curves.len()
        {
            return Err(Error::Invalid);
        }
        let mut neighbors = vec![Vec::new(); f];
        for pair in uses {
            work.step()?;
            let [(a, sa), (b, sb)] = pair.as_slice() else {
                return Err(Error::Invalid);
            };
            if a == b || sa == sb {
                return Err(Error::Invalid);
            }
            neighbors[*a].push(*b);
            neighbors[*b].push(*a);
        }
        connected(&neighbors, &mut work)?;
        for (vi, edges) in incident.iter().enumerate() {
            work.step()?;
            if edges.len() < 2 {
                return Err(Error::Invalid);
            }
            let edges = edges.iter().copied().collect::<Vec<_>>();
            let mut fan = vec![Vec::new(); edges.len()];
            for &(a, b) in &links[vi] {
                let a = edges.iter().position(|x| *x == a).ok_or(Error::Invalid)?;
                let b = edges.iter().position(|x| *x == b).ok_or(Error::Invalid)?;
                fan[a].push(b);
                fan[b].push(a);
            }
            if fan.iter().any(|n| n.len() != 2) {
                return Err(Error::Invalid);
            }
            connected(&fan, &mut work)?;
        }
        Ok(Self { input })
    }
}
struct Work<'a> {
    remaining: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl Work<'_> {
    fn step(&mut self) -> Result<(), Error> {
        self.charge(1)
    }
    fn charge(&mut self, n: usize) -> Result<(), Error> {
        if self.interrupt.is_some_and(|f| f()) {
            return Err(Error::Cancelled);
        }
        self.remaining = self.remaining.checked_sub(n).ok_or(Error::Limit)?;
        Ok(())
    }
}
fn connected(neighbors: &[Vec<usize>], work: &mut Work<'_>) -> Result<(), Error> {
    let mut seen = vec![false; neighbors.len()];
    let mut pending = vec![0];
    while let Some(i) = pending.pop() {
        work.step()?;
        if seen[i] {
            continue;
        }
        seen[i] = true;
        pending.extend(neighbors[i].iter().copied().filter(|&j| !seen[j]));
    }
    if seen.into_iter().all(|x| x) {
        Ok(())
    } else {
        Err(Error::Invalid)
    }
}
/// Construct one fixed-frame external NPT head from an already admitted actual
/// thread boundary. The qualified boundary is retained separately by the caller;
/// the shell alone earns ownership and never inherits its material admission.
pub fn external_npt_head(
    thread: &AuthoredThread,
    boundary: &RadialBoundary,
    across: f64,
    height: f64,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<Ownership, Error> {
    if !across.is_finite() || across <= 0. || !height.is_finite() || height <= 0. {
        return Err(Error::Invalid);
    }
    if across > 1e7 || height > 1e7 {
        return Err(Error::Limit);
    }
    if thread.spec.radial_taper <= 0.
        || thread.spec.crest_radius_limit_mm.is_some()
        || !matches!(thread.spec.section, AuthoredSection::NptFlat { .. })
    {
        return Err(Error::Unsupported);
    }
    let axes = boundary.cap_axes();
    let l = thread.spec.length_mm;
    if axes[0].origin != [0., 0., 0.]
        || axes[1].origin != [0., 0., l]
        || axes
            .iter()
            .any(|a| a.normal != [0., 0., 1.] && a.normal != [0., 0., -1.])
    {
        return Err(Error::Unsupported);
    }
    let original = boundary.input();
    let expected = thread.raw_boundary_ownership(thread.tolerance, interrupt)?;
    if expected.ownership.input() != original {
        return Err(Error::Unsupported);
    }
    let mut input = Input {
        vertices: original.vertices.clone(),
        curves: original
            .curves
            .iter()
            .cloned()
            .map(Curve::Original)
            .collect(),
        edges: original.edges.clone(),
        surfaces: original.surfaces.clone(),
        faces: original
            .faces
            .iter()
            .map(|f| Face {
                surface: f.surface,
                outer: f.outer.clone(),
                inner: vec![],
                forward: f.forward,
            })
            .collect(),
    };
    let cap = input
        .faces
        .iter()
        .position(
            |f| matches!(&input.surfaces[f.surface],Surface::Plane(a) if a.origin==axes[1].origin),
        )
        .ok_or(Error::Invalid)?;
    // Original top-cap edge uses already oppose the unchanged thread sides.
    let inner = input.faces[cap].outer.clone();
    let start = input.vertices.len();
    let half = across / 2.;
    let h = half / 3f64.sqrt();
    let polygon = [
        [2. * h, 0.],
        [h, half],
        [-h, half],
        [-2. * h, 0.],
        [-h, -half],
        [h, -half],
    ];
    let ztop = l + height;
    if !ztop.is_finite() || ztop > 1e7 {
        return Err(Error::Limit);
    }
    for z in [l, ztop] {
        for p in polygon {
            input.vertices.push([p[0], p[1], z]);
        }
    }
    let mut edge = |a, b| {
        let curve = input.curves.len();
        input.curves.push(Curve::Line([a, b]));
        let index = input.edges.len();
        input.edges.push(radial::Edge {
            vertices: [a, b],
            curve,
            forward: true,
        });
        index
    };
    let lower: [usize; 6] = std::array::from_fn(|i| edge(start + i, start + (i + 1) % 6));
    let upper: [usize; 6] = std::array::from_fn(|i| edge(start + 6 + i, start + 6 + (i + 1) % 6));
    let vertical: [usize; 6] = std::array::from_fn(|i| edge(start + i, start + 6 + i));
    // Viewed from +z the shoulder outer cycle is clockwise and its hole CCW.
    input.faces[cap] = Face {
        surface: input.faces[cap].surface,
        outer: lower.iter().rev().map(|&e| (e, false)).collect(),
        inner: vec![inner],
        forward: false,
    };
    for i in 0..6 {
        let a = input.vertices[start + i];
        let b = input.vertices[start + (i + 1) % 6];
        let dx = b[0] - a[0];
        let dy = b[1] - a[1];
        let norm = dx.portable_hypot(dy);
        let surface = input.surfaces.len();
        input.surfaces.push(Surface::Plane(Axis {
            origin: a,
            normal: [dy / norm, -dx / norm, 0.],
            reference: [0., 0., 1.],
        }));
        input.faces.push(Face {
            surface,
            outer: vec![
                (lower[i], true),
                (vertical[(i + 1) % 6], true),
                (upper[i], false),
                (vertical[i], false),
            ],
            inner: vec![],
            forward: true,
        });
    }
    let surface = input.surfaces.len();
    input.surfaces.push(Surface::Plane(Axis {
        origin: [0., 0., ztop],
        normal: [0., 0., 1.],
        reference: [1., 0., 0.],
    }));
    input.faces.push(Face {
        surface,
        outer: upper.into_iter().map(|e| (e, true)).collect(),
        inner: vec![],
        forward: true,
    });
    Ownership::new(input, interrupt)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn fixture() -> (AuthoredThread, RadialBoundary) {
        let t = AuthoredThread::new(
            AuthoredSpec {
                major_radius_mm: 5.,
                pitch_mm: 1.,
                length_mm: 6.,
                radial_taper: 1. / 32.,
                section: AuthoredSection::NptFlat {
                    depth_mm: 0.6,
                    flat_width_mm: 0.2,
                },
                crest_radius_limit_mm: None,
            },
            1e-7,
            1e-6,
        )
        .unwrap();
        let boundary = t
            .raw_boundary_ownership(1e-7, None)
            .unwrap()
            .ownership
            .bind_surfaces(1e-7, None)
            .unwrap()
            .certify_axial(1e-7, 1e-6, None)
            .unwrap()
            .certify_caps(None)
            .unwrap()
            .certify_ordering(None)
            .unwrap()
            .certify_material(1e-7, None)
            .unwrap();
        (t, boundary)
    }
    #[test]
    fn original_npt_rims_earn_annular_closed_head_ownership_without_meshes() {
        let (t, b) = fixture();
        let shell = external_npt_head(&t, &b, 16., 3., None).unwrap();
        let s = shell.input();
        assert_eq!(s.faces.len(), 13);
        assert_eq!(&s.vertices[..b.input().vertices.len()], b.input().vertices);
        assert_eq!(&s.edges[..b.input().edges.len()], b.input().edges);
        for (i, c) in b.input().curves.iter().enumerate() {
            assert_eq!(s.curves[i], Curve::Original(c.clone()));
        }
        for (i, surface) in b.input().surfaces.iter().enumerate() {
            assert_eq!(&s.surfaces[i], surface);
        }
        let shoulder = s.faces.iter().find(|f| !f.inner.is_empty()).unwrap();
        assert_eq!(shoulder.inner[0], b.input().faces[1].outer);
        assert_eq!(
            s.vertices.len() as isize - s.edges.len() as isize + s.faces.len() as isize,
            3
        );
        assert!(Ownership::new(s.clone(), None).is_ok());
        let mut bad = s.clone();
        bad.faces[1].inner.clear();
        assert!(matches!(Ownership::new(bad, None), Err(Error::Invalid)));
        let mut bad = s.clone();
        bad.faces[1].inner[0][0].1 ^= true;
        assert!(matches!(Ownership::new(bad, None), Err(Error::Invalid)));
        let mut bad = s.clone();
        bad.faces[6].outer[0].1 ^= true;
        assert!(matches!(Ownership::new(bad, None), Err(Error::Invalid)));
        let mut bad = s.clone();
        bad.vertices.push([0., 0., 0.]);
        assert!(matches!(Ownership::new(bad, None), Err(Error::Invalid)));
        let mut bad = s.clone();
        bad.edges.push(bad.edges[0]);
        assert!(matches!(Ownership::new(bad, None), Err(Error::Invalid)));
        let mut bad = s.clone();
        bad.vertices[0][0] = f64::NAN;
        assert!(matches!(Ownership::new(bad, None), Err(Error::Invalid)));
        assert!(matches!(
            Ownership::new(s.clone(), Some(&|| true)),
            Err(Error::Cancelled)
        ));
        assert!(matches!(
            external_npt_head(&t, &b, 16., 3., Some(&|| true)),
            Err(Error::Cancelled)
        ));
        assert!(external_npt_head(&t, &b, 16., 3., None).is_ok());
        let occupied =
            occupied::certify(Ownership::new(s.clone(), None).unwrap(), 1e-7, 1e-6, None).unwrap();
        assert_eq!(occupied.head_height_mm(), [6., 9.]);
        assert!(occupied.planar_sewing_mm() < 1e-7 / 16.);
        assert!(
            occupied.join_sewing_mm()
                >= occupied.planar_sewing_mm() + occupied.capped_thread().boundary_error_mm()
        );
        assert_eq!(occupied.capped_thread().input().faces.len(), 6);
        let receipt = occupied.receipt();
        assert_eq!(receipt.delegated_work_ceiling_reserved, 75_000_000);
        assert_eq!(
            receipt.aggregate_work_reserved,
            75_000_000 + receipt.local_work_used
        );
        assert!(receipt.aggregate_work_reserved <= 96_000_000);
        assert!(receipt.retained_source_bytes <= receipt.retained_result_bytes_upper_bound);
        assert!(receipt.retained_result_bytes_upper_bound < receipt.phase_peak_bytes_reserved);
        assert!(receipt.phase_peak_bytes_reserved <= 64 * 1024 * 1024);
        for options in [
            occupied::Options {
                maximum_work: 74_999_999,
                ..Default::default()
            },
            occupied::Options {
                maximum_work: 75_000_001,
                ..Default::default()
            },
            occupied::Options {
                maximum_work: receipt.aggregate_work_reserved - 1,
                ..Default::default()
            },
            occupied::Options {
                maximum_bytes: receipt.phase_peak_bytes_reserved - 1,
                ..Default::default()
            },
            occupied::Options {
                additional_resident_bytes: usize::MAX,
                ..Default::default()
            },
        ] {
            assert!(matches!(
                occupied::certify_with_options(
                    Ownership::new(s.clone(), None).unwrap(),
                    1e-7,
                    1e-6,
                    options,
                    None
                ),
                Err(Error::Limit)
            ));
        }
        let mut excess = s.clone();
        let spline = excess
            .curves
            .iter_mut()
            .find_map(|c| match c {
                Curve::Original(radial::Curve::Spline { poles, .. }) => Some(poles),
                _ => None,
            })
            .unwrap();
        let original_len = spline.len();
        spline.reserve_exact(3_000_000);
        assert_eq!(spline.len(), original_len);
        assert!(matches!(
            occupied::certify(Ownership::new(excess, None).unwrap(), 1e-7, 1e-6, None),
            Err(Error::Limit)
        ));
        let exact = occupied::Options {
            maximum_work: receipt.aggregate_work_reserved,
            maximum_bytes: receipt.phase_peak_bytes_reserved,
            additional_resident_bytes: 0,
        };
        let calls = std::cell::Cell::new(0usize);
        let count = || {
            calls.set(calls.get() + 1);
            false
        };
        let retry = occupied::certify_with_options(
            Ownership::new(s.clone(), None).unwrap(),
            1e-7,
            1e-6,
            exact,
            Some(&count),
        )
        .unwrap();
        assert_eq!(retry.receipt().aggregate_work_reserved, exact.maximum_work);
        assert_eq!(
            retry.receipt().phase_peak_bytes_reserved,
            exact.maximum_bytes
        );
        assert_eq!(retry.shell().input(), s);
        assert!(calls.get() > 100);
        for stop in [calls.get() / 2, calls.get() - 1] {
            let visited = std::cell::Cell::new(0usize);
            let cancel = || {
                let n = visited.get() + 1;
                visited.set(n);
                n >= stop
            };
            assert!(matches!(
                occupied::certify_with_options(
                    Ownership::new(s.clone(), None).unwrap(),
                    1e-7,
                    1e-6,
                    exact,
                    Some(&cancel)
                ),
                Err(Error::Cancelled)
            ));
            assert_eq!(visited.get(), stop);
        }
        // Unchanged corners/ownership do not certify the support between them.
        // A bowed original spline must never acquire straight prism admission.
        for face in [1, 6, 12] {
            let mut bad = s.clone();
            let edge = bad.edges[bad.faces[face].outer[0].0];
            let [a, b] = edge.vertices.map(|v| bad.vertices[v]);
            let mut poles = (0..10)
                .map(|i| std::array::from_fn(|k| a[k] + (b[k] - a[k]) * i as f64 / 9.))
                .collect::<Vec<_>>();
            poles[0] = a;
            poles[9] = b;
            poles[4][2] += 1.;
            bad.curves[edge.curve] = Curve::Original(radial::Curve::Spline { spans: 1, poles });
            let ownership = Ownership::new(bad, None).unwrap();
            assert!(matches!(
                occupied::certify(ownership, 1e-7, 1e-6, None),
                Err(Error::Unsupported)
            ));
        }
        // This ordering has six distinct corners and all local left turns,
        // but winds a star. Even genuine straight edges and recomputed planar
        // supports cannot turn its topological sphere into embedded material.
        let mut star = s.clone();
        let base = b.input().vertices.len();
        for (i, angle) in [0., 130., 260., 30., 160., 290.].into_iter().enumerate() {
            let angle: f64 = angle;
            let (sin, cos) = angle.to_radians().portable_sin_cos();
            for z in [0, 6] {
                star.vertices[base + z + i][0] = 20. * cos;
                star.vertices[base + z + i][1] = 20. * sin;
            }
        }
        for i in 0..6 {
            let a = star.vertices[base + i];
            let q = star.vertices[base + (i + 1) % 6];
            let dx = q[0] - a[0];
            let dy = q[1] - a[1];
            let norm = dx.portable_hypot(dy);
            let surface = star.faces[6 + i].surface;
            star.surfaces[surface] = Surface::Plane(Axis {
                origin: a,
                normal: [dy / norm, -dx / norm, 0.],
                reference: [0., 0., 1.],
            });
            // The stronger star's kernel really contains the entire actual
            // thread, so a later containment test cannot hide the winding bug.
            let normal = [dy / norm, -dx / norm, 0.];
            let upper = b.support_bounds(normal, 0.01, None).unwrap()[1];
            let limit = normal[0] * a[0] + normal[1] * a[1];
            assert!(limit - upper > 1.);
        }
        let ownership = Ownership::new(star, None).unwrap();
        assert!(matches!(
            occupied::certify(ownership, 1e-7, 1e-6, None),
            Err(Error::Unresolved)
        ));
        for mode in 0..3 {
            let mut bad = s.clone();
            let side = bad.faces[6].surface;
            let Surface::Plane(ref mut a) = bad.surfaces[side] else {
                panic!("actual side plane")
            };
            match mode {
                0 => a.normal[0] = f64::NAN,
                1 => a.normal = [0.; 3],
                _ => a.reference = a.normal,
            };
            assert!(matches!(
                occupied::certify(Ownership::new(bad, None).unwrap(), 1e-7, 1e-6, None),
                Err(Error::Invalid)
            ));
        }
        let mut bad = s.clone();
        let side = bad.faces[6].surface;
        let Surface::Plane(ref mut a) = bad.surfaces[side] else {
            panic!("expected actual side plane")
        };
        a.origin[0] += 0.01;
        assert!(matches!(
            occupied::certify(Ownership::new(bad, None).unwrap(), 1e-7, 1e-6, None),
            Err(Error::Unresolved)
        ));
        assert!(matches!(
            occupied::certify(
                external_npt_head(&t, &b, 8., 3., None).unwrap(),
                1e-7,
                1e-6,
                None
            ),
            Err(Error::Unresolved)
        ));
        assert!(matches!(
            occupied::certify(
                external_npt_head(&t, &b, 16., 1e-8, None).unwrap(),
                1e-7,
                1e-6,
                None
            ),
            Err(Error::Unresolved)
        ));
        let contact_thin = 32e-7 + occupied.join_sewing_mm();
        assert!(matches!(
            occupied::certify(
                external_npt_head(&t, &b, 16., contact_thin, None).unwrap(),
                1e-7,
                1e-6,
                None
            ),
            Err(Error::Unresolved)
        ));
        assert!(matches!(
            occupied::certify(
                Ownership::new(s.clone(), None).unwrap(),
                1e-7,
                1e-6,
                Some(&|| true)
            ),
            Err(Error::Cancelled)
        ));
    }
}
