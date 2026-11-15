//! Derived stock candidate. Original owner bijection is retained; no raw STEP
//! occupied set is inferred from coordinate projection or binary64 caches.
use super::*;
#[path = "fields.rs"]
mod fields;
use crate::{
    axial_pipe_thread::authored::StockClipPlacement,
    orthogonal::OrthogonalSolid,
    orthogonal_boundary::{self, point_key},
    tapered_thread::interval::I,
};
use std::rc::Rc;
struct Stock {
    geometry: Rc<OrthogonalSolid>,
    placement: StockClipPlacement,
    raw_faces: Vec<usize>,
    raw_vertices: BTreeMap<[u64; 3], u32>,
    maximum_vertex_error_mm: f64,
}
fn scalar(value: f64, spelling: &str) -> I {
    let mantissa = spelling.split_once('E').map_or(spelling, |(m, _)| m);
    let digits = mantissa
        .bytes()
        .filter(u8::is_ascii_digit)
        .collect::<Vec<_>>();
    if !digits.is_empty() && digits.iter().all(|b| *b == b'0') {
        I::point(0.)
    } else {
        I::new(value.next_down(), value.next_up())
    }
}
fn raw_point(
    boundary: &Boundary<'_, '_>,
    vertex: u32,
    budget: &mut Budget<'_>,
) -> std::result::Result<[I; 3], Error> {
    budget.charge(128)?;
    let [_, point] = boundary.shell.args(vertex, "VERTEX_POINT")? else {
        return Err(Error::Invalid);
    };
    let point = reference(point)?;
    let entity = boundary.shell.simple(point)?;
    let [Value::String(_), values] = boundary.shell.args(point, "CARTESIAN_POINT")? else {
        return Err(Error::Invalid);
    };
    let values = numbers(values)?;
    let source = entity.source.ok_or(Error::Invalid)?;
    budget.charge(source.len())?;
    let spellings = part21::numeric_spellings(entity).ok_or(Error::Invalid)?;
    let [x, y, z] = spellings.as_slice() else {
        return Err(Error::Invalid);
    };
    let result = std::array::from_fn(|i| scalar(values[i], [x, y, z][i]));
    if result.iter().any(|i| !i.valid()) {
        return Err(Error::Unsupported);
    }
    Ok(result)
}
fn difference(a: [I; 3], b: [I; 3]) -> f64 {
    (0..3)
        .fold(I::point(0.), |sum, i| {
            sum + I::point((a[i] - b[i]).magnitude())
        })
        .hi
}
fn world(placement: StockClipPlacement, p: Point) -> std::result::Result<[I; 3], Error> {
    placement
        .point_bounds(p.map(|v| [v, v]))
        .map(|p| p.map(|[lo, hi]| I::new(lo, hi)))
        .map_err(|_| Error::Unsupported)
}
fn stock(
    boundary: &Boundary<'_, '_>,
    position_limit: f64,
    budget: &mut Budget<'_>,
) -> std::result::Result<Stock, Error> {
    if !position_limit.is_finite() || position_limit <= 0. {
        return Err(Error::Invalid);
    }
    let planar = boundary
        .faces
        .iter()
        .enumerate()
        .filter_map(|(i, f)| matches!(f, Face::Planar { .. }).then_some(i))
        .collect::<Vec<_>>();
    // Current material qualification targets the actual small Angle source.
    // This bounds the existing recovery's grid by16^3 cells before it runs;
    // larger independent analytic boundary evidence remains available.
    if planar.len() > 8 {
        return Err(Error::Unsupported);
    }
    let mut indices = BTreeMap::new();
    let mut faces = Vec::new();
    for &original in &planar {
        budget.charge(32)?;
        indices.insert(original, faces.len());
        let Face::Planar {
            origin,
            normal,
            area,
        } = boundary.faces[original]
        else {
            return Err(Error::Invalid);
        };
        faces.push(Face::Planar {
            origin,
            normal,
            area,
        });
    }
    let mut shell = Shell::new(boundary.document.records(), boundary.shell.tolerance, true);
    for (&id, edge) in &boundary.shell.edges {
        budget.charge(64)?;
        if !matches!(edge.curve, Curve::Line)
            || !edge.uses.iter().all(|u| indices.contains_key(&u.face))
        {
            continue;
        }
        if shell.edges.len() == 24 {
            return Err(Error::Unsupported);
        }
        for &v in &edge.vertices {
            shell.points.insert(v, boundary.shell.points[&v]);
        }
        shell.edges.insert(
            id,
            Edge {
                vertices: edge.vertices,
                curve: edge.curve,
                uses: edge
                    .uses
                    .iter()
                    .map(|u| EdgeUse {
                        face: indices[&u.face],
                        forward: u.forward,
                    })
                    .collect(),
            },
        );
    }
    if shell.points.len() > 16 {
        return Err(Error::Unsupported);
    }
    // Every original vertex must meet all three plane families. Therefore
    // local_points' face joins give each coordinate an actual source face
    // plane, and grid dimensions are bounded by these three face counts.
    let mut representatives: Vec<Point> = Vec::new();
    let mut families = Vec::new();
    let mut counts = [0usize; 3];
    for face in &faces {
        budget.charge(256)?;
        let Face::Planar { normal, .. } = face else {
            return Err(Error::Invalid);
        };
        let family = representatives
            .iter()
            .position(|r| dot(*r, *normal).abs() > 1. - 1e-10);
        let family = if let Some(family) = family {
            family
        } else {
            if representatives.len() == 3
                || representatives
                    .iter()
                    .any(|r| dot(*r, *normal).abs() > 1e-10)
            {
                return Err(Error::Unsupported);
            }
            representatives.push(*normal);
            representatives.len() - 1
        };
        families.push(family);
        counts[family] += 1;
    }
    if representatives.len() != 3 || counts.iter().any(|n| *n < 2) {
        return Err(Error::Unsupported);
    }
    for vertex in shell.points.keys() {
        let mut incident = [false; 3];
        for edge in shell.edges.values() {
            budget.charge(64)?;
            if edge.vertices.contains(vertex) {
                for u in &edge.uses {
                    incident[families[u.face]] = true;
                }
            }
        }
        if incident.iter().any(|v| !*v) {
            return Err(Error::Unsupported);
        }
    }
    let grid = counts.into_iter().product::<usize>();
    let v = shell.points.len();
    let e = shell.edges.len();
    let f = faces.len();
    // Finite grid painting/component/boundary primitives plus bounded frame,
    // sorting and original vertex/edge/face bookkeeping. No stage renewal.
    budget.charge(8192 * grid + 1024 * (v * v + e * e + f * f))?;

    let recovered = shell.orthogonal(&faces)?;
    let ImportProfile::Orthogonal { solid } = recovered.profile else {
        return Err(Error::Unsupported);
    };
    let w = cross(recovered.u, recovered.v);
    let matrix = std::array::from_fn(|i| [recovered.u[i], recovered.v[i], w[i]]);
    let placement =
        StockClipPlacement::new(matrix, recovered.origin_mm).map_err(|_| Error::Unsupported)?;
    let geometry = Rc::new(solid);
    if geometry.quads().len() > 256 || geometry.cells().len() > 64 {
        return Err(Error::Unsupported);
    }
    budget.charge(128 * geometry.quads().len() * geometry.quads().len())?;
    let recovered_boundary =
        orthogonal_boundary::build(geometry.quads()).map_err(|_| Error::Unsupported)?;
    if recovered_boundary.faces.len() != planar.len() {
        return Err(Error::Unsupported);
    }
    let mut raw_vertices = BTreeMap::new();
    let mut used = BTreeSet::new();
    let mut maximum_vertex_error_mm = 0f64;
    for p in recovered_boundary
        .faces
        .iter()
        .flat_map(|f| &f.rings)
        .flatten()
    {
        let key = point_key(*p);
        if raw_vertices.contains_key(&key) {
            continue;
        }
        let expected = world(placement, *p)?;
        let mut found = None;
        for &vertex in shell.points.keys() {
            let actual = raw_point(boundary, vertex, budget)?;
            let error = difference(expected, actual);
            if error <= position_limit && found.replace((vertex, error)).is_some() {
                return Err(Error::Unsupported);
            }
        }
        let (vertex, error) = found.ok_or(Error::Unsupported)?;
        if !used.insert(vertex) {
            return Err(Error::Unsupported);
        }
        raw_vertices.insert(key, vertex);
        maximum_vertex_error_mm = maximum_vertex_error_mm.max(error);
    }
    if used.len() != shell.points.len() {
        return Err(Error::Unsupported);
    }
    let mut raw_faces = Vec::new();
    let mut used_faces = BTreeSet::new();
    for face in &recovered_boundary.faces {
        let rings = face
            .rings
            .iter()
            .map(|r| {
                r.iter()
                    .map(|p| raw_vertices[&point_key(*p)])
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let mut found = None;
        for &original in &planar {
            let mut successors = BTreeMap::new();
            for edge in shell.edges.values() {
                budget.charge(32)?;
                if let Some(u) = edge.uses.iter().find(|u| u.face == indices[&original]) {
                    let [a, b] = if u.forward {
                        edge.vertices
                    } else {
                        [edge.vertices[1], edge.vertices[0]]
                    };
                    if successors.insert(a, b).is_some() {
                        return Err(Error::Invalid);
                    }
                }
            }
            let count = rings.iter().map(Vec::len).sum::<usize>();
            let matches = successors.len() == count
                && rings.iter().all(|r| {
                    r.iter()
                        .zip(r.iter().cycle().skip(1))
                        .all(|(a, b)| successors.get(a) == Some(b))
                });
            if matches && found.replace(original).is_some() {
                return Err(Error::Unsupported);
            }
        }
        let original = found.ok_or(Error::Unsupported)?;
        if !used_faces.insert(original) {
            return Err(Error::Unsupported);
        }
        raw_faces.push(original);
    }
    if used_faces.len() != planar.len() {
        return Err(Error::Unsupported);
    }
    Ok(Stock {
        geometry,
        placement,
        raw_faces,
        raw_vertices,
        maximum_vertex_error_mm,
    })
}

#[test]
fn actual_cardinal_and_posed_stock_recovery_preserve_original_owner_cycles() {
    for placement in [
        StockClipPlacement::new([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]], [0.; 3]).unwrap(),
        StockClipPlacement::new(
            [[0.8, -0.6, 0.], [0.6, 0.8, 0.], [0., 0., 1.]],
            [11., -7., 5.],
        )
        .unwrap(),
    ] {
        let text = super::tests::text_placed(placement);
        let document =
            CompoundDocument::parse(&text, owned_insert::DocumentOptions::default(), None).unwrap();
        let boundary = decode(
            &document,
            document.represented_members()[0].item_entity(),
            Options::default(),
            &|| false,
        )
        .unwrap();
        let mut budget = Budget {
            work: 0,
            maximum: 16_000_000,
            interrupt: &|| false,
        };
        let source = stock(&boundary, 1e-7 / 16., &mut budget).unwrap();
        assert_eq!(source.raw_faces.len(), 8);
        assert_eq!(source.raw_vertices.len(), 12);
        assert!((source.geometry.volume_mm3() - 15000.).abs() < 1e-8);
        assert!(source.maximum_vertex_error_mm <= 1e-7 / 16.);
        let determinant = source.placement.determinant_bounds();
        assert!(determinant[0] <= 1. && determinant[1] >= 1.);
    }
}

struct Geometry {
    source: crate::orthogonal_bore_seats::Seats,
    stock: Stock,
    /// Original finite source owners in derived bore order. No index equality
    /// with an authored source is inferred from this reconstruction.
    raw_chains: Vec<usize>,
}
fn local(placement: StockClipPlacement, p: Point) -> std::result::Result<Point, Error> {
    placement
        .inverse_point_bounds(p.map(|v| [v, v]))
        .map(|p| p.map(|[lo, hi]| lo.midpoint(hi)))
        .map_err(|_| Error::Unsupported)
}
#[cfg(test)]
fn geometry(
    chains: &Chains<'_, '_, '_>,
    position_limit: f64,
    budget: &mut Budget<'_>,
) -> std::result::Result<Geometry, Error> {
    geometry_options(chains, position_limit, budget, MaterialOptions::default())
}
fn geometry_options(
    chains: &Chains<'_, '_, '_>,
    position_limit: f64,
    budget: &mut Budget<'_>,
    options: MaterialOptions,
) -> std::result::Result<Geometry, Error> {
    use crate::{
        orthogonal_bore_seats::{Error as SourceError, Options as SourceOptions, Seat, Seats},
        orthogonal_pockets::OrthogonalThroughBore,
    };
    let boundary = chains.boundary();
    let stock = stock(boundary, position_limit, budget)?;
    let source_boundary =
        orthogonal_boundary::build(stock.geometry.quads()).map_err(|_| Error::Unsupported)?;
    let mut bores = Vec::new();
    let mut seats = Vec::new();
    let mut raw_chains = Vec::new();
    for (index, chain) in chains.bores().iter().enumerate() {
        budget.charge(256)?;
        let caps = chain.cap_face_entities().map(|id| {
            stock
                .raw_faces
                .iter()
                .position(|raw| boundary.face_entities[*raw] == id)
                .ok_or(Error::Unsupported)
        });
        let [a, b] = caps;
        let caps = [a?, b?];
        let normal = source_boundary.faces[caps[0]].normal;
        let axis = normal
            .iter()
            .position(|n| n.abs() == 1.)
            .ok_or(Error::Unsupported)?;
        let centers = chain
            .rim_circles()
            .iter()
            .map(|circle| local(stock.placement, circle.center_mm))
            .collect::<std::result::Result<Vec<_>, _>>()?;
        let first = centers[0];
        let last = *centers.last().ok_or(Error::Unsupported)?;
        let transverse = std::array::from_fn(|i| first[i].midpoint(last[i]));
        let endpoints = caps.map(|face| {
            let mut p = transverse;
            p[axis] = source_boundary.faces[face].rings[0][0][axis];
            p
        });
        let cylinder = chain.cylinder_index();
        let radius = chain.rim_circles()[cylinder].radius_mm;
        let bore = bores.len();
        bores.push(OrthogonalThroughBore {
            faces: caps,
            centers_mm: endpoints,
            radius_mm: radius,
        });
        if cylinder == 1 {
            seats.push(Seat {
                bore,
                end: false,
                radius_mm: chain.rim_circles()[0].radius_mm,
                depth_mm: (endpoints[0][axis] - centers[1][axis]).abs(),
            });
        }
        if cylinder + 1 < chain.wall_face_entities().len() {
            seats.push(Seat {
                bore,
                end: true,
                radius_mm: chain
                    .rim_circles()
                    .last()
                    .ok_or(Error::Unsupported)?
                    .radius_mm,
                depth_mm: (endpoints[1][axis] - centers[cylinder + 1][axis]).abs(),
            });
        }
        raw_chains.push(index);
    }
    let remaining = budget.maximum - budget.work;
    if remaining == 0 {
        return Err(Error::Limit);
    }
    let source = Seats::new(
        Rc::clone(&stock.geometry),
        Rc::from(bores),
        &seats,
        SourceOptions {
            clearance_mm: options.clearance_mm,
            minimum_feature_mm: options.minimum_feature_mm,
            maximum_work: remaining,
            maximum_bytes: options
                .transfer
                .maximum_bytes
                .checked_sub(MATERIAL_RESERVE)
                .filter(|n| *n > 0)
                .ok_or(Error::Limit)?,
        },
        budget.interrupt,
    )
    .map_err(|e| match e {
        SourceError::Invalid => Error::Invalid,
        SourceError::Unsupported => Error::Unsupported,
        SourceError::Limit => Error::Limit,
        SourceError::Cancelled => Error::Cancelled,
    })?;
    budget.charge(source.work_used)?;
    Ok(Geometry {
        source,
        stock,
        raw_chains,
    })
}
#[test]
fn derive_perpendicular_receiving_geometry_from_actual_finite_chains() {
    for placement in [
        StockClipPlacement::new([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]], [0.; 3]).unwrap(),
        StockClipPlacement::new(
            [[0.8, -0.6, 0.], [0.6, 0.8, 0.], [0., 0., 1.]],
            [11., -7., 5.],
        )
        .unwrap(),
    ] {
        let text = super::tests::text_placed(placement);
        let document =
            CompoundDocument::parse(&text, owned_insert::DocumentOptions::default(), None).unwrap();
        let boundary = decode(
            &document,
            document.represented_members()[0].item_entity(),
            Options::default(),
            &|| false,
        )
        .unwrap();
        let chains = boundary
            .through_chains(Options::default(), &|| false)
            .unwrap();
        let mut budget = Budget {
            work: 0,
            maximum: 16_000_000,
            interrupt: &|| false,
        };
        let derived = geometry(&chains, 1e-7 / 16., &mut budget).unwrap();
        assert_eq!(derived.source.bores().len(), 2);
        assert_eq!(derived.source.seats().len(), 2);
        assert_eq!(derived.raw_chains, [0, 1]);
        assert_eq!(derived.stock.raw_faces.len(), 8);
        let actual = super::super::read_analytic_solid(&text, 1e-7).unwrap();
        let local = derived.source.volume_bounds_mm3();
        let d = derived.stock.placement.determinant_bounds();
        let physical = I::new(local[0], local[1]) * I::new(d[0], d[1]);
        assert!((actual.volume_mm3 - physical.midpoint()).abs() < 1e-7);
    }
}

const MATERIAL_RESERVE: usize = 1024 * 1024;
#[derive(Clone, Copy)]
pub struct TransferOptions {
    pub maximum_position_error_mm: f64,
    /// Circles use theta/(2pi), walls their finite normalized axial span,
    /// lines their endpoint span, and planes common owner-ID barycentric patches.
    pub maximum_normalized_derivative_error_mm: f64,
    pub maximum_sewing_error_mm: f64,
    pub maximum_work: usize,
    /// Additional query storage; immutable document/boundary prerequisites
    /// retain their separate ledgers and remain live during this operation.
    pub maximum_bytes: usize,
}
impl Default for TransferOptions {
    fn default() -> Self {
        Self {
            maximum_position_error_mm: 1e-7,
            maximum_normalized_derivative_error_mm: 1e-7,
            maximum_sewing_error_mm: 1e-7,
            maximum_work: 16_000_000,
            maximum_bytes: 64 * 1024 * 1024,
        }
    }
}
impl TransferOptions {
    pub(crate) fn validate(self) -> std::result::Result<(), Error> {
        if [
            self.maximum_position_error_mm,
            self.maximum_normalized_derivative_error_mm,
            self.maximum_sewing_error_mm,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v <= 0.)
            || self.maximum_work == 0
            || self.maximum_work > 16_000_000
            || self.maximum_bytes == 0
            || self.maximum_bytes > 64 * 1024 * 1024
        {
            return Err(Error::Invalid);
        }
        Ok(())
    }
    fn check(self, report: &fields::Report) -> std::result::Result<(), Error> {
        if report.position > self.maximum_position_error_mm
            || report.derivative > self.maximum_normalized_derivative_error_mm
            || report.sewing > self.maximum_sewing_error_mm
        {
            return Err(Error::Unsupported);
        }
        Ok(())
    }
}
#[derive(Clone, Copy)]
pub struct MaterialOptions {
    pub transfer: TransferOptions,
    pub clearance_mm: f64,
    pub minimum_feature_mm: f64,
}
impl MaterialOptions {
    pub(crate) fn validate(self) -> std::result::Result<(), Error> {
        self.transfer.validate()?;
        if !self.clearance_mm.is_finite()
            || self.clearance_mm < 0.
            || !self.minimum_feature_mm.is_finite()
            || self.minimum_feature_mm <= 0.
        {
            return Err(Error::Invalid);
        }
        Ok(())
    }
}
impl Default for MaterialOptions {
    fn default() -> Self {
        Self {
            transfer: TransferOptions::default(),
            clearance_mm: 1e-7,
            minimum_feature_mm: 1e-5,
        }
    }
}
#[derive(Clone, Copy)]
pub struct Correspondence {
    position: f64,
    derivative: f64,
    sewing: f64,
    faces: usize,
    edges: usize,
    vertices: usize,
}
impl Correspondence {
    pub fn maximum_position_error_mm(&self) -> f64 {
        self.position
    }
    pub fn maximum_normalized_derivative_error_mm(&self) -> f64 {
        self.derivative
    }
    pub fn maximum_sewing_error_mm(&self) -> f64 {
        self.sewing
    }
    pub fn original_face_count(&self) -> usize {
        self.faces
    }
    pub fn original_edge_count(&self) -> usize {
        self.edges
    }
    pub fn original_vertex_count(&self) -> usize {
        self.vertices
    }
}
impl From<fields::Report> for Correspondence {
    fn from(r: fields::Report) -> Self {
        Self {
            position: r.position,
            derivative: r.derivative,
            sewing: r.sewing,
            faces: r.faces,
            edges: r.edges,
            vertices: r.vertices,
        }
    }
}
/// Closed, embedded material derived from original decoded stock and finite
/// rim records. It is an explicit bounded healing of raw supports/trims, not
/// the literal raw STEP occupied set or its mass.
pub struct HealedMaterial<'boundary, 'document, 'input> {
    chains: Chains<'boundary, 'document, 'input>,
    geometry: Geometry,
    correspondence: Correspondence,
    volume: [f64; 2],
    work: usize,
    bytes: usize,
}
impl<'boundary, 'document, 'input> Boundary<'document, 'input> {
    pub fn heal_material(
        &'boundary self,
        options: MaterialOptions,
        interrupt: &dyn Fn() -> bool,
    ) -> std::result::Result<HealedMaterial<'boundary, 'document, 'input>, Error> {
        let mut budget = Budget {
            work: 0,
            maximum: options.transfer.maximum_work,
            interrupt,
        };
        budget.charge(0)?;
        options.validate()?;
        if options.transfer.maximum_bytes <= MATERIAL_RESERVE {
            return Err(Error::Limit);
        }
        let chains = self.through_chains(
            Options {
                maximum_work: options.transfer.maximum_work,
                maximum_bytes: options.transfer.maximum_bytes - MATERIAL_RESERVE,
            },
            interrupt,
        )?;
        budget.charge(chains.work_used())?;
        let position = (options.transfer.maximum_position_error_mm / 16.).min(1e-7 / 16.);
        let geometry = geometry_options(&chains, position, &mut budget, options)?;
        let report =
            fields::correspondence(&fields::View::derived(&geometry), &chains, &mut budget)?;
        options.transfer.check(&report)?;
        let local = geometry.source.volume_bounds_mm3();
        let determinant = geometry.stock.placement.determinant_bounds();
        let volume = I::new(local[0], local[1]) * I::new(determinant[0], determinant[1]);
        if !volume.valid() || volume.lo <= 0. {
            return Err(Error::Unsupported);
        }
        let bytes = geometry
            .source
            .accounted_peak_bytes
            .checked_add(MATERIAL_RESERVE)
            .ok_or(Error::Limit)?;
        if bytes > options.transfer.maximum_bytes {
            return Err(Error::Limit);
        }
        budget.charge(0)?;
        Ok(HealedMaterial {
            chains,
            geometry,
            correspondence: report.into(),
            volume: [volume.lo, volume.hi],
            work: budget.work,
            bytes,
        })
    }
}
pub struct SourceMatch {
    correspondence: Correspondence,
    original_volume: [f64; 2],
    composed_position: f64,
    composed_derivative: f64,
    work: usize,
    bytes: usize,
    known: usize,
}
impl SourceMatch {
    pub fn source_to_raw_correspondence(&self) -> Correspondence {
        self.correspondence
    }
    pub fn original_source_volume_bounds_mm3(&self) -> [f64; 2] {
        self.original_volume
    }
    pub fn source_to_healed_position_error_mm(&self) -> f64 {
        self.composed_position
    }
    pub fn source_to_healed_normalized_derivative_error_mm(&self) -> f64 {
        self.composed_derivative
    }
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn query_budget_bytes(&self) -> usize {
        self.bytes
    }
    pub fn aggregate_known_work(&self) -> usize {
        self.known
    }
}
impl<'boundary, 'document, 'input> HealedMaterial<'boundary, 'document, 'input> {
    pub fn boundary(&self) -> &'boundary Boundary<'document, 'input> {
        self.chains.boundary()
    }
    pub fn chains(&self) -> &Chains<'boundary, 'document, 'input> {
        &self.chains
    }
    pub fn healed_source(&self) -> &crate::orthogonal_bore_seats::Seats {
        &self.geometry.source
    }
    pub fn healed_placement(&self) -> StockClipPlacement {
        self.geometry.stock.placement
    }
    pub fn healed_volume_bounds_mm3(&self) -> [f64; 2] {
        self.volume
    }
    pub fn raw_to_healed_correspondence(&self) -> Correspondence {
        self.correspondence
    }
    pub fn stock_vertex_correspondence_mm(&self) -> f64 {
        self.geometry.stock.maximum_vertex_error_mm
    }
    /// Includes finite-chain work once, excluding the borrowed boundary/document.
    pub fn work_used(&self) -> usize {
        self.work
    }
    pub fn query_budget_bytes(&self) -> usize {
        self.bytes
    }
    pub fn aggregate_known_work(&self) -> usize {
        self.boundary().document.work_used() + self.boundary().work_used() + self.work
    }
    /// Original source coefficients and its literal affine map are checked
    /// separately against the same immutable raw finite fields as the healing.
    /// expected_faces is the producer's original source-order face ownership.
    pub fn check_original(
        &self,
        source: &crate::orthogonal_bore_seats::Seats,
        placement: StockClipPlacement,
        expected_faces: &[u32],
        options: TransferOptions,
        interrupt: &dyn Fn() -> bool,
    ) -> std::result::Result<SourceMatch, Error> {
        let mut budget = Budget {
            work: 0,
            maximum: options.maximum_work,
            interrupt,
        };
        budget.charge(0)?;
        options.validate()?;
        if options.maximum_bytes < MATERIAL_RESERVE {
            return Err(Error::Limit);
        }
        let owners = fields::original_owners(
            source,
            placement,
            expected_faces,
            &self.geometry,
            &self.chains,
            (options.maximum_position_error_mm / 16.).min(1e-7 / 16.),
            &mut budget,
        )?;
        let view = fields::View {
            source,
            placement,
            raw_faces: &owners.faces,
            raw_vertices: &owners.vertices,
            raw_chains: &owners.chains,
            reversed: &owners.reversed,
        };
        let report = fields::correspondence(&view, &self.chains, &mut budget)?;
        options.check(&report)?;
        let position = (I::point(report.position) + I::point(self.correspondence.position)).hi;
        let derivative =
            (I::point(report.derivative) + I::point(self.correspondence.derivative)).hi;
        if position > options.maximum_position_error_mm
            || derivative > options.maximum_normalized_derivative_error_mm
        {
            return Err(Error::Unsupported);
        }
        let local = source.volume_bounds_mm3();
        let determinant = placement.determinant_bounds();
        let volume = I::new(local[0], local[1]) * I::new(determinant[0], determinant[1]);
        if !volume.valid() || volume.lo <= 0. {
            return Err(Error::Unsupported);
        }
        budget.charge(0)?;
        let original_prerequisite = if std::ptr::eq(source, &self.geometry.source) {
            0
        } else {
            source.work_used
        };
        Ok(SourceMatch {
            correspondence: report.into(),
            original_volume: [volume.lo, volume.hi],
            composed_position: position,
            composed_derivative: derivative,
            work: budget.work,
            bytes: MATERIAL_RESERVE,
            known: self
                .aggregate_known_work()
                .checked_add(original_prerequisite)
                .and_then(|n| n.checked_add(budget.work))
                .ok_or(Error::Limit)?,
        })
    }
}

#[test]
fn public_healed_material_and_original_gate_preserve_separate_mass_and_limits() {
    use crate::step::orthogonal_seats;
    let original = super::tests::source();
    let placement = StockClipPlacement::new(
        [[0.8, -0.6, 0.], [0.6, 0.8, 0.], [0., 0., 1.]],
        [11., -7., 5.],
    )
    .unwrap();
    let candidate = orthogonal_seats::candidate(
        &original,
        placement,
        orthogonal_seats::Options::default(),
        &|| false,
    )
    .unwrap();
    let document = CompoundDocument::parse(
        candidate.as_str(),
        owned_insert::DocumentOptions::default(),
        None,
    )
    .unwrap();
    let boundary = decode(
        &document,
        document.represented_members()[0].item_entity(),
        Options::default(),
        &|| false,
    )
    .unwrap();
    let material = boundary
        .heal_material(MaterialOptions::default(), &|| false)
        .unwrap();
    let matched = material
        .check_original(
            &original,
            placement,
            candidate.face_entities(),
            TransferOptions::default(),
            &|| false,
        )
        .unwrap();
    assert!(std::ptr::eq(material.boundary(), &boundary));
    assert!(std::ptr::eq(material.chains().boundary(), &boundary));
    assert_eq!(material.healed_source().bores().len(), 2);
    let c = material.raw_to_healed_correspondence();
    assert_eq!(c.original_face_count(), 12);
    assert_eq!(c.original_edge_count(), boundary.edge_count());
    assert_eq!(c.original_vertex_count(), boundary.vertex_count());
    assert!(
        c.maximum_position_error_mm() < 1e-7
            && c.maximum_normalized_derivative_error_mm() < 1e-7
            && c.maximum_sewing_error_mm() < 1e-7
    );
    assert!(material.stock_vertex_correspondence_mm() < 1e-7 / 16.);
    assert!(
        matched.source_to_healed_position_error_mm() < 1e-7
            && matched.source_to_healed_normalized_derivative_error_mm() < 1e-7
    );
    let original_volume = original.volume_bounds_mm3();
    let d = placement.determinant_bounds();
    let v = I::new(original_volume[0], original_volume[1]) * I::new(d[0], d[1]);
    assert_eq!(matched.original_source_volume_bounds_mm3(), [v.lo, v.hi]);
    let healed = material.healed_volume_bounds_mm3();
    assert!(healed[0] > 0. && healed[1] >= healed[0]);
    assert!(material.healed_placement().determinant_bounds()[0] > 0.);
    assert_eq!(
        material.aggregate_known_work(),
        document.work_used() + boundary.work_used() + material.work_used()
    );
    assert_eq!(
        matched.aggregate_known_work(),
        material.aggregate_known_work() + original.work_used + matched.work_used()
    );
    let exact = MaterialOptions {
        transfer: TransferOptions {
            maximum_work: material.work_used(),
            maximum_bytes: material.query_budget_bytes(),
            ..TransferOptions::default()
        },
        ..MaterialOptions::default()
    };
    assert!(boundary.heal_material(exact, &|| false).is_ok());
    assert!(matches!(
        boundary.heal_material(
            MaterialOptions {
                transfer: TransferOptions {
                    maximum_work: exact.transfer.maximum_work - 1,
                    ..exact.transfer
                },
                ..exact
            },
            &|| false
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        boundary.heal_material(
            MaterialOptions {
                transfer: TransferOptions {
                    maximum_bytes: exact.transfer.maximum_bytes - 1,
                    ..exact.transfer
                },
                ..exact
            },
            &|| false
        ),
        Err(Error::Limit)
    ));
    let match_exact = TransferOptions {
        maximum_work: matched.work_used(),
        maximum_bytes: matched.query_budget_bytes(),
        ..TransferOptions::default()
    };
    assert!(
        material
            .check_original(
                &original,
                placement,
                candidate.face_entities(),
                match_exact,
                &|| false
            )
            .is_ok()
    );
    assert!(matches!(
        material.check_original(
            &original,
            placement,
            candidate.face_entities(),
            TransferOptions {
                maximum_work: match_exact.maximum_work - 1,
                ..match_exact
            },
            &|| false
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        material.check_original(
            &original,
            placement,
            candidate.face_entities(),
            TransferOptions {
                maximum_bytes: match_exact.maximum_bytes - 1,
                ..match_exact
            },
            &|| false
        ),
        Err(Error::Limit)
    ));
    assert!(matches!(
        material.check_original(
            &original,
            placement,
            candidate.face_entities(),
            TransferOptions {
                maximum_normalized_derivative_error_mm: matched
                    .source_to_healed_normalized_derivative_error_mm()
                    .next_down(),
                ..TransferOptions::default()
            },
            &|| false
        ),
        Err(Error::Unsupported)
    ));
    let mut wrong = candidate.face_entities().to_vec();
    let offset = original.boundary().faces.len();
    wrong.swap(offset, offset + 1);
    assert!(matches!(
        material.check_original(
            &original,
            placement,
            &wrong,
            TransferOptions::default(),
            &|| false
        ),
        Err(Error::Unsupported)
    ));
    wrong = candidate.face_entities().to_vec();
    wrong.swap(
        offset + original.bores().len(),
        offset + original.bores().len() + 1,
    );
    assert!(matches!(
        material.check_original(
            &original,
            placement,
            &wrong,
            TransferOptions::default(),
            &|| false
        ),
        Err(Error::Unsupported)
    ));
    assert!(matches!(
        boundary.heal_material(
            MaterialOptions {
                minimum_feature_mm: 0.,
                ..MaterialOptions::default()
            },
            &|| false
        ),
        Err(Error::Invalid)
    ));
    assert!(matches!(
        boundary.heal_material(
            MaterialOptions {
                transfer: TransferOptions {
                    maximum_position_error_mm: f64::NAN,
                    ..TransferOptions::default()
                },
                ..MaterialOptions::default()
            },
            &|| false
        ),
        Err(Error::Invalid)
    ));
    use std::cell::Cell;
    let calls = Cell::new(0);
    boundary
        .heal_material(MaterialOptions::default(), &|| {
            calls.set(calls.get() + 1);
            false
        })
        .unwrap();
    let stop = calls.get() - 3;
    calls.set(0);
    assert!(matches!(
        boundary.heal_material(MaterialOptions::default(), &|| {
            calls.set(calls.get() + 1);
            calls.get() > stop
        }),
        Err(Error::Cancelled)
    ));
    assert!(boundary.heal_material(exact, &|| false).is_ok());
    assert_eq!(document.literal_input(), candidate.as_str());
}
#[test]
fn actual_cone_semiaxis_change_cannot_hide_behind_unchanged_rims() {
    let original = super::tests::text_placed(
        StockClipPlacement::new([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]], [0.; 3]).unwrap(),
    );
    let line = original
        .lines()
        .find(|l| l.contains("CONICAL_SURFACE("))
        .unwrap();
    let (prefix, angle) = line.rsplit_once(',').unwrap();
    let angle = angle.strip_suffix(");").unwrap().parse::<f64>().unwrap();
    let changed_line = format!("{prefix},{:.16E});", angle + 4e-9);
    let changed = original.replacen(line, &changed_line, 1);
    assert_ne!(original, changed);
    let document =
        CompoundDocument::parse(&changed, owned_insert::DocumentOptions::default(), None).unwrap();
    let boundary = decode(
        &document,
        document.represented_members()[0].item_entity(),
        Options::default(),
        &|| false,
    )
    .unwrap();
    assert!(matches!(
        boundary.heal_material(MaterialOptions::default(), &|| false),
        Err(Error::Unsupported)
    ));
    let material = boundary
        .heal_material(
            MaterialOptions {
                transfer: TransferOptions {
                    maximum_normalized_derivative_error_mm: 1e-6,
                    ..TransferOptions::default()
                },
                ..MaterialOptions::default()
            },
            &|| false,
        )
        .unwrap();
    assert!(
        material
            .raw_to_healed_correspondence()
            .maximum_normalized_derivative_error_mm()
            > 1e-7
    );
    assert!(
        material
            .raw_to_healed_correspondence()
            .maximum_position_error_mm()
            < 1e-7
    );
    assert!(material.healed_volume_bounds_mm3()[0] > 0.);
}
