//! Rust transfer data and explicit refusals for the removed comparison provider.
use std::{convert::Infallible, fmt, marker::PhantomData, ops::Deref, pin::Pin};

#[derive(Debug)]
pub(crate) struct Unavailable;

impl fmt::Display for Unavailable {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str("OpenCascade is unavailable in this build")
    }
}

impl std::error::Error for Unavailable {}

type Result<T> = std::result::Result<T, Unavailable>;

/// Geometry cannot create a native handle in this configuration.
pub(crate) struct UniquePtr<T>(Infallible, PhantomData<T>);

impl<T> UniquePtr<T> {
    pub(crate) fn is_null(&self) -> bool {
        match self.0 {}
    }

    pub(crate) fn pin_mut(&mut self) -> Pin<&mut T> {
        match self.0 {}
    }
}

impl<T> Deref for UniquePtr<T> {
    type Target = T;

    fn deref(&self) -> &T {
        match self.0 {}
    }
}

pub struct Mesh {
    pub vertices: Vec<f64>,
    pub indices: Vec<u32>,
    #[doc = " One evaluation-local B-rep face index per display triangle."]
    pub face_indices: Vec<u32>,
}
pub struct ManufacturingMesh {
    pub mesh: Mesh,
    pub maximum_displacement_mm: f64,
    pub edge_pairs: u32,
    pub vertex_links: u32,
    pub removed_native_pole_triangles: u32,
    pub removed_exact_zero_triangles: u32,
    pub exact_vertex_pairs: u32,
    pub internal_constraint_edges: u32,
    pub removed_meshing_internal_edges: u32,
    pub removed_meshing_constraint_wire_containers: u32,
    pub removed_meshing_internal_vertices: u32,
    pub remeshed_without_internal_constraints: bool,
}
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct MeshCacheStats {
    pub requests: u64,
    pub hits: u64,
    pub builds: u64,
    pub entries: u64,
    pub bytes: u64,
    pub evictions: u64,
    pub unretained: u64,
}
#[derive(Clone, Copy, Debug)]
pub struct Point3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}
pub struct Frame3 {
    pub origin: Point3,
    pub u: Point3,
    pub v: Point3,
}
pub struct CurveEdge {
    pub kind: u8,
    pub a: Point3,
    pub mid: Point3,
    pub b: Point3,
    pub radius: f64,
    pub source_index: u32,
}
pub struct Bounds3 {
    pub min: Point3,
    pub max: Point3,
}
#[doc = " One solid read from a STEP file: its representation item name, if"]
#[doc = " any, and its exact located B-rep as OCCT BRep text."]
pub struct StepSolidData {
    pub name: String,
    pub brep: String,
}
pub struct SweepLimits {
    pub tolerance: f64,
    pub minimum_feature: f64,
    pub maximum_length: f64,
}
#[derive(Clone, Copy, Debug)]
pub struct ThreadData {
    pub radius: f64,
    pub pitch: f64,
    pub height: f64,
    pub depth: f64,
    pub crest_width: f64,
    pub rounding: f64,
    pub lower_slope: f64,
    pub upper_slope: f64,
    pub taper: f64,
    pub overlap: f64,
    pub tolerance: f64,
    pub crest_limit: f64,
    pub root_rounding: f64,
    pub knuckle_root_radius: f64,
    pub starts: u32,
    pub left: bool,
    pub axial_flat: bool,
}
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct ThreadFactoryStats {
    pub requests: u64,
    pub hits: u64,
    pub constructions: u64,
    pub continuous_attempts: u64,
    pub bounded_attempts: u64,
    pub entries: u64,
    pub serialized_bytes: u64,
    pub evictions: u64,
    pub unretained: u64,
}
pub struct WormWheelData {
    pub teeth: u32,
    pub module: f64,
    pub width: f64,
    pub bore: f64,
    pub thinning: f64,
    pub tolerance: f64,
    pub minimum: f64,
}
pub struct BevelData {
    pub teeth: u32,
    pub pitch: f64,
    pub base: f64,
    pub tip: f64,
    pub root: f64,
    pub distance: f64,
    pub width: f64,
    pub bore: f64,
    pub half_base: f64,
    pub fillet: f64,
    pub tolerance: f64,
    pub minimum: f64,
}
pub struct GearData {
    pub teeth: u32,
    pub pitch_radius: f64,
    pub base_radius: f64,
    pub tip_radius: f64,
    pub root_radius: f64,
    pub fillet_radius: f64,
    pub half_base: f64,
    pub height: f64,
    pub bore: f64,
    pub outer: f64,
    pub hub_radius: f64,
    pub hub_start: f64,
    pub hub_end: f64,
    pub twist: f64,
    pub gap: f64,
    pub style: u8,
    pub tolerance: f64,
}
pub struct FaceInfo {
    pub index: u32,
    pub planar: bool,
    pub origin: Point3,
    pub normal: Point3,
    pub center: Point3,
    pub area: f64,
}
pub struct TopologyNode {
    pub kind: u8,
    pub index: u32,
    pub incident_faces: Vec<u32>,
}
pub struct Ancestry {
    pub kind: u8,
    pub operand: u8,
    pub output: u32,
    pub input: u32,
}
pub struct Seed {
    pub face: u32,
    pub source_index: u32,
    pub role: u8,
}
pub struct EdgeInfo {
    pub index: u32,
    pub points: Vec<f64>,
}
pub struct VertexInfo {
    pub index: u32,
    pub point: Point3,
}
pub struct AnalyticEdgeData {
    pub index: u32,
    pub kind: u8,
    pub a: Point3,
    pub b: Point3,
    pub frame: Frame3,
    pub radius: f64,
    pub start: f64,
    pub sweep: f64,
}
pub struct SplineData {
    pub degree: u32,
    pub poles: Vec<f64>,
    pub knots: Vec<f64>,
    pub multiplicities: Vec<u32>,
    pub weights: Vec<f64>,
}
pub struct DistanceOptionsData {
    pub a_kind: u8,
    pub a_index: u32,
    pub b_kind: u8,
    pub b_index: u32,
    pub tolerance: f64,
    pub maximum: u32,
}
pub struct DistanceWitnessData {
    pub a: Point3,
    pub b: Point3,
}
pub struct DistanceData {
    pub distance: f64,
    pub total_solutions: u32,
    pub inner: bool,
    pub complete: bool,
    pub witnesses: Vec<DistanceWitnessData>,
}
pub struct DrawingOptionsData {
    pub frame: Frame3,
    pub section: bool,
    pub section_offset: f64,
    pub section_only: bool,
    pub depth: bool,
    pub depth_mm: f64,
    pub window_u: bool,
    pub u_min: f64,
    pub u_max: f64,
    pub window_v: bool,
    pub v_min: f64,
    pub v_max: f64,
    pub detail: bool,
    pub detail_center: Point3,
    pub detail_radius: f64,
    pub tolerance: f64,
    pub hidden: bool,
    pub tangents: bool,
    pub max_curves: u32,
    pub max_points: u32,
}
pub struct DrawingCurveData {
    pub source_index: u32,
    pub role: u8,
    pub kind: u8,
    pub degree: u32,
    pub periodic: bool,
    pub data: Vec<f64>,
    pub knots: Vec<f64>,
    pub multiplicities: Vec<u32>,
    pub preview: Vec<f64>,
}
pub struct PlanarSegmentData {
    pub edge_index: u32,
    pub kind: u8,
    pub a: Point3,
    pub b: Point3,
    pub center: Point3,
    pub sweep: f64,
}
pub struct PlanarWireData {
    pub outer: bool,
    pub segments: Vec<PlanarSegmentData>,
}
pub struct HoleAxisData {
    pub origin: Point3,
    pub axis: Point3,
    pub radius: f64,
}
pub struct CurvedEdgeData {
    pub kind: u8,
    pub a: Point3,
    pub b: Point3,
    pub center: Point3,
    pub radius: f64,
    pub clockwise: bool,
    pub spline: SplineData,
    pub source_index: u32,
}
pub struct CurveSampleData {
    pub point: Point3,
    pub derivative: Point3,
}
pub struct CurveProjectionData {
    pub point: Point3,
    pub tangent: Point3,
    pub parameter: f64,
    pub distance: f64,
}
pub struct SketchSupportData {
    pub kind: u8,
    pub index: u32,
}
pub struct SketchReferenceData {
    pub geometry: u8,
    pub point: Point3,
    pub curve: CurvedEdgeData,
    pub support: Vec<SketchSupportData>,
    pub has_parameter: bool,
    pub parameter: f64,
    pub has_interval: bool,
    pub first: f64,
    pub last: f64,
    pub tangent: bool,
    pub crossing: i8,
}
pub struct SketchReferenceCatalogData {
    pub results: Vec<SketchReferenceData>,
    pub coincident: bool,
}
pub struct LoftSectionData {
    pub frame: Frame3,
    pub edges: Vec<CurvedEdgeData>,
}
pub struct CurveIntersectionData {
    pub point: Point3,
    pub first: f64,
    pub second: f64,
}
pub struct CurveQueryData {
    pub points: Vec<CurveIntersectionData>,
    pub overlaps: u32,
}
pub(crate) fn helical_thread(spec: &ThreadData) -> Result<UniquePtr<Shape>> {
    Err(Unavailable)
}
pub(crate) fn helical_thread_fresh(spec: &ThreadData) -> Result<UniquePtr<Shape>> {
    Err(Unavailable)
}
pub(crate) fn cached_helical_thread(spec: &ThreadData) -> Result<UniquePtr<Shape>> {
    Err(Unavailable)
}
pub(crate) fn helical_thread_after_miss(spec: &ThreadData) -> Result<UniquePtr<Shape>> {
    Err(Unavailable)
}
pub(crate) fn spherical_bevel_gear(spec: &BevelData) -> Result<UniquePtr<Shape>> {
    Err(Unavailable)
}
pub(crate) fn za_worm(spec: &ThreadData, bore: f64) -> Result<UniquePtr<Shape>> {
    Err(Unavailable)
}
pub(crate) fn za_worm_wheel(spec: &WormWheelData) -> Result<UniquePtr<Shape>> {
    Err(Unavailable)
}
pub(crate) fn involute_gear(spec: &GearData) -> Result<UniquePtr<Shape>> {
    Err(Unavailable)
}
pub(crate) fn sketch_reference_geometry(
    shape: &Shape,
    kind: u8,
    index: u32,
    frame: &Frame3,
    method: u8,
    tolerance: f64,
) -> Result<SketchReferenceCatalogData> {
    Err(Unavailable)
}
pub(crate) enum CompoundJob {}
pub(crate) fn compound_job() -> Result<UniquePtr<CompoundJob>> {
    Err(Unavailable)
}
pub(crate) fn compound_add(job: Pin<&mut CompoundJob>, shape: &Shape) -> Result<()> {
    Err(Unavailable)
}
pub(crate) fn compound_finish(job: Pin<&mut CompoundJob>) -> Result<UniquePtr<Shape>> {
    Err(Unavailable)
}
pub(crate) fn compound_property_fallback_probe(shape: &Shape) -> Result<bool> {
    Err(Unavailable)
}
pub(crate) fn solid_component_at(shape: &Shape, index: u32) -> Result<UniquePtr<Shape>> {
    Err(Unavailable)
}
pub(crate) fn shell_solid(
    shape: &Shape,
    faces: &[u32],
    thickness: f64,
    tolerance: f64,
) -> Result<UniquePtr<Shape>> {
    Err(Unavailable)
}
pub(crate) fn draft_solid(
    shape: &Shape,
    faces: &[u32],
    direction: &Point3,
    neutral: &Frame3,
    angle: f64,
) -> Result<UniquePtr<Shape>> {
    Err(Unavailable)
}
pub(crate) fn minimum_distance(
    a: &Shape,
    b: &Shape,
    options: &DistanceOptionsData,
) -> Result<DistanceData> {
    Err(Unavailable)
}
pub(crate) enum Shape {}
pub(crate) enum DrawingJob {}
pub(crate) fn drawing_job(options: &DrawingOptionsData) -> Result<UniquePtr<DrawingJob>> {
    Err(Unavailable)
}
pub(crate) fn drawing_add(
    job: Pin<&mut DrawingJob>,
    shape: &Shape,
    pattern: u8,
    angle: f64,
    spacing: f64,
    line_width: f64,
) -> Result<()> {
    Err(Unavailable)
}
pub(crate) fn planar_face_contours(
    shape: &Shape,
    index: u32,
    frame: &Frame3,
    tolerance: f64,
) -> Result<Vec<PlanarWireData>> {
    Err(Unavailable)
}
pub(crate) fn transfer_hole_axis(
    shape: &Shape,
    kind: u8,
    index: u32,
    tolerance: f64,
) -> Result<HoleAxisData> {
    Err(Unavailable)
}
pub(crate) fn drawing_run(job: Pin<&mut DrawingJob>) -> Result<Vec<DrawingCurveData>> {
    Err(Unavailable)
}
pub(crate) fn extrude(width: f64, depth: f64, height: f64) -> Result<UniquePtr<Shape>> {
    Err(Unavailable)
}
pub(crate) fn cut(
    shape: &Shape,
    x: f64,
    width: f64,
    depth: f64,
    height: f64,
) -> Result<UniquePtr<Shape>> {
    Err(Unavailable)
}
pub(crate) fn mesh(shape: &Shape, deflection: f64) -> Result<Mesh> {
    Err(Unavailable)
}
pub(crate) fn volume(shape: &Shape) -> Result<f64> {
    Err(Unavailable)
}
pub(crate) fn volume_integration_count(shape: &Shape) -> Result<u64> {
    Err(Unavailable)
}
pub(crate) fn shares_volume_cache(a: &Shape, b: &Shape) -> bool {
    false
}
pub(crate) fn volume_cache_concurrent_probe(shape: &Shape) -> Result<bool> {
    Err(Unavailable)
}
pub(crate) fn face_integration_count(shape: &Shape) -> Result<u64> {
    Err(Unavailable)
}
pub(crate) fn shares_face_cache(a: &Shape, b: &Shape) -> bool {
    false
}
pub(crate) fn face_cache_concurrent_probe(shape: &Shape) -> Result<bool> {
    Err(Unavailable)
}
pub(crate) fn thread_factory_stats() -> Result<ThreadFactoryStats> {
    Err(Unavailable)
}
pub(crate) fn thread_meridian_frame_probe() -> Result<bool> {
    Err(Unavailable)
}
pub(crate) fn remember_validated_thread(shape: &Shape) -> Result<bool> {
    Err(Unavailable)
}
pub(crate) fn remember_thread_prototype(shape: &Shape, spec: &ThreadData) -> Result<bool> {
    Err(Unavailable)
}
pub(crate) fn reset_thread_factory_cache() -> Result<()> {
    Err(Unavailable)
}
pub(crate) fn thread_factory_test_limits(entries: u32, bytes: u64) -> Result<()> {
    Err(Unavailable)
}
pub(crate) fn thread_factory_contains(data: &ThreadData) -> Result<bool> {
    Err(Unavailable)
}
pub(crate) fn valid(shape: &Shape) -> Result<bool> {
    Err(Unavailable)
}
pub(crate) fn validity_check_count(shape: &Shape) -> Result<u64> {
    Err(Unavailable)
}
pub(crate) fn shares_validity_cache(a: &Shape, b: &Shape) -> bool {
    false
}
pub(crate) fn validity_cache_concurrent_probe(shape: &Shape) -> Result<bool> {
    Err(Unavailable)
}
pub(crate) fn validity_cache_rejection_probe() -> Result<bool> {
    Err(Unavailable)
}
pub(crate) fn write_step(shape: &Shape, path: &str) -> Result<()> {
    Err(Unavailable)
}
pub(crate) fn read_step(path: &str) -> Result<UniquePtr<Shape>> {
    Err(Unavailable)
}
pub(crate) fn step_solids(path: &str) -> Result<Vec<StepSolidData>> {
    Err(Unavailable)
}
pub(crate) fn shape_to_brep(shape: &Shape) -> Result<String> {
    Err(Unavailable)
}
pub(crate) fn shape_binary_geometry(shape: &Shape, maximum_bytes: u64) -> Result<Vec<u8>> {
    Err(Unavailable)
}
pub(crate) fn shape_from_brep(text: &str) -> Result<UniquePtr<Shape>> {
    Err(Unavailable)
}
pub(crate) fn shape_from_thread_prototype(text: &str, seeds: &[Seed]) -> Result<UniquePtr<Shape>> {
    Err(Unavailable)
}
pub(crate) fn route_native_output_to_stderr() {}
pub(crate) fn cylinder(frame: &Frame3, radius: f64, height: f64) -> Result<UniquePtr<Shape>> {
    Err(Unavailable)
}
pub(crate) fn profile_extrude(
    frame: &Frame3,
    edges: &[CurveEdge],
    offsets: &[u32],
    height: f64,
) -> Result<UniquePtr<Shape>> {
    Err(Unavailable)
}
pub(crate) fn boolean_shapes(
    a: &Shape,
    b: &Shape,
    operation: u8,
    fuzzy: f64,
) -> Result<UniquePtr<Shape>> {
    Err(Unavailable)
}
pub(crate) fn intersect_parallel(a: &Shape, b: &Shape, fuzzy: f64) -> Result<UniquePtr<Shape>> {
    Err(Unavailable)
}
pub(crate) fn clip_half_space(
    shape: &Shape,
    plane: &Frame3,
    positive: bool,
    fuzzy: f64,
) -> Result<UniquePtr<Shape>> {
    Err(Unavailable)
}
pub(crate) fn inspection_union(shape: &Shape) -> Result<UniquePtr<Shape>> {
    Err(Unavailable)
}
pub(crate) fn share_shape(shape: &Shape) -> Result<UniquePtr<Shape>> {
    Err(Unavailable)
}
pub(crate) fn transform_shape(
    shape: &Shape,
    rotation: &[f64],
    translation: &Point3,
) -> Result<UniquePtr<Shape>> {
    Err(Unavailable)
}
pub(crate) fn affine_transform_shape(
    shape: &Shape,
    matrix: &[f64],
    translation: &Point3,
) -> Result<UniquePtr<Shape>> {
    Err(Unavailable)
}
pub(crate) fn mirror_shape(shape: &Shape, plane: &Frame3) -> Result<UniquePtr<Shape>> {
    Err(Unavailable)
}
pub(crate) fn treat_edges(
    shape: &Shape,
    edges: &[u32],
    mode: u8,
    size: f64,
) -> Result<UniquePtr<Shape>> {
    Err(Unavailable)
}
pub(crate) fn treat_edges_sided(
    shape: &Shape,
    edges: &[u32],
    sides: &[u32],
    mode: u8,
    first: f64,
    second: f64,
) -> Result<UniquePtr<Shape>> {
    Err(Unavailable)
}
pub(crate) fn bounds(shape: &Shape) -> Result<Bounds3> {
    Err(Unavailable)
}
pub(crate) fn bounds_cache_concurrent_probe(shape: &Shape) -> Result<bool> {
    Err(Unavailable)
}
pub(crate) fn faces(shape: &Shape) -> Result<Vec<FaceInfo>> {
    Err(Unavailable)
}
pub(crate) fn topology(shape: &Shape) -> Result<Vec<TopologyNode>> {
    Err(Unavailable)
}
pub(crate) fn ancestry(shape: &Shape) -> Result<Vec<Ancestry>> {
    Err(Unavailable)
}
pub(crate) fn seeds(shape: &Shape) -> Result<Vec<Seed>> {
    Err(Unavailable)
}
pub(crate) fn display_edges(shape: &Shape, deflection: f64) -> Result<Vec<EdgeInfo>> {
    Err(Unavailable)
}
pub(crate) fn vertices(shape: &Shape) -> Result<Vec<VertexInfo>> {
    Err(Unavailable)
}
pub(crate) fn analytic_edges(shape: &Shape) -> Result<Vec<AnalyticEdgeData>> {
    Err(Unavailable)
}
pub(crate) fn solid_count(shape: &Shape) -> Result<u32> {
    Err(Unavailable)
}
pub(crate) fn curve_evaluate(curve: &CurvedEdgeData, parameter: f64) -> Result<CurveSampleData> {
    Err(Unavailable)
}
pub(crate) fn curve_project(
    curve: &CurvedEdgeData,
    point: &Point3,
    tolerance: f64,
) -> Result<Vec<CurveProjectionData>> {
    Err(Unavailable)
}
pub(crate) fn spline_segment(curve: &CurvedEdgeData, from: f64, to: f64) -> Result<SplineData> {
    Err(Unavailable)
}
pub(crate) fn curve_intersections(
    first: &CurvedEdgeData,
    second: &CurvedEdgeData,
    tolerance: f64,
    self_query: bool,
) -> Result<CurveQueryData> {
    Err(Unavailable)
}
pub(crate) fn curved_profile_area(
    frame: &Frame3,
    edges: &[CurvedEdgeData],
    offsets: &[u32],
    tolerance: f64,
) -> Result<f64> {
    Err(Unavailable)
}
pub(crate) fn curved_profile_signed_area(
    frame: &Frame3,
    edges: &[CurvedEdgeData],
    offsets: &[u32],
    tolerance: f64,
) -> Result<f64> {
    Err(Unavailable)
}
pub(crate) fn curved_profile_classify(
    frame: &Frame3,
    edges: &[CurvedEdgeData],
    offsets: &[u32],
    tolerance: f64,
    point: &Point3,
) -> Result<u8> {
    Err(Unavailable)
}
pub(crate) fn curved_profile_extrude_vector(
    frame: &Frame3,
    edges: &[CurvedEdgeData],
    offsets: &[u32],
    tolerance: f64,
    displacement: &Point3,
    developed_sheet: bool,
) -> Result<UniquePtr<Shape>> {
    Err(Unavailable)
}
pub(crate) fn curved_profile_revolve(
    frame: &Frame3,
    edges: &[CurvedEdgeData],
    offsets: &[u32],
    tolerance: f64,
    axis_origin: &Point3,
    axis_direction: &Point3,
    angle: f64,
) -> Result<UniquePtr<Shape>> {
    Err(Unavailable)
}
pub(crate) fn curved_profile_to_plane(
    frame: &Frame3,
    edges: &[CurvedEdgeData],
    offsets: &[u32],
    limits: &SweepLimits,
    direction: &Point3,
    termination: &Frame3,
) -> Result<UniquePtr<Shape>> {
    Err(Unavailable)
}
pub(crate) fn curved_profile_ruled(
    start: &Frame3,
    end: &Frame3,
    start_edges: &[CurvedEdgeData],
    end_edges: &[CurvedEdgeData],
    offsets: &[u32],
    tolerance: f64,
) -> Result<UniquePtr<Shape>> {
    Err(Unavailable)
}
pub(crate) fn curved_profile_pipe(
    frame: &Frame3,
    edges: &[CurvedEdgeData],
    path: &[Point3],
    tolerance: f64,
) -> Result<UniquePtr<Shape>> {
    Err(Unavailable)
}
pub(crate) fn curved_profile_loft(
    sections: &[LoftSectionData],
    ruled: bool,
    tolerance: f64,
) -> Result<UniquePtr<Shape>> {
    Err(Unavailable)
}
pub(crate) fn remember_validated_mesh(
    shape: &Shape,
    deflection: f64,
    angle: f64,
    max_triangles: u32,
    mesh: &Mesh,
) -> Result<bool> {
    Err(Unavailable)
}
pub(crate) fn tessellate_manufacturing(
    shape: &Shape,
    deflection: f64,
    angle: f64,
    max_triangles: u32,
    seam_tolerance: f64,
) -> Result<ManufacturingMesh> {
    Err(Unavailable)
}
pub(crate) fn tessellate_fresh(
    shape: &Shape,
    deflection: f64,
    angle: f64,
    max_triangles: u32,
) -> Result<Mesh> {
    Err(Unavailable)
}
pub(crate) fn mesh_cache_stats() -> Result<MeshCacheStats> {
    Err(Unavailable)
}
pub(crate) fn reset_mesh_cache() -> Result<()> {
    Err(Unavailable)
}
pub(crate) fn mesh_cache_test_limits(
    entries: u32,
    total_bytes: u64,
    receipt_bytes: u64,
) -> Result<()> {
    Err(Unavailable)
}
pub(crate) fn shares_mesh_identity(a: &Shape, b: &Shape) -> bool {
    false
}
pub(crate) fn tessellate(
    shape: &Shape,
    deflection: f64,
    angle: f64,
    max_triangles: u32,
) -> Result<Mesh> {
    Err(Unavailable)
}
