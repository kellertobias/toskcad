//! Bounded, transactional proposals. The output remains a candidate until the
//! actual independent decoder and complete source correspondence admit it.
use super::numbers;
use crate::{
    axial_pipe_thread::{
        Error,
        authored::hardware::shell::{self, occupied::OccupiedHead},
    },
    radial_tensor::{self as radial, Axis, Surface},
    step::{Entities, Entity, StepError, StepFile, Writer},
};
use std::{fmt, fmt::Write as _, mem::size_of};
type Result<T> = std::result::Result<T, Error>;
#[derive(Clone, Copy, Debug)]
pub struct Options {
    pub maximum_work: usize,
    pub maximum_bytes: usize,
    pub additional_resident_bytes: usize,
}
impl Default for Options {
    fn default() -> Self {
        Self {
            maximum_work: 32_000_000,
            maximum_bytes: 64 * 1024 * 1024,
            additional_resident_bytes: 0,
        }
    }
}
#[derive(Clone, Copy, Debug)]
pub struct Receipt {
    pub work_used: usize,
    pub output_bytes_upper_bound: usize,
    pub phase_peak_bytes_reserved: usize,
}
pub struct CandidateFile {
    text: String,
    receipt: Receipt,
    pub(super) trace: Vec<Record>,
    pub(super) root: u32,
}
pub(super) struct Record {
    id: u32,
    body: std::ops::Range<usize>,
}
impl CandidateFile {
    pub fn as_str(&self) -> &str {
        &self.text
    }
    pub fn into_string(self) -> String {
        self.text
    }
    pub fn receipt(&self) -> Receipt {
        self.receipt
    }
    pub(super) fn retained_bytes(&self) -> Result<usize> {
        self.text
            .capacity()
            .checked_add(
                self.trace
                    .capacity()
                    .checked_mul(size_of::<Record>())
                    .ok_or(Error::Limit)?,
            )
            .and_then(|n| n.checked_add(size_of::<Self>()))
            .ok_or(Error::Limit)
    }
}
struct AppendError(Error);
impl From<StepError> for AppendError {
    fn from(e: StepError) -> Self {
        Self(step_error(e))
    }
}
fn step_error(e: StepError) -> Error {
    match e {
        StepError::Limit => Error::Limit,
        StepError::InvalidGeometry | StepError::InvalidPlacement => Error::Invalid,
        StepError::Unsupported => Error::Unsupported,
    }
}
pub fn propose(
    source: &OccupiedHead,
    options: Options,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<CandidateFile> {
    let mut writer = Writer::new();
    let (receipt, trace, root) = append_candidate(&mut writer, source, options, interrupt)?;
    let text = writer.finish().map_err(step_error)?;
    if text.len() > receipt.output_bytes_upper_bound {
        return Err(Error::Limit);
    }
    Ok(CandidateFile {
        text,
        receipt,
        trace,
        root,
    })
}
impl crate::step_read::literal_member_binding::RetainedRecords for CandidateFile {
    fn record_count(&self) -> usize {
        self.trace.len()
    }
    fn record_index(&self, id: u32) -> Option<usize> {
        self.trace.binary_search_by_key(&id, |r| r.id).ok()
    }
    fn component(&self, record: usize, index: usize) -> Option<&str> {
        if index != 0 {
            return None;
        }
        self.text.get(self.trace.get(record)?.body.clone())
    }
}
// This private transaction is deliberately used only with the fresh Writer
// owned by propose. It makes no allocation claim about arbitrary old point
// caches left by other Writer operations/rollback histories.
fn append_candidate(
    writer: &mut Writer,
    source: &OccupiedHead,
    options: Options,
    interrupt: Option<&dyn Fn() -> bool>,
) -> Result<(Receipt, Vec<Record>, u32)> {
    if options.maximum_work == 0
        || options.maximum_work > 32_000_000
        || options.maximum_bytes == 0
        || options.maximum_bytes > 64 * 1024 * 1024
    {
        return Err(Error::Invalid);
    }
    let mut receipt = None;
    writer
        .append_typed(|file| {
            let result = (|| {
                let mut measure = Measure {
                    work: Work {
                        remaining: options.maximum_work,
                        interrupt,
                    },
                    bytes: 4096,
                    points: 0,
                };
                measure.input(source.shell().input())?;
                let output = measure
                    .bytes
                    .checked_add(file.data.len())
                    .ok_or(Error::Limit)?;
                let records = measure.points.checked_add(4096).ok_or(Error::Limit)?;
                let map = file
                    .points
                    .capacity()
                    .checked_add(measure.points)
                    .and_then(|n| n.checked_add(1))
                    .and_then(|n| n.checked_mul(8 * (size_of::<([u64; 3], Entity)>() + 1)))
                    .ok_or(Error::Limit)?;
                let peak = source
                    .receipt()
                    .retained_result_bytes_upper_bound
                    .checked_add(options.additional_resident_bytes)
                    .and_then(|n| n.checked_add(file.data.capacity().checked_mul(2)?))
                    .and_then(|n| n.checked_add(output.checked_mul(3)?))
                    .and_then(|n| n.checked_add(map))
                    .and_then(|n| n.checked_add(records.checked_mul(size_of::<Record>())?))
                    .and_then(|n| n.checked_add(1024 * 1024))
                    .ok_or(Error::Limit)?;
                if peak > options.maximum_bytes {
                    return Err(Error::Limit);
                }
                measure
                    .work
                    .charge(output.div_ceil(size_of::<usize>()) + measure.points)?;
                file.data
                    .try_reserve_exact(output)
                    .map_err(|_| Error::Limit)?;
                file.points
                    .try_reserve(measure.points)
                    .map_err(|_| Error::Limit)?;
                let mut trace = Vec::new();
                trace.try_reserve_exact(records).map_err(|_| Error::Limit)?;
                let mut encoder = Encoder {
                    file,
                    work: measure.work,
                    trace,
                    output_limit: output,
                };
                let solid = encoder.body(source.shell().input())?;
                let result = Receipt {
                    work_used: options.maximum_work - encoder.work.remaining,
                    output_bytes_upper_bound: output,
                    phase_peak_bytes_reserved: peak,
                };
                receipt = Some((
                    result,
                    encoder.trace,
                    u32::try_from(solid.0).map_err(|_| Error::Limit)?,
                ));
                Ok(solid)
            })();
            result.map_err(AppendError)
        })
        .map_err(|e| e.0)?;
    receipt.ok_or(Error::Invalid)
}
struct Work<'a> {
    remaining: usize,
    interrupt: Option<&'a dyn Fn() -> bool>,
}
impl Work<'_> {
    fn charge(&mut self, n: usize) -> Result<()> {
        if self.interrupt.is_some_and(|f| f()) {
            return Err(Error::Cancelled);
        }
        self.remaining = self.remaining.checked_sub(n).ok_or(Error::Limit)?;
        Ok(())
    }
}
struct Measure<'a> {
    work: Work<'a>,
    bytes: usize,
    points: usize,
}
impl Measure<'_> {
    fn add(&mut self, n: usize) -> Result<()> {
        self.bytes = self.bytes.checked_add(n).ok_or(Error::Limit)?;
        self.work.charge(1)
    }
    fn scalar(&mut self, x: f64) -> Result<()> {
        let (w, b) = numbers::cost(x)?;
        self.work.charge(w)?;
        self.add(b)
    }
    fn tuple(&mut self, p: [f64; 3]) -> Result<()> {
        for x in p {
            self.scalar(x)?;
        }
        Ok(())
    }
    fn point(&mut self, p: [f64; 3]) -> Result<()> {
        self.points = self.points.checked_add(1).ok_or(Error::Limit)?;
        self.tuple(p)?;
        self.add(80)
    }
    fn axis(&mut self, a: Axis) -> Result<()> {
        self.point(a.origin)?;
        self.tuple(a.normal)?;
        self.tuple(a.reference)?;
        self.add(240)
    }
    fn canonical(&mut self, a: Axis) -> Result<()> {
        if a.origin[0] != 0.
            || a.origin[1] != 0.
            || (a.normal != [0., 0., 1.] && a.normal != [0., 0., -1.])
            || a.reference != [1., 0., 0.]
        {
            return Err(Error::Unsupported);
        }
        Ok(())
    }
    fn input(&mut self, s: &shell::Input) -> Result<()> {
        for &p in &s.vertices {
            self.point(p)?;
            self.add(80)?;
        }
        for curve in &s.curves {
            match curve {
                shell::Curve::Line(ends) => {
                    let [a, b] = ends.map(|i| s.vertices[i]);
                    let delta = std::array::from_fn(|k| b[k] - a[k]);
                    for k in 0..3 {
                        if !numbers::exact_difference(b[k], a[k], delta[k])? {
                            return Err(Error::Unsupported);
                        }
                    }
                    self.tuple(delta)?;
                    self.add(240)?;
                }
                shell::Curve::Original(radial::Curve::Spline { spans, poles }) => {
                    for &p in poles {
                        self.point(p)?;
                    }
                    self.add(
                        poles
                            .len()
                            .checked_mul(12)
                            .and_then(|n| n.checked_add((spans + 1) * 12 + 256))
                            .ok_or(Error::Limit)?,
                    )?;
                }
                shell::Curve::Original(radial::Curve::Circle { axis, radius_mm }) => {
                    self.canonical(*axis)?;
                    self.axis(*axis)?;
                    self.scalar(*radius_mm)?;
                    self.add(80)?;
                }
            }
        }
        for surface in &s.surfaces {
            match surface {
                Surface::Tensor(t) => {
                    for &p in t.poles() {
                        self.point(p)?;
                    }
                    let [a, b] = t.spans();
                    self.add(
                        t.poles()
                            .len()
                            .checked_mul(12)
                            .and_then(|n| {
                                n.checked_add((a + b + 2) * 12 + t.dimensions()[0] * 2 + 256)
                            })
                            .ok_or(Error::Limit)?,
                    )?;
                }
                Surface::Plane(a) => {
                    self.axis(*a)?;
                    self.add(80)?;
                }
                Surface::Cone {
                    axis,
                    radius_mm,
                    semi_angle_rad,
                } => {
                    self.canonical(*axis)?;
                    self.axis(*axis)?;
                    self.scalar(*radius_mm)?;
                    self.scalar(*semi_angle_rad)?;
                    self.add(80)?;
                }
            }
        }
        self.add(s.edges.len() * 128)?;
        for face in &s.faces {
            self.add(128)?;
            for ring in std::iter::once(&face.outer).chain(&face.inner) {
                self.add(ring.len() * 128 + 160)?;
            }
        }
        Ok(())
    }
}
struct Counter(usize);
impl fmt::Write for Counter {
    fn write_str(&mut self, s: &str) -> fmt::Result {
        self.0 = self.0.checked_add(s.len()).ok_or(fmt::Error)?;
        Ok(())
    }
}
struct Encoder<'a, 'c> {
    file: &'a mut StepFile,
    work: Work<'c>,
    trace: Vec<Record>,
    output_limit: usize,
}
impl Encoder<'_, '_> {
    fn entity(&mut self, args: fmt::Arguments<'_>) -> Result<Entity> {
        let mut n = Counter(0);
        writeln!(&mut n, "{}={args};", Entity(self.file.next)).map_err(|_| Error::Limit)?;
        self.work.charge(n.0 + 1)?;
        if self.file.data.len().checked_add(n.0).ok_or(Error::Limit)? > self.output_limit
            || self.trace.len() == self.trace.capacity()
        {
            return Err(Error::Limit);
        }
        let before = self.file.data.len();
        let id = self.file.entity(args);
        let start = before + id.to_string().len() + 1;
        let end = self.file.data.len() - 2;
        self.trace.push(Record {
            id: u32::try_from(id.0).map_err(|_| Error::Limit)?,
            body: start..end,
        });
        Ok(id)
    }
    fn scalar(&mut self, n: f64) -> Result<String> {
        self.work.charge(numbers::cost(n)?.0)?;
        numbers::literal(n)
    }
    fn tuple(&mut self, p: [f64; 3]) -> Result<String> {
        let [x, y, z] = p;
        let x = self.scalar(x)?;
        let y = self.scalar(y)?;
        let z = self.scalar(z)?;
        Ok(format!("({x},{y},{z})"))
    }
    fn point(&mut self, p: [f64; 3]) -> Result<Entity> {
        self.work.charge(1)?;
        let p = p.map(|x| if x == 0. { 0. } else { x });
        let key = p.map(f64::to_bits);
        if let Some(&id) = self.file.points.get(&key) {
            return Ok(id);
        }
        let tuple = self.tuple(p)?;
        let id = self.entity(format_args!("CARTESIAN_POINT('',{tuple})"))?;
        self.file.points.insert(key, id);
        Ok(id)
    }
    fn direction(&mut self, p: [f64; 3]) -> Result<Entity> {
        let tuple = self.tuple(p)?;
        self.entity(format_args!("DIRECTION('',{tuple})"))
    }
    fn axis(&mut self, a: Axis) -> Result<Entity> {
        let p = self.point(a.origin)?;
        let n = self.direction(a.normal)?;
        let r = self.direction(a.reference)?;
        self.entity(format_args!("AXIS2_PLACEMENT_3D('',{p},{n},{r})"))
    }
    fn knots(n: usize) -> (String, String) {
        let mult = (0..=n)
            .map(|i| if i == 0 || i == n { "10" } else { "9" })
            .collect::<Vec<_>>()
            .join(",");
        let knots = (0..=n).map(|i| i.to_string()).collect::<Vec<_>>().join(",");
        (mult, knots)
    }
    fn body(&mut self, s: &shell::Input) -> Result<Entity> {
        let mut vertices = Vec::new();
        for &p in &s.vertices {
            let point = self.point(p)?;
            vertices.push(self.entity(format_args!("VERTEX_POINT('',{point})"))?);
        }
        let mut curves = Vec::new();
        for curve in &s.curves {
            let curve = match curve {
                shell::Curve::Line(ends) => {
                    let [a, b] = ends.map(|i| s.vertices[i]);
                    let delta = std::array::from_fn(|k| b[k] - a[k]);
                    let point = self.point(a)?;
                    let d = self.direction(delta)?;
                    let vector = self.entity(format_args!("VECTOR('',{d},1.)"))?;
                    self.entity(format_args!("LINE('',{point},{vector})"))?
                }
                shell::Curve::Original(radial::Curve::Spline { spans, poles }) => {
                    let mut points = Vec::new();
                    for &p in poles {
                        points.push(self.point(p)?);
                    }
                    let (m, k) = Self::knots(*spans);
                    self.entity(format_args!("B_SPLINE_CURVE_WITH_KNOTS('',9,{},.UNSPECIFIED.,.F.,.F.,({m}),({k}),.UNSPECIFIED.)",Entities(&points)))?
                }
                shell::Curve::Original(radial::Curve::Circle { axis, radius_mm }) => {
                    let a = self.axis(*axis)?;
                    let r = self.scalar(*radius_mm)?;
                    self.entity(format_args!("CIRCLE('',{a},{r})"))?
                }
            };
            curves.push(curve);
        }
        let mut surfaces = Vec::new();
        for surface in &s.surfaces {
            let surface = match surface {
                Surface::Plane(a) => {
                    let a = self.axis(*a)?;
                    self.entity(format_args!("PLANE('',{a})"))?
                }
                Surface::Cone {
                    axis,
                    radius_mm,
                    semi_angle_rad,
                } => {
                    let a = self.axis(*axis)?;
                    let r = self.scalar(*radius_mm)?;
                    if *semi_angle_rad == 0. {
                        self.entity(format_args!("CYLINDRICAL_SURFACE('',{a},{r})"))?
                    } else {
                        let angle = self.scalar(*semi_angle_rad)?;
                        self.entity(format_args!("CONICAL_SURFACE('',{a},{r},{angle})"))?
                    }
                }
                Surface::Tensor(t) => {
                    let mut points = Vec::new();
                    for &p in t.poles() {
                        points.push(self.point(p)?);
                    }
                    let rows = points
                        .chunks_exact(t.dimensions()[1])
                        .map(|row| Entities(row).to_string())
                        .collect::<Vec<_>>()
                        .join(",");
                    let [nu, nv] = t.spans();
                    let (um, uk) = Self::knots(nu);
                    let (vm, vk) = Self::knots(nv);
                    self.entity(format_args!("B_SPLINE_SURFACE_WITH_KNOTS('',9,9,({rows}),.UNSPECIFIED.,.F.,.F.,.F.,({um}),({vm}),({uk}),({vk}),.UNSPECIFIED.)"))?
                }
            };
            surfaces.push(surface);
        }
        let mut edges = Vec::new();
        for e in &s.edges {
            let [a, b] = e.vertices.map(|v| vertices[v]);
            let curve = curves[e.curve];
            let fwd = if e.forward { ".T." } else { ".F." };
            edges.push(self.entity(format_args!("EDGE_CURVE('',{a},{b},{curve},{fwd})"))?);
        }
        let mut faces = Vec::new();
        for face in &s.faces {
            let mut bounds = Vec::new();
            for (i, ring) in std::iter::once(&face.outer).chain(&face.inner).enumerate() {
                let mut wire = Vec::new();
                for &(edge, direction) in ring {
                    let edge = edges[edge];
                    let sense = if direction { ".T." } else { ".F." };
                    wire.push(self.entity(format_args!("ORIENTED_EDGE('',*,*,{edge},{sense})"))?);
                }
                let wire = self.entity(format_args!("EDGE_LOOP('',{})", Entities(&wire)))?;
                let name = if i == 0 {
                    "FACE_OUTER_BOUND"
                } else {
                    "FACE_BOUND"
                };
                bounds.push(self.entity(format_args!("{name}('',{wire},.T.)"))?);
            }
            let surface = surfaces[face.surface];
            let sense = if face.forward { ".T." } else { ".F." };
            faces.push(self.entity(format_args!(
                "ADVANCED_FACE('',{},{surface},{sense})",
                Entities(&bounds)
            ))?);
        }
        let shell = self.entity(format_args!("CLOSED_SHELL('',{})", Entities(&faces)))?;
        self.entity(format_args!("MANIFOLD_SOLID_BREP('',{shell})"))
    }
}
