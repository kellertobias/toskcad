//! Owning ACTUAL file geometry; no authored thread or head recipe is inferred.
use super::*;
use crate::step_read::literal_member_binding::RetainedRecords;
use crate::tapered_thread::interval::I;
use std::sync::Arc;
#[cfg(test)]
mod tests;
#[derive(Clone, Copy, Debug)]
pub struct OwnedHeadOptions {
    pub resources: Options,
    pub material_tolerance_mm: f64,
    pub minimum_wall_mm: f64,
}
impl Default for OwnedHeadOptions {
    fn default() -> Self {
        Self {
            resources: Options::default(),
            material_tolerance_mm: 1e-7,
            minimum_wall_mm: 1e-6,
        }
    }
}
#[derive(Clone, Copy, Debug)]
pub struct OwnedHeadReceipt {
    pub material: Receipt,
    pub retained_auxiliary_work_with_material_prefix: usize,
    pub fresh_retention_work: usize,
    pub aggregate_work: usize,
    pub accounted_peak_bytes: usize,
    pub retained_bytes: usize,
}
struct Row {
    id: u32,
    first: usize,
    count: usize,
}
struct ActualRecords {
    input: Arc<str>,
    rows: Vec<Row>,
    components: Vec<Range<usize>>,
}
impl ActualRecords {
    fn capture(
        input: Arc<str>,
        doc: &CompoundDocument<'_>,
        b: &mut Budget<'_>,
    ) -> Result<(Self, usize)> {
        let mut count = 0usize;
        for values in doc.records().values() {
            b.step(1)?;
            count = count.checked_add(values.len()).ok_or(Error::Limit)?;
        }
        let mut rows = b.allocate(doc.records().len())?;
        let mut components = b.allocate(count)?;
        for (&id, values) in doc.records() {
            b.step(32)?;
            let first = components.len();
            for e in values {
                let source = e.source.expect("actual neutral component");
                b.step(source.len() + 16)?;
                let start = (source.as_ptr() as usize)
                    .checked_sub(input.as_ptr() as usize)
                    .ok_or(Error::Invalid)?;
                let end = start.checked_add(source.len()).ok_or(Error::Limit)?;
                if input.get(start..end) != Some(source) {
                    return Err(Error::Invalid);
                }
                components.push(start..end);
            }
            rows.push(Row {
                id,
                first,
                count: values.len(),
            });
        }
        let heap =
            rows.capacity() * size_of::<Row>() + components.capacity() * size_of::<Range<usize>>();
        Ok((
            Self {
                input,
                rows,
                components,
            },
            heap,
        ))
    }
    fn heap_bytes(&self) -> usize {
        self.rows.capacity() * size_of::<Row>()
            + self.components.capacity() * size_of::<Range<usize>>()
    }
}
impl RetainedRecords for ActualRecords {
    fn record_count(&self) -> usize {
        self.rows.len()
    }
    fn record_index(&self, id: u32) -> Option<usize> {
        self.rows.binary_search_by_key(&id, |r| r.id).ok()
    }
    fn component(&self, record: usize, index: usize) -> Option<&str> {
        let r = self.rows.get(record)?;
        if index >= r.count {
            return None;
        }
        self.input
            .get(self.components.get(r.first + index)?.clone())
    }
}
/// SAME declared monotone axial transfer used by the occupied head proof.
/// Raw coefficients stay unchanged. Geometry consumers apply this map to the
/// auxiliary thread and remove only its derived upper verification cap.
#[derive(Clone, Copy, Debug)]
pub struct AxialSewingTransfer {
    lower: [f64; 2],
    upper: [f64; 2],
    target_upper: f64,
    factor: [f64; 2],
    displacement: f64,
}
impl AxialSewingTransfer {
    // Callers charge the fixed interval arithmetic before entering. Raw bounds
    // are the earned complete retained control hull, not a cap-plane proxy.
    pub(super) fn from_raw_domain(
        lower: [f64; 2],
        upper: [f64; 2],
        target_upper: f64,
        raw_bounds: [[f64; 2]; 3],
    ) -> Result<Self> {
        let lo = I::new(lower[0], lower[1]);
        let hi = I::new(upper[0], upper[1]);
        let factor = (I::point(target_upper) - lo) / (hi - lo);
        let displacement =
            (factor - I::point(1.)) * (I::new(raw_bounds[2][0], raw_bounds[2][1]) - lo);
        if !factor.valid() || factor.lo <= 0. || !displacement.valid() {
            return Err(Error::Unsupported);
        }
        Ok(Self {
            lower,
            upper,
            target_upper,
            factor: [factor.lo, factor.hi],
            displacement: displacement.magnitude(),
        })
    }
    pub fn lower_raw_height_mm(self) -> [f64; 2] {
        self.lower
    }
    pub fn upper_raw_height_mm(self) -> [f64; 2] {
        self.upper
    }
    pub fn target_upper_height_mm(self) -> f64 {
        self.target_upper
    }
    pub fn positive_axial_scale_bounds(self) -> [f64; 2] {
        self.factor
    }
    /// Bound on the owned raw thread image, not on arbitrary exterior points.
    pub fn domain_displacement_mm(self) -> f64 {
        self.displacement
    }
    pub fn map_raw_box_mm(self, input: [[f64; 2]; 3]) -> Result<[[f64; 2]; 3]> {
        let mut p = input.map(|b| I::new(b[0], b[1]));
        if p.iter().any(|v| !v.valid()) {
            return Err(Error::Invalid);
        }
        let lower = I::new(self.lower[0], self.lower[1]);
        p[2] = lower + (p[2] - lower) * I::new(self.factor[0], self.factor[1]);
        if p.iter().any(|v| !v.valid()) {
            return Err(Error::Unsupported);
        }
        Ok(p.map(|v| [v.lo, v.hi]))
    }
    pub fn map_raw_vector_bounds(self, input: [[f64; 2]; 3]) -> Result<[[f64; 2]; 3]> {
        let mut v = input.map(|b| I::new(b[0], b[1]));
        if v.iter().any(|v| !v.valid()) {
            return Err(Error::Invalid);
        }
        v[2] = v[2] * I::new(self.factor[0], self.factor[1]);
        if v.iter().any(|v| !v.valid()) {
            return Err(Error::Unsupported);
        }
        Ok(v.map(|v| [v.lo, v.hi]))
    }
    /// Exact positive diagonal-map cofactor, before any normal normalization.
    pub fn map_raw_normal_bounds(self, input: [[f64; 2]; 3]) -> Result<[[f64; 2]; 3]> {
        let mut v = input.map(|b| I::new(b[0], b[1]));
        if v.iter().any(|v| !v.valid()) {
            return Err(Error::Invalid);
        }
        let scale = I::new(self.factor[0], self.factor[1]);
        v[0] = v[0] * scale;
        v[1] = v[1] * scale;
        if v.iter().any(|v| !v.valid()) {
            return Err(Error::Unsupported);
        }
        Ok(v.map(|v| [v.lo, v.hi]))
    }
}
pub struct OwnedActualHead {
    records: ActualRecords,
    auxiliary: OwnedMaterial,
    head: HeadInput,
    chart: CapChart,
    roles: PhysicalRoles,
    polygon: [[f64; 2]; 6],
    stations: [f64; 2],
    planar_sewing: f64,
    join_sewing: f64,
    minimum_wall: f64,
    local_bounds: [[f64; 2]; 3],
    sewing: AxialSewingTransfer,
    volume: [f64; 2],
    receipt: OwnedHeadReceipt,
}
impl OwnedActualHead {
    /// Caller supplies the immutable actual bytes. Capturing its Arc does not
    /// copy the input allocation; the single shared input is counted once.
    pub fn read(
        input: Arc<str>,
        options: OwnedHeadOptions,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Self> {
        let r = options.resources;
        if !options.material_tolerance_mm.is_finite()
            || options.material_tolerance_mm <= 0.
            || options.material_tolerance_mm > 1e-4
            || !options.minimum_wall_mm.is_finite()
            || options.minimum_wall_mm <= 0.
        {
            return Err(Error::Invalid);
        }
        // Input Arc reference counts/alignment and the complete output header.
        // All variable buffers are separately charged by their actual capacity.
        let header = size_of::<Self>() + 2 * size_of::<usize>() + std::mem::align_of::<usize>();
        let additional = r
            .additional_resident_bytes
            .checked_add(header)
            .ok_or(Error::Limit)?;
        let (packet, records) = DecodedFile::parse_with_capture(
            &input,
            Options {
                additional_resident_bytes: additional,
                ..r
            },
            interrupt,
            |doc, b| {
                b.step(1)?;
                ActualRecords::capture(Arc::clone(&input), doc, b)
            },
        )?;
        let doc = CompoundDocument::parse(
            packet.derived_auxiliary_text(),
            DocumentOptions {
                maximum_work: packet.derived_parser_work_reserved,
                retained_budget_bytes: r
                    .maximum_bytes
                    .checked_sub(packet.retained_beyond_auxiliary_input)
                    .filter(|n| *n > 0)
                    .ok_or(Error::Limit)?,
            },
            interrupt,
        )
        .map_err(neutral_error)?;
        let parse_work = r
            .maximum_work
            .checked_sub(packet.receipt().aggregate_work)
            .filter(|n| *n > 0)
            .ok_or(Error::Limit)?;
        let parse_bytes = r
            .maximum_bytes
            .checked_sub(packet.retained_beyond_auxiliary_input)
            .and_then(|n| n.checked_sub(doc.retained_budget_bytes()))
            .filter(|n| *n > 0)
            .ok_or(Error::Limit)?;
        let selected = Selected::parse(
            &doc,
            doc.represented_members()[0].item_entity(),
            super::super::Options {
                maximum_work: parse_work.min(16_000_000),
                maximum_bytes: parse_bytes,
            },
            interrupt,
        )?;
        let proof = packet.certify(
            selected.parsed(),
            options.material_tolerance_mm,
            options.minimum_wall_mm,
            interrupt,
        )?;
        let material = proof.material();
        let material_receipt = proof.receipt();
        let old_head_bytes = head::storage(material.head())?;
        let outside = doc
            .retained_budget_bytes()
            .checked_add(packet.retained_beyond_auxiliary_input)
            .and_then(|n| n.checked_add(old_head_bytes))
            .ok_or(Error::Limit)?;
        let child_bytes = r
            .maximum_bytes
            .checked_sub(outside)
            .filter(|n| *n > 0)
            .ok_or(Error::Limit)?;
        let inherited = material.capped_thread().work_used();
        let child_work = r
            .maximum_work
            .checked_sub(material_receipt.aggregate_work)
            .and_then(|n| n.checked_add(inherited))
            .filter(|n| *n > 0)
            .ok_or(Error::Limit)?
            .min(64_000_000);
        let auxiliary = selected.retain_occupied(
            material.capped_thread(),
            child_work,
            child_bytes,
            interrupt,
        )?;
        let retention = auxiliary
            .work_used()
            .checked_sub(inherited)
            .ok_or(Error::Limit)?;
        let mut aggregate = material_receipt
            .aggregate_work
            .checked_add(retention)
            .ok_or(Error::Limit)?;
        let peak = material_receipt.accounted_peak_bytes.max(
            outside
                .checked_add(auxiliary.construction_query_peak_bytes())
                .ok_or(Error::Limit)?,
        );
        if peak > r.maximum_bytes {
            return Err(Error::Limit);
        }
        let mut work = Budget {
            used: 0,
            remaining: r.maximum_work.checked_sub(aggregate).ok_or(Error::Limit)?,
            bytes: peak,
            peak,
            limit: r.maximum_bytes,
            interrupt,
        };
        // Clone admission includes all live selected arrays, old head data and
        // the actual raw/derived buffers BEFORE allocating the new head arrays.
        let clone_bytes = outside
            .checked_add(selected.parsed().query_budget_bytes())
            .and_then(|n| n.checked_add(auxiliary.retained_budget_bytes()))
            .and_then(|n| n.checked_add(old_head_bytes))
            .ok_or(Error::Limit)?;
        if clone_bytes > r.maximum_bytes {
            return Err(Error::Limit);
        }
        work.step(512)?;
        let head = material.head().clone();
        let polygon = *material.head_polygon_mm();
        let stations = material.head_stations_mm();
        let planar_sewing = material.planar_sewing_mm();
        let join_sewing = material.join_sewing_mm();
        let minimum_wall = material.minimum_wall_mm();
        let mut cap_heights = [[0.; 2]; 2];
        let mut cap_count = 0;
        work.step(32)?;
        for face in auxiliary.faces() {
            if let Support::Plane { height } = face.support() {
                if cap_count == 2 {
                    return Err(Error::Unsupported);
                }
                cap_heights[cap_count] =
                    auxiliary.point(height).ok_or(Error::Invalid)?.coordinates()[2].bounds();
                cap_count += 1;
            }
        }
        if cap_count != 2 {
            return Err(Error::Unsupported);
        }
        cap_heights.sort_by(|a, b| a[0].total_cmp(&b[0]));
        let lower = I::new(cap_heights[0][0], cap_heights[0][1]);
        let upper = I::new(cap_heights[1][0], cap_heights[1][1]);
        let factor = (I::point(stations[0]) - lower) / (upper - lower);
        if !factor.valid() || factor.lo <= 0. {
            return Err(Error::Unsupported);
        }
        let mut twice = I::point(0.);
        for i in 0..6 {
            work.step(64)?;
            let a = polygon[i].map(I::point);
            let b = polygon[(i + 1) % 6].map(I::point);
            twice = twice + a[0] * b[1] - a[1] * b[0];
        }
        let area = twice / I::point(2.);
        let prism = area * (I::point(stations[1]) - I::point(stations[0]));
        let thread = auxiliary.local_volume_mm3();
        let volume = I::new(thread[0], thread[1]) * factor + prism;
        if !area.valid() || area.lo <= 0. || !volume.valid() || volume.lo <= 0. {
            return Err(Error::Unsupported);
        }
        // Complete actual control hull, including any cap-plane collars.
        work.step(64)?;
        let raw_bounds = auxiliary.local_enclosure_mm();
        let sewing = AxialSewingTransfer::from_raw_domain(
            cap_heights[0],
            cap_heights[1],
            stations[0],
            raw_bounds,
        )?;
        work.step(256)?;
        let mut bounds = sewing.map_raw_box_mm(raw_bounds)?;
        for vertex in &head.vertices {
            work.step(12)?;
            for k in 0..3 {
                bounds[k][0] = bounds[k][0].min(vertex[k]);
                bounds[k][1] = bounds[k][1].max(vertex[k]);
            }
        }
        for range in &mut bounds {
            work.step(8)?;
            *range = [
                (I::point(range[0]) - I::point(join_sewing)).lo,
                (I::point(range[1]) + I::point(join_sewing)).hi,
            ];
        }
        aggregate = aggregate.checked_add(work.used).ok_or(Error::Limit)?;
        let retained = header
            .checked_add(input.len())
            .and_then(|n| n.checked_add(records.heap_bytes()))
            .and_then(|n| n.checked_add(auxiliary.retained_budget_bytes()))
            .and_then(|n| {
                n.checked_add(
                    head::storage(&head)
                        .ok()?
                        .checked_sub(size_of::<HeadInput>())?,
                )
            })
            .ok_or(Error::Limit)?;
        let peak = peak.max(clone_bytes).max(
            retained
                .checked_add(r.additional_resident_bytes)
                .ok_or(Error::Limit)?,
        );
        if aggregate > r.maximum_work || peak > r.maximum_bytes {
            return Err(Error::Limit);
        }
        work.step(0)?;
        let chart = packet.chart;
        let roles = packet.roles;
        Ok(Self {
            records,
            auxiliary,
            head,
            chart,
            roles,
            polygon,
            stations,
            planar_sewing,
            join_sewing,
            minimum_wall,
            local_bounds: bounds,
            sewing,
            volume: [volume.lo, volume.hi],
            receipt: OwnedHeadReceipt {
                material: material_receipt,
                retained_auxiliary_work_with_material_prefix: inherited + retention,
                fresh_retention_work: retention,
                aggregate_work: aggregate,
                accounted_peak_bytes: peak,
                retained_bytes: retained,
            },
        })
    }
    pub fn actual_file(&self) -> &str {
        &self.records.input
    }
    pub fn original_component(&self, id: u32, component: usize) -> Option<&str> {
        self.records
            .component(self.records.record_index(id)?, component)
    }
    pub fn original_entity_at(&self, index: usize) -> Option<u32> {
        self.records.rows.get(index).map(|r| r.id)
    }
    pub fn original_record_count(&self) -> usize {
        self.records.record_count()
    }
    /// Independently earned auxiliary ONLY. Its upper verification cap is not
    /// a physical whole-head face and must be removed by the joining stage.
    pub fn auxiliary_thread(&self) -> &OwnedMaterial {
        &self.auxiliary
    }
    pub fn head(&self) -> &HeadInput {
        &self.head
    }
    pub fn cap_chart(&self) -> CapChart {
        self.chart
    }
    pub fn physical_roles(&self) -> PhysicalRoles {
        self.roles
    }
    pub fn head_polygon_mm(&self) -> &[[f64; 2]; 6] {
        &self.polygon
    }
    pub fn head_stations_mm(&self) -> [f64; 2] {
        self.stations
    }
    pub fn planar_sewing_mm(&self) -> f64 {
        self.planar_sewing
    }
    pub fn join_sewing_mm(&self) -> f64 {
        self.join_sewing
    }
    pub fn minimum_wall_mm(&self) -> f64 {
        self.minimum_wall
    }
    pub fn axial_sewing_transfer(&self) -> AxialSewingTransfer {
        self.sewing
    }
    pub fn local_sewn_enclosure_mm(&self) -> [[f64; 2]; 3] {
        self.local_bounds
    }
    /// Volume of the explicitly declared closed sewn reference, including the
    /// actual axial transfer and actual polygon prism, never an authored mass.
    pub fn local_sewn_volume_mm3(&self) -> [f64; 2] {
        self.volume
    }
    pub fn receipt(&self) -> OwnedHeadReceipt {
        self.receipt
    }
}
