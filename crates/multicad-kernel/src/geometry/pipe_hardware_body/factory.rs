//! Caller-owned bounded reuse of completed identity-frame thread preparation.
//! Stock, placement and manufacturing certification are never cached here.
use super::*;
use std::{mem::size_of, rc::Weak};
const MAX_BYTES: usize = 32 * 1024 * 1024;
const MAX_ENTRIES: usize = 16;
// An expired Weak keeps the inline Rc allocation until its candidate drops.
const WEAK_CELL_BYTES: usize =
    std::mem::size_of::<ThreadSource>() + 2 * std::mem::size_of::<usize>();
#[derive(Clone, Copy, PartialEq, Eq)]
struct Key {
    owner: EntityId,
    values: [u64; 15],
}
impl Key {
    fn new(owner: EntityId, spec: AuthoredSpec, policy: GeometryPolicy) -> Self {
        // Exhaustive destructuring forces new source/policy fields to be audited.
        let AuthoredSpec {
            major_radius_mm,
            pitch_mm,
            length_mm,
            radial_taper,
            section,
            crest_radius_limit_mm,
        } = spec;
        let GeometryPolicy {
            linear_tolerance_mm,
            boolean_fuzzy_mm,
            minimum_feature_mm,
            maximum_coordinate_mm,
        } = policy;
        let (tag, depth, flat) = match section {
            multicad_rust_kernel::axial_pipe_thread::authored::AuthoredSection::NptFlat {
                depth_mm,
                flat_width_mm,
            } => (0, depth_mm, flat_width_mm),
            multicad_rust_kernel::axial_pipe_thread::authored::AuthoredSection::BsptRounded {
                depth_mm,
            } => (1, depth_mm, 0.),
        };
        Self {
            owner,
            values: [
                major_radius_mm.to_bits(),
                pitch_mm.to_bits(),
                length_mm.to_bits(),
                radial_taper.to_bits(),
                tag,
                depth.to_bits(),
                flat.to_bits(),
                u64::from(crest_radius_limit_mm.is_some()),
                crest_radius_limit_mm.unwrap_or(0.).to_bits(),
                linear_tolerance_mm.to_bits(),
                boolean_fuzzy_mm.to_bits(),
                minimum_feature_mm.to_bits(),
                maximum_coordinate_mm.to_bits(),
                source::MEASURE_ACCURACY.to_bits(),
                1,
            ],
        }
    }
}
struct Entry {
    key: Key,
    source: Rc<ThreadSource>,
    bytes: usize,
    used: u64,
}
struct Pending {
    key: Key,
    source: Weak<ThreadSource>,
    used: u64,
}
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct PipeThreadFactoryStats {
    pub requests: usize,
    pub hits: usize,
    pub prepared: usize,
    pub entries: usize,
    pub pending: usize,
    /// Conservatively accounted actual retained source + fixed metadata bound.
    pub retained_bytes: usize,
    pub evicted: usize,
}
/// Owned by a compiler/evaluation cache. There is no global/thread-local factory.
/// Pending weak candidates retain no extra source allocation. Only a successful
/// authoritative compile publishes them; dropping/aborting preserves old entries.
pub struct PipeThreadFactory {
    entries: Vec<Entry>,
    pending: Vec<Pending>,
    maximum_bytes: usize,
    maximum_entries: usize,
    clock: u64,
    active: bool,
    stats: PipeThreadFactoryStats,
}
impl PipeThreadFactory {
    pub fn new(maximum_bytes: usize, maximum_entries: usize) -> Self {
        let maximum_bytes = maximum_bytes.min(MAX_BYTES);
        let maximum_entries = maximum_entries.min(MAX_ENTRIES).min(
            maximum_bytes / (size_of::<Entry>() + size_of::<Pending>() + WEAK_CELL_BYTES).max(1),
        );
        let entries = Vec::new();
        let pending = Vec::new();
        Self {
            entries,
            pending,
            maximum_bytes,
            maximum_entries,
            clock: 0,
            active: false,
            stats: Default::default(),
        }
    }
    fn metadata_budget(&self) -> usize {
        self.maximum_entries * (size_of::<Entry>() + size_of::<Pending>() + WEAK_CELL_BYTES)
    }
    fn metadata_bytes(&self) -> usize {
        self.entries.capacity() * size_of::<Entry>()
            + self.pending.capacity() * (size_of::<Pending>() + WEAK_CELL_BYTES)
    }
    pub fn retained_bytes(&self) -> usize {
        self.metadata_bytes() + self.entries.iter().map(|e| e.bytes).sum::<usize>()
    }
    pub fn stats(&self) -> PipeThreadFactoryStats {
        PipeThreadFactoryStats {
            entries: self.entries.len(),
            pending: self.pending.len(),
            retained_bytes: self.retained_bytes(),
            ..self.stats
        }
    }
    pub fn clear(&mut self) {
        self.entries = Vec::new();
        self.pending = Vec::new();
        self.active = false;
        self.stats = Default::default();
    }
    pub fn begin_compile(&mut self) {
        self.pending.clear();
        self.active = true;
        self.stats = Default::default();
    }
    pub fn finish_compile(&mut self, success: bool) {
        self.active = false;
        if !success {
            self.pending.clear();
            return;
        }
        while let Some(candidate) = self.pending.pop() {
            let Some(source) = candidate.source.upgrade() else {
                continue;
            };
            let bytes = source.retained_bytes();
            if self.maximum_entries == 0
                || bytes.saturating_add(self.metadata_budget()) > self.maximum_bytes
            {
                continue;
            }
            if self.entries.iter().any(|entry| entry.key == candidate.key) {
                continue;
            }
            while self.entries.len() >= self.maximum_entries
                || self.retained_bytes().saturating_add(bytes) > self.maximum_bytes
            {
                let Some((index, _)) = self
                    .entries
                    .iter()
                    .enumerate()
                    .min_by_key(|(_, entry)| entry.used)
                else {
                    break;
                };
                self.entries.swap_remove(index);
                self.stats.evicted += 1;
            }
            if self.entries.capacity() == 0 {
                self.entries = Vec::with_capacity(self.maximum_entries);
            }
            self.entries.push(Entry {
                key: candidate.key,
                source,
                bytes,
                used: candidate.used,
            });
        }
    }
    #[allow(clippy::too_many_arguments)]
    pub fn build(
        &mut self,
        spec: AuthoredSpec,
        kind: PipeHardwareKind,
        minimum_wall_mm: f64,
        feature: EntityId,
        frame: PlaneFrame,
        policy: GeometryPolicy,
        interrupt: Option<&dyn Fn() -> bool>,
    ) -> Result<Solid, GeometryError> {
        checkpoint(interrupt)?;
        let key = Key::new(feature, spec, policy);
        self.clock = self.clock.wrapping_add(1);
        self.stats.requests += 1;
        let prepared = if let Some(entry) = self.entries.iter().find(|e| e.key == key) {
            Some(Rc::clone(&entry.source))
        } else {
            self.pending
                .iter()
                .find(|e| e.key == key)
                .and_then(|e| e.source.upgrade())
        };
        let hit = prepared.is_some();
        // Complete fresh intent, clearance and placement validation surrounds
        // reuse. Only the unchanged identity-frame source preparation is shared.
        let (solid, source) = construct(
            spec,
            kind,
            minimum_wall_mm,
            feature,
            frame,
            policy,
            interrupt,
            prepared,
        )?;
        checkpoint(interrupt)?;
        if hit {
            if let Some(entry) = self.entries.iter_mut().find(|e| e.key == key) {
                entry.used = self.clock;
            }
            self.stats.hits += 1;
        } else {
            self.stats.prepared += 1;
        }
        if self.active
            && !hit
            && self.pending.len() < self.maximum_entries
            && source
                .retained_bytes()
                .saturating_add(self.metadata_budget())
                <= self.maximum_bytes
        {
            if self.pending.capacity() == 0 {
                self.pending = Vec::with_capacity(self.maximum_entries);
            }
            self.pending.push(Pending {
                key,
                source: Rc::downgrade(&source),
                used: self.clock,
            });
        }
        Ok(solid)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn parameters(bspt: bool) -> AuthoredSpec {
        use multicad_core::{
            document::Document,
            model::{
                NptGaugeConvention, NptThreadForm, NptThreadIntent, NptThreadSize,
                NptTruncationTarget, PipeGaugeConvention, PipeThreadForm, PipeThreadIntent,
                PipeThreadSize, ThreadMode,
            },
            parameters,
        };
        let doc = Document::new("Prepared actual pipe").unwrap();
        let eval = parameters::evaluate(&doc).unwrap();
        let t = if bspt {
            PipeThreadIntent {
                form: PipeThreadForm::Rc,
                size: PipeThreadSize::OneEighth,
                gauge: PipeGaugeConvention::Iso7Nominal1994,
            }
            .parameters()
        } else {
            NptThreadIntent {
                form: NptThreadForm::External,
                size: NptThreadSize::OneEighth,
                gauge: NptGaugeConvention::Nbs1959PlainFace,
                truncation: NptTruncationTarget::MidpointConeParallel,
            }
            .parameters()
        };
        let r = t
            .resolve(
                &eval,
                if bspt {
                    ThreadMode::CutInternal
                } else {
                    ThreadMode::NewExternal
                },
            )
            .unwrap();
        AuthoredSpec {
            major_radius_mm: r.major_radius_mm,
            pitch_mm: r.pitch_mm,
            length_mm: r.length_mm,
            radial_taper: r.profile.radial_taper,
            section: if bspt {
                multicad_rust_kernel::axial_pipe_thread::authored::AuthoredSection::BsptRounded {
                    depth_mm: r.profile.depth_mm,
                }
            } else {
                let n = r.npt.as_ref().unwrap();
                multicad_rust_kernel::axial_pipe_thread::authored::AuthoredSection::NptFlat {
                    depth_mm: n.depth_mm,
                    flat_width_mm: n.flat_width_mm,
                }
            },
            crest_radius_limit_mm: r.npt.as_ref().and_then(|n| n.stock_crest_radius_limit_mm),
        }
    }
    fn kind(bspt: bool, af: f64) -> PipeHardwareKind {
        if bspt {
            PipeHardwareKind::HexNut {
                across_flats_mm: af,
            }
        } else {
            PipeHardwareKind::HexHead {
                across_flats_mm: af,
                height_mm: 3.,
            }
        }
    }
    fn build(
        factory: &mut PipeThreadFactory,
        spec: AuthoredSpec,
        owner: EntityId,
        af: f64,
        frame: PlaneFrame,
    ) -> Solid {
        factory
            .build(
                spec,
                kind(spec.radial_taper < 0., af),
                1.,
                owner,
                frame,
                GeometryPolicy::default(),
                None,
            )
            .unwrap()
    }
    #[test]
    fn key_covers_every_spec_and_policy_field_and_variant() {
        let spec = parameters(false);
        let owner = EntityId::new();
        let policy = GeometryPolicy::default();
        let original = Key::new(owner, spec, policy);
        let mut variants = Vec::new();
        let mut value = spec;
        value.major_radius_mm += 0.01;
        variants.push(value);
        value = spec;
        value.pitch_mm += 0.01;
        variants.push(value);
        value = spec;
        value.length_mm += 0.01;
        variants.push(value);
        value = spec;
        value.radial_taper += 0.01;
        variants.push(value);
        value = spec;
        value.crest_radius_limit_mm = None;
        variants.push(value);
        value = spec;
        value.crest_radius_limit_mm = Some(20.);
        variants.push(value);
        value = spec;
        if let multicad_rust_kernel::axial_pipe_thread::authored::AuthoredSection::NptFlat {
            depth_mm,
            flat_width_mm,
        } = spec.section
        {
            value.section =
                multicad_rust_kernel::axial_pipe_thread::authored::AuthoredSection::NptFlat {
                    depth_mm: depth_mm + 0.01,
                    flat_width_mm,
                };
            variants.push(value);
            value.section =
                multicad_rust_kernel::axial_pipe_thread::authored::AuthoredSection::NptFlat {
                    depth_mm,
                    flat_width_mm: flat_width_mm + 0.01,
                };
            variants.push(value);
            value.section =
                multicad_rust_kernel::axial_pipe_thread::authored::AuthoredSection::BsptRounded {
                    depth_mm,
                };
            variants.push(value);
        }
        for changed in variants {
            assert!(Key::new(owner, changed, policy) != original);
        }
        assert!(Key::new(EntityId::new(), spec, policy) != original);
        for field in 0..4 {
            let mut changed = policy;
            match field {
                0 => changed.linear_tolerance_mm *= 2.,
                1 => changed.boolean_fuzzy_mm = 0.01,
                2 => changed.minimum_feature_mm *= 2.,
                _ => changed.maximum_coordinate_mm *= 2.,
            }
            assert!(Key::new(owner, spec, changed) != original);
        }
    }
    #[test]
    fn actual_bspt_af_and_frame_edits_reuse_original_source_and_manufacture_fresh() {
        let spec = parameters(true);
        let owner = EntityId::new();
        let mut factory = PipeThreadFactory::new(MAX_BYTES, 16);
        factory.begin_compile();
        let first = build(&mut factory, spec, owner, 18., PlaneFrame::default());
        assert_eq!(factory.stats().entries, 0);
        assert_eq!(factory.stats().pending, 1);
        factory.finish_compile(true);
        assert_eq!(factory.stats().entries, 1);
        assert!(factory.retained_bytes() < MAX_BYTES);
        eprintln!(
            "PREPARED_BSPT bytes={} limit={MAX_BYTES}",
            factory.retained_bytes()
        );
        let frame = PlaneFrame {
            origin: [11., -7., 5.],
            u: [0., 1., 0.],
            v: [0., 0., 1.],
        };
        factory.begin_compile();
        let hit_at = std::time::Instant::now();
        let edited = build(&mut factory, spec, owner, 20., frame);
        eprintln!(
            "PREPARED_BSPT_AF_FRAME_HIT ms={:.6}",
            hit_at.elapsed().as_secs_f64() * 1000.
        );
        factory.finish_compile(true);
        assert_eq!(factory.stats().hits, 1);
        assert_eq!(factory.stats().prepared, 0);
        assert!(Rc::ptr_eq(
            &first.pipe_hardware_body().unwrap().source._prepared,
            &edited.pipe_hardware_body().unwrap().source._prepared
        ));
        let direct = Solid::rust_authored_axial_pipe_hardware_with_interrupt(
            spec,
            kind(true, 20.),
            1.,
            owner,
            frame,
            GeometryPolicy::default(),
            None,
        )
        .unwrap();
        assert_eq!(edited.bounds().unwrap(), direct.bounds().unwrap());
        assert_eq!(edited.volume_mm3().unwrap(), direct.volume_mm3().unwrap());
        assert_eq!(edited.faces().unwrap(), direct.faces().unwrap());
        assert_eq!(
            edited.topology_snapshot(owner, 3).unwrap(),
            direct.topology_snapshot(owner, 3).unwrap()
        );
        let q = MeshOptions {
            linear_deflection_mm: 0.1,
            angular_deflection_rad: 0.3,
            maximum_triangles: 500_000,
        };
        let a = edited.manufacturing_mesh_with_receipt(q).unwrap();
        let b = edited.manufacturing_mesh_with_receipt(q).unwrap();
        let c = direct.manufacturing_mesh_with_receipt(q).unwrap();
        assert_eq!(a.mesh, b.mesh);
        assert_eq!(a.mesh, c.mesh);
        assert_eq!(a.edge_pairs, c.edge_pairs);
        assert_eq!(a.vertex_links, c.vertex_links);
        assert_eq!(edited.rust_pipe_hardware_mesh_build_count(), Some(2));
        factory.clear();
        assert_eq!(factory.retained_bytes(), 0);
        assert_eq!(edited.volume_mm3().unwrap(), direct.volume_mm3().unwrap());
    }
    #[test]
    fn actual_owner_spec_policy_misses_and_success_eviction_are_bounded() {
        let spec = parameters(false);
        let owner = EntityId::new();
        let mut factory = PipeThreadFactory::new(MAX_BYTES, 1);
        factory.begin_compile();
        let first = build(&mut factory, spec, owner, 18., PlaneFrame::default());
        factory.finish_compile(true);
        for (other, changed, policy) in [
            (EntityId::new(), spec, GeometryPolicy::default()),
            (
                owner,
                AuthoredSpec {
                    major_radius_mm: spec.major_radius_mm + 0.001,
                    ..spec
                },
                GeometryPolicy::default(),
            ),
            (
                owner,
                spec,
                GeometryPolicy {
                    maximum_coordinate_mm: 1e6,
                    ..GeometryPolicy::default()
                },
            ),
        ] {
            factory.begin_compile();
            let solid = factory
                .build(
                    changed,
                    kind(false, 18.),
                    1.,
                    other,
                    PlaneFrame::default(),
                    policy,
                    None,
                )
                .unwrap();
            assert_eq!(factory.stats().hits, 0);
            factory.finish_compile(true);
            assert_eq!(factory.stats().entries, 1);
            assert_eq!(factory.stats().evicted, 1);
            assert!(factory.retained_bytes() <= MAX_BYTES);
            assert!(solid.volume_mm3().unwrap() > 0.);
        }
        let mut disabled = PipeThreadFactory::new(1, 1);
        disabled.begin_compile();
        let solid = build(&mut disabled, spec, owner, 18., PlaneFrame::default());
        disabled.finish_compile(true);
        assert_eq!(disabled.stats().entries, 0);
        assert_eq!(disabled.retained_bytes(), 0);
        assert_eq!(solid.volume_mm3().unwrap(), first.volume_mm3().unwrap());
    }
    #[test]
    fn actual_late_cancel_error_and_abort_never_publish_or_remove_prior_source() {
        use std::cell::Cell;
        let spec = parameters(false);
        let owner = EntityId::new();
        let mut factory = PipeThreadFactory::new(MAX_BYTES, 2);
        factory.begin_compile();
        let first = build(&mut factory, spec, owner, 18., PlaneFrame::default());
        factory.finish_compile(true);
        let before = factory.entries[0].source.clone();
        factory.begin_compile();
        assert!(
            factory
                .build(
                    spec,
                    kind(false, 1.),
                    1.,
                    owner,
                    PlaneFrame::default(),
                    GeometryPolicy::default(),
                    None
                )
                .is_err()
        );
        factory.finish_compile(false);
        assert!(Rc::ptr_eq(&before, &factory.entries[0].source));
        let fresh = EntityId::new();
        let ticks = Cell::new(0);
        let observe = || {
            ticks.set(ticks.get() + 1);
            false
        };
        factory.begin_compile();
        let candidate = factory
            .build(
                spec,
                kind(false, 18.),
                1.,
                fresh,
                PlaneFrame::default(),
                GeometryPolicy::default(),
                Some(&observe),
            )
            .unwrap();
        let total = ticks.get();
        assert_eq!(factory.stats().pending, 1);
        factory.finish_compile(false);
        assert_eq!(factory.stats().entries, 1);
        assert_eq!(factory.stats().pending, 0);
        ticks.set(0);
        let cancel = || {
            ticks.set(ticks.get() + 1);
            ticks.get() == total
        };
        factory.begin_compile();
        assert!(matches!(
            factory.build(
                spec,
                kind(false, 18.),
                1.,
                fresh,
                PlaneFrame::default(),
                GeometryPolicy::default(),
                Some(&cancel)
            ),
            Err(GeometryError {
                kind: GeometryErrorKind::Cancelled,
                ..
            })
        ));
        factory.finish_compile(false);
        assert_eq!(factory.stats().entries, 1);
        assert_eq!(factory.stats().pending, 0);
        assert!(Rc::ptr_eq(&before, &factory.entries[0].source));
        factory.begin_compile();
        let retry = build(&mut factory, spec, fresh, 18., PlaneFrame::default());
        factory.finish_compile(true);
        assert_eq!(retry.volume_mm3().unwrap(), candidate.volume_mm3().unwrap());
        assert_eq!(factory.stats().entries, 2);
        assert!(first.volume_mm3().unwrap() > 0.);
    }
}
