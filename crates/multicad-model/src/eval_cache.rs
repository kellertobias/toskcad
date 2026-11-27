//! Cross-compile feature result cache (CAD-150).
//!
//! A full compile still walks the complete canonical document in dependency
//! order, but a Feature whose result key is unchanged since an earlier
//! compile reuses that compile's exact solid (shared, immutable B-rep handle)
//! and its display/topology snapshot instead of re-running the native
//! provider, the topology export and the tessellation.
//!
//! # Result key
//!
//! The key of a feature is a SHA-256 over:
//! - a cache format tag, solid backend, geometry policy and display options;
//! - the feature's own identity;
//! - the set of every entity in the feature's *input closure*, each as its
//!   identity plus a hash of its complete canonical source (name, links and
//!   definition) and, for Parameters, the evaluated value.
//!
//! The input closure starts at the feature and follows (a) its evaluation
//! dependencies (`Document::evaluation_dependencies`, the same graph that
//! orders the compile) and (b) every entity identity mentioned anywhere in
//! its source, including expression bindings and identities embedded in
//! strings. Both kinds of edge are followed transitively. A Part reached only
//! through a mention (for example the feature's own output Part) contributes
//! its own source plus the identity, type and suppression of its listed
//! features, not their closures: geometry can only flow along evaluation
//! dependencies, because nothing else orders a producer before its consumer.
//!
//! The evaluation revision is not part of the key. It only stamps topology
//! snapshots, so a reused snapshot is re-stamped with the current revision.
//!
//! # Purity
//!
//! Only provider types whose result is a function of that closure are
//! cached ([`cacheable_type`]). Weldment members find their connections by
//! scanning the document; those connections and the connected partner
//! members are declared as hidden inputs and join the closure. Providers that
//! otherwise look entities up by scanning the document (walls/openings,
//! levels), resolve
//! assembly poses (transfer, placement-based scale) or are unknown to this
//! module are always evaluated, and so is every feature whose closure
//! contains such a feature, an Instance, an Assembly, or a sketch targeted by
//! an occurrence reference. A new provider is therefore uncached until it is
//! audited and added here. Original Rust joining-plate outputs additionally use
//! the shared bounded occurrence-pose projection: absolute driver/structure
//! fingerprints replace only their explicit pose links, while every geometry,
//! parameter and reference supplier remains in the ordinary audited closure.
//! Generic Instance/Assembly consumers remain uncacheable.
//!
//! Fence geometry reads its declared run/type/path/level/parameter closure.
//! A miss validates the complete run before returning any member solid, even
//! for source-scoped compilation. The compiler checks complete member identity
//! and material assignments before every hit; those checks cannot be skipped
//! by a cached shape. Failed partial run generation returns no cacheable member.
//!
//! Native symbol-template validation additionally retains up to 64 completed
//! receipts (also bounded by `max_entries`). Only audited self-contained source
//! providers qualify. Receipts cover the entire pin and compiler configuration,
//! and are published only after the enclosing compile completes successfully.
//!
//! # Assembly and occurrence composition (CAD-396)
//!
//! Reusable Assembly compounds and placed Instance occurrences are not
//! feature results, but a full compile still rebuilds their exact solids,
//! topology, display payloads and canonical leaves from the member solids.
//! When every input of such a node is the *same* immutable B-rep as in an
//! earlier successful full compile ([`Solid::shares_geometry_with`], pointer
//! identity, never geometric equality) and every other input recorded in its
//! key (member list, component, producer, resolved pose bits, leaf metadata
//! and compiler configuration) is equal, the node's earlier result is the
//! same function of the same inputs and is reused, re-stamped with the
//! current revision. Entries retain their input solids, so an input address
//! cannot be reused by different geometry while its entry exists. Poses are
//! still resolved and refined every compile; source-material checks and
//! display budgets are still applied to reused nodes. Entries are staged and
//! promoted only after the enclosing full compile succeeds.
//!
//! # Lifetime, memory and threads
//!
//! Entries are only stored for successfully evaluated features, so a
//! superseded or failed compile can leave valid entries behind but never a
//! partial or stale one: a key describes all inputs of its result. The cache
//! holds kernel solids and is therefore `!Send`/`!Sync`; it lives on the
//! thread that compiles (the model runtime worker). After every compile the
//! entries used by that compile are kept (their geometry is shared with the
//! returned model, bounded by [`crate::MAX_MODEL_FEATURES`]); older entries
//! are evicted least-recently-used first until their estimated size is at
//! most `retained_bytes` and the entry count at most `max_entries`.
use crate::{Compiler, GeometryBackend, ShapeSnapshot, SheetLayout};
use multicad_core::{
    document::{Document, EntityData, EntityId},
    model,
    parameters::EvaluatedParameters,
    project_file::ContentHash,
    sketch_external,
};
use multicad_kernel::{PlaneFrame, Solid};
use std::{
    collections::{BTreeMap, BTreeSet},
    sync::Arc,
};

mod profile_connections;

/// Default bound on the estimated size of entries not used by the latest
/// compile (kept for undo/redo and toggled edits).
pub const DEFAULT_RETAINED_BYTES: usize = 256 * 1024 * 1024;
/// Default bound on the total number of entries.
pub const DEFAULT_MAX_ENTRIES: usize = 4 * crate::MAX_MODEL_FEATURES;
const FORMAT: &[u8] = b"multicad.evaluation-cache/1";

/// Feature provider types whose evaluation reads only the input closure
/// (plus the [`hidden_inputs`] declared for its type).
pub fn cacheable_type(type_id: &str) -> bool {
    [
        model::EXTRUDE_TYPE,
        model::PRIMITIVE_TYPE,
        model::THREAD_TYPE,
        model::GEAR_TYPE,
        model::MECHANICAL_PART_TYPE,
        model::REVOLVE_TYPE,
        model::BOOLEAN_TYPE,
        model::SPLIT_TYPE,
        model::TRANSFORM_TYPE,
        model::EDGE_TREATMENT_TYPE,
        model::HELICAL_TUBE_CUT_TYPE,
        model::ADVANCED_TYPE,
        model::ALIGN_TYPE,
        model::PANEL_TYPE,
        // Every miss proves the complete run before returning a member. Typed
        // run/type/path/level/parameter links cover geometry; the compiler also
        // checks the full member/material set before accepting a cache hit.
        model::FENCE_MEMBER_TYPE,
        model::PANEL_JOINT_CUT_TYPE,
        model::LOOSE_TENON_TYPE,
        model::MACHINING_CUT_TYPE,
        model::SHEET_TYPE,
        model::SHEET_PROFILE_CUT_TYPE,
        model::MEMBER_TYPE,
        model::HARDWARE_TYPE,
        model::PROFILE_CONNECTION_OUTPUT_TYPE,
        model::AID_TYPE,
        model::imports::STEP_BODY_TYPE,
        model::imports::BOARD_TYPE,
        model::imports::COMPONENT_TYPE,
    ]
    .contains(&type_id)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct CacheKey([u8; 32]);

/// Side results a provider records next to its solid.
#[derive(Clone, Debug, Default)]
pub(crate) struct Extras {
    pub panel_frame: Option<PlaneFrame>,
    pub sheet_layout: Option<SheetLayout>,
    pub member_path: Option<model::PathSegment>,
    pub consumed: Vec<(EntityId, EntityId)>,
}

struct Entry {
    solid: Solid,
    shape: Arc<ShapeSnapshot>,
    extras: Extras,
    bytes: usize,
    used: u64,
}

/// What one compile did with the cache.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct CompileStats {
    /// Features whose geometry provider ran (cache misses and uncacheable).
    pub evaluated: BTreeSet<EntityId>,
    /// Features reused from an earlier compile.
    pub reused: BTreeSet<EntityId>,
    /// The identical full Rust model was reused without scheduling features.
    pub reused_full_model: bool,
    /// Evaluated features that are not cacheable (see the module docs).
    pub uncacheable: BTreeSet<EntityId>,
    /// Same-owner audited literal symbol nodes restored this compile.
    pub symbol_reused: BTreeSet<EntityId>,
    /// Audited literal symbol nodes constructed this compile.
    pub symbol_evaluated: BTreeSet<EntityId>,
    /// Assembly/Instance compositions restored from identical inputs (CAD-396).
    pub composition_reused: BTreeSet<EntityId>,
    /// Assembly/Instance compositions constructed this compile.
    pub composition_evaluated: BTreeSet<EntityId>,
    /// Entries evicted after this compile.
    pub evicted: usize,
    /// Entries and their estimated size after eviction.
    pub entries: usize,
    pub estimated_bytes: usize,
}

pub(crate) struct SymbolNode {
    pub solid: Solid,
    pub shape: Option<Arc<ShapeSnapshot>>,
    pub occurrence: Option<Arc<crate::OccurrenceSnapshot>>,
    pub leaves: crate::assemblies::CachedLeaves,
}
impl SymbolNode {
    fn shared(&self, revision: u64) -> Option<Self> {
        let shape = self.shape.as_ref().map(|shape| {
            let mut copy = (**shape).clone();
            copy.topology.revision = revision;
            Arc::new(copy)
        });
        let occurrence = self.occurrence.as_ref().map(|occurrence| {
            let mut copy = (**occurrence).clone();
            copy.topology.revision = revision;
            Arc::new(copy)
        });
        Some(Self {
            solid: self.solid.share().ok()?,
            shape,
            occurrence,
            leaves: self.leaves.shared(revision)?,
        })
    }
    fn bytes(&self) -> usize {
        let primary = self.shape.as_ref().map_or_else(
            || {
                self.occurrence
                    .as_ref()
                    .map_or(0, |o| o.topology.elements.len().saturating_mul(1280))
            },
            |s| estimate(s),
        );
        primary.saturating_add(self.leaves.estimated_bytes())
    }
}
struct SymbolEntry {
    node: SymbolNode,
    bytes: usize,
    used: u64,
}
/// A composed Assembly/Instance node and the exact solids it was built from.
struct CompositionEntry {
    inputs: Vec<Solid>,
    node: SymbolNode,
    bytes: usize,
    used: u64,
}
impl CompositionEntry {
    fn matches(&self, inputs: &[Solid]) -> bool {
        self.inputs.len() == inputs.len()
            && self
                .inputs
                .iter()
                .zip(inputs)
                .all(|(a, b)| a.shares_geometry_with(b))
    }
}
/// Bounds on retained compositions (alternatives per key, entries). Their
/// estimated size is bounded separately by the owner's `retained_bytes`.
const COMPOSITION_ALTERNATIVES: usize = 4;
const COMPOSITION_ENTRIES: usize = 4096;
#[derive(Clone)]
struct LiteralSymbol {
    entities: Vec<EntityId>,
    feature: EntityId,
}
fn literal_symbols(document: &Document) -> BTreeMap<EntityId, LiteralSymbol> {
    use multicad_core::symbols;
    let mut output = BTreeMap::new();
    // The ordinary compiler validates ownership before constructing its cache
    // session. Keep this additional eligibility check fail-closed as well.
    if symbols::validate(document).is_err() {
        return output;
    }
    for entity in document.entities() {
        let EntityData::Extension { definition, .. } = &entity.data else {
            continue;
        };
        if definition.type_id != symbols::INSTANCE_TYPE {
            continue;
        }
        let Ok(instance) = symbols::read_instance(document, entity.id) else {
            continue;
        };
        let Ok(definition) = symbols::read_definition(document, instance.definition) else {
            continue;
        };
        let mut feature = None;
        let mut parts = 0;
        let mut eligible = true;
        for item in definition.template.entities() {
            match &item.data {
                EntityData::Parameter { .. } | EntityData::Material { .. } => {}
                EntityData::Part { .. } => {
                    parts += 1;
                    if item.id != definition.root {
                        eligible = false;
                    }
                }
                EntityData::Feature {
                    definition: source,
                    suppressed: false,
                    ..
                } if source.type_id == model::MECHANICAL_PART_TYPE => {
                    if feature.replace(item.id).is_some() {
                        eligible = false;
                    }
                    if !model::read_mechanical_part(&definition.template, item.id).is_ok_and(
                        |source| {
                            matches!(
                                source.parameters,
                                model::MechanicalPartParameters::MachineScrew { .. }
                            )
                        },
                    ) {
                        eligible = false;
                    }
                }
                _ => eligible = false,
            }
        }
        let Some(feature) = feature.and_then(|id| instance.ids.get(&id).copied()) else {
            continue;
        };
        let Some(part) = instance.ids.get(&definition.root).copied() else {
            continue;
        };
        if !eligible
            || parts != 1
            || !matches!(document.get(part).map(|e|&e.data),Some(EntityData::Part {features,..}) if features==&[feature])
        {
            continue;
        }
        if !matches!(document.get(instance.instance).map(|e|&e.data),Some(EntityData::Instance {component,placement_source:None,..}) if *component==instance.wrapper)
            || !matches!(document.get(instance.member).map(|e|&e.data),Some(EntityData::Instance {component,placement_source:None,..}) if *component==part)
            || !matches!(document.get(instance.wrapper).map(|e|&e.data),Some(EntityData::Assembly {instances}) if instances==&[instance.member])
        {
            continue;
        }
        let mut entities = instance.ids.values().copied().collect::<Vec<_>>();
        entities.extend([
            entity.id,
            instance.definition,
            instance.instance,
            instance.wrapper,
            instance.wrapper_metadata,
            instance.member,
        ]);
        // Conservatively fingerprint canonical material definitions/assignments
        // as well as each private Part/material entity in the ownership map.
        entities.extend(document.entities().filter(|e|matches!(&e.data,EntityData::Material {..})||matches!(&e.data,EntityData::Extension {definition,..} if definition.type_id.contains("material"))).map(|e|e.id));
        entities.sort();
        entities.dedup();
        let graph = LiteralSymbol { entities, feature };
        for id in [instance.member, instance.wrapper, instance.instance] {
            output.insert(id, graph.clone());
        }
    }
    output
}

struct FullRustModel {
    compiler: Compiler,
    model: crate::Model,
    features: BTreeSet<EntityId>,
    /// Audited literal symbol nodes contained in the retained model.
    symbols: BTreeSet<EntityId>,
    bytes: usize,
}

/// Exact results reused across compiles on one model worker.
pub struct EvaluationCache {
    entries: BTreeMap<CacheKey, Entry>,
    clock: u64,
    retained_bytes: usize,
    max_entries: usize,
    last: CompileStats,
    pub(crate) gear_pairs: crate::PairCycleProofStore,
    pub(crate) native_intersections: crate::native_intersection_proofs::ProofStore,
    symbol_proofs: BTreeSet<ContentHash>,
    pending_symbol_proofs: BTreeSet<ContentHash>,
    symbol_proof_hits: usize,
    symbol_native_checks: usize,
    symbol_nodes: BTreeMap<CacheKey, SymbolEntry>,
    pending_symbol_nodes: BTreeMap<CacheKey, SymbolEntry>,
    compositions: BTreeMap<CacheKey, Vec<CompositionEntry>>,
    composition_bytes: usize,
    pending_compositions: Vec<(CacheKey, CompositionEntry)>,
    full_rust: Option<FullRustModel>,
    pub(crate) pipe_threads: multicad_kernel::PipeThreadFactory,
}
impl Default for EvaluationCache {
    fn default() -> Self {
        Self::new(DEFAULT_RETAINED_BYTES, DEFAULT_MAX_ENTRIES)
    }
}
impl EvaluationCache {
    pub fn new(retained_bytes: usize, max_entries: usize) -> Self {
        Self {
            entries: BTreeMap::new(),
            clock: 0,
            retained_bytes,
            max_entries,
            last: CompileStats::default(),
            gear_pairs: Default::default(),
            native_intersections: crate::native_intersection_proofs::ProofStore::new(
                retained_bytes,
                max_entries,
            ),
            symbol_proofs: BTreeSet::new(),
            pending_symbol_proofs: BTreeSet::new(),
            symbol_proof_hits: 0,
            symbol_native_checks: 0,
            symbol_nodes: BTreeMap::new(),
            pending_symbol_nodes: BTreeMap::new(),
            compositions: BTreeMap::new(),
            composition_bytes: 0,
            pending_compositions: Vec::new(),
            full_rust: None,
            pipe_threads: multicad_kernel::PipeThreadFactory::new(retained_bytes, max_entries),
        }
    }
    pub fn pipe_thread_factory_stats(&self) -> multicad_kernel::PipeThreadFactoryStats {
        self.pipe_threads.stats()
    }
    /// Statistics of the most recent cached compile (also when it failed or
    /// was superseded).
    pub fn last_compile(&self) -> &CompileStats {
        &self.last
    }
    pub(crate) fn reset_compile_stats(&mut self) {
        self.last = CompileStats {
            entries: self.entries.len(),
            estimated_bytes: self.last.estimated_bytes,
            ..CompileStats::default()
        };
    }
    pub fn len(&self) -> usize {
        self.entries.len()
    }
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
    pub fn estimated_bytes(&self) -> usize {
        self.entries
            .values()
            .map(|e| e.bytes)
            .sum::<usize>()
            .saturating_add(self.full_rust.as_ref().map_or(0, |full| full.bytes))
            .saturating_add(self.pipe_threads.retained_bytes())
    }
    pub fn clear(&mut self) {
        self.pipe_threads.clear();
        self.entries.clear();
        self.full_rust = None;
        self.gear_pairs.clear();
        self.native_intersections.clear();
        self.symbol_proofs.clear();
        self.pending_symbol_proofs.clear();
        self.symbol_nodes.clear();
        self.pending_symbol_nodes.clear();
        self.compositions.clear();
        self.composition_bytes = 0;
        self.pending_compositions.clear();
        self.last = CompileStats::default();
    }

    pub(crate) fn reuse_full_rust(
        &mut self,
        compiler: &Compiler,
        document: &Document,
        revision: u64,
    ) -> Option<crate::Model> {
        let full = self.full_rust.as_ref()?;
        let stamp = full.model.snapshot().stamp;
        if full.compiler != *compiler
            || stamp.document != document.id()
            || stamp.revision != revision
            || !ContentHash::document(document).is_ok_and(|hash| hash == stamp.source_hash)
        {
            self.full_rust = None;
            return None;
        }
        let model = full.model.share_completed();
        self.last.reused.clone_from(&full.features);
        self.last.symbol_reused.clone_from(&full.symbols);
        self.last.reused_full_model = true;
        Some(model)
    }

    pub(crate) fn store_full_rust(&mut self, compiler: &Compiler, model: &crate::Model) {
        if self.max_entries == 0 {
            return;
        }
        let snapshot = model.snapshot();
        // Conservatively charge shared display data against both caches.
        let shapes = snapshot.shapes.values().map(|shape| estimate(shape));
        let occurrences = snapshot
            .occurrences
            .values()
            .map(|occurrence| occurrence.topology.elements.len().saturating_mul(1280));
        let leaves = snapshot
            .leaf_occurrences
            .values()
            .map(|leaf| leaf.topology.elements.len().saturating_mul(1280));
        let bytes = shapes
            .chain(occurrences)
            .chain(leaves)
            .fold(0, usize::saturating_add);
        if bytes > self.retained_bytes.saturating_sub(self.estimated_bytes()) {
            return;
        }
        self.full_rust = Some(FullRustModel {
            compiler: *compiler,
            model: model.share_completed(),
            features: self
                .last
                .evaluated
                .union(&self.last.reused)
                .copied()
                .collect(),
            symbols: self
                .last
                .symbol_evaluated
                .union(&self.last.symbol_reused)
                .copied()
                .collect(),
            bytes,
        });
        self.last.estimated_bytes = self.estimated_bytes();
    }
    /// Same-owner, exact native Intersect(None) receipts. These never replace
    /// current source, phase, engagement, leaf enumeration or budget checks.
    pub fn native_intersection_proof_stats(
        &self,
    ) -> crate::native_intersection_proofs::NativeIntersectionProofStats {
        self.native_intersections.stats()
    }
    /// Completed, bounded, backend-specific template receipts retained by this owner.
    pub fn symbol_template_proofs(&self) -> usize {
        self.symbol_proofs.len()
    }
    /// Cumulative successful template-proof hits (independent of feature stats).
    pub fn symbol_template_proof_hits(&self) -> usize {
        self.symbol_proof_hits
    }
    /// Cumulative template geometry checks, including failed checks.
    pub fn symbol_template_native_checks(&self) -> usize {
        self.symbol_native_checks
    }
    /// Successfully committed same-owner literal symbol node results.
    pub fn symbol_node_entries(&self) -> usize {
        self.symbol_nodes.len()
    }
    /// Retained Assembly/Instance compositions (CAD-396).
    pub fn composition_entries(&self) -> usize {
        self.compositions.values().map(Vec::len).sum()
    }
    /// Estimated size of the retained compositions.
    pub fn composition_bytes(&self) -> usize {
        self.composition_bytes
    }
    pub(crate) fn begin_symbol_proofs(&mut self) {
        self.pending_compositions.clear();
        self.pending_symbol_nodes.clear();
        self.pending_symbol_proofs.clear();
    }
    pub(crate) fn finish_symbol_proofs(&mut self, success: bool) {
        let pending = std::mem::take(&mut self.pending_symbol_proofs);
        let nodes = std::mem::take(&mut self.pending_symbol_nodes);
        let compositions = std::mem::take(&mut self.pending_compositions);
        if success {
            self.promote_compositions(compositions);
            let limit = self.max_entries.min(768);
            let budget = self.retained_bytes.min(64 * 1024 * 1024);
            for (key, entry) in nodes {
                if limit == 0 || entry.bytes > budget {
                    continue;
                }
                self.symbol_nodes.remove(&key);
                while self.symbol_nodes.len() >= limit
                    || self.symbol_nodes.values().map(|e| e.bytes).sum::<usize>()
                        > budget - entry.bytes
                {
                    let Some(oldest) = self
                        .symbol_nodes
                        .iter()
                        .min_by_key(|(_, e)| e.used)
                        .map(|(key, _)| *key)
                    else {
                        break;
                    };
                    self.symbol_nodes.remove(&oldest);
                }
                self.symbol_nodes.insert(key, entry);
            }
            let limit = self.max_entries.min(64);
            for key in pending {
                if limit == 0 {
                    break;
                }
                if self.symbol_proofs.len() >= limit {
                    self.symbol_proofs.pop_first();
                }
                self.symbol_proofs.insert(key);
            }
        }
    }
    fn promote_compositions(&mut self, staged: Vec<(CacheKey, CompositionEntry)>) {
        let limit = self.max_entries.min(COMPOSITION_ENTRIES);
        let budget = self.retained_bytes;
        for (key, entry) in staged {
            if limit == 0 || entry.bytes > budget {
                continue;
            }
            let alternatives = self.compositions.entry(key).or_default();
            let mut freed = 0;
            alternatives.retain(|old| {
                let keep = !old.matches(&entry.inputs);
                if !keep {
                    freed += old.bytes;
                }
                keep
            });
            if alternatives.len() >= COMPOSITION_ALTERNATIVES
                && let Some(oldest) = (0..alternatives.len()).min_by_key(|&i| alternatives[i].used)
            {
                freed += alternatives.remove(oldest).bytes;
            }
            self.composition_bytes -= freed;
            while self.composition_entries() >= limit
                || self.composition_bytes > budget - entry.bytes
            {
                let Some((oldest, index)) = self
                    .compositions
                    .iter()
                    .flat_map(|(key, list)| {
                        list.iter()
                            .enumerate()
                            .map(move |(index, e)| (e.used, *key, index))
                    })
                    .min()
                    .map(|(_, key, index)| (key, index))
                else {
                    break;
                };
                let list = self.compositions.get_mut(&oldest).expect("listed key");
                self.composition_bytes -= list.remove(index).bytes;
                if list.is_empty() {
                    self.compositions.remove(&oldest);
                }
            }
            self.composition_bytes += entry.bytes;
            self.compositions.entry(key).or_default().push(entry);
        }
    }
    pub(crate) fn symbol_proof_hit(&mut self, key: ContentHash) -> bool {
        let hit = self.symbol_proofs.contains(&key);
        self.symbol_proof_hits += usize::from(hit);
        hit
    }
    pub(crate) fn symbol_native_check(&mut self) {
        self.symbol_native_checks += 1;
    }
    pub(crate) fn stage_symbol_proof(&mut self, key: ContentHash) {
        // Current compile also has a 64-definition source budget.
        self.pending_symbol_proofs.insert(key);
    }
    fn evict(&mut self) -> usize {
        let current = self.clock;
        let mut stale: Vec<_> = self
            .entries
            .iter()
            .filter(|(_, e)| e.used < current)
            .map(|(k, e)| (e.used, *k, e.bytes))
            .collect();
        stale.sort();
        let mut stale_bytes: usize = stale
            .iter()
            .map(|s| s.2)
            .sum::<usize>()
            .saturating_add(self.pipe_threads.retained_bytes());
        let mut evicted = 0;
        for (_, key, bytes) in stale {
            if stale_bytes <= self.retained_bytes && self.entries.len() <= self.max_entries {
                break;
            }
            self.entries.remove(&key);
            stale_bytes -= bytes;
            evicted += 1;
        }
        evicted
    }
}

fn estimate(shape: &ShapeSnapshot) -> usize {
    let mesh = shape.mesh.vertices.len() * 8
        + shape.mesh.indices.len() * 4
        + shape.mesh.face_indices.len() * 4;
    let edges: usize = shape
        .edges
        .iter()
        .map(|e| 32 + e.points_mm.len() * 24)
        .sum();
    // Topology records and the native B-rep are estimated per element.
    let topology = shape.topology.elements.len() * (256 + 1024);
    std::mem::size_of::<ShapeSnapshot>() + mesh + edges + shape.vertices.len() * 32 + topology
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum Visit {
    Full,
    /// A Part reached only through a mention.
    PartSource,
}

/// Evaluation dependencies and other mentioned entities of one entity.
type Edges = (Vec<EntityId>, Vec<EntityId>);

/// One compile's view of the cache.
pub(crate) struct Session<'a> {
    cache: &'a mut EvaluationCache,
    backend: GeometryBackend,
    profile_poses: BTreeMap<EntityId, Option<Arc<profile_connections::PoseInputs>>>,
    document: &'a Document,
    parameters: &'a EvaluatedParameters,
    salt: Vec<u8>,
    hashes: BTreeMap<(EntityId, Visit), Option<[u8; 32]>>,
    edges: BTreeMap<EntityId, Option<Edges>>,
    sources: BTreeMap<EntityId, Option<std::rc::Rc<Source>>>,
    volatile_sketches: BTreeSet<EntityId>,
    stats: CompileStats,
    literal_symbols: BTreeMap<EntityId, LiteralSymbol>,
}
impl<'a> Session<'a> {
    pub(crate) fn new(
        cache: &'a mut EvaluationCache,
        compiler: &Compiler,
        backend: GeometryBackend,
        document: &'a Document,
        parameters: &'a EvaluatedParameters,
    ) -> Self {
        cache.clock += 1;
        let mut salt = FORMAT.to_vec();
        salt.extend(
            format!(
                "|{backend:?}|{:?}|{:?}|{}|",
                compiler.policy, compiler.mesh, compiler.omit_display
            )
            .bytes(),
        );
        let volatile_sketches = document
            .entities()
            .filter(|e| sketch_external::is_occurrence_reference(document, e.id))
            .filter_map(|e| sketch_external::read_occurrence(document, e.id).ok())
            .map(|request| request.query.sketch)
            .collect();
        Self {
            cache,
            backend,
            profile_poses: BTreeMap::new(),
            document,
            parameters,
            salt,
            hashes: BTreeMap::new(),
            edges: BTreeMap::new(),
            sources: BTreeMap::new(),
            volatile_sketches,
            stats: CompileStats::default(),
            literal_symbols: literal_symbols(document),
        }
    }
    pub(crate) fn pipe_thread_factory(&mut self) -> &mut multicad_kernel::PipeThreadFactory {
        &mut self.cache.pipe_threads
    }
    fn profile_pose(&mut self, id: EntityId) -> Option<Arc<profile_connections::PoseInputs>> {
        if let Some(found) = self.profile_poses.get(&id) {
            return found.clone();
        }
        let found = profile_connections::project(self.document, id).map(Arc::new);
        self.profile_poses.insert(id, found.clone());
        found
    }
    /// The canonical serialization of an entity and the identities it
    /// mentions, computed once per compile for both its hash and its edges.
    fn source(&mut self, id: EntityId) -> Option<std::rc::Rc<Source>> {
        if let Some(found) = self.sources.get(&id) {
            return found.clone();
        }
        let computed = self
            .document
            .get(id)
            .and_then(Source::of)
            .map(std::rc::Rc::new);
        self.sources.insert(id, computed.clone());
        computed
    }
    /// Evaluation-dependency and mention edges of an entity.
    fn edges(&mut self, id: EntityId) -> Option<Edges> {
        if let Some(found) = self.edges.get(&id) {
            return found.clone();
        }
        let computed = (|| {
            let projection = if profile_connections::source(self.document, id) {
                Some(self.profile_pose(id)?)
            } else {
                None
            };
            self.document.get(id)?;
            let mut evaluation: BTreeSet<_> = self.document.evaluation_dependencies(id).ok()?;
            let mut mentioned = self.source(id)?.mentioned.clone();
            mentioned.remove(&id);
            if let Some(projection) = projection {
                for id in &projection.omitted {
                    evaluation.remove(id);
                    mentioned.remove(id);
                }
                evaluation.extend(projection.geometry.iter().copied());
            }
            let mentioned = mentioned
                .into_iter()
                .filter(|m| !evaluation.contains(m) && self.document.get(*m).is_some())
                .collect();
            Some((evaluation.into_iter().collect(), mentioned))
        })();
        self.edges.insert(id, computed.clone());
        computed
    }
    fn hash(&mut self, id: EntityId, visit: Visit) -> Option<[u8; 32]> {
        if let Some(found) = self.hashes.get(&(id, visit)) {
            return *found;
        }
        let computed = (|| {
            let projection = if profile_connections::source(self.document, id) {
                Some(self.profile_pose(id)?)
            } else {
                None
            };
            let entity = self.document.get(id)?;
            let mut bytes = self.source(id)?.bytes.clone();
            if let Some(projection) = projection {
                bytes.extend(b"|absolute-profile-poses|");
                bytes.extend(projection.fingerprint.0);
            }
            match (&entity.data, visit) {
                (EntityData::Parameter { .. }, _) => {
                    bytes.extend(b"|value|");
                    bytes.extend(serde_json::to_vec(self.parameters.values.get(&id)?).ok()?);
                }
                (EntityData::Part { features, .. }, Visit::PartSource) => {
                    for feature in features {
                        let Some(EntityData::Feature {
                            definition,
                            suppressed,
                            ..
                        }) = self.document.get(*feature).map(|e| &e.data)
                        else {
                            return None;
                        };
                        bytes.extend(
                            format!(
                                "|{feature}|{}|{}|{suppressed}",
                                definition.type_id, definition.schema_version
                            )
                            .bytes(),
                        );
                    }
                }
                _ => {}
            }
            Some(ContentHash::of(&bytes).0)
        })();
        self.hashes.insert((id, visit), computed);
        computed
    }
    /// The result key of a feature, or None when it must be evaluated.
    pub(crate) fn key(&mut self, feature: EntityId) -> Option<CacheKey> {
        let mut visited: BTreeMap<(EntityId, Visit), [u8; 32]> = BTreeMap::new();
        let mut pending = vec![(feature, Visit::Full)];
        while let Some((id, visit)) = pending.pop() {
            if visited.contains_key(&(id, visit)) {
                continue;
            }
            let entity = self.document.get(id)?;
            visited.insert((id, visit), self.hash(id, visit)?);
            if visit == Visit::PartSource {
                continue;
            }
            match &entity.data {
                EntityData::Feature { definition, .. }
                    if definition.type_id == model::PROFILE_CONNECTION_OUTPUT_TYPE
                        && (self.backend != GeometryBackend::Rust
                            || !profile_connections::output(self.document, id)) =>
                {
                    return None;
                }
                EntityData::Instance { .. } | EntityData::Assembly { .. } => return None,
                EntityData::Feature { definition, .. } if !cacheable_type(&definition.type_id) => {
                    return None;
                }
                EntityData::Sketch { .. } if self.volatile_sketches.contains(&id) => return None,
                _ => {}
            }
            if let EntityData::Feature { definition, .. } = &entity.data {
                pending.extend(
                    hidden_inputs(self.document, id, &definition.type_id)?
                        .into_iter()
                        .map(|d| (d, Visit::Full)),
                );
            }
            let (evaluation, mentioned) = self.edges(id)?;
            pending.extend(evaluation.into_iter().map(|d| (d, Visit::Full)));
            for m in mentioned {
                let part = matches!(
                    self.document.get(m).map(|e| &e.data),
                    Some(EntityData::Part { .. })
                );
                pending.push((m, if part { Visit::PartSource } else { Visit::Full }));
            }
        }
        // The closure is a set; two features can share one (for example
        // connected weldment members), so the key also names its root.
        let mut bytes = self.salt.clone();
        bytes.extend(format!("root|{feature}|").bytes());
        for ((id, visit), hash) in &visited {
            bytes.extend(id.to_string().bytes());
            bytes.push(*visit as u8);
            bytes.extend(hash);
        }
        Some(CacheKey(ContentHash::of(&bytes).0))
    }
    pub(crate) fn symbol_key(&mut self, id: EntityId) -> Option<CacheKey> {
        let graph = self.literal_symbols.get(&id)?.clone();
        // Keep the ordinary audited feature dependency key as an additional
        // guard: an unsupported or contextual producer cannot enter this path.
        let feature = self.key(graph.feature)?;
        let mut bytes = self.salt.clone();
        bytes.extend(b"literal-machine-screw-symbol/1");
        bytes.extend(id.to_string().bytes());
        bytes.extend(feature.0);
        for entity in graph.entities {
            bytes.extend(entity.to_string().bytes());
            bytes.extend(self.hash(entity, Visit::Full)?);
        }
        Some(CacheKey(ContentHash::of(&bytes).0))
    }
    pub(crate) fn lookup_symbol(
        &mut self,
        key: &CacheKey,
        id: EntityId,
        revision: u64,
    ) -> Option<SymbolNode> {
        let entry = self.cache.symbol_nodes.get_mut(key)?;
        let node = entry.node.shared(revision)?;
        entry.used = self.cache.clock;
        self.stats.symbol_reused.insert(id);
        Some(node)
    }
    pub(crate) fn evaluated_symbol(&mut self, id: EntityId) {
        self.stats.symbol_evaluated.insert(id);
    }
    pub(crate) fn stage_symbol(&mut self, key: CacheKey, node: SymbolNode) {
        let bytes = node.bytes();
        let budget = self.cache.retained_bytes.min(64 * 1024 * 1024);
        if self.cache.max_entries == 0 || bytes > budget {
            return;
        }
        // Stage only complete node+leaf results; promotion happens after the
        // outer full compile passes its final reference/cancellation checks.
        let used = self.cache.clock;
        let limit = self.cache.max_entries.min(768);
        self.cache.pending_symbol_nodes.remove(&key);
        while self.cache.pending_symbol_nodes.len() >= limit
            || self
                .cache
                .pending_symbol_nodes
                .values()
                .map(|entry| entry.bytes)
                .sum::<usize>()
                > budget - bytes
        {
            let Some(oldest) = self
                .cache
                .pending_symbol_nodes
                .iter()
                .min_by_key(|(_, entry)| entry.used)
                .map(|(key, _)| *key)
            else {
                break;
            };
            self.cache.pending_symbol_nodes.remove(&oldest);
        }
        self.cache
            .pending_symbol_nodes
            .insert(key, SymbolEntry { node, bytes, used });
    }
    /// The key of an Assembly/Instance composition. `content` must record
    /// every non-solid input the node reads; solids are matched separately.
    pub(crate) fn composition_key(&self, id: EntityId, content: &[u8]) -> CacheKey {
        let mut bytes = self.salt.clone();
        bytes.extend(b"composition/1|");
        bytes.extend(id.to_string().bytes());
        bytes.push(b'|');
        bytes.extend(content);
        CacheKey(ContentHash::of(&bytes).0)
    }
    /// An earlier composition built from exactly these input solids.
    pub(crate) fn lookup_composition(
        &mut self,
        key: &CacheKey,
        inputs: &[Solid],
        id: EntityId,
        revision: u64,
    ) -> Option<SymbolNode> {
        let entry = self
            .cache
            .compositions
            .get_mut(key)?
            .iter_mut()
            .find(|entry| entry.matches(inputs))?;
        let node = entry.node.shared(revision)?;
        entry.used = self.cache.clock;
        self.stats.composition_reused.insert(id);
        Some(node)
    }
    pub(crate) fn evaluated_composition(&mut self, id: EntityId) {
        self.stats.composition_evaluated.insert(id);
    }
    /// Stage a complete composition; it is promoted after a successful compile.
    pub(crate) fn stage_composition(
        &mut self,
        key: CacheKey,
        inputs: Vec<Solid>,
        node: SymbolNode,
    ) {
        let bytes = node.bytes();
        let used = self.cache.clock;
        self.cache.pending_compositions.push((
            key,
            CompositionEntry {
                inputs,
                node,
                bytes,
                used,
            },
        ));
    }
    /// Reuse an earlier exact result, re-stamped for this revision.
    pub(crate) fn lookup(
        &mut self,
        key: &CacheKey,
        revision: u64,
    ) -> Option<(Solid, Arc<ShapeSnapshot>, Extras)> {
        let clock = self.cache.clock;
        let entry = self.cache.entries.get_mut(key)?;
        let solid = entry.solid.share().ok()?;
        if entry.shape.topology.revision != revision {
            let mut shape = ShapeSnapshot::clone(&entry.shape);
            shape.topology.revision = revision;
            entry.shape = Arc::new(shape);
        }
        entry.used = clock;
        self.stats.reused.insert(entry.shape.producer);
        Some((solid, entry.shape.clone(), entry.extras.clone()))
    }
    pub(crate) fn evaluated(&mut self, feature: EntityId, cacheable: bool) {
        self.stats.evaluated.insert(feature);
        if !cacheable {
            self.stats.uncacheable.insert(feature);
        }
    }
    pub(crate) fn store(
        &mut self,
        key: CacheKey,
        solid: &Solid,
        shape: &Arc<ShapeSnapshot>,
        extras: Extras,
    ) {
        if let Ok(solid) = solid.share() {
            self.cache.entries.insert(
                key,
                Entry {
                    solid,
                    shape: shape.clone(),
                    extras,
                    bytes: estimate(shape),
                    used: self.cache.clock,
                },
            );
        }
    }
}
impl Drop for Session<'_> {
    fn drop(&mut self) {
        if std::thread::panicking() {
            self.cache.pipe_threads.finish_compile(false);
        }
        let mut stats = std::mem::take(&mut self.stats);
        stats.evicted = self.cache.evict();
        if self.cache.full_rust.is_some()
            && self.cache.estimated_bytes() > self.cache.retained_bytes
        {
            self.cache.full_rust = None;
        }
        stats.entries = self.cache.entries.len();
        stats.estimated_bytes = self.cache.estimated_bytes();
        self.cache.last = stats;
    }
}

/// Entities a provider finds by scanning the document rather than through
/// its own source. None when they cannot be determined.
fn hidden_inputs(document: &Document, feature: EntityId, type_id: &str) -> Option<Vec<EntityId>> {
    if type_id == model::LOOSE_TENON_TYPE {
        let source = model::read_loose_tenon(document, feature).ok()?;
        let joint = model::read_panel_joint(document, source.joint).ok()?;
        let pair = joint.pairs.get(source.pair as usize)?;
        // The provider resolves each Part's current base panel from its
        // feature list. Include those definitions even for older joint links.
        let mut inputs = vec![
            model::panel_feature(document, pair.a)?,
            model::panel_feature(document, pair.b)?,
        ];
        inputs.extend(model::panel_placements(document, pair.a));
        inputs.extend(model::panel_placements(document, pair.b));
        return Some(inputs);
    }
    if type_id == model::PANEL_JOINT_CUT_TYPE {
        // The cut is resolved where both panels are now, so every Move /
        // rotate of either panel changes it (CAD-431).
        let cut = model::read_panel_joint_cut(document, feature).ok()?;
        let joint = model::read_panel_joint(document, cut.joint).ok()?;
        let pair = joint.pairs.get(cut.pair as usize)?;
        let mut inputs = model::panel_placements(document, pair.a);
        inputs.extend(model::panel_placements(document, pair.b));
        return Some(inputs);
    }
    if type_id != model::MEMBER_TYPE {
        return Some(vec![]);
    }
    // A member's end cuts come from every connection naming its Part and
    // from each connected partner member's resolved path and profile.
    let part = model::read_member(document, feature).ok()?.output.body;
    let mut inputs = vec![];
    for (connection, definition) in model::member_connections(document, part).ok()? {
        inputs.push(connection);
        for member in definition.members {
            if member != part {
                inputs.push(model::member_feature(document, member)?);
            }
        }
    }
    Some(inputs)
}

/// An entity's `serde_json::to_vec` bytes and every identity [`mentions`]
/// finds in its `serde_json::to_value` form, from one serialization pass.
struct Source {
    bytes: Vec<u8>,
    mentioned: BTreeSet<EntityId>,
}
impl Source {
    fn of(entity: &multicad_core::document::Entity) -> Option<Self> {
        use serde::Serialize as _;
        let mut bytes = Vec::new();
        let mut formatter = Scanning::default();
        let serialized = entity
            .serialize(&mut serde_json::Serializer::with_formatter(
                &mut bytes,
                &mut formatter,
            ))
            .is_ok();
        if serialized && !formatter.unusual_key {
            return Some(Self {
                bytes,
                mentioned: formatter.found,
            });
        }
        // A non-string map key renders differently in a Value; keep the
        // original two-pass derivation for that rare shape (and its errors).
        let bytes = serde_json::to_vec(entity).ok()?;
        let mut mentioned = BTreeSet::new();
        mentions(&serde_json::to_value(entity).ok()?, &mut mentioned);
        Some(Self { bytes, mentioned })
    }
}
/// The default compact JSON output that also scans every string, object key
/// included, exactly as decoded text, the way [`mentions`] scans a Value.
#[derive(Default)]
struct Scanning {
    string: Option<String>,
    found: BTreeSet<EntityId>,
    /// A number or literal was written inside a string (a non-string key).
    unusual_key: bool,
}
macro_rules! scanning_numbers {
    ($($name:ident: $ty:ty),*) => {$(
        fn $name<W: ?Sized + std::io::Write>(&mut self, writer: &mut W, value: $ty) -> std::io::Result<()> {
            self.unusual_key |= self.string.is_some();
            serde_json::ser::CompactFormatter.$name(writer, value)
        }
    )*};
}
impl serde_json::ser::Formatter for &mut Scanning {
    scanning_numbers!(
        write_bool: bool, write_i8: i8, write_i16: i16, write_i32: i32, write_i64: i64,
        write_i128: i128, write_u8: u8, write_u16: u16, write_u32: u32, write_u64: u64,
        write_u128: u128, write_f32: f32, write_f64: f64
    );
    fn write_null<W: ?Sized + std::io::Write>(&mut self, writer: &mut W) -> std::io::Result<()> {
        self.unusual_key |= self.string.is_some();
        serde_json::ser::CompactFormatter.write_null(writer)
    }
    fn write_number_str<W: ?Sized + std::io::Write>(
        &mut self,
        writer: &mut W,
        value: &str,
    ) -> std::io::Result<()> {
        self.unusual_key |= self.string.is_some();
        serde_json::ser::CompactFormatter.write_number_str(writer, value)
    }
    fn begin_string<W: ?Sized + std::io::Write>(&mut self, writer: &mut W) -> std::io::Result<()> {
        self.string = Some(String::new());
        serde_json::ser::CompactFormatter.begin_string(writer)
    }
    fn end_string<W: ?Sized + std::io::Write>(&mut self, writer: &mut W) -> std::io::Result<()> {
        if let Some(text) = self.string.take() {
            scan(&text, &mut self.found);
        }
        serde_json::ser::CompactFormatter.end_string(writer)
    }
    fn write_string_fragment<W: ?Sized + std::io::Write>(
        &mut self,
        writer: &mut W,
        fragment: &str,
    ) -> std::io::Result<()> {
        if let Some(text) = &mut self.string {
            text.push_str(fragment);
        }
        serde_json::ser::CompactFormatter.write_string_fragment(writer, fragment)
    }
    fn write_char_escape<W: ?Sized + std::io::Write>(
        &mut self,
        writer: &mut W,
        escape: serde_json::ser::CharEscape,
    ) -> std::io::Result<()> {
        use serde_json::ser::CharEscape;
        if let Some(text) = &mut self.string {
            text.push(match escape {
                CharEscape::Quote => '"',
                CharEscape::ReverseSolidus => '\\',
                CharEscape::Solidus => '/',
                CharEscape::Backspace => '\u{8}',
                CharEscape::FormFeed => '\u{c}',
                CharEscape::LineFeed => '\n',
                CharEscape::CarriageReturn => '\r',
                CharEscape::Tab => '\t',
                CharEscape::AsciiControl(byte) => char::from(byte),
            });
        }
        serde_json::ser::CompactFormatter.write_char_escape(writer, escape)
    }
}

/// Every existing-format entity identity mentioned in a source value:
/// identity strings, identities embedded in longer strings and map keys.
fn mentions(value: &serde_json::Value, found: &mut BTreeSet<EntityId>) {
    match value {
        serde_json::Value::String(text) => scan(text, found),
        serde_json::Value::Array(items) => items.iter().for_each(|v| mentions(v, found)),
        serde_json::Value::Object(map) => {
            for (key, value) in map {
                scan(key, found);
                mentions(value, found);
            }
        }
        _ => {}
    }
}
fn scan(text: &str, found: &mut BTreeSet<EntityId>) {
    if let Ok(id) = EntityId::parse(text) {
        found.insert(id);
        return;
    }
    let bytes = text.as_bytes();
    if bytes.len() < 36 {
        return;
    }
    for start in 0..=bytes.len() - 36 {
        let window = &bytes[start..start + 36];
        // Reject on the four separator positions before the hex digits.
        if [8, 13, 18, 23].into_iter().any(|i| window[i] != b'-') {
            continue;
        }
        let shaped = window.iter().enumerate().all(|(i, b)| match i {
            8 | 13 | 18 | 23 => *b == b'-',
            _ => b.is_ascii_hexdigit(),
        });
        if shaped
            && let Ok(text) = std::str::from_utf8(window)
            && let Ok(id) = EntityId::parse(text)
        {
            found.insert(id);
        }
    }
}

/// Side-result maps a provider may write next to its solid.
pub(crate) struct Sides<'a> {
    pub panel_frames: &'a mut BTreeMap<EntityId, PlaneFrame>,
    pub sheet_layouts: &'a mut BTreeMap<EntityId, SheetLayout>,
    pub member_paths: &'a mut BTreeMap<EntityId, model::PathSegment>,
    pub consumed: &'a mut BTreeSet<(EntityId, EntityId)>,
}
impl Extras {
    pub(crate) fn apply(self, feature: EntityId, sides: Sides) {
        if let Some(frame) = self.panel_frame {
            sides.panel_frames.insert(feature, frame);
        }
        if let Some(layout) = self.sheet_layout {
            sides.sheet_layouts.insert(feature, layout);
        }
        if let Some(path) = self.member_path {
            sides.member_paths.insert(feature, path);
        }
        sides.consumed.extend(self.consumed);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mentions_find_identities_in_values_strings_and_keys() {
        let [a, b, c, d] = std::array::from_fn(|_| EntityId::new());
        let value = serde_json::json!({
            "target": a.to_string(),
            "role": format!("side/{b}/panel.v_start"),
            "by_id": { c.to_string(): [1, 2, format!("x{d}")] },
            "noise": "0123456789abcdef-not-an-identity",
        });
        let mut found = BTreeSet::new();
        mentions(&value, &mut found);
        assert_eq!(found, BTreeSet::from([a, b, c, d]));
    }

    #[test]
    fn one_pass_sources_equal_the_value_derivation() {
        let check = |entity: &multicad_core::document::Entity| {
            let source = Source::of(entity).unwrap();
            assert_eq!(source.bytes, serde_json::to_vec(entity).unwrap());
            let mut expected = BTreeSet::new();
            mentions(&serde_json::to_value(entity).unwrap(), &mut expected);
            assert_eq!(source.mentioned, expected);
        };
        let document = multicad_core::project_file::read(
            &mut include_bytes!("../../../examples/toolcase.tcd").as_slice(),
            Default::default(),
            &multicad_core::project_file::supported_capabilities(),
        )
        .unwrap()
        .project
        .document;
        let mut mentioning = 0;
        for entity in document.entities() {
            check(entity);
            mentioning += usize::from(!Source::of(entity).unwrap().mentioned.is_empty());
        }
        assert!(mentioning > 100, "{mentioning}");
        let [a, b, c] = std::array::from_fn(|_| EntityId::new());
        for name in [
            format!("{a}"),
            format!("x\"{a}\"\n\t{b}\\/\u{1}{c}"),
            a.to_string().replace('-', ""),
            format!("\u{1}{a}")[..20].to_string() + &b.to_string(),
            format!("{a}{b}"),
        ] {
            check(&multicad_core::document::Entity::new(
                name,
                EntityData::Material {
                    density_kg_m3: Some(-0.),
                    color_rgba: [1, 2, 3, 4],
                    section_hatch: None,
                },
            ));
        }
    }

    #[test]
    fn cacheable_types_exclude_scanning_and_unknown_providers() {
        assert!(cacheable_type(model::EXTRUDE_TYPE));
        assert!(cacheable_type(model::MEMBER_TYPE));
        assert!(cacheable_type(model::FENCE_MEMBER_TYPE));
        for excluded in [
            model::WALL_TYPE,
            model::SLAB_TYPE,
            model::STAIR_TYPE,
            model::TRANSFER_TYPE,
            model::SCALE_TYPE,
            "vendor.unknown",
        ] {
            assert!(!cacheable_type(excluded), "{excluded}");
        }
    }
}
