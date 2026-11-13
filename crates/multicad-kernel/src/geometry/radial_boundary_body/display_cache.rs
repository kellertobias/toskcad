//! Bounded detached display data. Occupied geometry and manufacturing stay fresh.
use super::*;
use multicad_rust_kernel::display::Tessellation;
use std::{collections::VecDeque, rc::Weak};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Quality {
    deflection: u64,
    angle: u64,
    triangles: usize,
}
impl From<Tessellation> for Quality {
    fn from(options: Tessellation) -> Self {
        Self {
            deflection: options.deflection.to_bits(),
            angle: options.angle.to_bits(),
            triangles: options.triangle_limit,
        }
    }
}
// A Weak keeps its Rc allocation alive after the value is dropped. Reserve
// the inline Placement and both counters, with conservative alignment slack.
const WEAK_OWNER_BYTES: usize =
    size_of::<Placement>() + 2 * size_of::<usize>() + 2 * align_of::<Placement>();
struct Entry {
    owner: Weak<Placement>,
    quality: Quality,
    mesh: Mesh,
    bytes: usize,
}
struct Cache {
    entries: VecDeque<Entry>,
    bytes: usize,
    byte_limit: usize,
    entry_limit: usize,
    receipt_limit: usize,
}
impl Default for Cache {
    fn default() -> Self {
        Self::new(
            128 * 1024 * 1024 - size_of::<RefCell<Cache>>(),
            1024,
            32 * 1024 * 1024,
        )
    }
}
impl Cache {
    fn new(byte_limit: usize, entry_limit: usize, receipt_limit: usize) -> Self {
        // Reserve the bounded descriptor queue once and charge its actual
        // capacity, including every key, Weak owner and detached Vec header.
        let entries = VecDeque::with_capacity(entry_limit.min(byte_limit / size_of::<Entry>()));
        let bytes = entries.capacity() * size_of::<Entry>();
        let entry_limit = entry_limit.min(entries.capacity());
        Self {
            entries,
            bytes,
            byte_limit,
            entry_limit,
            receipt_limit,
        }
    }
    fn remove(&mut self, index: usize) -> Entry {
        let entry = self.entries.remove(index).expect("live cache index");
        self.bytes -= entry.bytes;
        entry
    }
    fn expire(&mut self) {
        let mut index = 0;
        while index < self.entries.len() {
            if self.entries[index].owner.upgrade().is_none() {
                self.remove(index);
            } else {
                index += 1;
            }
        }
    }
    fn owner_index(&self, owner: &Rc<Placement>) -> Option<usize> {
        self.entries.iter().position(|entry| {
            entry
                .owner
                .upgrade()
                .is_some_and(|cached| Rc::ptr_eq(&cached, owner))
        })
    }
    fn get(&mut self, owner: &Rc<Placement>, quality: Quality) -> Option<Mesh> {
        self.expire();
        let index = self.owner_index(owner)?;
        if self.entries[index].quality != quality {
            return None;
        }
        let entry = self.remove(index);
        let mesh = entry.mesh.clone();
        self.bytes += entry.bytes;
        self.entries.push_back(entry);
        Some(mesh)
    }
    fn remember(&mut self, owner: &Rc<Placement>, quality: Quality, mesh: &Mesh) {
        self.expire();
        // Preflight the source capacities before cloning rejected large data.
        if payload_bytes(mesh).is_some_and(|bytes| self.admits(bytes)) {
            self.remember_owned(owner, quality, mesh.clone());
        }
    }
    fn remember_owned(&mut self, owner: &Rc<Placement>, quality: Quality, mesh: Mesh) {
        if let Some(bytes) = payload_bytes(&mesh).filter(|&bytes| self.admits(bytes)) {
            if let Some(index) = self.owner_index(owner) {
                self.remove(index);
            }
            while self.entries.len() >= self.entry_limit || self.bytes > self.byte_limit - bytes {
                self.remove(0);
            }
            self.bytes += bytes;
            self.entries.push_back(Entry {
                owner: Rc::downgrade(owner),
                quality,
                mesh,
                bytes,
            });
        }
    }
    fn admits(&self, bytes: usize) -> bool {
        bytes <= self.receipt_limit
            && bytes
                <= self
                    .byte_limit
                    .saturating_sub(self.entries.capacity() * size_of::<Entry>())
            && self.entry_limit > 0
    }
}
fn payload_bytes(mesh: &Mesh) -> Option<usize> {
    payload_capacities([
        mesh.vertices.capacity(),
        mesh.indices.capacity(),
        mesh.face_indices.capacity(),
    ])?
    .checked_add(WEAK_OWNER_BYTES)
}
fn payload_capacities([vertices, indices, faces]: [usize; 3]) -> Option<usize> {
    vertices
        .checked_mul(size_of::<f64>())?
        .checked_add(indices.checked_mul(size_of::<u32>())?)?
        .checked_add(faces.checked_mul(size_of::<u32>())?)
}

thread_local! {
    static DISPLAY: RefCell<Cache> = RefCell::new(Cache::default());
}
pub(super) fn mesh(
    owner: &Rc<Placement>,
    options: Tessellation,
    usage: MeshUse,
    interrupt: Option<&dyn Fn() -> bool>,
    mut build: impl FnMut() -> Result<Mesh, GeometryError>,
) -> Result<Mesh, GeometryError> {
    // Keep one cache/publication path regardless of the source builder's type.
    mesh_erased(owner, options, usage, interrupt, &mut build)
}
fn mesh_erased(
    owner: &Rc<Placement>,
    options: Tessellation,
    usage: MeshUse,
    interrupt: Option<&dyn Fn() -> bool>,
    build: &mut dyn FnMut() -> Result<Mesh, GeometryError>,
) -> Result<Mesh, GeometryError> {
    checkpoint(interrupt)?;
    let quality = Quality::from(options);
    if usage == MeshUse::Display
        && let Some(mesh) = DISPLAY.with(|cache| cache.borrow_mut().get(owner, quality))
    {
        checkpoint(interrupt)?;
        return Ok(mesh);
    }
    let mesh = build()?;
    checkpoint(interrupt)?;
    if !valid_parts(&mesh.vertices, &mesh.indices, &mesh.face_indices) {
        return Err(error(
            "retained display mesh",
            GeometryErrorKind::InvalidResult,
            "Invalid detached mesh arrays",
        ));
    }
    if usage == MeshUse::Display {
        DISPLAY.with(|cache| cache.borrow_mut().remember(owner, quality, &mesh));
    }
    Ok(mesh)
}

#[cfg(test)]
mod tests;
