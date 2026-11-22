//! Request-local compression of successful, exact literal-cell certificates.
//! Grid and mesh coordinates are shared metadata, never a geometry substitute.
use super::*;
use std::mem::size_of;

pub(super) struct Snapshot {
    u: Box<[Box<[f64]>]>,
    v: Box<[f64]>,
    vertices: Box<[f64]>,
    indices: Box<[u32]>,
    matched: Box<[bool]>,
}
impl Snapshot {
    fn bytes(&self) -> usize {
        size_of::<Self>()
            + self.u.len() * size_of::<Box<[f64]>>()
            + self.u.iter().map(|u| size_of_val(&**u)).sum::<usize>()
            + size_of_val(&*self.v)
            + size_of_val(&*self.vertices)
            + size_of_val(&*self.indices)
            + size_of_val(&*self.matched)
    }
    fn matches(&self, key: &[u64; 18], work: &mut Work<'_>) -> Result<bool, Error> {
        // Binary searches, reconstruction and the full 18-word comparison are
        // charged even on a miss. Source and options are immutable Prepared fields.
        work.charge(64)?;
        let side = usize::try_from(key[0]).map_err(|_| Error::Invalid)?;
        let Some(nodes) = self.u.get(side) else {
            return Ok(false);
        };
        let locate = |nodes: &[f64], lo: u64, hi: u64| {
            nodes
                .binary_search_by(|x| x.total_cmp(&f64::from_bits(lo)))
                .ok()
                .filter(|&i| nodes.get(i + 1).is_some_and(|x| x.to_bits() == hi))
        };
        let Some(u) = locate(nodes, key[1], key[2]) else {
            return Ok(false);
        };
        let Some(row) = locate(&self.v, key[3], key[4]) else {
            return Ok(false);
        };
        work.charge(self.u.len().checked_mul(2).ok_or(Error::Limit)?)?;
        let n: usize = self.u.iter().map(|u| u.len() - 1).sum();
        let offset: usize = self.u[..side].iter().map(|u| u.len() - 1).sum();
        let i = offset + u;
        let cell = row * n + i;
        if !self.matched[cell] {
            return Ok(false);
        }
        let ids = [
            row * n + i,
            row * n + (i + 1) % n,
            (row + 1) * n + (i + 1) % n,
            (row + 1) * n + i,
        ];
        work.charge(128)?;
        let mut reconstructed = [0; 18];
        reconstructed[..5].copy_from_slice(&[
            side as u64,
            nodes[u].to_bits(),
            nodes[u + 1].to_bits(),
            self.v[row].to_bits(),
            self.v[row + 1].to_bits(),
        ]);
        for (j, &index) in self.indices[6 * cell..6 * cell + 6].iter().enumerate() {
            let Some(local) = ids.iter().position(|&id| id == index as usize) else {
                return Ok(false);
            };
            reconstructed[5] |= (local as u64) << (2 * j);
        }
        for (slot, value) in reconstructed[6..].iter_mut().zip(
            ids.into_iter()
                .flat_map(|i| &self.vertices[3 * i..3 * i + 3]),
        ) {
            *slot = value.to_bits();
        }
        Ok(reconstructed == *key)
    }
}
impl Prepared<'_> {
    fn cache_bytes(&self) -> usize {
        // Existing certificate allowance, counting only the old key payload:
        // this deliberately does not borrow uncounted BTree node overhead.
        self.limits
            .certificates
            .saturating_mul(size_of::<[u64; 18]>())
    }
    pub(super) fn begin_cell_cache(
        &mut self,
        grid: &Grid,
        ring: &Ring,
        work: &mut Work<'_>,
    ) -> Result<(), Error> {
        self.cell_matches_pending = None;
        if !self.limits.compact_certificates {
            return Ok(());
        }
        let count = (grid.v.len() - 1)
            .checked_mul(ring.entries.len())
            .ok_or(Error::Limit)?;
        let resident = self.cell_cache.as_ref().map_or(0, Snapshot::bytes);
        let bytes = count
            .checked_mul(size_of::<bool>())
            .and_then(|b| b.checked_add(size_of::<Box<[bool]>>()))
            .ok_or(Error::Limit)?;
        if resident
            .checked_add(bytes)
            .is_some_and(|b| b <= self.cache_bytes())
        {
            work.charge(count)?;
            self.cell_matches_pending = Some(vec![false; count].into_boxed_slice());
        }
        Ok(())
    }
    pub(super) fn cached_cell(&self, key: &[u64; 18], work: &mut Work<'_>) -> Result<bool, Error> {
        if self.limits.compact_certificates {
            return self
                .cell_cache
                .as_ref()
                .map_or(Ok(false), |s| s.matches(key, work));
        }
        Ok(self.certificates.contains(key))
    }
    pub(super) fn mark_cell_match(&mut self, index: usize) {
        if let Some(flags) = &mut self.cell_matches_pending {
            flags[index] = true;
        }
    }
    pub(super) fn finish_cell_cache(
        &mut self,
        mesh: &DisplayMesh,
        grid: &Grid,
        work: &mut Work<'_>,
    ) -> Result<(), Error> {
        let Some(matched) = self.cell_matches_pending.take() else {
            return Ok(());
        };
        let bytes = size_of::<Snapshot>()
            .checked_add(grid.u.len() * size_of::<Box<[f64]>>())
            .and_then(|b| b.checked_add(grid.u.iter().map(|u| size_of_val(&**u)).sum::<usize>()))
            .and_then(|b| b.checked_add(size_of_val(&*grid.v)))
            .and_then(|b| b.checked_add(size_of_val(&*mesh.vertices)))
            .and_then(|b| b.checked_add(size_of_val(&*mesh.indices)))
            .and_then(|b| b.checked_add(size_of_val(&*matched)))
            .ok_or(Error::Limit)?;
        if bytes > self.cache_bytes() {
            return Ok(());
        }
        // No old/new metadata overlap: only the bounded success flags coexist
        // with the old snapshot. Failed cells never become cache hits.
        self.cell_cache = None;
        let elements = grid.u.iter().map(Vec::len).sum::<usize>()
            + grid.v.len()
            + mesh.vertices.len()
            + mesh.indices.len();
        work.charge(elements)?;
        let snapshot = Snapshot {
            u: grid
                .u
                .iter()
                .map(|u| u.clone().into_boxed_slice())
                .collect(),
            v: grid.v.clone().into_boxed_slice(),
            vertices: mesh.vertices.clone().into_boxed_slice(),
            indices: mesh.indices.clone().into_boxed_slice(),
            matched,
        };
        debug_assert_eq!(snapshot.bytes(), bytes);
        self.cell_cache = Some(snapshot);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn fixture() -> (Snapshot, [u64; 18]) {
        let vertices = vec![0., 0., 0., 1., 0., 0., 0., 1., 0., 1., 1., 0.];
        let ids = [0, 1, 3, 2];
        let mut key = [0; 18];
        key[..5].copy_from_slice(&[
            0,
            0f64.to_bits(),
            0.5f64.to_bits(),
            0f64.to_bits(),
            1f64.to_bits(),
        ]);
        for (i, local) in [0, 1, 2, 0, 2, 3].into_iter().enumerate() {
            key[5] |= local << (2 * i);
        }
        for (slot, value) in key[6..]
            .iter_mut()
            .zip(ids.into_iter().flat_map(|i| &vertices[3 * i..3 * i + 3]))
        {
            *slot = f64::to_bits(*value);
        }
        (
            Snapshot {
                u: vec![vec![0., 0.5, 1.].into_boxed_slice()].into_boxed_slice(),
                v: vec![0., 1.].into_boxed_slice(),
                vertices: vertices.into_boxed_slice(),
                indices: vec![0, 1, 3, 0, 3, 2, 1, 0, 2, 1, 2, 3].into_boxed_slice(),
                matched: vec![true, false].into_boxed_slice(),
            },
            key,
        )
    }
    #[test]
    fn exact_shared_metadata_key_and_success_only_controls() {
        let (mut snapshot, key) = fixture();
        let mut work = Work {
            remaining: 100_000,
            interrupt: None,
        };
        assert!(snapshot.matches(&key, &mut work).unwrap());
        for slot in 0..18 {
            let mut changed = key;
            changed[slot] ^= 1;
            assert!(
                !snapshot.matches(&changed, &mut work).unwrap(),
                "slot {slot}"
            );
        }
        snapshot.matched[0] = false;
        assert!(!snapshot.matches(&key, &mut work).unwrap());
        snapshot.matched[0] = true;
        snapshot.indices[0] = 1;
        assert!(!snapshot.matches(&key, &mut work).unwrap());
    }
    #[test]
    fn shared_metadata_lookup_charges_and_cancels_then_retries() {
        let (snapshot, key) = fixture();
        let mut work = Work {
            remaining: 0,
            interrupt: None,
        };
        assert_eq!(snapshot.matches(&key, &mut work), Err(Error::Limit));
        let interrupt = || true;
        let mut work = Work {
            remaining: 100_000,
            interrupt: Some(&interrupt),
        };
        assert_eq!(snapshot.matches(&key, &mut work), Err(Error::Cancelled));
        let mut work = Work {
            remaining: 100_000,
            interrupt: None,
        };
        assert!(snapshot.matches(&key, &mut work).unwrap());
        assert!(work.remaining < 100_000);
    }
}
