//! Bounded Gauss2 with the full source fourth derivative and exact dyadic domains.
use super::*;
struct Cell {
    patch: usize,
    ends: [f64; 2],
    depth: u8,
    value: Values,
    slot: usize,
}
fn zero() -> Values {
    [I::point(0.); 4]
}

fn cell(
    atlas: &PatchAtlas<'_>,
    patch: usize,
    ends: [f64; 2],
    depth: u8,
    work: &mut Work<'_>,
) -> Result<Cell, Error> {
    work.step()?;
    let width = I::point(ends[1]) - I::point(ends[0]);
    let mid = (I::point(ends[0]) + I::point(ends[1])) / I::point(2.);
    let offset = width / I::point(12.).sqrt();
    let a = density::evaluate(atlas, patch, mid - offset)?;
    let b = density::evaluate(atlas, patch, mid + offset)?;
    let d = density::evaluate(atlas, patch, I::new(ends[0], ends[1]))?;
    let mut value = zero();
    for i in 0..4 {
        let error =
            I::point(d[i].c[4].magnitude()) * width.square().square() * width / I::point(180.);
        value[i] = ((a[i].c[0] + b[i].c[0]) * width / I::point(2.) + I::new(-error.hi, error.hi))
            .checked()
            .map_err(|_| Error::Unresolved)?;
    }
    Ok(Cell {
        patch,
        ends,
        depth,
        value,
        slot: 0,
    })
}
fn halves(ends: [f64; 2], depth: u8) -> Result<[[f64; 2]; 2], Error> {
    // All domains descend from [0,1] by exact binary bisection. At most48
    // halvings keep every endpoint and midpoint exactly representable.
    if depth >= 48 {
        return Err(Error::Unresolved);
    }
    let mid = ends[0].midpoint(ends[1]);
    // These private domains are exclusively the depth<=47 dyadic children
    // of [0,1]; their midpoint is strictly interior and exactly representable.
    Ok([[ends[0], mid], [mid, ends[1]]])
}
fn ready(
    atlas: &PatchAtlas<'_>,
    patch: usize,
    ends: [f64; 2],
    depth: u8,
    work: &mut Work<'_>,
) -> Result<Vec<Cell>, Error> {
    let mut pending = vec![(ends, depth)];
    let mut result = Vec::new();
    while let Some((ends, depth)) = pending.pop() {
        match cell(atlas, patch, ends, depth, work) {
            Ok(cell) => result.push(cell),
            Err(Error::Unresolved) => {
                for child in halves(ends, depth)? {
                    pending.push((child, depth + 1));
                }
            }
            Err(e) => return Err(e),
        }
    }
    Ok(result)
}
fn insert(heap: &mut BTreeMap<u64, Vec<Cell>>, cell: Cell, length: f64) {
    // Priority is a heuristic only; acceptance uses final dimensioned widths.
    let score = cell
        .value
        .into_iter()
        .enumerate()
        .map(|(i, v)| v.width() / if i == 0 { 1. } else { length })
        .fold(0_f64, f64::max);
    heap.entry(score.to_bits()).or_default().push(cell);
}
// Each retained cell occupies one fixed tree leaf. Replacing a cell recomputes
// only its O(log maximum_cells) ancestors, never subtracting rounded partial
// sums or repeatedly rebuilding every live leaf. All leaf slots are bounded
// by cumulative successful cell proofs; the hard65536 ceiling bounds this
// complete tree to8MiB plus the separately bounded source/heap descriptors.
struct SumTree {
    nodes: Vec<Values>,
    capacity: usize,
    next: usize,
}
impl SumTree {
    fn new(limit: usize) -> Self {
        let capacity = limit.next_power_of_two();
        Self {
            nodes: vec![zero(); 2 * capacity],
            capacity,
            next: 0,
        }
    }
    fn set(&mut self, slot: usize, value: Values) {
        let mut at = self.capacity + slot;
        self.nodes[at] = value;
        while at > 1 {
            at /= 2;
            self.nodes[at] =
                std::array::from_fn(|i| self.nodes[2 * at][i] + self.nodes[2 * at + 1][i]);
        }
    }
    fn insert(&mut self, value: Values) -> usize {
        let slot = self.next;
        self.next += 1;
        // Every insertion follows a charged successful cell proof; there
        // cannot be more insertions than the original cumulative proof quota.
        assert!(
            slot < self.capacity,
            "cell proof quota bounds sum leaf slots"
        );
        self.set(slot, value);
        slot
    }
    fn total(&self) -> Values {
        self.nodes[1]
    }
}
pub(super) fn face(
    atlas: &PatchAtlas<'_>,
    face: usize,
    key: FaceKey,
    options: FacePropertyOptions,
    work: &mut Work<'_>,
) -> Result<FaceProperties, Error> {
    let mut heap = BTreeMap::<u64, Vec<Cell>>::new();
    let length = atlas
        .topology
        .events()
        .source()
        .bounding_box_mm()
        .into_iter()
        .flatten()
        .map(f64::abs)
        .fold(1., f64::max);
    let mut total = SumTree::new(options.maximum_cells);
    for (index, _) in atlas
        .patches
        .iter()
        .enumerate()
        .filter(|(_, p)| p.face == face)
    {
        for mut c in ready(atlas, index, [0., 1.], 0, work)? {
            c.slot = total.insert(c.value);
            insert(&mut heap, c, length);
        }
    }
    loop {
        work.check()?;
        if let Some(result) = receipt(key, total.total(), options) {
            return Ok(result);
        }
        let (score, mut bucket) = heap.pop_last().ok_or(Error::Unresolved)?;
        let old = bucket
            .pop()
            .expect("only nonempty score buckets are inserted");
        if !bucket.is_empty() {
            heap.insert(score, bucket);
        }
        total.set(old.slot, zero());
        for ends in halves(old.ends, old.depth)? {
            for mut c in ready(atlas, old.patch, ends, old.depth + 1, work)? {
                c.slot = total.insert(c.value);
                insert(&mut heap, c, length);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn exact_dyadic_query_domains_stop_at_the_declared_derivative_resolution() {
        let mut ends = [0., 1.];
        for depth in 0..48 {
            let [left, right] = halves(ends, depth).unwrap();
            assert_eq!(left[0], ends[0]);
            assert_eq!(right[1], ends[1]);
            assert_eq!(left[1], right[0]);
            assert!(left[1] > ends[0] && left[1] < ends[1]);
            assert_eq!(left[1], 2_f64.powi(-i32::from(depth) - 1));
            ends = left;
        }
        assert!(matches!(halves(ends, 48), Err(Error::Unresolved)));
    }
}
