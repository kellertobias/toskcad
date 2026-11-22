//! Sufficient original-position C0 and same-source Q7+R mean-value gates.
//! The prefix marks EVERY undecided leaf for the unchanged B7 fallback.
use super::*;
pub(super) struct Needs {
    prefix: Vec<usize>,
    divisions: [usize; 2],
}
impl Needs {
    pub(super) fn bytes(&self) -> usize {
        self.prefix.capacity() * std::mem::size_of::<usize>()
    }
    pub(super) fn contains(&self, domain: [[f64; 2]; 2]) -> bool {
        let [u, v] = std::array::from_fn::<_, 2, _>(|k| {
            domain[k].map(|x| (x * self.divisions[k] as f64) as usize)
        });
        let width = self.divisions[1] + 1;
        self.prefix[u[1] * width + v[1]] + self.prefix[u[0] * width + v[0]]
            > self.prefix[u[0] * width + v[1]] + self.prefix[u[1] * width + v[0]]
    }
    pub(super) fn discover(
        grid: &mut Grid,
        root: &position::Root,
        lattice: &position::Lattice,
        normal: &ReducedNormalCell11,
        patch: Patch,
        quality: Tessellation,
        work: &mut Work<'_>,
    ) -> R<Self> {
        let [nu, nv] = grid.divisions;
        let count = (nu + 1)
            .checked_mul(nv + 1)
            .ok_or(PhillipsReadError::Limit)?;
        work.storage(count.saturating_mul(std::mem::size_of::<usize>()) + SCRATCH)?;
        work.charge(count)?;
        let mut result = Self {
            prefix: vec![0; count],
            divisions: grid.divisions,
        };
        work.storage(result.bytes() + SCRATCH)?;
        work.held += result.bytes();
        // Fixed26 slots cover the <=24 directed U subdivisions. Full stack
        // metadata/cells and the row constructor fit existing2MiB scratch.
        work.charge(2048 + ReducedNormal11::CELL_ROOT_WORK)?;
        let mut stack: [Option<(ReducedNormalCell11, u32)>; 26] = std::array::from_fn(|_| None);
        stack[0] = Some((normal.clone(), nu.trailing_zeros()));
        let mut top = 1;
        let mut u = 0;
        while top > 0 {
            work.charge(64)?;
            top -= 1;
            let (row, depth) = stack[top].take().expect("nonempty DFS prefix");
            if depth > 0 {
                work.charge(ReducedNormal11::CELL_SPLIT_WORK)?;
                let [left, right] = row.subdivide(0).map_err(tensor_error)?;
                stack[top] = Some((right, depth - 1));
                stack[top + 1] = Some((left, depth - 1));
                top += 2;
                continue;
            }
            work.charge(ReducedNormal11::MEAN_STRIP_WORK)?;
            // The strip constructor32KiB is inside the existing2MiB transient
            // reservation; the dynamic prefix is separately resident above.
            const {
                assert!(ReducedNormal11::MEAN_STRIP_BYTES < SCRATCH);
            }
            let ud = row.domain()[0];
            let strip = match row.mean_value_strip() {
                Ok(strip) => Some(strip),
                Err(crate::tapered_thread::Error::Unresolved) => None,
                Err(e) => return Err(tensor_error(e)),
            };
            for v in 0..nv {
                if strip.is_none() {
                    work.charge(16)?;
                    result.prefix[(u + 1) * (nv + 1) + v + 1] = 1;
                    continue;
                }
                // Combined executed-path audit: <=420 interval operations;
                // canonical lookups/midpoints/copies, BOTH facet checks,
                // validation, emission and loops fit the remaining860.
                // 512 point/MVT +256 C0 +512 canonical/facet/emission =1280.
                work.charge(ReducedNormal11::MEAN_POINT_WORK + 256 + 512)?;
                let vd = [v as f64 / nv as f64, (v + 1) as f64 / nv as f64];
                let domain = [ud, vd];
                let corners = lattice.corners(domain)?;
                let ids = [[u, v], [u + 1, v], [u + 1, v + 1], [u, v + 1]].map(|[u, v]| {
                    (((patch.offset[0] + u) % grid.shape[0]) * (grid.shape[1] + 1)
                        + patch.offset[1]
                        + v) as u32
                });
                let points: [P; 4] = std::array::from_fn(|j| {
                    *grid.vertices[ids[j] as usize]
                        .get_or_insert_with(|| corners[j].map(I::midpoint))
                });
                let literal = (0..4)
                    .map(|j| {
                        sum_norm(std::array::from_fn(|k| {
                            corners[j][k] - I::point(points[j][k])
                        }))
                        .hi
                    })
                    .fold(0., f64::max);
                let error = (I::point(root.interpolation_error(domain)) + I::point(literal)).hi;
                let mut easy = error <= quality.deflection;
                let mut triangles = [[0, 1, 2], [0, 2, 3]];
                if easy {
                    let bounds = match strip.as_ref().expect("validated numeric strip").bounds(vd) {
                        Ok(bounds) => Some(bounds),
                        Err(crate::tapered_thread::Error::Unresolved) => None,
                        Err(e) => return Err(tensor_error(e)),
                    };
                    let Some(bounds) = bounds else {
                        result.prefix[(u + 1) * (nv + 1) + v + 1] = 1;
                        continue;
                    };
                    let mut bounds = bounds.map(|[a, b]| I::new(a, b));
                    if !patch.forward {
                        bounds = bounds.map(|v| -v);
                    }
                    for tri in &mut triangles {
                        if !patch.forward {
                            tri.swap(1, 2);
                        }
                        easy &= normal_match(tri.map(|j| points[j]), bounds, quality.angle);
                    }
                }
                if easy {
                    for tri in triangles {
                        grid.indices.extend(tri.map(|j| ids[j]));
                        grid.faces.push(patch.face);
                    }
                    grid.error = grid.error.max(error);
                } else {
                    result.prefix[(u + 1) * (nv + 1) + v + 1] = 1;
                }
            }
            u += 1;
        }
        // Exact integer integral image: every undecided leaf contributes1.
        // A zero rectangle can omit only leaves already completely proved.
        for u in 1..=nu {
            for v in 1..=nv {
                work.charge(32)?;
                let i = u * (nv + 1) + v;
                result.prefix[i] += result.prefix[i - (nv + 1)] + result.prefix[i - 1]
                    - result.prefix[i - (nv + 2)];
            }
        }
        Ok(result)
    }
}
#[cfg(test)]
#[test]
fn actual_two_late_source_normal_violations_cannot_become_easy_mvt_facets() {
    let source = super::super::tests::actual_owned();
    let captures = [
        (
            2usize,
            [0, 3, 2],
            [
                [13829193282708167439, 4605821245853391634, 0],
                [13827324337644860486, 4607423297782739558, 0],
                [
                    13827972181160491941,
                    4607210246063216462,
                    4573055141623720871,
                ],
                [
                    13829708017205777490,
                    4605269260684184137,
                    4573055141623720871,
                ],
            ],
        ),
        (
            3usize,
            [0, 2, 1],
            [
                [4608619248475498913, 4603885847775034730, 0],
                [4608901185673223724, 4602365225963095902, 0],
                [
                    4608734419475705069,
                    4603379534665025238,
                    4573055141623720871,
                ],
                [
                    4608405601244196287,
                    4604700703773068354,
                    4573055141623720871,
                ],
            ],
        ),
    ];
    for (face, tri, bits) in captures {
        let field = Arc::new(
            ReducedNormal11::new(source.thread_normal(face, [0, 0]).unwrap().clone()).unwrap(),
        );
        assert!(!source.thread_faces()[face].forward());
        let mut row = field.root_cell();
        for _ in 0..3 {
            row = row.subdivide(0).unwrap()[0].clone();
        }
        let strip = row.mean_value_strip().unwrap();
        let bounds = strip
            .bounds([0., 0.03125])
            .unwrap()
            .map(|[a, b]| -I::new(a, b));
        let points = bits.map(|p| p.map(f64::from_bits));
        assert!(!normal_match(tri.map(|j| points[j]), bounds, 0.3));
    }
}
