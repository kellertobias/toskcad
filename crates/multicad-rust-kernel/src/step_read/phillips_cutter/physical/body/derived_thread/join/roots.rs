//! Exact cell sign domains; isolators name algebraic roots, never trim away their interiors.
use super::*;
pub(super) struct Budget {
    pub(super) remaining_nodes: usize,
}
impl Budget {
    pub(super) fn new() -> Self {
        Self {
            remaining_nodes: 8192,
        }
    }
    fn node(&mut self, a: &mut Arithmetic<'_, '_>, work: usize) -> Result<()> {
        a.work.step(work)?;
        self.remaining_nodes = self
            .remaining_nodes
            .checked_sub(1)
            .ok_or(ReadError::Limit)?;
        Ok(())
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DyadicParameter {
    pub(super) numerator: u64,
    pub(super) depth: u32,
}
impl DyadicParameter {
    pub fn numerator(self) -> u64 {
        self.numerator
    }
    pub fn denominator_exponent(self) -> u32 {
        self.depth
    }
    pub(super) const ZERO: Self = Self {
        numerator: 0,
        depth: 0,
    };
    pub(super) const ONE: Self = Self {
        numerator: 1,
        depth: 0,
    };
    fn middle(a: Self, b: Self) -> Result<Self> {
        let depth = a.depth.max(b.depth);
        if depth >= 60 {
            return Err(ReadError::Limit);
        }
        Ok(Self {
            numerator: (a.numerator << (depth - a.depth)) + (b.numerator << (depth - b.depth)),
            depth: depth + 1,
        })
    }
    pub(super) fn cmp(self, other: Self) -> std::cmp::Ordering {
        let depth = self.depth.max(other.depth);
        (u128::from(self.numerator) << (depth - self.depth))
            .cmp(&(u128::from(other.numerator) << (depth - other.depth)))
    }
    fn narrow(lo: Self, hi: Self) -> bool {
        let depth = lo.depth.max(hi.depth).max(40);
        let lo = u128::from(lo.numerator) << (depth - lo.depth);
        let hi = u128::from(hi.numerator) << (depth - hi.depth);
        hi - lo <= 1u128 << (depth - 40)
    }
}
#[derive(Debug)]
pub(super) struct Root {
    pub(super) bracket: [DyadicParameter; 2],
}
#[derive(Debug)]
pub(super) struct Domain {
    pub(super) ends: [usize; 2],
    pub(super) sign: i8,
}
#[derive(Debug)]
pub(super) struct Certificate {
    pub(super) roots: Vec<Root>,
    pub(super) domains: Vec<Domain>,
}
fn variation(
    p: &[Poly],
    x: DyadicParameter,
    side: i8,
    a: &mut Arithmetic<'_, '_>,
) -> Result<usize> {
    let mut previous = 0;
    let mut count = 0;
    for polynomial in p {
        let mut derivative = polynomial.clone();
        let mut order = 0;
        let sign = loop {
            let value = a.evaluate(&derivative, x)?;
            if !value.is_zero() {
                let mut s = if value.is_positive() { 1 } else { -1 };
                if side < 0 && order % 2 == 1 {
                    s = -s;
                }
                break s;
            }
            // Every retained Sturm polynomial is nonzero. Repeated exact
            // differentiation therefore reaches a nonzero leading constant.
            derivative = a.derivative(&derivative)?;
            order += 1;
        };
        if previous != 0 && sign != previous {
            count += 1;
        }
        previous = sign;
    }
    Ok(count)
}
pub(super) fn certify_sturm(
    p: Poly,
    a: &mut Arithmetic<'_, '_>,
    budget: &mut Budget,
) -> Result<Certificate> {
    // Only the nonconstant Bernstein fallback enters here; it already owns
    // the finite original polynomial. No constant/zero Sturm chain is formed.
    a.work.step(1 + 19 * 128)?; // Bound/poll before the full coefficient clone.
    let derivative = a.derivative(&p)?;
    let gcd = a.gcd(p.clone(), derivative)?;
    let squarefree = a.divide(p.clone(), &gcd)?;
    let sturm = a.sturm(squarefree.clone())?;
    let mut roots = Vec::new();
    if a.evaluate(&squarefree, DyadicParameter::ZERO)?.is_zero() {
        roots.push(Root {
            bracket: [DyadicParameter::ZERO; 2],
        });
    }
    let count = variation(&sturm, DyadicParameter::ZERO, 1, a)?
        - variation(&sturm, DyadicParameter::ONE, -1, a)?;
    let mut pending = vec![(DyadicParameter::ZERO, DyadicParameter::ONE, count)];
    while let Some((lo, hi, count)) = pending.pop() {
        budget.node(a, 1)?;
        if count == 0 {
            continue;
        }
        if count == 1
            && DyadicParameter::narrow(lo, hi)
            && !a.evaluate(&squarefree, lo)?.is_zero()
            && !a.evaluate(&squarefree, hi)?.is_zero()
        {
            roots.push(Root { bracket: [lo, hi] });
            continue;
        }
        let middle = DyadicParameter::middle(lo, hi)?;
        let at = a.evaluate(&squarefree, middle)?.is_zero();
        let left = variation(&sturm, lo, 1, a)? - variation(&sturm, middle, -1, a)?;
        let right = variation(&sturm, middle, 1, a)? - variation(&sturm, hi, -1, a)?;
        // Exact one-sided Sturm counts partition the parent into these two
        // open children plus the retained midpoint root, without rounding.
        pending.push((middle, hi, right));
        if at {
            roots.push(Root {
                bracket: [middle; 2],
            });
        }
        pending.push((lo, middle, left));
    }
    if a.evaluate(&squarefree, DyadicParameter::ONE)?.is_zero() {
        roots.push(Root {
            bracket: [DyadicParameter::ONE; 2],
        });
    }
    finish(p, roots, a)
}

fn finish(p: Poly, mut roots: Vec<Root>, a: &mut Arithmetic<'_, '_>) -> Result<Certificate> {
    a.work.step(1 + roots.len() * roots.len())?;
    roots.sort_by(|a, b| a.bracket[0].cmp(b.bracket[0]));
    let mut domains = Vec::new();
    // Boundary0/last are literal endpoints; middle IDs own unique algebraic
    // roots of this ORIGINAL coefficient polynomial, not rounded parameters.
    let bounds = std::iter::once([DyadicParameter::ZERO; 2])
        .chain(roots.iter().map(|r| r.bracket))
        .chain(std::iter::once([DyadicParameter::ONE; 2]))
        .collect::<Vec<_>>();
    for (i, pair) in bounds.windows(2).enumerate() {
        let left = pair[0][1];
        let right = pair[1][0];
        // Completed exact bisection supplies ordered disjoint isolators;
        // accepted closed intervals have nonzero endpoints. Only duplicated
        // literal endpoint contacts can make a zero-length boundary pair.
        if left.cmp(right).is_eq() && a.evaluate(&p, left)?.is_zero() {
            continue;
        }
        let witness = if left.cmp(right).is_eq() {
            left
        } else {
            DyadicParameter::middle(left, right)?
        };
        let value = a.evaluate(&p, witness)?;
        // This witness lies between consecutive exhaustive algebraic roots,
        // or at a nonroot isolator endpoint, so its exact sign is nonzero.
        domains.push(Domain {
            ends: [i, i + 1],
            sign: if value.is_positive() { 1 } else { -1 },
        });
    }
    Ok(Certificate { roots, domains })
}
// Exact Bernstein-Descartes subdivision: zero variations excludes every
// interior root; one variation proves one simple interior root. Exact midpoint
// and endpoint zeros are retained separately, including their multiplicity.
// Non-dyadic multiple contacts use the exact square-free Sturm fallback;
// resource exhaustion refuses the request, never discards a contact.
pub(super) fn certify(
    p: Poly,
    a: &mut Arithmetic<'_, '_>,
    budget: &mut Budget,
) -> Result<Certificate> {
    if p.len() == 1 {
        budget.node(a, 1)?;
        // An identically zero comparison has infinitely many contacts. This
        // finite-event receipt refuses it explicitly; it never invents a
        // sampled representative or drops a coincident source interval.
        if p[0].is_zero() {
            return Err(ReadError::Unsupported);
        }
        return Ok(Certificate {
            roots: Vec::new(),
            domains: vec![Domain {
                ends: [0, 1],
                sign: if p[0].is_positive() { 1 } else { -1 },
            }],
        });
    }
    let coefficients = a.bernstein(&p)?;
    certify_coefficients(p, coefficients, a, budget)
}
/// Distinct fixed degree22 entry. Its producer accepts exactly12 actual poles;
/// it cannot silently widen the old degree18 scratch/work contract.
pub(super) fn certify22(
    p: exact::Polynomial22,
    a: &mut Arithmetic<'_, '_>,
    budget: &mut Budget,
) -> Result<Certificate> {
    // The common Sturm fallback already charges19*128 for source clones. Add
    // the extra four coefficients here BEFORE any possible degree22 fallback.
    a.work.step(4 * 128)?;
    if p.degree() == 0 {
        return certify(p.into_polynomial(), a, budget);
    }
    let coefficients = a.bernstein22(&p)?;
    certify_coefficients(p.into_polynomial(), coefficients, a, budget)
}
fn certify_coefficients(
    p: Poly,
    coefficients: Poly,
    a: &mut Arithmetic<'_, '_>,
    budget: &mut Budget,
) -> Result<Certificate> {
    let mut roots = Vec::new();
    if coefficients[0].is_zero() {
        roots.push(Root {
            bracket: [DyadicParameter::ZERO; 2],
        });
    }
    if coefficients.last().unwrap().is_zero() {
        roots.push(Root {
            bracket: [DyadicParameter::ONE; 2],
        });
    }
    let mut pending = vec![(DyadicParameter::ZERO, DyadicParameter::ONE, coefficients)];
    while let Some((lo, hi, coefficients)) = pending.pop() {
        budget.node(a, coefficients.len())?;
        let mut previous = 0;
        let mut variations = 0;
        for c in &coefficients {
            let sign = if c.is_zero() {
                0
            } else if c.is_positive() {
                1
            } else {
                -1
            };
            if sign != 0 {
                if previous != 0 && sign != previous {
                    variations += 1;
                }
                previous = sign;
            }
        }
        if variations == 0 {
            continue;
        }
        if variations == 1
            && DyadicParameter::narrow(lo, hi)
            && !coefficients[0].is_zero()
            && !coefficients.last().unwrap().is_zero()
        {
            roots.push(Root { bracket: [lo, hi] });
            continue;
        }
        if variations > 1 && lo.depth.max(hi.depth) >= 8 {
            // Restart only this exact polynomial with its square-free Sturm
            // chain, preserving all spent work/node charges. Multiple roots
            // are retained, not inferred absent from a failed sign split.
            drop(pending);
            drop(coefficients);
            drop(roots);
            return certify_sturm(p, a, budget);
        }
        let middle = DyadicParameter::middle(lo, hi)?;
        let [left, right] = a.split_bernstein(coefficients)?;
        if left.last().unwrap().is_zero() {
            roots.push(Root {
                bracket: [middle; 2],
            });
        }
        pending.push((middle, hi, right));
        pending.push((lo, middle, left));
    }
    finish(p, roots, a)
}
