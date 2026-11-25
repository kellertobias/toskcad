# Placed Hanger STEP: transport bound

Status: derivation for CAD-352 item 4, written before the implementation.

## Problem

A local Hanger export is qualified once in its own frame. The qualification yields:

- `e0`: a source-to-healed position bound;
- `e1 = [eu, ev]`: source support derivative bounds;
- `ec`: a coedge derivative bound;
- a healed volume interval.

A placed export writes the same local records inside one `REPRESENTATION_MAP` / `MAPPED_ITEM` with a `CARTESIAN_TRANSFORMATION_OPERATOR_3D`. Exact literal binding proves the mapped member's local records equal the qualified candidate. Only the operator remains to be proven.

The authored occurrence is `y = R x + t`, where R has the frame's u, v and u×v columns. A reader decodes the file operator as `y = A x + o`, where `A = s·[a1 a2 a3]`. A normalizes the written direction ratios and applies the uniform scale s, so A can differ from R by a few rounding units. The local proofs must be transported to world coordinates with a sound bound on that difference. No tolerance changes.

## Bounds

All quantities are outward-rounded intervals. Norms are upper bounds.

1. **Operator norm.** Let `G = AᵀA`. Gershgorin on G gives
   `‖A‖₂² ≤ max_i Σ_j |G_ij|`.
   For orthonormal-to-rounding axes this is `s²(1 + m)` with `m = max_i Σ_j |G_ij − s²δ_ij| / s²`.
   The implementation uses the first, unconditional form.

2. **Frame difference.** `D = A − R` uses decoded intervals minus exact frame values. It is bounded by its Frobenius norm `‖D‖_F` (≥ ‖D‖₂).

3. **Position.** Let H be the healed local surface, S the local source, and B the local source box. For any local point x of H, with `|x − S(x)| ≤ e0`:
   `|A x + o − (R S(x) + t)| ≤ ‖A‖₂ e0 + |o − t + D c| + ‖D‖_F r`.
   Here c is the centre of B and r is half its diagonal, so `|D(x − c)| ≤ ‖D‖_F r`.
   The world bound must stay within the unchanged policy tolerance, otherwise the file is refused.

4. **Derivatives.** For a chart direction ∂ of a local support, with `|∂H − ∂S| ≤ e`:
   `|A ∂H − R ∂S| ≤ ‖A‖₂ e + ‖D‖_F · sup|∂S|`.
   `sup|∂S|` is bounded from the written records themselves, which literal binding proves identical to the source-qualified candidate:
   - **Polynomial B-spline, degree p, knots u.** The derivative curve has control points `p (P_{i+1} − P_i) / (u_{i+p+1} − u_{i+1})`. By the convex-hull property, `sup|∂C| ≤ max_i` of their norms. Zero-length knot intervals carry no basis support and are skipped. For a tensor surface, `∂/∂u` applies this to every row j. The v basis is a partition of unity, so the maximum over all (i, j) bounds `sup|∂S/∂u|`, and symmetrically for v.
   - **Rational B-spline with positive weights.** Write `S = N / W`, where `N = Σ w P B` and `W = Σ w B`. Then `∂S = (∂N − S ∂W) / W`. The bounds are:
     - `W ≥ w_min = min w_i > 0` and `|S| ≤ max_i |P_i|`, by the convex hull of positive-weight control points;
     - `|∂N|` from the hull of `p (w_{i+1}P_{i+1} − w_i P_i) / Δ`, and `|∂W|` from the hull of `p (w_{i+1} − w_i) / Δ`.

     This gives `sup|∂S| ≤ (sup|∂N| + max|P| sup|∂W|) / w_min`. A non-positive or non-finite weight refuses.
   - **Plane.** The chart derivatives are the two placement axes. Their norms are bounded by evaluating the written directions with outward rounding, not by assuming 1.
   - **Coedge curves (3D B-splines).** These use the curve bound above. 2D pcurves are unaffected by the operator.

5. **Volume.** Healed volume transports by `det A`. The bound is `det A = s³ det[a1 a2 a3]`, an interval product of the decoded axes, and must be positive. A reflection or singular map refuses. The source volume uses the authored frame, whose determinant is 1 for a proper frame.

## Consequences

- The work is linear in the number of control points. It is metered and refuses on limit or cancellation.
- No tolerance is widened. In practice the operator difference is rounding-level, so its contribution is far below 1e-7 mm. The position check still enforces the unchanged policy tolerance, and derivative bounds are reported exactly as the local receipt reports them.
- A mapped member is refused for an unplaced source, and an unmapped member is refused for a placed source.
- No product or tolerance decision is needed.
