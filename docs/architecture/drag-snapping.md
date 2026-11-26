# Dimension drag snapping

Status: implemented for primitive and panel corner drags, the primitive move gizmo and the Extrude distance arrow (CAD-401, CAD-402, CAD-403). Interactive macOS review pending.

## Policy

All dimension handles share `crates/multicad-ui/src/modeling/drag_snap.rs`.

- **Step.** The snap step is the floor-grid spacing the viewport displays at the current zoom (`viewport::grid_step`). It is the configured grid spacing coarsened by powers of ten, so lines stay at least 8 px apart and about 100 fit across the view. The floor grid and the drags use the same function, so the visible grid is the snap grid. Edge-on views hide the grid lines but keep the step of their zoom. The viewport reports the step to the modeling controller every frame (`ModelToolController::set_drag_grid`).
- **Bit-exact values.** A snapped value is `n × step`, formed from the integer `n` and the decimal digits of the step. A 10 mm grid gives exactly `90`, and a 0.1 mm grid gives the decimal `0.3`, never an accumulated `0.30000000000000004`.
- **Modifier.** Holding Shift during the drag releases the grid, as for Move and Rotate. It is read every frame, so it can be pressed or released mid-drag. A free drag still lands on 0.001 mm.
- **Typed values.** Values typed into corner inputs, the Inspector or the Extrude arrow field are never snapped.
- **Feedback.** Corner handles show the size readout, the active spacing and the modifier ("grid 10 mm · Shift for free"). The Extrude arrow label shows the same spacing text. The primitive Inspector names the spacing.

Sketch point snapping (CAD-136) is separate. It uses the configured sketch spacing and its own Settings toggles.

## Primitive and panel corners

- The side opposite the dragged corner never moves. When that fixed side lies exactly on a grid line of a world-aligned axis, the dragged side snaps to a grid line, so both the size and the position are grid multiples. Otherwise only the size snaps, for example on a rotated primitive.
- A size never drops below one step. Dragging through or across the opposite corner stops one step short of it.
- Drags are limited to 1,000,000 mm per size and coordinate.
- The exact Rust box stores extents, area and volume as bounded rationals. A free drag coarsens its resolution by decimal places until the exact volume fits. A grid drag that would not fit is refused.
- When the pointer gives no valid shape (off the drag plane, nonfinite, or refused), the last valid draft is kept and the readout says "outside the supported size range; keeping the last valid size".
- The typed API keeps its validation. A box with inadmissible or invalid extents is still refused with an explicit error.

## Extrude distance

Dragged distances are signed grid multiples along the extrusion direction. A drag never lands on zero: a value that would round to zero takes one step on the pointer's side. A typed `0` is still reported as invalid and never commits.

Existing geometry comes first. The viewport reports its current canonical scene (`ModelToolController::set_drag_scene`, cached by scene stamp). Every planar face perpendicular to the extrusion direction, such as the opposite face of a body, is a target at its signed offset from the start plane. A target within 8 px of the dragged value wins over the grid, and the nearest target wins. The arrow label then reads "face plane · Shift for free". Shift releases face targets as well as the grid. The start plane itself (offset zero) is never a target.

## Verification

- `modeling/drag_snap.rs` unit tests: bit-exact multiples, nonzero lengths and distances, the modifier, unusable grids, exact box admission, geometry targets before the grid.
- `viewport/grid.rs` unit test: the drag step equals the drawn line spacing at three zoom levels, plus edge-on views.
- `modeling/controller/drag_snap_tests.rs` drives real egui input through the production native worker:
  - every free drag state of CAD-401 (drag start, zero size, crossing the start point, far and extreme pointers, off-plane perspective) previews without an error;
  - a dragged box commits as one Undo step with the shown dimensions;
  - 10 mm and 100 mm zoom levels, Shift held and released, and a typed override for primitives and the Extrude distance;
  - an Extrude distance drag next to an existing box snaps to its 43 mm top face before the 10 mm grid, Shift releases it, and the commit stores exactly 43.
- `multicad-model/tests/primitive_box_admission.rs`: the typed API still refuses the inadmissible box sizes the UI used to send, and invalid boxes.

Not verified: interactive macOS feel and legibility, and Windows/Linux runtime behavior.
