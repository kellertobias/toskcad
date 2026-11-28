# Frame drawing: profile and tube frames drawn as lines (CAD-435 design)

Status: **design accepted for planning, 2026-11-27.** The product owner
answered the eleven CAD-435 questions (§10). Nothing in this document is
implemented yet except the CAD-434 stock picker work it builds on. Slice 1
(§9) is re-planned for the answers: it is a full 3D slice.

Related: [weldments.md](weldments.md) (what exists today),
[walls.md](walls.md) and CAD-427 (the wall drawing gesture this follows),
[shared-materials.md](shared-materials.md) (stock is authored in Project),
CAD-434 (the stock picker, which already follows the direction set here).

## 1. Goal

A user draws a metal frame the way they draw a wall: click, type a length,
turn, type the next length, close the loop. Every line becomes a member of a
chosen profile. They say where the profile sits on the line, and the joints
at the nodes appear on their own: mitres where two equal tubes meet at a
corner, a fish-mouth where a round tube meets another, a butt joint where a
small tube meets a bigger one. Clicking a node badge changes its joint. The
cut list, with lengths, both end angles and cope templates, is always
current.

What no other CAD does in one tool (see §8):

- **No skeleton step.** Frame generators make you draw a sketch or 3D
  sketch first and then "insert members" on it. Here the drawing *is* the
  frame.
- **Joints are node objects, decided automatically, visible and
  overridable in place.** Elsewhere the corner treatment is a feature
  option or an extra feature per pair of members.
- **The cut list is a live, exact output** (overall length, both cut angles,
  long-point side, cope pattern), not a table you export afterwards.
- **One canonical document for the window and MCP**: an agent can draw,
  inspect and change the same frame through the same typed transactions.

## 2. What exists today and how this relates to it

| Existing piece | Location | Decision |
| --- | --- | --- |
| Profile stock (`StockKind::Profile`, shapes rectangular, round, angle, item Line 8 40×40, sketch) | `crates/multicad-core/src/materials/source.rs` | **Reuse unchanged.** Stock stays authored in Project → Materials (CAD-224). Frame members select it. CAD-434 adds a corner radius and the visual picker. |
| `ProfileSection`, `ResolvedMember`, member frame, anchor, rotation, extensions | `crates/multicad-core/src/model/weldment/resolve.rs` | **Reuse.** A frame member resolves to exactly the same `ResolvedMember`, so the native solid, provenance roles and measurement code do not change. |
| `weldment.member` feature (one output Part per member), `MemberPath::{Points, SketchLine}` | `crates/multicad-core/src/model/weldment/mod.rs` | **Extend.** Add `MemberPath::FrameMember { frame, member }`. Frame members stay `weldment.member` features with their own Part, so BOM, drawings, materials and part properties keep working. |
| `weldment.connection` (pairwise Miter, Butt, Cope with gap) | `crates/multicad-core/src/model/weldment/connection.rs` | **Reuse the geometry, not the entity.** Node joints in a frame are stored on the frame node. The evaluator calls the same miter, butt and cope cut builders per member pair. Standalone connections remain for existing members. |
| Frame transaction: every solid line of a sketch becomes a member | `weldment` frame transaction | **Keep, and offer "Convert to frame"** later: an explicit, undoable command that makes a frame graph from a sketch-line frame and its connections. No automatic migration of files. |
| Member cut list, `group_cuts`, weld schedule, weld plans | `crates/multicad-model/src/weldment.rs`, `weld_schedule.rs` | **Reuse and extend** with cope patterns, bar nesting and mass (§6). |
| item Line 8 connectors (angle bracket, joining plate, universal end fasteners) on groove stations | `crates/multicad-core/src/model/profile_connections.rs` and `profile_connections/frames.rs` | **Reuse as the fastener joint for slot profiles.** A node joint of kind `SlotFastener` derives the two `ProfileAttachment` stations from the node geometry and creates the existing `model.profile_connection` recipe. The station, footprint and access-bore logic stays where it is. |
| Linked diagrams (CAD-395) | `crates/multicad-ui/src/modeling/diagram.rs` | **Reuse** for the anchor picker (§4.3), the joint editor (§5) and the stock picker (CAD-434). |
| Plan outline editor (walls, fences, slabs) | `crates/multicad-ui/src/modeling/outline.rs` | **Reuse the interaction, not the type.** Chaining, Backspace, Enter, closing on the first corner, live segment dimensions and checked clicks are the same gesture. Frames need a graph (branches, 3D), so the frame editor is a sibling that shares the snapping and value entry helpers. |

Why **extend** rather than migrate: the member solid, its end cuts, the cut
list and the weld schedule already compute exact results and are tested.
A frame is a better *source* for members and joints; it does not need a new
evaluator. Keeping `weldment.member` as the Part-producing feature means
every downstream consumer (BOM rows, drawing references, weld marks, MCP
tools) keeps working with frame-drawn members on day one.

Why **not reuse `weldment.connection` entities** for frame joints: a node is
one decision about several members (T, X, three-way corners). Storing it as
N pairwise connection entities would duplicate state that the node already
holds, and the two could disagree. The single canonical document rule says
the node is the source; the pairwise cuts are derived.

## 3. Model

### 3.1 Entities

```text
weldment.frame  (extension entity, the canonical source)
├── plane            base plane: sketch plane, work plane, level or face frame
├── nodes[]          stable ids, each with a position definition
├── members[]        stable ids, start node, end node, placement
└── joints[]         per-node overrides only (auto joints are derived)

weldment.member      (one per frame member, output Part)
└── path = FrameMember { frame, member }
```

The frame entity holds intent only. Each member feature evaluates from the
frame plus its stock, like today's sketch-line members evaluate from the
sketch. Member feature and Part identities are created once, when a member
is drawn, and stay through every later edit. Deleting a member deletes its
feature and Part (CAD-426 cascade rules apply).

### 3.2 Nodes

A node is a point in the frame's plane coordinates `(u, v, w)`, `w` along
the plane normal (0 for planar frames). Each node has one definition:

| Definition | Meaning | Typical source |
| --- | --- | --- |
| `Absolute { u, v, w }` | coordinates as `Dimension`s (literals or `=$param`) | first click, typed coordinates |
| `Polar { from, length, angle, rise }` | from another node: length along an in-plane angle, optional rise along the normal | chained drawing with typed length and angle |
| `Along { member, distance, from_end }` | a station on an existing member | snapping onto a member to start a T branch |
| `Bound { reference }` *(later)* | a model vertex or a point of another frame | snapping to model geometry |

`Polar` keeps what the user typed: change `$rail` and every node after it
moves. Dragging a node rewrites its own definition (a literal `Polar`
length or absolute coordinate) and never the definitions it depends on;
nodes that depend on it follow. Dependency cycles are refused at the
transaction. This is deliberately not a general constraint solver; see
question Q1.

### 3.3 Members

```text
Member {
    id, start: NodeId, end: NodeId,
    stock: EntityId,              // Project profile stock
    placement: Placement,
    extension: [Dimension; 2],    // start / end, negative trims (as today)
}
Placement {
    anchor: SectionAnchor,        // which point of the section is on the line
    offset: [Dimension; 2],       // section x / y shift after anchoring
    roll: Dimension,              // degrees about the member axis
    up: UpRule,                   // section height reference
}
SectionAnchor = Grid([AxisAnchor; 2])      // today's 3×3 bounding-box points
              | Centroid                   // area centroid (angles, custom)
              | Side(FrameSide)            // inner / outer / left / right edge
```

- **Grid** is today's anchor and stays the stored form for the nine
  bounding-box points.
- **Side** is resolved per member at evaluation: for a closed planar loop,
  *inner* is the side toward the loop interior; for an open chain, *left* and
  *right* look along the drawing direction (the wall rule). The section face
  on that side lies on the line. This is what makes "draw the outside
  dimensions of a table top, profiles inside" one choice instead of
  per-member anchors. Reversing a member keeps its side.
- **UpRule** gives the section height direction: the plane normal for
  in-plane members (today's sketch-line rule); for off-plane members the
  in-plane direction toward the frame centroid, overridable. Roll turns the
  section from there.
- **Defaults** come from the frame: the draw tool sets a frame default
  placement and stock; a member stores only what differs (serde
  `skip_serializing_if` default), so changing the frame default changes all
  members that never overrode it.

### 3.4 Joints

A node with two or more member ends or stations has a joint. The evaluator
derives an **automatic joint** from the node geometry and sections (§5.4).
An **override** is stored only when the user changes it:

```text
Joint { node, kind: JointKind, through: Option<MemberId>,
        gap: Dimension, prep: WeldPrep, fasteners: Option<SlotFastener> }
```

The derived joint is not stored, so it follows every edit; the override
says "auto" until the user picks something. An override that no longer
applies (a member was deleted, an angle became impossible) is reported on
the node, drawn red, and never silently replaced.

### 3.5 Capability and files

Files with frames require a new capability, `weldment.frame.1`. Existing
weldment files are unchanged. Stock gains only what CAD-434 adds (corner
radius with its own capability).

## 4. Interaction

The tool is **Model → Metal → Frame** (Q7: a new Model ribbon group
"Metal" with a "Frame" tool; Project stays ribbon-free).
It opens a form on the right, like walls: profile, placement, plane, then
the node and member lists. The viewport is where most of the work happens.

### 4.1 Drawing

```text
  click ──► node 1 ──► move: rubber-band member with live length + angle
                     ──► type "600" Tab "90" Enter  (or click)
                     ──► node 2 … repeat
  click node 1 again ──► loop closed, tool keeps chaining from a new click
  Enter ──► finish chain     Esc ──► cancel the pending segment
  Backspace ──► remove the last node
```

- **Chaining.** Each click ends a member and starts the next one at the
  same node. The pending member shows its real section (not a line), placed
  with the current anchor, so the user sees which side the tube goes.
- **Typed entry.** While moving, typing starts a length field; Tab moves
  to the angle field. Units follow the shared value entry (`mm`, `cm`, `m`,
  `in`, `ft`, expressions, `=$param`). A typed value locks that quantity;
  the mouse only chooses the rest. The values are stored as typed
  (`Polar` definition).
- **Angles.** Relative to the previous member by default (so "90" means a
  right-angle turn); holding Shift or typing `@` switches to absolute plane
  angle. Snaps at 15° steps and to parallel / perpendicular with existing
  members, with the same guide lines as sketch snapping (CAD-422).
- **Closing loops.** Clicking the first node of the chain, or pressing `C`,
  closes the loop. A closed loop gets *inner/outer* side meaning (§3.3).
- **Snapping.**
  - to nodes (joins the members at that node);
  - to a member's line or a station on it (creates an `Along` node: the
    member keeps running through, a T joint is derived);
  - to member ends extended (inference lines), mid points and to model
    geometry (vertices, edges, faces) using the shared drag-snapping
    (CAD-401/403) infrastructure.
- **Branching.** Starting a chain on an existing node or member continues
  the same frame. Starting on empty space starts a new chain in the same
  frame. A frame is a graph, not a polyline.

### 4.2 3D frames

The frame plane is chosen at start: a sketch plane, a work plane, a level,
or a clicked planar face (as Add Sketch does).

- **Off-plane members.** From a node, pressing `Z` (or the "Up" button)
  locks the pending member to the plane normal; typing a length gives a
  vertical leg. `Polar.rise` stores it.
- **Other planes.** Picking two existing members that share a node sets a
  temporary drawing plane through them, so the side of a cube can be drawn
  after its base. Nodes always live in the frame coordinates.
- **Copy up.** A selected closed loop can be copied along the normal by a
  typed distance, with connecting legs at every node: the common "rack /
  table" frame in two steps. *(Slice 3.)*
- Clicks never project silently: a click that hits no valid plane is
  refused with the reason, like the wall editor's checked clicks.

### 4.3 Placement

The form shows a small section diagram (the CAD-434 linked diagram) with
the nine grid anchors, centroid and the four side choices as clickable
points on the drawn section, plus offset and roll fields linked to their
dimensions. The viewport previews the placed section on the pending and
selected members. Changing the frame default re-places every member that
did not override it; a member's own placement is in its row.

### 4.4 Editing an existing frame

Opening a frame (from the tree or by double-clicking a member) reopens the
same tool:

- **Drag a node**: the node moves with live lengths of every attached
  member; snapping applies; the joints recompute live.
- **Click a member**: select it; its row shows stock, placement, length
  (read-only when the length comes from its nodes); Delete removes it.
- **Click a node badge**: opens the joint editor (§5).
- **Split / insert node** on a member, **merge** two collinear members,
  **reverse** a member.
- All edits are drafts until Apply; Discard Edit, Undo and Redo work as for
  walls. MCP has the same edits as typed transactions.

## 5. Connections at nodes

### 5.1 Joint kinds

| Kind | Members | What is cut | Exists today |
| --- | --- | --- | --- |
| **Mitre** | 2 ends | both on the bisector plane | yes, mirrored sections only |
| **Butt, A through** / **B through** | 2 ends or end + station | the terminating member square (or angled, later) to the through member's face; the through member's end covers it | square, rectangular through only |
| **Overlap** | 2 ends | none; members run past each other by their section size, side by side (flats, angles, bolted frames) | no |
| **Cope / saddle** | branch end + through | branch fitted to the through member's actual outside (fish-mouth for round) | yes |
| **Notch** | angle or channel end | one leg removed over the partner's width so the other leg runs on (bent-angle corner) | no |
| **Tab and slot** | 2 | a tab on the branch end, a matching slot through the receiving wall (tube-laser style) | no |
| **Gusset** | any, added to a joint | a plate Part from sheet stock, triangle with legs along both members, optional chamfer | no |
| **End cap** | 1 end | a plate Part closing an open tube end | no |
| **Slot fastener** | slot profiles | butt cut plus an item Line 8 recipe (universal end fastener, angle bracket, joining plate) | recipes yes, auto placement no |

**T joint:** one member runs through the node (an `Along` station or two
collinear members merged), the branch is butted or coped onto it.
**X joint:** one member runs through, the two halves of the other are each
butted or coped on it (they stay two members, two Parts); for slot profiles
an X is usually two fasteners.

**Three or more members at a corner** (cube corner, table leg under two
rails): the joint picks an order: one member through (the leg, by default),
the others cut against it and against each other (two mitres plus a butt
onto the leg top, or two butts). The editor lists the order explicitly.

**Leg corners** (one leg along the plane normal and two members in the
plane, nothing else) are built in one of three ways, stored as the node's
joint override `corner` (product owner, 2026-11-27: "ask: show a dropdown
at the edge location and let the user cycle with arrow keys and then with
enter to select. while an entry is selected, preview the selected entry"):

| `corner` | Construction | Table of SHS 40×40×2, inner rectangle 1000 × 600 |
| --- | --- | --- |
| `leg_through` (automatic) | The leg runs through to the rail tops; both rails butt onto its faces. | legs 720 mm, rails 920 / 520 mm, outside 1000 × 600 |
| `mitre_on_leg` | The rails mitre; the leg stops under them, butted square to a rail's underside. | rails 1000 / 600 mm mitred, legs 680 mm |
| `mitre_around_leg` | The leg is centred on the node and runs to the rail tops; the rails mitre and are also cut to the leg faces. | legs 720 mm standing 20 mm outside, rails cut back to u = max(v, 20 mm) |

A corner construction excludes a joint kind and a through member at the
same node; one stored at a node that is no longer a leg corner is reported,
never dropped. Two plane cuts on one member end are clipped in turn. In the
Frame tool, clicking a leg corner opens a dropdown at the corner: ↑ and ↓
cycle the entries and the highlighted one previews in the viewport, Enter
applies it as one transaction, Esc restores the previous joint. Choosing a
joint kind or through member in the joint editor replaces the corner
construction. MCP `edit_frame` takes the same `corner` in `set_joint`.

### 5.2 Gap and weld preparation

- **Gap** (Q9: "make configurable weld gap, defaults to 0") is a value on
  each joint, with a frame-level default of 0 mm. A non-zero gap moves each
  cut back by the gap (half per side for a mitre), as weldment connections
  do today, so the cut list shows real saw lengths (confirmed by the
  product owner: "Shortens members").
- **Preparation**: none, bevel (angle, land) or V on the cut edge. Modelled
  as an exact chamfer on the cut faces in a later slice; until then it is
  carried to the weld schedule and the cutting sheet as an annotation.
- A joint can create its **weld** record (the existing `weldment.weld`)
  with a default symbol from the joint kind: fillet for butt and cope,
  V-groove for mitre with V prep.

### 5.3 Fasteners for slot profiles

For item Line 8 (and later other slot systems) the joint uses
`profile_connections.rs`: the node gives each member its groove side and
station (end, or `Along` distance); the recipe validates the footprint
(`validate_station`), and the existing `model.profile_connection` creates
the connector, screw and nut Parts with their exact mounting frames and
access bores. No welding records are made for slot profiles.

### 5.4 Automatic defaults

Evaluated per node from current geometry; shown as a badge in the viewport
(`M` mitre, `B→` butt with the arrow on the through member, `C` cope,
`F` fastener, red `!` refused). Rules, in order:

1. Members meeting below 15°, or not meeting (paths pass each other):
   **refused**, shown red with the reason.
2. Collinear pair at a node with no other member: **one member through**
   (Q5); the node becomes a station.
3. Slot profiles: **butt, the longer member through, plus universal end
   fastener** if the footprint fits, else angle bracket (Q10).
4. Two round tubes, or a round branch onto anything: **cope**, larger
   diameter through; equal diameters at an L corner: mitre.
5. Two ends with mirroring sections at a corner (same stock, mirrored
   placement): **mitre** (Q4).
6. Two ends with different sections: **butt, larger section through**.
7. End on a station (T): **butt** for flat faces, **cope** otherwise.
8. Three or more: the member along the plane normal (a leg) runs through;
   in-plane members mitre or butt with each other by rules 4–6.

The user changes any node in the joint editor: a linked diagram of the two
(or more) sections at the node, with the kind as picture tiles, the through
member as a toggle, and gap and prep fields linked to the diagram.

### 5.5 Joints on channel, T and I/H sections

CAD-434 adds parametric channel, T and I/H stock (Q11). Until a later
slice adds their joint rules, members of these sections can be placed and
mitred when mirrored, and they butt square onto a flat face only as the
terminating member. Copes against them and butt joints with them running
through are refused with the reason, for example "Copes against a channel
are not supported yet; use a butt joint onto a flat face".

## 6. Outputs

### 6.1 Cut list (exact, live)

Built on the existing `member_cut_list` and `group_cuts`, which already
measure the final solid and give overall length, both end treatments with
cut angle and long-point face, quantity per group, and whether a member
fits one stock bar. Additions:

- **Cope pattern** per coped end: the intersection curve unrolled around
  the branch tube, sampled exactly from the cut solid, as a 1:1 printable
  wrap template (PDF) and as a table of cut depth per angle. This is what a
  workshop without a tube laser actually needs.
- **Bar nesting**: a 1D cutting-stock plan per stock with saw kerf
  (default 3 mm, a Project stock setting), so the list says "3 bars of
  6 m: bar 1 = A×2 + B, offcut 412 mm".
- **Mass and material**: section area (exact, including CAD-434 corner
  radii) × length × Material density per member and per group, totals per
  stock and per frame.

### 6.2 Cutting sheet

A drawing page generated from the cut list: one row per cut group with a
small dimensioned elevation of the member showing both end angles, the
long-point side, the length, quantity and mark; cope templates on following
pages. It uses the existing drawing tables and annotations and follows
geometry edits like other generated drawing content.

### 6.3 Tube laser export (feasibility)

- **Per-member STEP**: feasible now; tube laser CAM (for example
  TruTops Tube, Lantek Flex3d Tubes, BLM ArtTube) imports a single straight
  member solid and derives its own cut paths. The export puts each member
  in its member frame (axis along +X), named with its cut-group mark.
- **DSTV NC1** (steel structures standard): feasible for rectangular,
  round and angle sections with plane cuts and holes; copes are not
  representable in plain NC1 and need the STEP route. Later slice.
- Machine-specific formats are out of scope.

## 7. Diagrams

Placement choices on a rectangular tube, looking along the member:

```text
  Grid anchors (today)          Side anchors (new, per loop)
  ┌───────────────┐                outer
  │●      ●      ●│             ┌─────────┐
  │               │       ──────┤ line on ├────── outer edge on line
  │●      ●      ●│             │  outer  │
  │               │             └─────────┘
  │●      ●      ●│                inner
  └───────────────┘
```

A rectangle drawn by its outside dimensions with profiles inside:

```text
   node 1 ●━━━━━━━━━━━━━━━━━━━━━━━━━● node 2      ━ drawn lines
          ┃▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓┃           ▓ profile, side = inner
          ┃▓                       ▓┃           M mitre badges at nodes
          ┃▓                       ▓┃
   node 4 ●━━━━━━━━━━━━━━━━━━━━━━━━━● node 3
```

## 8. Prior art

| Tool | How frames are made | Joints | Cut list | Weak points this design fixes |
| --- | --- | --- | --- | --- |
| SOLIDWORKS Weldments | Structural Member feature on 2D/3D sketch segments; profile library sketches; pierce point and rotation per group | Corner treatment per group (end miter, end butt 1/2), Trim/Extend feature, gussets and end caps as features | Cut list folder with length and angle properties | Skeleton sketch first; corner treatment is per group, not per node; T joints need a separate Trim/Extend; one multibody part |
| Autodesk Inventor Frame Generator | Insert members on skeleton edges or points; Content Center profiles; offsets and rotation | Miter, Notch, Trim/Extend, Lengthen/Shorten as separate commands | Frame member BOM; cut lengths | Skeleton part required; many files; edits re-run commands |
| Onshape Frames | Frame feature on sketch/3D paths; standard content profiles; trim per corner override | Per-vertex corner overrides (miter, butt, coped); Trim feature; Gusset and End cap features | Cut list table with angles | Still skeleton first; no drawing gesture; no cope templates |
| Fusion | No native frame tool (add-ins, sweeps) | manual | manual | — |
| FreeCAD | Frame workbench / BIM structure on sketch edges | limited, add-on dependent | add-on | fragmented |
| Tekla Structures | Steel members placed by points; connection components with rules | Rich automatic connections (structural steel) | NC/DSTV, drawings | Building-steel scale and cost, not for workshop frames |
| SketchUp + plugins | draw lines, extrude profiles | none or plugin | plugin | not parametric |

How this is clearly better:

1. **Draw, do not prepare.** The wall-style gesture with typed length and
   angle creates nodes, members and joints in one pass.
2. **Per-node joints with automatic, visible defaults**, overridable with a
   click on the badge; Tekla-like intelligence at workshop scale.
3. **Side alignment relative to a loop** (inner/outer), so outside
   dimensions are what you draw.
4. **Exact, live cut list with cope wrap templates and bar nesting**, from
   the final solids, not from nominal lengths.
5. **Slot profiles and welded tubes in one tool**: the same frame can be
   item Line 8 with real fasteners or welded steel.
6. **Separate Parts per member** with stable identity, so drawings, BOM
   and MCP keep their references through edits.

## 9. Phasing

Slice 1 is a **full 3D slice** (Q3), tracked as CAD-439 (model, joints,
MCP), CAD-440 (planar tool), CAD-441 (3D, T/X, multi-member nodes) and
CAD-442 (slot-profile fasteners). Later slices: CAD-443 (outputs),
CAD-444 (more joints), CAD-445 (section tables). It is ordered so a planar frame works
first, but it is only done when off-plane members, temporary planes,
multi-member nodes and T/X joints work too.

### Slice 1: 3D frames with automatic joints

Work order inside the slice (each step is usable on its own):

1. **Core model and transactions.** `weldment.frame` with nodes
   (`Absolute`, `Polar` with `rise`, `Along`), members (stock, placement
   with grid and side anchors, offset, roll, up rule, extensions) and joint
   overrides (kind, through member, gap with frame default 0, order);
   `MemberPath::FrameMember`; capability `weldment.frame.1`; typed
   transactions shared by UI and MCP (create frame, add chain, move node,
   split/merge member, set placement, set joint, delete member/node).
2. **Planar evaluation.** Two-member nodes derive cuts with the existing
   miter, butt and cope builders; default rules 1, 4, 5, 6 (§5.4); per-node
   refusals in plain language.
3. **Frame tool, planar.** Model → Metal → Frame: plane pick, chained
   drawing with live section preview, typed length and relative angle with
   units, closing loops, snapping to nodes and angle steps, Backspace,
   Enter, Esc, Apply, Discard; node dragging; joint badges and the joint
   editor (kind, through member, gap). Placement: frame default stock and
   placement, per-member override, inner/outer for loops (Q6) and
   left/right for open chains.
4. **T and X joints.** `Along` nodes by snapping onto a member; collinear
   chains run through as one member (Q5); split and merge; X joints as one
   through member and two butted or coped halves.
5. **3D.** Off-plane members (`Z` lock and typed rise), temporary drawing
   planes through two members, off-plane up rule, snapping in 3D.
6. **Multi-member nodes.** Three or more members at a node with an explicit
   order (leg through, rule 8), angled butt joints and butt against round
   members.
7. **Outputs in the tool.** Live cut list with both end angles, long-point
   faces, groups, mass per member and total.
8. **Slot profiles.** Automatic, overridable fasteners at slot-profile
   nodes (Q10) through `profile_connections.rs`.

Acceptance criteria for slice 1:

- **Planar rectangle.** A user draws a closed 4-node rectangle
  1000 × 600 mm of SHS 40×40×2 with side = inner by typing lengths and
  angles, and gets four member Parts whose outside is exactly
  1000 × 600 mm and four mitres (Q4); `member_cut_list` reports
  2 × 1000 mm and 2 × 600 mm, 45° at both ends, long points on the outer
  face.
- **Parametric.** Changing a `Polar` length to `=$width` and editing
  `$width` resizes the frame; member feature and Part ids are unchanged.
- **Editing.** Dragging a node updates every attached member and joint
  live; Discard restores the source; Undo/Redo restore exactly.
- **Defaults and overrides.** A round branch onto a larger round tube is
  coped; a 30×30 end on a 40×40 is butted with the 40×40 through;
  overriding to "A through" is stored and survives edits; deleting a member
  reports a stale override instead of replacing it.
- **Gap.** A frame gap default of 0 builds members that meet exactly; a
  2 mm gap on one mitre shortens both members by their share and the cut
  list shows the saw lengths.
- **T and X.** Starting a chain on a member's midpoint makes a T: the
  member runs through and the branch is butted or coped. Drawing straight
  through a node keeps one member. Crossing a member at a node makes an X
  with one through member and two halves.
- **3D table.** From the closed rectangle, four 700 mm legs drawn along
  the plane normal give a table frame whose corner nodes each hold two
  mitred rails and a leg running through to the top; all eight members are
  separate Parts and the cut list lists 4 legs with square ends.
- **Refusals.** Angles below 15°, members that do not meet and joints a
  section cannot take (for example a cope against a channel, §5.5) are
  refused at the node with the reason, never built.
- **Slot profiles.** An item Line 8 rectangle gets universal end fasteners
  at its corners by default; removing one at a node is stored.
- **MCP parity.** MCP builds the same frames with the same transactions and
  returns the same cut list.
- **Files.** Frames save, reopen and round-trip in files requiring
  `weldment.frame.1`; files without frames are unchanged.
- **Tests.** UI tests drive the tool by real input (clicks, typed entry)
  and wait on readiness conditions; model tests measure the native solids.

### Slice 2: fabrication outputs (Q8)

Cope wrap templates (1:1 PDF), a generated cutting-sheet drawing page and
bar nesting with kerf; then per-member STEP for tube lasers; DSTV NC1
later.

### Slice 3: more joints and parts

Overlap, notch, tab and slot, gussets and end caps as Parts from Project
sheet stock, weld preparation as exact chamfers, joint-created weld
records, joints against channel, T and I/H sections.

### Slice 4: standard section tables

EN 10365 IPE, HEA/HEB, UPE/UPN size tables with root radii for the
parametric shapes CAD-434 adds.

## 10. Decisions (answered 2026-11-27)

| # | Question | Answer |
| --- | --- | --- |
| Q1 | Source of truth for node positions | Own node graph (recommended) |
| Q2 | Relationship to existing members and connections | No preference; the team applies the recommendation: extend, frame-owned `weldment.member` features, joints on nodes |
| Q3 | Scope of the first slice | **Full 3D first** (not the recommendation); see §9 |
| Q4 | Default at an L corner of equal profiles | Mitre (recommended) |
| Q5 | Drawing straight through a node | One member runs through (recommended) |
| Q6 | Meaning of "inner/outer" alignment | Relative to the loop (recommended) |
| Q7 | Where the tool lives | Free text: "Model -> Metal -> Frame". Profile picker and inline New profile stock as proposed |
| Q8 | First fabrication output after slice 1 | Cope templates and cutting sheet (recommended) |
| Q9 | Weld gap | Free text: "make configurable weld gap, defaults to 0". A per-joint value with a frame default of 0; a non-zero gap shortens members (confirmed) |
| Q10 | Slot profile fasteners at nodes | Automatic, overridable (recommended) |
| Q11 | Channels, T and I/H beams | **Add now in CAD-434** (not the recommendation): parametric shapes in the picker; joints frames cannot yet handle are refused in plain language; standard size tables later |

## 11. Spike

None was needed: every geometric operation of slice 1 (sweeping a placed
section, plane and cope cuts, measuring cut lists) already exists and is
tested in the weldment evaluator. The new work is the source model, the
derivation of joints from nodes and the drawing tool.
