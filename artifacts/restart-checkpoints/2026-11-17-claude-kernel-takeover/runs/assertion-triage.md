# Assertion / wrong-result triage (workspace run 96cbc320)

Branch `claude/assertion-triage-56e2d5`, based on 3c58812e (96cbc320 + 5 STEP
commits). Scope: the `c_rust_wrong_result_or_assertion` group of
`runs/workspace-96cbc320/summary.json`, minus the gear tests (bevel ×2,
gear_pairs, internal_gears, shared_pair_cycles ×3, worm_gears; owned by the
gear session), plus the three refusal-classified single failures the task
named (occurrence_external_review, weldment copes).

Categories: **(a)** code regressed → code fixed; **(b)** test encoded an
OpenCascade-only premise or wording → test updated to the equally strict Rust
contract; **(c)** missing Rust capability → left failing.

Build env: pinned 1.99.1, `--locked`, `DEVELOPER_DIR=/Applications/Xcode.app/Contents/Developer SDKROOT=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk`, worktree-local target dir. macOS only.

## Per-test table

| Test | Cat. | Evidence | Action (commit) | Result |
|---|---|---|---|---|
| model `symbol_occurrence_cache::second_actual_catalog_placement_reuses_only_old_owned_nodes_with_complete_current_metadata` | a | `literal_symbols` was built only for `GeometryBackend::OpenCascade` since e78ea17f, so Rust compiles never reused symbol nodes | enable symbol-node reuse for Rust (a5d874e2) | pass |
| model `symbol_occurrence_cache::source_override_pose_pin_policy_mesh_and_material_changes_cannot_reuse_old_node_proofs` | a | as above; additionally an unchanged recompile hits `reuse_full_rust`, which reported features but not the symbol nodes it restored | same commit: full Rust reuse reports its symbol nodes in `symbol_reused` | pass |
| model `symbol_occurrence_cache::superseded_and_failed_candidates_publish_no_new_symbol_nodes_and_preserve_old_receipts` | a | as first row | a5d874e2 | pass |
| app `equipment_library::included_equipment_files_match_native_generation_and_graphical_tool_pins` | b (intended output change) | regenerated vs shipped bytes: manifest and every document field identical except `definition.preview.data` (PNG rendered from the display mesh, now Rust tessellation). Renders compared side by side: same parts, different faceting | regenerate 12 equipment files + library.json (fb0866d9) | pass |
| app `plant_library::included_plant_projects_match_native_generation_and_the_creation_tool_pin` | b (intended output change) | same: only the preview PNG differs | regenerate 2 plant files (fb0866d9) | pass |
| app `mechanical_library::shipped_bearing_defaults_reproduce_exactly_and_match_independent_native_dimensions` | b (intended output change) | same: only the preview PNG differs; bearing dimensions/pins unchanged | regenerate 4 bearing files (fb0866d9) | pass |
| app `wood_screw_catalog::exact_44_row_matrix_retains_nonmetric_source_limits_and_all_128_old_project_bytes` | b | legacy list still held pre-60096b7b hashes of kitchen sink and cooktop (that commit regenerated them and updated library.json only) | update the 128-entry legacy list for 60096b7b and the 18 re-pinned files (fb0866d9) | pass |
| kernel `advanced_basics::asymmetric_chamfer_uses_the_explicit_face_and_distance_angle` | a | Rust refusal "must touch the selected edge" dropped the adjacency wording of the OCCT refusal both kernel and UI tests pin | refusal text now "is not adjacent to the selected edge" (9447926e) | pass |
| ui `modeling::controller::chamfer_tests::native_sided_chamfer_mode_switch_retains_source_and_failed_nonadjacent_face_is_repairable` | a | same refusal surfaced in the UI | 9447926e | pass |
| ui `modeling::controller::edge_treatment::tests::failed_native_fillet_receipt_is_visible_repairable_and_cannot_commit` | a | impossible 2 mm fillet on 1 mm stock surfaced the generic policy text "Length is below minimum feature size" | Rust box fillet now reports "Edge treatment failed: the fillet radius leaves a straight box face below the minimum feature size", same error kind (4758b02b, f473fcd9) | pass |
| kernel `planar_contours::circular_and_semicircular_faces_remain_exact_and_curved_or_elliptic_targets_fail` | a | oblique-bore body has no `planar_loops` branch; the cylinder fallback read face 5 as edge 5 (a vertical side line) and failed with Input "Boundary is outside the contour plane tolerance" | refuse oblique-bore faces and non-line/arc edges as Unsupported ellipse/spline (4265c80b) | pass |
| model `buttress_threads::buttress_wall_and_cosmetic_manufacturing_failures_remain_explicit` | a | refusal is correct (wall < `minimum_wall`) but said only "insufficient wall" | text keeps "insufficient wall" and names the minimum wall (989a3ddd, f473fcd9) | pass |
| kernel `rust_manufacturing_mesh::rust_batched_cuts_keep_empty_and_refused_outcomes_distinct` | b | `Solid::box_solid` is Rust now (OCCT_LINKED = false); the "native" tool became a corner-touching Rust box → InvalidResult | Unsupported case uses a disjoint exact cylinder the rectilinear batch cannot represent (e9f2927a) | pass |
| kernel `interference::rust_pair_verifies_gap_contact_interference_and_containment` | b | `box_solid` "mixed" operand is Rust; actual result Overlapping, 1000 mm³, distance 0 — correct for identical boxes | assert that result; keep an Unsupported check with an oblique cylinder (e9f2927a) | pass |
| kernel `rust_distance::whole_compound_containment_classifies_members_independently_and_refuses_curves` | b | `stock(false, …)` = `box_solid` is Rust; containment distance 0 is correct | assert distance 0, inner solution and the witness contract; curve refusals unchanged (e9f2927a) | pass |
| kernel `interference::rotated_intersecting_solids_keep_only_consistent_native_witnesses` | b | required OCCT's filtered set (`total > retained`, incomplete). Rust planar distance enumerates every distinct witness (`complete ⇔ total == retained`, as `rust_distance::check_witnesses`) | assert the Rust contract (e9f2927a). See note 1 | pass |
| kernel `edge_treatment::named_cylinder_rims_treat_exactly_but_seam_and_impossible_radius_fail_without_mutation` | b | Rust cylinder topology: 2 rim arcs, no seam line (probe) | assert exactly the two rim edges exist; impossible-radius and immutability checks kept; renamed `…_without_a_seam_and_impossible_radius_fails_without_mutation` (e9f2927a) | pass |
| model `rust_rounded_threads::unqualified_partial_and_external_receiving_tools_refuse_real_stock_atomically` | b | partial CutInternal now builds a `RoundedBlindBore` (blind rounded receivers added 2026-11-17 in d5e2c4c3/622286ca, after the test) | assert blind receiver, major radius, stock height, floor > 0, volume < stock; external refusal kept; renamed `partial_internal_receivers_are_blind_and_external_…` (197e4700) | pass |
| kernel `sketch_references::rational_and_nonrational_spline_contacts_are_native_parameterized_and_tangent` | c | spline side edges have no analytic geometry; every query refuses "Rust sketch references support finite vertices, lines and conics…" (probe), so no B-spline edge is found | none | still fails |
| kernel `sketch_references::singular_conic_collapsed_and_closed_retraced_spline_projections_fail_but_open_line_locus_survives` | c | circle edge is found; the spline edge lookup fails for the same reason | none | still fails |
| model `three_mf_indexed_manifold::native_faceted_box_exports_an_indexed_shell_and_exact_world_build_transforms` | b | Rust box tessellation already shares its 8 corners | fixture splits every triangle into its own vertices; writer must still recover the 8-vertex shell (af2871f2) | pass |
| model `three_mf_indexed_manifold::transformed_fused_plate_and_blind_cross_holes_use_proven_native_seams_without_moving_the_solid` | b | raw Rust mesh: 114 vertices, manifold, writer accepts; receipt displacement 0 | assert raw is a manifold accepted shell and displacement == 0; corrected-mesh, volume, immutability checks kept (af2871f2) | pass |
| model `developed_sheet_budget::larger_developed_profile_keeps_exact_overlap_validation_and_the_query_limit` | c | duplicate hole is refused (Input) but as "Invalid analytic normal extrusion: InvalidProfile"; the Rust `normal_sweep` validator reports all crossing/nesting/degenerate loops as one `InvalidProfile` variant (its own unit tests pin that), so no overlap diagnosis exists | none (needs an overlap-specific validation error in the Rust kernel) | still fails |
| model `mesh_exports::project_and_tag_scopes_keep_glb_membership_and_part_appearances` | b | 18781.42 vs 18785.40: boxes exact, pin is a 36-facet inscribed cylinder, which is exactly what the recipe's 10° angular deviation allows; ±2 was met only by OCCT's finer-than-requested mesh | per part: boxes 6000 to f32 storage; pin between the inscribed 36-gon and the true cylinder (af2871f2). See note 2 | pass |
| model `custom_connection_hardware::all_three_engineered_hardware_forms_have_actual_helices_independent_mass_material_and_step` | c | case 0 passes; case 1 (wood insert) names thread faces `wood-insert/…` without a `thread/` role, and its STEP readback is refused "Imported boundary is outside the retained Rust families" | none (STEP readback of the wood-insert helix is a Rust gap; role naming left unchanged to keep persisted identities) | still fails |
| ui `mechanical_import::tests::cancel_and_concurrent_edits_never_commit_a_stale_candidate` | a | `read_step_file` still called the OpenCascade reader ("OpenCascade is unavailable"); `read_kicad_file` already used the default backend | `read_step_file` follows `GeometryBackend::default()` (e7a5ec65) | pass |
| ui `mechanical_import::tests::update_target_follows_the_selected_import_part` | a | same | e7a5ec65 | pass |
| ui `mechanical_import::tests::minimum_hole_size_reads_the_board_again_off_the_window_thread` | b | 21 holes at 0.9 mm; OCCT refused via the 16-hole interactive profile budget, the Rust board is one exact normal sweep without it | assert 21 holes read and checked OK (re-read and 12-hole checks kept) (4db699f4) | pass |
| model `occurrence_external_review::current_native_spline_basis_regenerates_without_rewriting_authored_output_or_pole_ids` | c | "This rational spline profile cannot yet be extruded" | none | still fails |
| model `occurrence_external_review::editable_native_review_generator_roundtrips_twelve_projects_and_refuses_overwrite` | c | "Rust sketch references support finite vertices, lines and conics…" | none | still fails |
| model `occurrence_external_review::finite_curved_face_section_drives_a_real_closed_solid_and_empty_section_cannot_reuse_guesses` | c | same sketch-reference refusal | none | still fails |
| model `weldment::copes_fit_the_actual_receiving_member_with_clearance` | c | "Rust coped ends require straight-sided member sections" (CHS round tube copes) | none | still fails |

## Notes

1. The Rust planar distance marks a witness list complete when it holds every
   distinct witness it enumerated. For overlapping solids the true
   zero-distance set is a continuum, so "complete" there means complete within
   that enumeration, not every closest pair. Not changed here; worth a product
   decision on the `witnesses_complete` meaning for overlap/contact.
2. The pin bound is wider than the old fixed ±2 mm³ in the deficit direction
   (allowed ≈ 3.98 mm³) and every box is now checked per part to f32 storage
   instead of only the total. The bound follows from the recipe's declared
   angular deviation; a mesh that ignored it would still fail.
3. Re-pinned library previews are produced on macOS. Determinism of the
   preview renderer across Linux/Windows CI was not checked here.
4. The new chamfer/fillet/wall wording is user-visible refusal text only;
   refusal kinds and conditions are unchanged.

## Verification

- `cargo fmt --all -- --check`: clean.
- `cargo clippy --workspace --all-targets --locked -- -D warnings`: clean.
- Affected tests rerun individually after each change (results above).
- Full `cargo test --workspace --locked --no-fail-fast` at 4db699f4..fb0866d9
  (all fixes except f473fcd9): 7714 passed, 259 failed, 14 ignored, against
  284 failures at 96cbc320. 25 baseline failures now pass: 22 from this
  table plus three STEP trim-validation tests (`mcp.rs` ×2,
  `rust_helical_tube_cut.rs`) fixed by the 5 STEP commits this branch is
  based on. The 2 renamed tests pass under their new names. Two new failures came from my first refusal
  rewording (`multicad-model/tests/lead_screw_threads.rs` pinned
  "insufficient wall", `multicad-model/tests/edge_treatment.rs` box fillet
  pinned "minimum feature size"); f473fcd9 keeps both terms, and those two
  targets plus `buttress_threads` and the UI fillet receipt were rerun and
  pass (`buttress_threads::named_buttress_native_pairs…` is a baseline
  failure). No other test changed state.
- Of the 28 assertion-group tests in scope, 24 now pass and 4 remain failing
  as Rust capability gaps (sketch_references ×2, developed_sheet_budget,
  custom_connection_hardware). The 4 refusal-classified tests the task also
  named (occurrence_external_review ×3, weldment copes) are category c and
  still fail.
