# Cross-CAD STEP exchange corpus

This folder holds the independent STEP corpus for [CAD-345](https://plans.tokenet.de/CAD/345). The `manifest.json` file is authoritative. It pins every fixture with:

- producing CAD, version and `FILE_NAME` strings;
- units, which are millimetres for every entry;
- source, provenance and redistribution basis;
- SHA-256 and geometry tags;
- the Rust import outcome and verdict;
- for every refusal, the PLAINER item that owns the gap.

The runner is `crates/multicad-app/tests/exchange_corpus.rs`. The results are recorded in the [ledger](../../../../../../docs/architecture/exchange-corpus-ledger.md).

## Locations

| Location | Meaning |
| --- | --- |
| `vendored` | The file is committed in this folder, or in `../`, and always runs. |
| `repository` | The file is committed elsewhere in this repository and always runs. |
| `kicad-library` | The file is referenced by path and SHA-256 inside the installed KiCad 10.0.0 3D-model library (`/Applications/KiCad/KiCad.app/Contents/SharedSupport/3dmodels`, or `MULTICAD_KICAD_3DMODELS`). It runs only when the installed bytes match. Otherwise the ledger reports `unavailable`, never a pass. |

## Producers and redistribution

- **FreeCAD 1.1.1**: `../freecad-drilled-box.step` was authored for this repository and is covered by its GPL-3.0-or-later licence.
- **OpenCascade 7.9.3**: the `occt-*.step` files were written by the repository's optional native reference build. To regenerate them:

  ```sh
  CARGO_TARGET_DIR=target/occt cargo test --locked -p multicad-app --features occt \
    --test exchange_corpus -- --ignored regenerate_opencascade_fixtures
  ```

  The geometry was authored for this repository and is covered by GPL-3.0-or-later. Regeneration changes the header timestamps, so the manifest digests must be renewed after it.
- **KiCad 10.0.0 (`kicad-cli`, OpenCascade 7.8.1 STEP processor)**: the `kicad-*.step` board exports come from the boards in `kicad-boards/`, which were authored for this repository and are covered by GPL-3.0-or-later. Each export was produced from the folder that holds the board and its model:

  ```sh
  kicad-cli pcb export step --board-only -f -o kicad-plate-holes.step kicad-boards/plate-holes.kicad_pcb
  ```

  `board-with-part` was exported without `--board-only`, with `occt-cylinder.step` beside the board as its component model.
- **KiCad 3D-model library (kicad-packages3D)**: the five files in `kicad-library/` are byte-for-byte copies, with their licence headers intact. Their licence is CC-BY-SA 4.0 with the KiCad library exception.
  - Copyright the KiCad Libraries team.
  - Source: <https://gitlab.com/kicad/libraries/kicad-packages3D>.
  - Licence terms: <https://gitlab.com/kicad/libraries/kicad-packages3D/-/blob/master/LICENSE.md>.

  These individual test fixtures are shared under the same licence; they are not a redistributed library collection. Other library entries are referenced only, never copied. Library files without a licence header are never vendored.

## Commands

```sh
# Pinned ledger, manifest provenance and the real-executable MCP workflow.
cargo test --locked -p multicad-app --test exchange_corpus

# The same ledger with OpenCascade 7.9.3 source and re-export readback.
CARGO_TARGET_DIR=target/occt cargo test --locked -p multicad-app --features occt \
  --test exchange_corpus corpus_matches_the_pinned_ledger

# Aggregate sweep over the whole installed KiCad library (about 12 minutes, debug build).
cargo test --locked -p multicad-app --test exchange_corpus -- --ignored kicad_library_sweep --nocapture
```

The ledger test writes `exchange-corpus-ledger.json` to the Cargo test temporary directory, or to the path in `MULTICAD_CORPUS_LEDGER`.
