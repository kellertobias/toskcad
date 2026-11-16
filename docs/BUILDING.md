# Building the native development host

## Commands

From any working directory, invoke the repository's `build` executable (Windows: `build.ps1`):

```sh
./build                 # fetch locked Cargo dependencies; build every component
./build open            # rebuild; launch another native application window
./build check           # formatting, Clippy and tests
./build archive         # runtime ZIP and matching source ZIP
./build archive install # explicit macOS install/start into /Applications/MultiCAD.app
```

ZIP output supports the current host: macOS 13+ ARM64, Linux x64 or Windows x64. Other operating systems build on their own GitHub runners; the wrapper does not dispatch remote workflows. On PowerShell use `./build.ps1` with the same arguments. Install is not performed by any other command and has not been authorized/tested in this session.

## Prerequisites and bootstrap

Install Node.js 22.14+ for developer orchestration and Rust through rustup. macOS needs Apple Command Line Tools (`xcode-select --install`) for the AppKit/Foundation adapter. Windows needs the MSVC linker and Windows SDK supplied by Visual Studio Build Tools. Linux needs its system linker and Wayland/xkbcommon development packages (Ubuntu: `sudo apt-get install build-essential libwayland-dev libxkbcommon-dev`). Platform SDK installation can require administrator interaction; the wrapper reports missing bootstrap tools without installing system packages.

The wrapper selects the Rust version in `rust-toolchain.toml` without changing the global default. It fetches locked Cargo dependencies and builds the complete workspace. Geometry uses the owned Rust kernel. The macOS window adapter is the remaining application CXX boundary. Node.js is developer tooling; the running application requires no Node.js or Python interpreter.

## Cargo output and maintenance

Repository Cargo commands default to `target/cargo/`. `CARGO_TARGET_DIR` can select a separate output directory. The development wrapper and release packager obtain the executable path from Cargo metadata, so they package the selected build rather than an older default-path binary.

Before clearing Rust products, finish or cancel builds and tests using that output directory. `cargo clean` removes the configured Cargo output; application bundles and project files live elsewhere. For a separate build, choose a fresh `CARGO_TARGET_DIR`. Do not delete or move an active build directory. No automatic cleanup runs while building.

## Outputs and runtime checks

Default builds assemble a new development layout under `artifacts/development/build-*/`; `latest.json` records its location. macOS includes `MultiCAD.app`; Windows/Linux include the native executable and resources. `open` starts a new host after a successful build and records its process. Existing editable windows keep their documents and close prompts; the wrapper does not terminate them or unrelated applications.

Archive builds require clean committed source, stamp the requested version/source SHA, run workspace tests, build release products and verify the packaged executable's version/identity, dependency closure and real STEP smoke result. `MULTICAD_RELEASE_VERSION` overrides the development archive version. ZIPs, manifests and SHA-256 checksums are under `artifacts/release/`:

- `multicad-<version>-<target>.zip`: native application layout and third-party notices.
- `multicad-<version>-<target>-sources.zip`: exact application source, lockfile, vendored Cargo sources and offline rebuild instructions.

Source archives include the dependency source needed for offline rebuilding. Package verification rejects retired geometry archives, native libraries and prototype assets, and requires Rust kernel diagnostics in every runtime ZIP. macOS packages have local ad-hoc signatures and are not notarized; Windows/Linux packages are unsigned portable development layouts, not publisher-signed installers. Linux targets an Ubuntu 24.04+ desktop/system graphics stack; Windows requires a supported system graphics driver. The MSVC target uses the static CRT and packaging rejects imports of an unbundled VC runtime. No Windows/Linux binary is claimed tested until its runner actually executes.

`archive install` requires an explicit invocation. It refuses while an installed application is running, stages the new macOS bundle, retains the previous installation as a sibling backup, restores it on replacement failure and opens the new copy. Existing application configuration is not removed. `/Applications` must be writable through normal OS permissions; the script does not hide privilege escalation. Install/start is not part of the checks performed here.

Useful headless checks:

```sh
cargo run --locked -- --version
cargo run --locked -- --diagnostics
cargo run --locked -- --self-test artifacts/native-smoke
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --locked -- -D warnings
cargo test --workspace --locked
```

With no arguments `multicad` opens the native host. `--self-test` checks a Rust solid, Boolean cut, tessellation and STEP round-trip, and opens the bundled examples. It is a packaged runtime smoke check; complete CAD workflows require their own tests.

## Components and limits

| Component | Responsibility |
| --- | --- |
| `multicad-core` | Canonical document, typed commands, parameters and project containers. |
| `multicad-precision` | Owned exact arithmetic and bounded numerical primitives. |
| `multicad-rust-kernel` | Owned geometry sources and independent manufacturing/exchange qualification. |
| `multicad-kernel` | Native facade for retained Rust geometry, workers and geometry policy. |
| `multicad-solver` | Rust constraint solver for sketch geometry. |
| `multicad-model` | Parametric document evaluation and manufacturing workflows. |
| `multicad-ui` | Native application shell, viewport and macOS window adapter. |
| `multicad-mcp` | Typed external application control using the shared document API. |
| `multicad-app` | Entry point, diagnostics and packaged runtime smoke checks. |

Supported geometry and exports remain subject to the [Rust kernel qualification contract](architecture/rust-kernel-qualification.md). The [foundation study](architecture/foundation/verification.md) and [kernel replacement baseline](architecture/kernel-replacement-baseline.md) retain historical evidence about the retired geometry provider. Their build commands are obsolete; their results do not establish current Rust parity or coverage.

GitHub Actions builds three targets using the same package code. Read-only build/check jobs are separate from the gated tag-only publisher. Hosted execution, physical Windows/Linux behavior, signing, low-power performance and human acceptance remain separate evidence. GitHub is primary; see [RELEASING.md](RELEASING.md). Nothing is pushed automatically.
