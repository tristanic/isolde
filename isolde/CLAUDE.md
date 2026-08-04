# ISOLDE — Claude guidance

## Project overview

ISOLDE is a [ChimeraX](https://www.cgl.ucsf.edu/chimerax/) plugin for interactive molecular dynamics-based model building into low-resolution crystallographic and cryo-EM density maps. It combines real-time GPU-accelerated simulation (OpenMM) with ChimeraX's molecular graphics and Qt-based UI.

**Language mix:** Python (plugin logic, UI, command interface) + C++ (performance-critical geometry, validation, OpenMM bindings via pybind11).

**Sister package required:** ChimeraX-Clipper (map handling).

---

## Build & install

All changes — Python or C++ — require a full wheel rebuild and ChimeraX relaunch. There is no editable/live-reload install.

**Windows:**
```bat
make_win.bat clean app-install
```

Note the absence of a `release` token: ISOLDE requires `ChimeraX-Core ==1.13.*`,
which is currently the **daily** build (`make_win.bat`'s default). Passing
`release` targets the stable install and the build is rejected outright.

**Linux / macOS:**
```sh
make wheel
make install
```

**Prerequisites:** ChimeraX with BundleBuilder ≥ 1.6.0, ChimeraX-Core ==1.13.\*,
ChimeraX-Atomic ~= 1.67, ChimeraX-AtomicLibrary ~= 14.4, ChimeraX-Arrays ~= 1.1,
ChimeraX-Clipper ==0.28.\*, and a C++ compiler (MSVC / GCC / Xcode).

**OpenMM is no longer a separate install** — ChimeraX 1.13 bundles it (8.4 as of
the 2026-08 daily). `prep_toml.py` locates its libs/headers at build time via
`openmm.version.openmm_library_path` and substitutes them into `pyproject.toml`.

**Windows toolchain:** ChimeraX 1.13 ships Python 3.14 built with VS2022
(MSC v.1944), and its prebuilt C++ libraries use the v143 toolset. MSVC only
guarantees binary compatibility when the *linking* toolset is at least as new as
every input, so `make_win.bat` takes vcvars' default (newest installed) and then
verifies it against the ChimeraX being targeted, aborting if it is older. It
deliberately does **not** pin a toolset: a pin cannot express "match ChimeraX"
(UCSF does not publish theirs), and the previous `-vcvars_ver=14.1` pin silently
outlived ChimeraX's move off Python 3.11.

The build preprocesses `pyproject.toml.in` via `prep_toml.py` before invoking the
ChimeraX bundle builder (PEP 517). `pyproject.toml` is generated and git-ignored.

---

## Testing

There is no pytest runner or CI pipeline. Tests live in `isolde/src/tests/` and
run inside ChimeraX. Two kinds:

- **Self-running headless tests** follow a convention: a
  `if session is not None: run(session)` footer (ChimeraX injects `session`
  under `--script`), `run()` prints `ALL PASS` on success, and failures call a
  `_fail()` helper that raises `SystemExit`. Run one directly:

  ```bat
  run_chimerax.bat --nogui --exit --script src/tests/test_session_roundtrip.py
  ```

  `run_tests.bat` (repo root) discovers every `src/tests/test_*.py`, runs each
  self-running test headlessly, and reports pass/fail per file — it judges
  success by the `ALL PASS` sentinel in the output (a `SystemExit` from a
  ChimeraX `--script` does not reliably surface as a process exit code). Files
  with no `ALL PASS` sentinel are skipped, not failed. Examples:
  `test_session_roundtrip.py`, `test_cmd_surface.py`.

- **Manual/interactive harness:** `test_simulation.py` (`SimTester` class) is
  driven by hand inside a GUI ChimeraX session; it has no headless footer and is
  skipped by `run_tests.bat`.

When making changes, also load a structure in ChimeraX and exercise the affected
feature interactively.

**Do NOT try to test ISOLDE (or Clipper map setup) headless via `--nogui` or
`--offscreen`.** It does not work and is not worth attempting:
- `isolde start` explicitly refuses non-GUI mode (`cmd.py`: *"ISOLDE currently
  requires ChimeraX to be in GUI mode"*), so `session.isolde` is never created.
- Even constructing the map layer directly, Clipper's `SymmetryManager` enables
  spotlight mode on model-add, which calls `camera.camera(session, 'ortho')` →
  `view.render.opengl_context` is `None` with no live GL context (true for both
  `--nogui` and `--offscreen` on this platform) → `AttributeError`. Map display,
  contouring and surface upload all assume a real GL context too.

Verification is **interactive in the GUI**. For scripted/agentic control of an
*already-running* GUI ISOLDE, the MCP control surface is an option
(`isolde/src/mcp/README.md`). Otherwise, ask the user to run the manual steps.

---

## Code style

**Python formatting:** YAPF with the configuration in `setup.cfg` (95-character line limit). Run YAPF before committing. flake8 is also configured in `setup.cfg`; E501 (line length) is ignored there.

**Naming:** `snake_case` for functions and variables; `CamelCase` for classes.

**File headers:** Every new Python file should carry the standard author/license block:
```python
# @Author: Tristan Croll
# @Date:   DD-Mon-YYYY
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: DD-Mon-YYYY
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: YYYY Tristan Croll
```

---

## High-care areas

Do **not** modify the following without an explicit instruction to do so. They contain subtle physics assumptions, tight interdependencies, or require a full build to verify safely:

| File / directory | Why |
|---|---|
| `isolde/src/openmm/openmm_interface.py` | Core simulation engine — physics assumptions baked in throughout |
| `isolde/src/openmm/custom_forces.py` | Custom OpenMM force expressions — math errors break silently |
| `isolde/src/molobject.py` | Central molecular object layer — heavy interdependencies across the codebase |
| `isolde/src_cpp/` | C++ extensions — changes require a rebuild and are harder to spot-test |

---

## Common task patterns

**New CLI command:**
1. Implement in `isolde/src/cmd/` (or a relevant submodule)
2. Register the command descriptor in `isolde/src/__init__.py` under `bundle_api.register_command()`

**New UI widget or panel:**
- Add to `isolde/src/ui/` for standalone widgets, or extend `isolde/src/tool.py` (`ISOLDE_ToolUI`) for changes to the main panel

**Documentation:**
- User-facing docs live in `docs/source/` as Sphinx RST files
- Build with `make docs` (requires LaTeX for PDF output)
- Docstrings use triple-quoted format with parameter descriptions

---

## Key file map

| File | Role |
|---|---|
| `isolde/src/__init__.py` | Bundle API entry point; command and tool registration |
| `isolde/src/isolde.py` | Session-level ISOLDE singleton (1,500 LOC) |
| `isolde/src/molobject.py` | Restraint managers, validation objects (5,600 LOC) — the core data model |
| `isolde/src/molarray.py` | Bulk array operations on molecular collections |
| `isolde/src/openmm/openmm_interface.py` | OpenMM simulation engine wrapper |
| `isolde/src/openmm/custom_forces.py` | Custom force field terms |
| `isolde/src/validation/` | Ramachandran, rotamer, clash analysis |
| `isolde/src/cmd/` | CLI command implementations |
| `isolde/src/ui/` | Qt widget components |
| `isolde/src/tool.py` | Main GUI tool (`ISOLDE_ToolUI`) |
| `isolde/src_cpp/` | C++ geometry, validation, OpenMM bindings |
| `isolde/extern/pybind11` | Git submodule — Python ↔ C++ bindings |
| `pyproject.toml.in` | ChimeraX bundle manifest template — tool/command declarations, dependencies, C++ library/extension definitions. Processed by `prep_toml.py` into the generated (git-ignored) `pyproject.toml`. Replaced the old `bundle_info.xml`. |
