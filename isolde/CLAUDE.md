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
make_win.bat release clean app-install
```

**Linux / macOS:**
```sh
make wheel
make install
```

**Prerequisites:** ChimeraX with BundleBuilder ≥ 1.4.0, ChimeraX-Core ~= 1.11, ChimeraX-Atomic ~= 1.61, OpenMM headers at the paths referenced in `bundle_info.xml`, a compatible C++ compiler (MSVC 2015+ / GCC 4.9+ / Xcode).

The build preprocesses `pyproject.toml.in` via `prep_toml.py` before invoking the ChimeraX bundle builder (PEP 517).

---

## Testing

There is no pytest runner or CI pipeline. Tests live in `isolde/src/tests/` and
run inside ChimeraX. Two kinds:

- **Self-running headless tests** follow a convention: a
  `if session is not None: run(session)` footer (ChimeraX injects `session`
  under `--script`), `run()` prints `ALL PASS` on success, and failures call a
  `_fail()` helper that raises `SystemExit`. Run one directly:

  ```bat
  run_chimerax.bat --nogui --exit --script src/tests/test_mmff_parameterisation.py
  ```

  `run_tests.bat` (repo root) discovers every `src/tests/test_*.py`, runs each
  self-running test headlessly, and reports pass/fail per file — it judges
  success by the `ALL PASS` sentinel in the output (a `SystemExit` from a
  ChimeraX `--script` does not reliably surface as a process exit code). Files
  with no `ALL PASS` sentinel are skipped, not failed. Examples:
  `test_mmff_parameterisation.py`, `test_atom_deletion_robustness.py`.

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
| `bundle_info.xml` | ChimeraX bundle manifest (tool/command declarations, dependencies) |
| `pyproject.toml.in` | Build config template (processed by `prep_toml.py`) |
