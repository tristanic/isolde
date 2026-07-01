# ISOLDE — Claude guidance

## Project overview

ISOLDE is a [ChimeraX](https://www.cgl.ucsf.edu/chimerax/) plugin for interactive molecular dynamics-based model building into low-resolution crystallographic and cryo-EM density maps. It combines real-time GPU-accelerated simulation (OpenMM) with ChimeraX's molecular graphics and Qt-based UI.

**Language mix:** Python (plugin logic, UI, command interface) + C++ (performance-critical geometry, validation, OpenMM bindings via pybind11).

**Sister package required:** ChimeraX-Clipper (map handling).

**ISOLDE assumes a GUI session — it is not designed to run fully headlessly.**
Code reachable from normal interactive use (camera/mouse-mode changes, Qt
dialogs, spotlight/OpenGL-touching commands) is free to assume
`session.ui.is_gui` is true. An `AttributeError` from `session.ui.mouse_modes`,
`view.render`, etc. under `--nogui` is not, by itself, a bug — do not treat
"doesn't work headlessly" as something that generally needs fixing.

The one deliberate exception is the MCP/agent remote-control surface
(`isolde/src/remote_control/`), which does need to drive ISOLDE — including
starting/stopping simulations via commands like `isolde sim start` — from a
session that may have no GUI. Guard a GUI-only side effect with the existing
`self.gui_mode` (`session.ui.is_gui`) pattern (see `isolde.py`'s
`_set_right_mouse_mode_tug_atom`/`_sim_end_cb`) only when the code path is
actually reachable from that surface or from a command, not preemptively
elsewhere.

Most of `isolde/src/tests/` runs headlessly (`run_chimerax.bat --nogui --exit
--script ...`). When a test needs to exercise logic that incidentally sits
behind a GUI-only side effect the test doesn't care about (e.g. camera
repositioning in `navigate.py`), stub that side effect out in the test itself
rather than adding headless support to the underlying code.

---

## Build & install

All changes — Python or C++ — require a full wheel rebuild and ChimeraX relaunch. There is no editable/live-reload install.

**Windows:**
```bat
make_win.bat clean app-install          # daily/nightly ChimeraX
make_win.bat release clean app-install  # stable-release ChimeraX
```
The `release` token selects which ChimeraX installation to target (stable vs. daily); omit it when building against the daily build.

**Linux / macOS:**
```sh
make wheel
make install
```

**Prerequisites:** ChimeraX with BundleBuilder ≥ 1.4.0, ChimeraX-Core ~= 1.11, ChimeraX-Atomic ~= 1.61, OpenMM headers at the paths referenced in `bundle_info.xml`, a compatible C++ compiler (MSVC 2015+ / GCC 4.9+ / Xcode).

The build preprocesses `pyproject.toml.in` via `prep_toml.py` before invoking the ChimeraX bundle builder (PEP 517).

---

## Testing

There is one test file: `isolde/src/tests/test_simulation.py` (`SimTester` class). Tests are run manually inside ChimeraX — there is no pytest runner or CI pipeline. When making changes, load a structure in ChimeraX and exercise the affected feature interactively.

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

## Live editing: atoms and residues are never stable

The user can add or delete atoms/residues (and bonds) at essentially any point during an ISOLDE session — via ChimeraX's own editing commands, ISOLDE's model-building tools, or interactively mid-simulation. `Atom`/`Residue`/`Bond` objects are thin Python wrappers over C++ objects; deleting the underlying structure does not (by itself) stop a stale Python reference from existing. Dereferencing a deleted atomic object does not reliably raise a catchable Python exception — it can crash ChimeraX outright (null/dangling pointer read), which is far worse than an unhandled exception.

Guidance for any code (existing or new) that touches atomic objects:
- Don't hold onto `Atom`/`Residue`/`Bond` objects, or index positions into a `Collection`, across any interval where user interaction or Qt event-loop processing could occur (a callback, a trigger handler, a stored instance attribute meant to be read later). Re-fetch fresh collections at the point of use instead.
- Before dereferencing anything held across time (a model, manager, or residue stashed on `self` or in a closure), check `.deleted` first — this is already the established pattern throughout the codebase (see `isolde.py`, `openmm_interface.py`, `navigate.py`, `molobject.py`).
- Prefer ChimeraX's own attribute-registration mechanism (`register_attr` + get/set on the object itself) for any per-atom/per-residue data you need to persist, rather than a separate Python-side `{atom: value}` dict keyed by object identity — the former is tied to the object's own lifecycle and cannot go stale into a dangling reference; the latter can silently hold a reference to something already deleted.

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

## Change summaries

When summarising a set of file changes at the end of a task, use a three-column table:

| File | Status | >5 |
|---|---|---|
| `src/foo.py` | New | ✓ |
| `src/bar.py` | 1-line fix | |

- **File** — path relative to the working directory, as a markdown link.
- **Status** — plain-English description of what changed (e.g. "New", "Boot provider swapped", "CMAP fix").
- **>5** — contains a ✓ if more than 5 lines were substantively added or changed. Purely **moved** code counts as zero change regardless of line count.

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
