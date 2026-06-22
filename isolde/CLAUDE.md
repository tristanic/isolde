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
Python-only changes can skip `clean`: `make_win.bat release app-install`.

**Invoking the Windows build from Git Bash / an agent shell (non-cmd):**
`make_win.bat` uses paths relative to its own directory (`prep_toml.py`,
`devel install .`), and Git Bash/MSYS mangles `cmd.exe` flags (`/c`, `/d`) and
Windows paths. Two gotchas compound: MSYS rewrites `/d`, and `cmd.exe` won't
search the current directory for the batch file (NoDefaultCurrentDirectoryInExePath).
The invocation that works — **`cd` into the bundle dir first** (its internal
relative paths resolve against the inherited cwd), set `MSYS_NO_PATHCONV=1`, use a
single `/c`, and call the batch file by **absolute path**:
```sh
cd /c/Users/tcroll/my_gits/isolde/isolde && \
MSYS_NO_PATHCONV=1 cmd.exe /c "C:\Users\tcroll\my_gits\isolde\isolde\make_win.bat release app-install"
```
(`//c` only works *without* `MSYS_NO_PATHCONV`, but then `/d`/paths get mangled —
so prefer single `/c` + `MSYS_NO_PATHCONV=1` + absolute bat path.)

**Gotcha — silent build failure on wrong cwd:** if cwd is *not* the bundle dir
(e.g. you `cd`'d to the repo root for git first), `prep_toml.py` / `devel install .`
aren't found, the bat logs the error mid-stream **but still prints a stale
"Installed …" line and exits 0** — so the install silently keeps the previous
build. Always run the build from the bundle dir, and grep the log for
`error`/`prep_toml`/`does not exist`, not just "Installed".

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

---

## TODO / future work

### Generalise & expose the pyranose "chair" torsion restraints
ISOLDE ships a hand-curated set of torsion-restraint targets that pull pyranose
sugars toward the canonical **chair** conformation — overwhelmingly favoured in
reality, but in simulation it's easy to get trapped in a higher-energy
conformation and not notice. Two problems to solve (this mirrors the chirality
work — see the `chemcomp`-backed chiral generator in `src/atomic/chirality.py`):

1. **No user control (real problem).** There's currently no user-facing way to
   turn these restraints off or dial in a *different* target conformation. When
   non-chair conformations genuinely occur they're usually functionally
   important, and the restraints would silently fight them — so users need a way
   to disable / re-target per residue (or per selection).
2. **Narrow coverage.** Like the old curated chiral set, the targets cover only
   the small set of sugars curated so far. A CCD/`chemcomp`-driven approach
   (analogous to how `chiral_definitions_from_ccd` now generates chiral defs from
   the local CCD store) could generalise them to arbitrary pyranoses.

### Per-validator control over markup visibility for hidden atoms
ISOLDE's validation glyphs (Ramachandran, rotamer, chirality, …) only draw on
atoms that are actually visible (`display & !hide`) — so hiding part of the model
hides its markup. That's the right default, but for outlier markup there's a case
for the opposite: keep the marker visible even when its atoms are hidden, so a
serious error in a hidden region can't be missed. Add a user-facing, per-validator
option to choose this behaviour (e.g. a "show outliers even when atoms hidden"
toggle per validator type). The chirality annotator (`ChiralAnnotator` in
`src/validation/chiral_annotation.py`) currently hard-codes the visible-atoms-only
convention via `chirals.chiral_atoms.visibles`; the Rama annotator already has a
related `ignore_ribbon_hides` knob to model the API on.
