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

### Show R/S absolute configuration on selected chiral centres — DONE
Opt-in in-model markup that labels chosen chiral centres with their as-modelled
absolute configuration ("R"/"S"), complementing the always-on outlier glyph.
As-built design (the original "session-level `...Mgr`" idea was dropped — a non-
`Model` singleton can't draw cleanly):
1. **State on the object, ChimeraX-style.** A per-`ChiralCenter` boolean `label`
   (C++ `chiral.h/.cpp`, bound through `chiral_ext.h` → `ChiralCenter.label` /
   vectorised `ChiralCenters.labels`, cf. `Atoms.displays`). Off by default —
   protein backbone centres stay unlabelled. NB centres are transient (not
   serialized), so persistence lives on the **annotator's** snapshot, which records
   the labelled centres' atoms and re-applies the flag on restore (version 4
   unchanged). Required fixing a pre-existing gap: `ChiralAnnotator` was missing
   from the bundle `get_class` table and `validation/__init__.py` exports.
2. **Rendered by the existing per-structure `ChiralAnnotator`** (not a new
   manager): it draws `chirals[chirals.labels]` on two child drawings (one "R",
   one "S"), reusing its placement/visibility/camera-scaling machinery.
3. **Real 3D letters, camera-facing, fixed-size.** Precomputed extruded "R"/"S"
   meshes (`src/validation/_rs_glyphs.py`, generated offline by
   `tools/gen_rs_glyphs.py` — ChimeraX's Python has no font→outline lib, so no
   runtime dep). They face the camera (re-oriented on rotation via the annotator's
   `'graphics update'` gate) yet stay depth-tested/occluded, unlike a 2D label.
   Unlike the view-width-scaled outlier glyph these are purely informational, so
   they sit at a fixed `LABEL_HEIGHT_A` (0.5 Å) beside the atom (no zoom ramp).
   Colour via `set_label_color`: `'auto'` (contrasts with the background, tracked
   live), `'fromatoms'` (per-instance = each chiral atom's colour), or a custom RGBA.
4. **As-modelled config, cheaply.** `chirality.as_modelled_configs` = the CCD
   reference config (`reference_cip_codes`: RDKit CIP on the ideal geometry,
   cached per id) flipped iff `sign(expected_volume)·true_chiral_volume < 0` — no
   per-residue RDKit at *draw* time (recomputed only when the chiral set changes).
5. **UI/command.** `chiral [<atoms>] label true|false [labelColor auto|fromAtoms|<color>]`
   toggles the labels (no atom spec → current selection) and sets their colour.
   "Label R/S (sel.)" / "Clear R/S (sel.)" buttons in the Validation-tab Chiral
   widget emit `chiral sel label true|false` (GUI actions log their command
   equivalent, per ChimeraX convention). Mainly for ligands.

### cis/trans markup + restraints + flip for ALL flippable double bonds
ISOLDE marks up cis/trans peptide bonds and restrains `omega` to stop accidental
flipping. The same need applies — arguably more so — to **every "flippable"
double bond** (one not locked by a ring or other geometry). Rationale: a cis
peptide is at least eyeball-detectable (its HN and carbonyl O are displayed by
default), but a *trans* fatty-acid double bond is essentially invisible — with
nonpolar H hidden there's no visual cue it's even a double bond, let alone the
wrong isomer. Three pieces, mirroring the chiral work just completed:
1. **Restraints** — analogous to the `omega` peptide restraints, applied to the
   relevant dihedral of each flippable double bond (skip ring/locked ones).
2. **Markup** — needs fresh design. The peptide markup (a pseudo-trapezoid filling
   the "cup" between CA atoms, shown only when cis/twisted) won't transfer: for a
   natural fatty acid *cis* is correct and *trans* is wrong — the opposite default
   from peptides. Decide what "wrong/strained" looks like here and when it shows.
3. **Manual flip command** — a user-facing flip (cf. `isolde cisflip` /
   `isolde chiralflip`) for when a non-default isomer is genuinely intended.
Detecting "flippable" double bonds is naturally an RDKit-layer job (bond order +
ring membership), so it dovetails with the `rdkit_bridge`/`chemcomp` chemistry.
