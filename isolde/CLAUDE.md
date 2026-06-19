# ISOLDE — Claude guidance

> ⚠️ **Read this first — and treat it as binding.** ISOLDE's stylistic, structural, and architectural
> conventions are deliberate and load-bearing: code, the Python↔C++ bindings, session save/restore, and
> the UI all depend on them holding consistently. **When you extend or modify this codebase, follow the
> patterns documented here exactly — do not introduce a parallel style, "cleaner" abstraction, or
> one-off shortcut.** If a change seems to require breaking one of these conventions, stop and ask the
> maintainer first rather than diverging. See **[Conventions to preserve](#conventions-to-preserve)** for
> the specifics; the rest of this file explains the layout and build that make them work.

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
The three trailing words are independent flags parsed by `make_win.bat`: `release` selects the stable
ChimeraX install (omit for the daily build), `clean` runs `devel clean .`, and `app-install` runs
`devel install .` (which compiles + installs). **Without `app-install` nothing is actually built or
installed** — the script just regenerates `pyproject.toml` and optionally cleans. `make_win.bat`
auto-detects Visual Studio 2022 via `vswhere` and initialises the **v141 (VS2017) toolset** to match the
ABI ChimeraX was compiled with; a pre-opened vcvars terminal is no longer required.

**Linux / macOS:**
```sh
make wheel
make install
```

**Prerequisites:** ChimeraX with BundleBuilder ≥ 1.4.0 (see versions below), a compatible C++ compiler
(MSVC v141 / GCC 4.9+ / Xcode), and an OpenMM importable inside ChimeraX's Python.

**Build chain (how it actually works).** Both the Windows and POSIX entry points do the same two things:
1. Run `prep_toml.py` **automatically** — it imports `openmm`, reads `openmm.version.openmm_library_path`,
   and substitutes the OpenMM lib/include dirs into `pyproject.toml.in` to produce the generated
   `pyproject.toml`. This is why OpenMM must be importable in ChimeraX's Python *at build time*.
2. Invoke the ChimeraX BundleBuilder (PEP 517 backend):
   `ChimeraX --nogui --safemode --cmd "devel build/install ."`.

There is **no `bundle_info.xml`** — the bundle migrated to `pyproject.toml.in` (PEP 517). That file is the
single source of truth for C++ build targets, command/tool declarations, and dependency pins.

**Dependency pins** (from `pyproject.toml.in` — treat *that file* as authoritative, not this list):
ChimeraX-Core `>=1.12rc0, ==1.12.*`, ChimeraX-Atomic `~=1.67`, ChimeraX-AtomicLibrary `~=14.4`,
ChimeraX-Arrays `~=1.0`, ChimeraX-Clipper `~=0.28.0`. ISOLDE itself is currently `1.13.0`
(`src/__init__.py`). Note Clipper appears only in the runtime `project.dependencies`, not in
`build-system.requires`.

---

## Testing

There is one test file: `isolde/src/tests/test_simulation.py` (`SimTester` class). It is a **manual
integration harness meant to be instantiated inside ChimeraX**, not a pytest module — there is no pytest
runner or CI pipeline. When making changes, load a structure in ChimeraX and exercise the affected
feature interactively.

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

## Architecture notes

A few non-obvious things that save search time:

- **Three-layer data model.** `isolde.py` is the session singleton and owns the restraint *managers*.
  `molobject.py` defines those managers (`_RestraintMgr` subclasses such as `MDFFMgr`,
  `PositionRestraintMgr`, `DistanceRestraintMgr`, `ProperDihedralRestraintMgr`, …) plus the individual
  per-object restraint classes. `molarray.py` wraps collections of those objects in `Collection`
  subclasses (`DistanceRestraints`, `ProperDihedrals`, `Ramas`, …) for **vectorized bulk property
  access** — it is the array counterpart to `molobject.py`'s per-object layer.

- **C++ binding split — know which side a symbol lives on.** ISOLDE uses *two* interop mechanisms:
  - Most subsystems (`atomic/`, `geometry/`, `interpolation/`, `validation/`, `restraints/`) use raw
    `extern "C"` functions exposed via `*_ext.h` headers, compiled together through `src_cpp/molc.cpp`
    into the **`molc` library** and called from Python with **ctypes**. Grep these headers for
    `extern "C" EXPORT`.
  - Only the OpenMM extensions (`_openmm_async`, `_openmm_force_ext` under `src_cpp/openmm/`) use
    **pybind11**. Grep for `PYBIND11_MODULE`.

- **OpenMM layer is a thin wrapper.** `openmm/openmm_interface.py` (`SimConstruct`, `SimManager`,
  `OpenmmThreadHandler`) builds the system + restraints and delegates GPU stepping to OpenMM on a
  background thread; coordinates flow back via the handler. Don't expect the physics loop to live in
  Python.

---

## Conventions to preserve

**These are rules, not suggestions.** The patterns below are load-bearing and consistent across the
entire codebase, and keeping them that way is a hard requirement — the maintainer considers consistency
here more important than any individual local improvement. **Match them exactly when extending existing
subsystems; do not introduce a parallel style, rename established symbols, or "modernise" a pattern in
passing.** If you believe a convention genuinely needs to change, raise it with the maintainer and change
it only on explicit instruction — never unilaterally.

### C-backed molecular objects (the central idiom)

Python classes that wrap a C++ object follow a fixed shape — see `PositionRestraint` /
`PositionRestraintMgr` / `PositionRestraints` in `molobject.py` + `molarray.py` as the canonical trio:

- **Class name ↔ C function prefix is mechanical.** A class maps to its C symbols by snake-casing its
  name (`_as_snake_case`): `PositionRestraintMgr` → `position_restraint_mgr_*`. New C functions must
  follow suit: `{prefix}_new`, `{prefix}_delete`, `set_{prefix}_py_instance`, and one accessor per
  property `{prefix}_<thing>`. Override the derived name only via the explicit `c_class_name=` arg (as
  `MDFFMgr` does with `"mdff_mgr"`).
- **Constructors take only `c_pointer`** and call `set_c_pointer(self, c_pointer)`. Per-object classes
  subclass `chimerax.core.state.State`; managers subclass `_RestraintMgr` (itself a ChimeraX `Model`),
  which handles the `{prefix}_new` / `set_{prefix}_py_instance` wiring for you.
- **Declare data as class-level `c_property(...)`**, never hand-written getters/setters. The per-object
  class in `molobject.py` uses singular names (`target`, `spring_constant`); the matching `Collection`
  in `molarray.py` uses `cvec_property(...)` with **plural** names (`targets`, `spring_constants`). Keep
  the two in sync — every per-object property should have its vectorized twin.
- **Don't loop in Python where a `cvec_property` exists** — bulk reads/writes go through the array layer
  for a reason (one C call, not N).
- **Pointer→Python converters** are module-level helpers named `_foo_or_none(p)` and `_foos(p)`, passed
  as `astype=`. They import the `molarray` `Collection` lazily *inside the function* to break the
  `molobject ↔ molarray` circular import.
- **Families of near-identical classes** define their properties once in a `_Base._init_methods()` and
  stamp each concrete subclass with `@delayed_class_init` plus `_C_FUNCTION_PREFIX` / `_MGR_GETTER`
  class attrs (see the `ProperDihedralRestraint` / `AdaptiveDihedralRestraint` pair).
- **Session save**: set `SESSION_SAVE=True`, implement `take_snapshot`/`restore_snapshot` that store the
  manager + atoms and *re-query the manager* on restore (don't serialize C state), and stamp
  `data['version'] = ISOLDE_STATE_VERSION`.

**Skeleton — adding a new C-backed restraint type.** Copy this trio; rename `Foo`/`foo` consistently so
the snake-cased class names line up with the C symbols (`FooRestraintMgr` ↔ `foo_restraint_mgr_*`). The
C side must export the matching `foo_restraint_mgr_new` / `_delete` / `set_..._py_instance` and the
per-property accessors.

```python
# --- molobject.py ---------------------------------------------------------

# pointer -> Python converters (used as astype=); import Collection lazily.
def _foo_restraint_or_none(p):
    return FooRestraint.c_ptr_to_py_inst(p) if p else None
def _foo_restraints(p):
    from .molarray import FooRestraints
    return FooRestraints(p)

class FooRestraintMgr(_RestraintMgr):
    '''Manages Foo restraints for a single atomic structure.'''
    SESSION_SAVE = True
    def __init__(self, model, c_pointer=None, auto_add_to_session=True):
        # _RestraintMgr derives the C prefix from the class name and calls
        # foo_restraint_mgr_new + set_foo_restraint_mgr_py_instance for you.
        super().__init__('Foo Restraints', model, c_pointer=c_pointer, allow_hydrogens='no')
        if auto_add_to_session:
            model.add([self])

    def get_restraints(self, atoms, create=False):
        f = c_function('foo_restraint_mgr_get_restraint',
            args=(ctypes.c_void_p, ctypes.c_void_p, ctypes.c_bool, ctypes.c_size_t, ctypes.c_void_p),
            ret=ctypes.c_size_t)
        n = len(atoms)
        ret = numpy.empty(n, cptr)
        num = f(self._c_pointer, atoms._c_pointers, create, n, pointer(ret))
        return _foo_restraints(ret[:num])

class FooRestraint(State):
    '''A single Foo restraint.'''
    SESSION_SAVE = True
    def __init__(self, c_pointer):
        set_c_pointer(self, c_pointer)

    mgr    = c_property('foo_restraint_get_manager', cptr, astype=_foo_restraint_mgr,
                        read_only=True, doc='The owning manager. Read only.')
    atom   = c_property('foo_restraint_atom', cptr, astype=convert.atom_or_none,
                        read_only=True, doc='The restrained atom. Read only.')
    target = c_property('foo_restraint_target', float64, doc='Target value. Writable.')
    spring_constant = c_property('foo_restraint_k', float64,
                        doc='Spring constant (kJ/mol/nm^2). Writable.')
    enabled = c_property('foo_restraint_enabled', bool, doc='Enable/disable. Writable.')

    def take_snapshot(self, session, flags):
        from . import ISOLDE_STATE_VERSION
        return {'restraint mgr': self.mgr, 'atom': self.atom, 'version': ISOLDE_STATE_VERSION}
    @staticmethod
    def restore_snapshot(session, data):
        mgr = data['restraint mgr']
        return mgr.get_restraint(data['atom']) if mgr is not None else None

# --- molarray.py ----------------------------------------------------------

class FooRestraints(Collection):
    def __init__(self, c_pointers=None):
        super().__init__(c_pointers, FooRestraint)
    # plural twins of every per-object c_property, via cvec_property:
    atoms   = cvec_property('foo_restraint_atom', cptr, astype=convert.atoms,
                            read_only=True, doc='Restrained atoms. Read only.')
    targets = cvec_property('foo_restraint_target', float64, doc='Target values. Writable.')
    spring_constants = cvec_property('foo_restraint_k', float64, doc='Spring constants. Writable.')
    enableds = cvec_property('foo_restraint_enabled', bool, doc='Enabled mask. Writable.')
```

### Events — ChimeraX `TriggerSet`

Declare trigger names as `UPPER_CASE` class constants bound to lowercase string literals, collect them
in a `trigger_names` tuple, `add_trigger` each in `__init__`, and fire with `activate_trigger(NAME,
data)`. **Always capture the handler returned by `add_handler` and remove it on teardown.** One-shot
callbacks return `chimerax.core.triggerset.DEREGISTER` to self-unregister.

### UI panels

New panels subclass `UI_Panel_Base` (`ui/ui_base.py`): register trigger handlers in `__init__` (honour
`sim_sensitive` and `expert_level`), call `gui.register_panel(self)`, and implement `cleanup()` to remove
every handler. Build Qt **programmatically** using `DefaultVLayout` / `DefaultHLayout` — do not add new
`.ui` Designer files (`IsoldeFrame.ui` is legacy).

### Parameters & settings

User-tunable parameters live in a `Param_Mgr` subclass decorated with `@param_properties` + `@autodoc`,
with a `_default_params` dict of `name: (default, unit_or_None)`. Read via attribute or `[]`, write via
`set_param`, and listen on the `PARAMETER_CHANGED` trigger. Persisted defaults go through a ChimeraX
`Settings` subclass with an `AUTO_SAVE` dict (`settings.py`).

### Imports, logging, errors

- **Import order**: stdlib → `chimerax.*` → `Qt.*` → relative (`from .module import …`). Use
  **in-function lazy imports** to break circular deps and for optional/conditional dependencies.
- **All user-facing output** goes through `session.logger.info/warning/error`. Raise
  `chimerax.core.errors.UserError` for bad user input — never a bare exception or `assert`. GUI
  warnings/confirmations use the helpers in `dialog.py`.

---

## Common task patterns

**New CLI command (worked example — `isolde brefine` / `isolde brsr`).** Registration is lazy and goes
through a three-hop chain, which is why a new command touches more than one file:
1. **Implement** the callback(s) in a subsystem module — e.g. `isolde_brefine_cmd` /
   `isolde_brsr_cmd` in `isolde/src/refine/bfactor_refine.py`. Custom argument parsing subclasses
   `AtomSpecArg` (e.g. `_StructuresOrSymmetryMgrsArg`).
2. **Write a `register_*_command(logger)`** in that module that builds a `CmdDesc` and calls
   `register('isolde <name>', desc, fn, logger=logger)` (see `register_bfactor_refine_commands` near the
   bottom of `bfactor_refine.py`).
3. **Wire it into `register_isolde(logger)`** in `isolde/src/cmd/cmd.py` — import and call your
   `register_*_command(logger)` alongside the others (the `brefine`/`brsr` wiring is there). That
   function is itself invoked lazily by `bundle_api.register_command()` in `isolde/src/__init__.py` the
   first time any `isolde …` command is referenced.

So: implement + per-module register function (steps 1–2), then one call added in `cmd/cmd.py` (step 3).
`__init__.py` rarely needs editing — it already routes the whole `isolde` command tree.

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
| `isolde/src/molobject.py` | Restraint *managers* + per-object restraint/validation classes (5,600 LOC) — the core data model |
| `isolde/src/molarray.py` | `Collection` subclasses for vectorized bulk access to restraint/geometry objects (array counterpart to `molobject.py`) |
| `isolde/src/openmm/openmm_interface.py` | OpenMM simulation engine wrapper (`SimConstruct`/`SimManager`/`OpenmmThreadHandler`) |
| `isolde/src/openmm/custom_forces.py` | Custom force field terms |
| `isolde/src/refine/bfactor_refine.py` | `isolde brefine`/`brsr` — ML isotropic B-factor/occupancy refinement; thin front-end over ChimeraX-Clipper, builds Barron robust-loss restraints |
| `isolde/src/validation/` | Ramachandran, rotamer, clash analysis |
| `isolde/src/cmd/` | CLI command implementations (`cmd.py` hosts `register_isolde`, the command-tree entry point) |
| `isolde/src/ui/` | Qt widget components |
| `isolde/src/tool.py` | Main GUI tool (`ISOLDE_ToolUI`) |
| `isolde/src_cpp/` | C++ geometry, validation, OpenMM bindings (see Architecture notes for the ctypes/pybind11 split) |
| `isolde/extern/pybind11` | Git submodule — used by the OpenMM extensions only |
| `prep_toml.py` | Pre-build step (runs automatically): injects OpenMM lib/include paths into `pyproject.toml` |
| `pyproject.toml.in` | Build config template + dependency/version source of truth (PEP 517; processed by `prep_toml.py`) |
