# Extension bundles for ChimeraX

*How to ship a bundle that installs **into another bundle's namespace** and registers
subcommands under its command prefix. ISOLDE is the host in the examples; the pattern
works for any host bundle.*

## The mechanism

A bundle named `ChimeraX-ISOLDE-Foo` installs as the Python package **`chimerax.isolde.foo`**
— nested inside the host's namespace, but built, versioned, licensed and distributed as its
**own** bundle (its own wheel and `dist-info`).

| Bundle name | Installed package |
|---|---|
| `ChimeraX-ISOLDE-Foo` | `chimerax.isolde.foo` |

ChimeraX's BundleBuilder derives that path straight from the `[project] name`: strip
`ChimeraX-`, lowercase, hyphens → dots. You do **not** need `module-name-override` — the name
alone places the package. The same rule nests under any host: `ChimeraX-<Host>-<Ext>` →
`chimerax.<host>.<ext>`.

## Why nest instead of fork or patch

- **Import the host freely.** Python has no privacy: `from chimerax.isolde.<anything> import …`
  just works, so the vast majority of an extension needs **zero** change to the host.
- **Independent release.** Your bundle has its own version, wheel and license, and never gates
  the host's releases — ideal for an add-on that needs separate sign-off.
- **A clean seam for the rare host change.** When you genuinely need a host internal that isn't
  public, add a small public API **to the host** and flag that one change; everything else stays
  in your bundle.

## The manifest (`pyproject.toml`)

```toml
[project]
name = "ChimeraX-ISOLDE-Foo"            # -> installs as chimerax.isolde.foo
dependencies = [
  # dev-inclusive: a bare ==1.13.* excludes the daily's 1.13.dev* build
  "ChimeraX-Core   >=1.13.0.dev0, ==1.13.*",
  "ChimeraX-ISOLDE >=1.13.dev0,    ==1.13.*",   # the host you extend
]
dynamic = ["classifiers", "requires-python", "version"]

[tool.chimerax]
min-session-version = 1
max-session-version = 1
categories = ["General"]

# one block per subcommand you register
[tool.chimerax.command."isolde foo bar"]
category = "General"
description = "..."

# optional: ship user docs (enables help/usage auto-linking)
[tool.chimerax.package-data]
"src/docs" = ["**"]
```

Depend on the host bundle explicitly and match its `ChimeraX-Core` series. If you build against
the **daily**, use a dev-inclusive lower bound (`>=X.dev0, ==X.*`) — a plain `==X.*` sorts the
daily out and the build is rejected with an incompatible-version error.

## Registering subcommands under the host's verb

A **second** bundle can register commands under the host's existing top-level verb (`isolde …`)
— this works; no standalone verb fallback is needed. `__init__.py` defines a `BundleAPI` subclass
whose `register_command(bi, ci, logger)` forwards to a dispatcher (`cmd.register_command(ci.name,
logger)`); ChimeraX calls it lazily the first time each declared command is used:

```python
# cmd.py
from chimerax.core.commands import CmdDesc, register

def register_command(command_name, logger):     # one entry per declared command
    if command_name == "isolde foo bar":
        _register_foo_bar(logger)

def _register_foo_bar(logger):
    # import Arg types lazily so a --nogui import stays light
    from chimerax.core.commands import StringArg, BoolArg
    from chimerax.atomic import ResiduesArg
    desc = CmdDesc(
        optional=[("residues", ResiduesArg)],        # defaults to the selection
        keyword=[("some_option", StringArg),         # CLI: someOption
                 ("dry_run", BoolArg)],              # CLI: dryRun
        required_arguments=["some_option"],
        synopsis="...",
    )
    register("isolde foo bar", desc, foo_bar_cmd, logger=logger)
```

Two rules keep it wired: the name in `register(...)` must match the
`[tool.chimerax.command."…"]` block exactly, and snake_case keyword names auto-expose as
camelCase on the command line (`some_option` → `someOption`).

## Coexistence is safe

ChimeraX installs are **file-based** — each bundle's `dist-info/RECORD` lists only its own files.
A nested install writes only `chimerax/isolde/foo/*` and leaves the host's `__init__.py` untouched.
Rebuilding or uninstalling either bundle never removes the other's files, so host and extension
update independently.

## Build & test

```bash
# build + install into the ChimeraX user tree; pure-Python bundles are fast
make_win.bat app-install          # or:  <chimerax> ... devel install .

# run a headless test against the INSTALLED package
run_chimerax.bat --nogui --silent --exit --script src/tests/test_foo.py
```

## Gotchas, in one place

- **Install before importing.** Relative imports (`from .cmd import …`) resolve only *after*
  install, so tests must import the installed package (`chimerax.isolde.foo.…`), not a
  source-dir path.
- **Dev-inclusive version pins.** Building against the daily needs `>=X.dev0, ==X.*`; a bare
  `==X.*` rejects the daily.
- **Lazy GUI imports.** Import Qt (and other heavy/GUI-only libs) *inside* functions, so a
  `--nogui` import of your module never requires a display. This is also what lets you unit-test
  GUI-coupled orchestration headlessly — inject the event-loop calls (e.g. `QTimer.singleShot`)
  as seams and drive them from a plain queue in tests.
- **Resist editing the host.** Import what you need; when a host internal must become public,
  make that a small, flagged change *to the host* and keep everything else in your bundle.
- **Register under the host verb.** A second bundle extending `isolde …` works — declare each
  command in `pyproject.toml` and `register()` the exact same name.
