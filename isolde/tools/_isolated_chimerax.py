# @Author: Tristan Croll
# @Date:   27-Jun-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 27-Jun-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
"""
Per-lane ChimeraX launch shim.

ChimeraX resolves its single per-user data/config/cache directory at startup in
``chimerax.core.__main__._set_app_dirs()`` via
``appdirs.AppDirs("ChimeraX", "UCSF", version=...)``, then assigns it to
``site.USER_BASE`` / ``site.USER_SITE`` (where user bundles install and load
from). There is no CLI flag or environment variable that redirects it, and the
``ChimeraX.exe`` launcher hardcodes ``python -I`` (isolated), so PYTHONPATH /
sitecustomize injection is impossible through the normal executable.

This shim is the escape hatch: it is run as the *main script* by the bundled
``bin/python.exe`` (``python.exe -I _isolated_chimerax.py <root> <args...>``).
``-I`` reproduces ChimeraX's own isolation (no stray env / OS user-site leakage);
the lane root arrives via ``argv`` so ``-I`` ignoring the environment does not
matter. We monkeypatch ``appdirs`` *before* importing ``chimerax.core`` so every
ChimeraX user directory lands under ``<root>`` instead of the shared per-user
location -- giving each development "lane" its own isolated install tree.

Driven by ``run_chimerax.bat`` / ``run_chimerax.sh``; not meant to be run by hand.
"""

import os
import sys
import runpy


def main():
    if len(sys.argv) < 2:
        sys.stderr.write(
            "usage: python -I _isolated_chimerax.py <lane-root> [chimerax args...]\n"
        )
        raise SystemExit(2)

    root = os.path.abspath(sys.argv[1])
    chimerax_argv = sys.argv[2:]

    os.makedirs(root, exist_ok=True)

    import appdirs

    if sys.platform == "win32":
        # On Windows appdirs derives every user_*/site_* directory from
        # _get_win_folder(<CSIDL>) (which normally calls SHGetFolderPathW and
        # ignores %LOCALAPPDATA%). Forcing it to return <root> for every folder
        # id re-roots the whole tree there while preserving appdirs' own
        # "<root>\UCSF\ChimeraX\<version>" layout below it.
        appdirs._get_win_folder = lambda csidl_name: root
    elif sys.platform == "darwin":
        # macOS appdirs hardcodes ~/Library/... and ignores XDG_*; left as a
        # stub -- fill in only if/when macOS testing is needed.
        sys.stderr.write(
            "_isolated_chimerax: macOS redirect not implemented; "
            "ChimeraX will use its normal shared user directory.\n"
        )
    else:
        # Linux appdirs honours XDG_*_HOME, read at call time, so setting these
        # before ChimeraX starts redirects data/config/cache/state.
        os.environ["XDG_DATA_HOME"] = os.path.join(root, "data")
        os.environ["XDG_CONFIG_HOME"] = os.path.join(root, "config")
        os.environ["XDG_CACHE_HOME"] = os.path.join(root, "cache")
        os.environ["XDG_STATE_HOME"] = os.path.join(root, "state")

    # Reproduce the argv[0] the ChimeraX executable would pass to chimerax.core.
    sys.argv = ["ChimeraX", *chimerax_argv]
    runpy.run_module("chimerax.core", run_name="__main__")


if __name__ == "__main__":
    main()
