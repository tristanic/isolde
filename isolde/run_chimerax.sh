#!/usr/bin/env bash
# ============================================================================
# run_chimerax.sh -- launch ChimeraX with a per-lane isolated user directory.
#
# POSIX counterpart of run_chimerax.bat. See that file (and
# tools/_isolated_chimerax.py) for the full rationale. In short: ChimeraX
# installs every user bundle into one shared per-user directory, so parallel
# worktrees clobber each other at "devel install" time. This wrapper launches
# the bundled python directly on the shim, re-rooting ChimeraX's whole user
# tree under a lane-specific directory.
#
# Usage:   ./run_chimerax.sh [release] [<ChimeraX args>...]
#   release   -> stable install; omitted -> the daily build.
#   CHIMERAX_APP=<dir>  (env) overrides install autodetection.
#
# Lane resolution: walk up from this script's dir for a ".chimerax-lane"
# marker; if found at <dir>, the isolated tree lives at <dir>/.chimerax-home.
# With no marker, fall back to a unique per-worktree home and warn.
#
# NOTE: the Windows .bat is the battle-tested path on the primary dev box;
# this .sh mirrors its logic but should be smoke-tested on Linux/macOS.
# ============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --- optional leading "release" token ---
RELEASE=0
if [ "${1:-}" = "release" ]; then
    RELEASE=1
    shift
fi

# --- locate the ChimeraX install and its bundled python ---
uname_s="$(uname -s)"
find_python() {
    local app="$1"
    local cand
    case "$uname_s" in
        Darwin) cand="$app/Contents/bin" ;;
        *)      cand="$app/bin" ;;
    esac
    # bundled interpreter is python3.<minor>; pick the highest match
    ls "$cand"/python3.* 2>/dev/null | sort -V | tail -n1
}

# Resolve both the bundled python (for isolated launches) and the normal
# ChimeraX launcher (for standard, un-laned launches).
PY=""
CXLAUNCH=""
resolve_install() {
    local app="$1"
    PY="$(find_python "$app")"
    case "$uname_s" in
        Darwin) CXLAUNCH="$app/Contents/bin/ChimeraX" ;;
        *)      CXLAUNCH="$app/bin/ChimeraX" ;;
    esac
}

if [ -n "${CHIMERAX_APP:-}" ]; then
    resolve_install "$CHIMERAX_APP"
else
    case "$uname_s" in
        Darwin)
            if [ "$RELEASE" = 1 ]; then resolve_install /Applications/ChimeraX.app
            else resolve_install /Applications/ChimeraX_Daily.app; fi
            ;;
        *)
            # Linux: the launcher is chimerax/chimerax-daily on PATH; the bundled
            # python lives under its install dir.
            if [ "$RELEASE" = 1 ]; then cxbin="$(command -v chimerax || true)"
            else cxbin="$(command -v chimerax-daily || command -v chimerax || true)"; fi
            if [ -n "$cxbin" ]; then
                CXLAUNCH="$cxbin"
                PY="$(find_python "$(dirname "$(dirname "$(readlink -f "$cxbin")")")")"
            fi
            ;;
    esac
fi

# --- resolve the lane: marker in this dir or up to 2 parents above (no further).
# ISOLDE's bundle is at <lane>/isolde/isolde (2 up); the others sit directly in
# the lane dir (1 up). No marker => launch ChimeraX normally (standard user dir).
LANE_ROOT=""
dir="$SCRIPT_DIR"
level=0
while :; do
    if [ -e "$dir/.chimerax-lane" ]; then
        LANE_ROOT="$dir/.chimerax-home"
        break
    fi
    [ "$level" -ge 2 ] && break
    parent="$(dirname "$dir")"
    [ "$parent" = "$dir" ] && break    # reached filesystem root
    dir="$parent"
    level=$((level + 1))
done

if [ -n "$LANE_ROOT" ]; then
    # Isolated lane.
    if [ -z "$PY" ] || [ ! -x "$PY" ]; then
        echo "ERROR: could not find the bundled ChimeraX python." >&2
        echo "       Set CHIMERAX_APP to the ChimeraX install dir and retry." >&2
        exit 1
    fi
    SHIM="$SCRIPT_DIR/tools/_isolated_chimerax.py"
    if [ ! -f "$SHIM" ]; then
        echo "ERROR: launch shim not found at $SHIM" >&2
        exit 1
    fi
    echo "[run_chimerax] isolated lane: $LANE_ROOT" >&2
    exec "$PY" -I "$SHIM" "$LANE_ROOT" "$@"
else
    # No marker -> standard, shared ChimeraX user directory (unwrapped launch).
    if [ -z "$CXLAUNCH" ] || [ ! -x "$CXLAUNCH" ]; then
        echo "ERROR: could not find the ChimeraX launcher." >&2
        echo "       Set CHIMERAX_APP to the ChimeraX install dir and retry." >&2
        exit 1
    fi
    echo "[run_chimerax] no .chimerax-lane marker; using the standard ChimeraX user directory." >&2
    exec "$CXLAUNCH" "$@"
fi
