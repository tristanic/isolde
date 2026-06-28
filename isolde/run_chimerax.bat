@echo off
setlocal EnableExtensions EnableDelayedExpansion
REM ===========================================================================
REM run_chimerax.bat -- launch ChimeraX with a per-lane isolated user directory.
REM
REM ChimeraX installs every user bundle into a single shared per-user directory,
REM so parallel development worktrees clobber each other at "devel install" time.
REM This wrapper bypasses the hardcoded "python -I" ChimeraX launcher by calling
REM the bundled bin\python.exe directly on tools\_isolated_chimerax.py, which
REM re-roots ChimeraX's whole user tree under a lane-specific directory.
REM
REM Usage:   run_chimerax.bat [release] [<ChimeraX args>...]
REM   release   -> use the stable install (C:\Program Files\ChimeraX);
REM               omitted -> the daily build (C:\Program Files\ChimeraX-Daily).
REM   All remaining args are forwarded verbatim to ChimeraX.
REM
REM Lane resolution: walks up from this script's directory for a ".chimerax-lane"
REM marker. If found at <dir>, the isolated tree lives at <dir>\.chimerax-home --
REM so every bundle worktree under that lane shares one ChimeraX user dir. With no
REM marker, falls back to a unique per-worktree home (collision-safe) and warns.
REM ===========================================================================

REM --- Split optional leading "release" from the rest, preserving quoting ---
REM Pure string substitution (no for/f): for/f's default eol=";" truncates
REM commands like --cmd "devel install . ; exit" and leaks the inner quotes.
set "RELEASE="
set "CXARGS=%*"
if /i "%~1"=="release" (
    set "RELEASE=1"
    set "CXARGS=!CXARGS:*release=!"
)

REM --- Locate the ChimeraX install (for its bundled python.exe) ---
if defined RELEASE (
    set "CX_APP=C:\Program Files\ChimeraX"
) else (
    set "CX_APP=C:\Program Files\ChimeraX-Daily"
)
set "PYEXE=%CX_APP%\bin\python.exe"
if not exist "%PYEXE%" (
    echo ERROR: ChimeraX python not found at "%PYEXE%". 1>&2
    echo        ^(pass "release" for the stable install, or fix CX_APP in this script^) 1>&2
    exit /b 1
)

REM --- Resolve the lane: look for a .chimerax-lane marker in this script's
REM directory and up to 2 parents above it. Our repo layout never needs more:
REM ISOLDE's bundle sits at <lane>\isolde\isolde (2 up); the others sit directly
REM in the lane dir (1 up). No marker => launch ChimeraX normally, into the
REM standard (shared) user directory -- exactly as an unwrapped ChimeraX would.
set "DIR=%~dp0"
if "%DIR:~-1%"=="\" set "DIR=%DIR:~0,-1%"
set "LANE_ROOT="
set /a LEVEL=0
:walkup
if exist "%DIR%\.chimerax-lane" (
    set "LANE_ROOT=%DIR%\.chimerax-home"
    goto isolated
)
if !LEVEL! geq 2 goto standard
REM step up one directory (stop early if we somehow reach a drive root: a bare
REM drive has no backslash -- substring test, not find, which is unreliable here)
set "NOBS=!DIR:\=!"
if "!NOBS!"=="!DIR!" goto standard
for %%I in ("%DIR%") do set "PARENT=%%~dpI"
if "%PARENT:~-1%"=="\" set "PARENT=%PARENT:~0,-1%"
set "DIR=%PARENT%"
set /a LEVEL+=1
goto walkup

:isolated
set "SHIM=%~dp0tools\_isolated_chimerax.py"
if not exist "%SHIM%" (
    echo ERROR: launch shim not found at "%SHIM%". 1>&2
    exit /b 1
)
echo [run_chimerax] isolated lane: %LANE_ROOT% 1>&2
"%PYEXE%" -I "%SHIM%" "%LANE_ROOT%" %CXARGS%
exit /b %ERRORLEVEL%

:standard
REM No lane marker -> standard shared user dir. Use the console launcher for
REM headless / "-m module" / "-c command" runs (so stdout shows in the calling
REM console, e.g. make_docs' "-m sphinx"); the GUI launcher only for an
REM interactive session. (-m/-c must lead, so a plain substring test is safe.)
set "STDEXE=%CX_APP%\bin\ChimeraX.exe"
if not "!CXARGS:--nogui=!"=="!CXARGS!" set "STDEXE=%CX_APP%\bin\ChimeraX-console.exe"
if not "!CXARGS:-m =!"=="!CXARGS!"     set "STDEXE=%CX_APP%\bin\ChimeraX-console.exe"
if not "!CXARGS:-c =!"=="!CXARGS!"     set "STDEXE=%CX_APP%\bin\ChimeraX-console.exe"
echo [run_chimerax] no .chimerax-lane marker; using the standard ChimeraX user directory. 1>&2
"%STDEXE%" %CXARGS%
exit /b %ERRORLEVEL%
