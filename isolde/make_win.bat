@echo off
setlocal EnableExtensions EnableDelayedExpansion

REM ChimeraX 1.13 ships Python 3.14, built with VS2022 (MSC v.1944), and its own
REM C++ libraries (atomstruct, arrays, ...) are built with the v143 toolset, so
REM build with v143 too.  We do that by NOT passing -vcvars_ver: vcvars64 then
REM selects its default, i.e. the newest toolset installed with VS2022.
REM
REM Rather than pin a toolset version, the "MSVC toolset check" block below
REM verifies whatever vcvars picked against the ChimeraX actually being targeted.
REM A pin cannot express "match ChimeraX": UCSF does not publish the toolset they
REM build with, so the only ground truth is that install, readable at build time.
REM The old pin here was -vcvars_ver=14.1 (v141/VS2017) "to match ChimeraX", and
REM being a hard-coded value is exactly why it silently outlived ChimeraX's move
REM from Python 3.11 to 3.14.
REM
REM Because a .bat always runs under cmd.exe, this works whether invoked from cmd
REM or PowerShell, so a pre-opened vcvars terminal is no longer required.  Guarded
REM on VCINSTALLDIR so running from an existing vcvars terminal is a harmless
REM no-op.
REM
REM Delayed expansion (!VSWHERE! / !VSPATH!) is used inside the if-block because
REM the resolved paths contain "(x86)" -- a literal ")" that would close the
REM block prematurely if expanded with %...% at parse time.
set "VSWHERE=%ProgramFiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe"
if not defined VCINSTALLDIR (
    if not exist "!VSWHERE!" (
        echo ERROR: vswhere.exe not found - is Visual Studio 2022 installed?
        exit /b 1
    )
    set "VSPATH="
    for /f "usebackq delims=" %%i in (`"!VSWHERE!" -latest -products * -property installationPath`) do set "VSPATH=%%i"
    if not defined VSPATH (
        echo ERROR: no Visual Studio installation found by vswhere.
        exit /b 1
    )
    call "!VSPATH!\VC\Auxiliary\Build\vcvars64.bat"
)
for /f "usebackq tokens=*" %%v in (`cl 2^>^&1 ^| findstr /C:"Optimizing Compiler"`) do echo Building with: %%v
set DISTUTILS_USE_SDK=1
set MSSdk=1

REM ChimeraX launches go through run_chimerax.bat so each worktree/lane installs
REM into its own isolated ChimeraX user directory (see run_chimerax.bat). The
REM "release" token (stable vs daily install) is forwarded to it; "clean" and
REM "app-install" are consumed here.
set "RELEASE_TOK="
set DO_CLEAN=
set DO_INSTALL=
for %%A in (%*) do (
    if /i "%%A"=="release"     set "RELEASE_TOK=release"
    if /i "%%A"=="clean"       set DO_CLEAN=1
    if /i "%%A"=="app-install" set DO_INSTALL=1
)

REM --- MSVC toolset check ----------------------------------------------------
REM MSVC guarantees binary compatibility across VS2015-VS2022 only when the
REM toolset doing the LINKING is at least as new as the newest toolset that built
REM any input.  ChimeraX ships prebuilt C++ libraries and a Python built with its
REM own toolset, so an older one here is unsupported even when it appears to work.
REM Reading ChimeraX's compiler costs ~0.3s: run_chimerax.bat routes "-c" to
REM ChimeraX-console.exe and keeps its own notices on stderr, so stdout is just
REM the number.  RELEASE_TOK must already be parsed so we query the right install.
set "CX_MSC="
for /f "usebackq tokens=*" %%v in (`call "%~dp0run_chimerax.bat" %RELEASE_TOK% -c "import platform; print(platform.python_compiler().split()[1][2:])" 2^>nul`) do set "CX_MSC=%%v"
set "VCMIN="
for /f "tokens=1,2 delims=." %%a in ("%VCToolsVersion%") do if /i "%%a"=="14" set "VCMIN=%%b"
if not defined CX_MSC (
    echo WARNING: could not read the compiler ChimeraX was built with; skipping toolset check.
) else if not defined VCMIN (
    echo WARNING: unrecognised VCToolsVersion "%VCToolsVersion%"; skipping toolset check.
) else (
    REM 1900+1<minor>-100, not 1900+<minor>: set /a reads a leading zero as octal.
    set /a ACTIVE_MSC=1900+1!VCMIN!-100
    if !ACTIVE_MSC! LSS !CX_MSC! (
        echo ERROR: the active MSVC toolset ^(VCToolsVersion %VCToolsVersion%, _MSC_VER !ACTIVE_MSC!^)
        echo        is OLDER than the one ChimeraX was built with ^(_MSC_VER !CX_MSC!^).
        echo        Linking against ChimeraX's libraries with it is unsupported.
        echo        Install a newer MSVC toolset, or drop any -vcvars_ver pin above.
        exit /b 1
    )
    if !ACTIVE_MSC! GTR !CX_MSC! (
        echo WARNING: MSVC toolset _MSC_VER !ACTIVE_MSC! is NEWER than ChimeraX's !CX_MSC!.
        echo          That is the safe direction within VS2015-VS2022, but rebuild with a
        echo          matching toolset if these ever straddle a major ABI change.
    )
    if !ACTIVE_MSC! EQU !CX_MSC! echo MSVC toolset _MSC_VER !ACTIVE_MSC! matches ChimeraX.
)

REM Regenerate pyproject.toml from pyproject.toml.in (substitutes OpenMM lib/include dirs).
call "%~dp0run_chimerax.bat" %RELEASE_TOK% --nogui --safemode --exit --script prep_toml.py

if defined DO_CLEAN   call "%~dp0run_chimerax.bat" %RELEASE_TOK% --nogui --safemode --exit --cmd "devel clean ."
if defined DO_INSTALL call "%~dp0run_chimerax.bat" %RELEASE_TOK% --nogui --safemode --exit --cmd "devel install ."

endlocal
