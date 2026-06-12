@echo off
setlocal EnableExtensions EnableDelayedExpansion

REM ChimeraX is built with the v141 (VS2017) toolset; initialise the matching
REM Visual Studio environment if it isn't already active.  Because a .bat always
REM runs under cmd.exe, this works whether invoked from cmd or PowerShell, so a
REM pre-opened vcvars terminal is no longer required.  Guarded on VCINSTALLDIR so
REM running from an existing vcvars terminal is a harmless no-op.
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
    call "!VSPATH!\VC\Auxiliary\Build\vcvars64.bat" -vcvars_ver=14.1
)
set DISTUTILS_USE_SDK=1
set MSSdk=1

REM Default to the ChimeraX daily build; "release" selects the stable install.
set CHIMERAX_EXE="C:\Program Files\ChimeraX-Daily\bin\ChimeraX-console.exe"
set DO_CLEAN=
set DO_INSTALL=
for %%A in (%*) do (
    if /i "%%A"=="release"     set CHIMERAX_EXE="C:\Program Files\ChimeraX\bin\ChimeraX-console.exe"
    if /i "%%A"=="clean"       set DO_CLEAN=1
    if /i "%%A"=="app-install" set DO_INSTALL=1
)

REM Regenerate pyproject.toml from pyproject.toml.in (substitutes OpenMM lib/include dirs).
%CHIMERAX_EXE% --nogui --safemode --exit --script prep_toml.py

if defined DO_CLEAN   %CHIMERAX_EXE% --nogui --safemode --exit --cmd "devel clean ."
if defined DO_INSTALL %CHIMERAX_EXE% --nogui --safemode --exit --cmd "devel install ."

endlocal
