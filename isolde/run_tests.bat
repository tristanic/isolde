@echo off
setlocal EnableExtensions EnableDelayedExpansion
REM ===========================================================================
REM run_tests.bat -- run every self-contained ISOLDE test under src\tests\
REM headlessly (via run_chimerax.bat) and report pass/fail per file.
REM
REM A "self-contained" test follows the repo convention: it runs itself when
REM ChimeraX injects `session` (the `if session is not None: run(session)`
REM footer) and prints "ALL PASS" on success / calls _fail() (which raises
REM SystemExit) otherwise. See test_mmff_parameterisation.py and
REM test_atom_deletion_robustness.py for the pattern.
REM
REM Files without an "ALL PASS" sentinel are SKIPPED, not failed -- e.g.
REM test_simulation.py is the interactive SimTester harness (CLAUDE.md: "run
REM manually inside ChimeraX"), not a headless pass/fail test. Any future test
REM that adopts the convention is picked up automatically.
REM
REM Success is judged by the presence of "ALL PASS" in the captured output, not
REM by exit code: a SystemExit raised inside a ChimeraX --script does not
REM reliably surface as a process exit code.
REM
REM Usage:  run_tests.bat [release]
REM   release -> forwarded to run_chimerax.bat to target the stable install
REM              (omitted -> the daily build), matching run_chimerax.bat's own
REM              argument convention.
REM ===========================================================================

set "HERE=%~dp0"
if "%HERE:~-1%"=="\" set "HERE=%HERE:~0,-1%"

set "RELEASE="
if /i "%~1"=="release" set "RELEASE=release"

set "TESTDIR=%HERE%\src\tests"
set "OUT=%TEMP%\isolde_test_out.txt"
set /a NPASS=0, NFAIL=0, NSKIP=0
set "FAILED="

REM run from the repo root so the forward-slash --script path resolves the same
REM way it does when invoked by hand.
pushd "%HERE%"

echo ===========================================================================
echo Running ISOLDE tests in %TESTDIR%
echo ===========================================================================

for %%F in ("%TESTDIR%\test_*.py") do (
    findstr /c:"ALL PASS" "%%F" >nul 2>&1
    if errorlevel 1 (
        echo [SKIP] %%~nxF  ^(no "ALL PASS" sentinel: not a self-running test^)
        set /a NSKIP+=1
    ) else (
        echo.
        echo --- %%~nxF ---
        call "%HERE%\run_chimerax.bat" !RELEASE! --nogui --exit --script "src/tests/%%~nxF" > "%OUT%" 2>&1
        type "%OUT%"
        findstr /c:"ALL PASS" "%OUT%" >nul 2>&1
        if errorlevel 1 (
            echo [FAIL] %%~nxF
            set /a NFAIL+=1
            set "FAILED=!FAILED! %%~nxF"
        ) else (
            echo [PASS] %%~nxF
            set /a NPASS+=1
        )
    )
)

echo.
echo ===========================================================================
echo Summary: !NPASS! passed, !NFAIL! failed, !NSKIP! skipped
if not "!FAILED!"=="" echo Failed:!FAILED!
echo ===========================================================================

popd
if !NFAIL! gtr 0 exit /b 1
exit /b 0
