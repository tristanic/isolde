@echo off
SET CHIMERAX_EXE="c:\Program Files\ChimeraX-Daily\bin\ChimeraX-console.exe"
SET DO_CLEAN=
FOR %%A in (%*) DO (
IF "%%A" == "release" set CHIMERAX_EXE="c:\Program Files\ChimeraX\bin\ChimeraX-console.exe"
IF "%%A" == "clean" set DO_CLEAN=1
)
@echo on
IF DEFINED DO_CLEAN (
%CHIMERAX_EXE% --nogui --safemode --exit --script clean_docs.py
) ELSE (
%CHIMERAX_EXE% -m sphinx docs/source src/docs/user
)
