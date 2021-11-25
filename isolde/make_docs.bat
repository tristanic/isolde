@echo off
SET CHIMERAX_EXE="c:\Program Files\ChimeraX-Daily\bin\ChimeraX-console.exe"
FOR %%A in (%*) DO (

IF "%%A" == "release" set CHIMERAX_EXE="c:\Program Files\ChimeraX\bin\ChimeraX-console.exe"
)
@echo on
%CHIMERAX_EXE% -m sphinx docs/source src/docs/user
