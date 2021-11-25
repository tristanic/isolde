
@echo off
SET CHIMERAX_EXE="c:\Program Files\ChimeraX-Daily\bin\ChimeraX-console.exe"
FOR %%A in (%*) DO (

IF "%%A" == "release" set CHIMERAX_EXE="c:\Program Files\ChimeraX\bin\ChimeraX-console.exe"
)

for %%A in (%*) DO (
IF "%%A" == "clean" (
	@echo on
	%CHIMERAX_EXE% --nogui --safemode --exit --cmd "devel clean ."
	BREAK
)
)

for %%A in (%*) DO (
IF "%%A" == "app-install" (
	@echo on
	%CHIMERAX_EXE% -m PyQt5.pyrcc_main -o src/resources/resources_rc.py src/resources/resources.qrc
	%CHIMERAX_EXE% --nogui --safemode --exit --cmd "devel install ."
	BREAK
)
)
