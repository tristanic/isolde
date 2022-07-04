REM Force use of the same compiler as used to build ChimeraX
call "%VS140COMNTOOLS%"\vsvars32.bat

REM Believe it or not, vcvarsall.bat fails when called with absolute path
REM if there is a space in the path...
for /f %%i in ('cd') do set CURDIR=%%i
cd "%VS140COMNTOOLS%"\..\..\VC\
call vcvarsall.bat amd64 8.1
cd %CURDIR%
set DISTUTILS_USE_SDK=1
set MSSdk=1
