@echo off
@set PATH=%PATH%;C:\msys64\mingw64\bin
@echo(
@echo **************** mingw64 in msys64 ****************
.\bin\Debug\arma.exe

if %1%.==. goto PAUSE
if %1%==nopause goto EXIT

:PAUSE
pause

:EXIT
