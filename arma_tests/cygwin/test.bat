@echo off
@set PATH=%PATH%;C:\cygwin64\bin
@echo(
@echo **************** cygwin64 ****************
.\bin\Debug\arma.exe


if %1%.==. goto PAUSE
if %1%==nopause goto EXIT

:PAUSE
pause

:EXIT
