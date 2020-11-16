@set PATH=%PATH%;C:\mingw\mingw64\bin
@echo(
@echo **************** mingw64 ****************
.\bin\Debug\arma.exe

if %1%.==. goto PAUSE
if %1%==nopause goto EXIT

:PAUSE
pause

:EXIT


rem C:\mingw\mingw64\lib\libopenblas.a