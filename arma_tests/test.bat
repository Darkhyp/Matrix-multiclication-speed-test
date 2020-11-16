@echo off
cd cygwin
call test.bat nopause
cd../mingw
call test.bat nopause
cd../msys
call test.bat nopause
cd ..

pause
