@echo off
cd D:\_install\_Programming\parallel_computing\_matrix_product_tests\arma_tests


set openblaspath=D:\_install\_Programming\parallel_computing\OpenBLAS\0.3.9 openmpcyg\64
set openblaslib=%openblaspath%\lib\cygopenblas_v0.3.9-gcc_9_3_0.a
rem set openblaspath=D:\_install\_Programming\parallel_computing\OpenBLAS\0.3.9 cyg\64
rem set openblaslib=%openblaspath%"\lib\cygopenblas_v0.3.9-gcc_9_3_0.a
rem set openblaspath=D:\_install\_Programming\parallel_computing\OpenBLAS\0.2.18
rem set openblaslib=%openblaspath%"\lib\libopenblas_haswellp-r0.2.18.a

set arma=D:\_install\_Programming\parallel_computing\armadillo-9.860.1
set comppath=C:\cygwin64
set comppath2=C:\cygwin64\lib\gcc\x86_64-pc-cygwin\9.3.0
set addlibs="%comppath2%\libgomp.a" "C:\Program Files\gnuplot\bin\lua53.dll"
set sigpack=D:\_install\_Programming\parallel_computing\sigpack-1.0.6


set pathold=%PATH%
set PATH="%comppath%\lib";"%comppath%\bin";"%comppath%\include";"%comppath2%"
echo %PATH%

echo compiler position:
which g++.exe
echo(
echo compiler version:
g++.exe -v

:compile
@echo on
@echo(
@echo compiling...
g++.exe -c main.cpp -o main_cyg.o
rem g++.exe -Wall -fexceptions -std=gnu++17 -m64 -g -O3 -static -static-libgcc -fpermissive -lopenblas -llapack -lgfortran -lquadmath -DUSE_OPENMP -I"%openblaspath%\include" -I"%arma%\include" -I"%comppath%\include" -I"%comppath%\usr\include" -I"%comppath2%\include" -c main.cpp -o main_cyg.o
rem g++.exe -Wall -fexceptions -std=gnu++17 -m64 -g -O3 -static -static-libgcc -fpermissive -lopenblas -llapack -lgfortran -lquadmath -DUSE_OPENMP -I"%openblaspath%\include" -I"%arma%\include" -I"%comppath%\include" -I"%comppath%\usr\include" -I"%comppath2%\include" -I"%sigpack%\source" -c main.cpp -o main_cyg.o

:link
@echo(
@echo linking...
g++.exe -L"%openblaspath%\lib" -L"%comppath%\lib" -L"%comppath2%" -L"C:\Program Files\gnuplot\bin" -o arma_cyg.exe main_cyg.o  -m64  "%openblaslib%" -lgfortran -lquadmath %addlibs%

set PATH=%pathold%
