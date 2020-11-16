/*

In addition. "Z boson", I have tested your C code on the laptop with intel i5 (2 physical cores or 4 logical ones). Unfortunately, the calculation speed is not very fast. For 2000x2000 random double matrices I obtained the following results (using VS 2010 with OpenMP 2.0):

Compiled for Win64: C = A*B, where A,B are matrices with the size (2000x2000):

max number of threads = 4
Create random matrices:      = 0.303555 s
no transpose  no openmp = 100.539924 s
no transpose, openmp = 47.876084 s
transpose, no openmp = 27.872169 s
transpose, openmp = 15.821010 s

Compiled for Win32: C = A*B, where A,B are matrices with the size (2000x2000):

max number of threads = 4
Create random matrices:      = 0.378804 s
no transpose  no openmp = 98.613992 s
no transpose, openmp = 48.233655 s
transpose, no openmp = 29.590350 s
transpose, openmp = 13.678097 s

Note that for the "Hynek Blaha" code the calculation time on my system is 739.208s (226.62s with openMP)!

Whereas in Matlab x64:

n = 2000;
A = rand(n); B = rand(n);

tic
C = A*B;
toc

the calculation time is 0.591440 seconds.

But using openBLAS package I reached a speed of 0.377814 seconds (using minGW with openMP 4.0). The Armadillo package provides a simple way (in my opinion) for connection of matrix operations with openBLAS (or other similar packages). In this case the code is

*/
#include <iostream>
#include <armadillo>
using namespace std;
using namespace arma;

int main(){
    int N,n = 2000;
    N = 10; // number of repetitions
    wall_clock timer;

    arma_rng::set_seed_random();

    mat A(n, n, fill::randu), B(n, n, fill::randu);

    timer.tic();
    // repeat simulation N times
    for(int n=1;n<N;n++){
      mat C = A*B;
    }
    cout << timer.toc()/double(N) << "s" << endl;

    return 0;
}

