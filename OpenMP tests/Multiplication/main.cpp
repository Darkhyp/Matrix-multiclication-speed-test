// g++-4.8 -ansi -O2 -Wall -pedantic -pthread -fopenmp -x c main.cpp && ./a.out

#include <SDKDDKVer.h>
#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <random>
//#include <functional>
#include <sstream>


using namespace std;

void transpose(double *A, double *B, int n) {
    int i,j;
	for(i=0; i<n; i++) {
		for(j=0; j<n; j++) {
			B[j*n+i] = A[i*n+j];
		}
	}
}
void transpose_omp(double *A, double *B, int n) {
	#pragma omp parallel
	{
        int i,j;
        #pragma omp for private(i, j)
        for(i=0; i<n; i++) {
            for(j=0; j<n; j++) {
                B[j*n+i] = A[i*n+j];
            }
        }
	}
}

void gemm(double *A, double *B, double *C, int n)
{
	int i, j, k;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			double dot  = 0;
			for (k = 0; k < n; k++) {
				dot += A[i*n+k]*B[k*n+j];
			}
			C[i*n+j ] = dot;
		}
	}
}

void gemm_omp(double *A, double *B, double *C, int n)
{
	#pragma omp parallel
	{
//        printf("number of threads = %i\n\n",omp_get_num_threads());
		int i, j, k;
//		#pragma omp for
        #pragma omp for private(i, j, k)
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				double dot  = 0;
				for (k = 0; k < n; k++) {
					dot += A[i*n+k]*B[k*n+j];
				}
				C[i*n+j ] = dot;
			}
		}

	}
}

void gemmT(double *A, double *B, double *C, int n)
{
	int i, j, k;
	double *B2;
	B2 = (double*)malloc(sizeof(double)*n*n);
    transpose(B,B2, n);
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			double dot  = 0;
			for (k = 0; k < n; k++) {
				dot += A[i*n+k]*B2[j*n+k];
			}
			C[i*n+j ] = dot;
		}
	}
	free(B2);
}

void gemmT_omp(double *A, double *B, double *C, int n)
{
	double *B2;
	B2 = (double*)malloc(sizeof(double)*n*n);
    transpose_omp(B,B2, n);

	#pragma omp parallel
	{
//        printf("number of threads = %i\n\n",omp_get_num_threads());
		int i, j, k;
//		#pragma omp for
        #pragma omp for private(i, j, k)
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				double dot  = 0;
				for (k = 0; k < n; k++) {
					dot += A[i*n+k]*B2[j*n+k];
				}
				C[i*n+j ] = dot;
			}
		}

	}
	free(B2);
}

int main(int argc, char *argv[])
{
    cout << "Win" << sizeof(void*)*8 << ":" << endl;
    printf("\nDate : %s",__DATE__);
    printf("\nTime : %s",__TIME__);
    printf("\n====================================================\n\n");

    int n = 2000;
    if (argc >= 2)
    {
        istringstream iss( argv[1] );

        if (iss >> n)
        {
            // Conversion successful
        }
        else
        {
            n = 512;
            n = 1000;
            //	n = 2000;
        }
    }

	double dtime[4];

	double *A, *B, *C;



	printf("C = A*B, where A,B are matrices with size (%ix%i): \n\n",n,n);

	printf("max number of threads = %i\n",omp_get_max_threads());
    omp_set_num_threads(omp_get_max_threads());


	dtime[0] = omp_get_wtime();
	A = (double*)malloc(sizeof(double)*n*n);
	B = (double*)malloc(sizeof(double)*n*n);
	C = (double*)malloc(sizeof(double)*n*n);

/*
    normal_distribution<> normal(10.0, 3.0); // mean 10, sigma 3
    random_device rd;
    mt19937 engine(rd()); // knuth_b fails in MSVC2010, but compiles in GCC
    function<double()> rnd = std::bind(normal, engine);
//  function rnd2 = std::bind(normal, engine);

    cout << rnd() << endl;
//  cout << rnd2() << endl;
*/

    long maxRandom = 100000000000000;
    random_device rd;   // non-deterministic generator
    mt19937 gen(rd());  // to seed mersenne twister. replace the call to rd() with a constant value to get repeatable results.
    uniform_int_distribution<> dist(0,maxRandom); // distribute results between 0 and maxRandom inclusive.
    /*
    default_random_engine The default engine.
    typedef mt19937 default_random_engine;

    knuth_b Knuth engine.
    typedef shuffle_order_engine<minstd_rand0, 256> knuth_b;

    minstd_rand0 1988 minimal standard engine (Lewis, Goodman, and Miller, 1969).
    typedef linear_congruential_engine<unsigned int, 16807, 0, 2147483647> minstd_rand0;

    minstd_rand Updated minimal standard engine minstd_rand0 (Park, Miller, and Stockmeyer, 1993).
    typedef linear_congruential_engine<unsigned int, 48271, 0, 2147483647> minstd_rand;

    mt19937 32-bit Mersenne twister engine (Matsumoto and Nishimura, 1998).
    typedef mersenne_twister_engine<unsigned int, 32, 624, 397, 31, 0x9908b0df, 11, 0xffffffff, 7, 0x9d2c5680, 15, 0xefc60000, 18, 1812433253> mt19937;

    mt19937_64 64-bit Mersenne twister engine (Matsumoto and Nishimura, 2000).
    typedef mersenne_twister_engine<unsigned long long, 64, 312, 156, 31, 0xb5026f5aa96619e9ULL, 29, 0x5555555555555555ULL, 17, 0x71d67fffeda60000ULL, 37, 0xfff7eee000000000ULL, 43, 6364136223846793005ULL> mt19937_64;

    ranlux24 24-bit RANLUX engine (Martin Lüscher and Fred James, 1994).
    typedef discard_block_engine<ranlux24_base, 223, 23> ranlux24;

    ranlux24_base Used as a base for ranlux24.
    typedef subtract_with_carry_engine<unsigned int, 24, 10, 24> ranlux24_base;

    ranlux48 48-bit RANLUX engine (Martin Lüscher and Fred James, 1994).
    typedef discard_block_engine<ranlux48_base, 389, 11> ranlux48;

    ranlux48_base Used as a base for ranlux48.
    typedef subtract_with_carry_engine<unsigned long long, 48, 5, 12> ranlux48_base;
    */

    #pragma omp parallel
    {
        int i;
        #pragma omp for private(i)

        for(i=0; i<n*n; i++) {
    /*
            A[i] = 1.0*rand()/RAND_MAX;
            B[i] = 1.0*rand()/RAND_MAX;
    */
            A[i] = (double)dist(gen)/(double)maxRandom;
            B[i] = (double)dist(gen)/(double)maxRandom;
        }
    }
	dtime[0] = omp_get_wtime() - dtime[0];
	printf("Create random matrices: \t = %f s\n", dtime[0]);

//	printf("RAND_MAX \t = %g, A[3] = %g \n", RAND_MAX,A[3]);
	printf("Check one of elements: A[3] = %g \n", A[3]);

	dtime[0] = omp_get_wtime();
	printf("1. no transpose  no openmp ...\n");
	gemm(A,B,C, n);
	dtime[0] = omp_get_wtime() - dtime[0];

	dtime[1] = omp_get_wtime();
	printf("2. no transpose, openmp ...\n");
	gemm_omp(A,B,C, n);
	dtime[1] = omp_get_wtime() - dtime[1];

	dtime[2] = omp_get_wtime();
	printf("3. transpose, no openmp ...\n");
	gemmT(A,B,C, n);
	dtime[2] = omp_get_wtime() - dtime[2];

	dtime[3] = omp_get_wtime();
	printf("4. transpose, openmp ...\n");
	gemmT_omp(A,B,C, n);
	dtime[3] = omp_get_wtime() - dtime[3];



	printf("\n\n");
	printf("no transpose  no openmp \t(0.94s) = %f s\n", dtime[0]);
	printf("no transpose, openmp \t\t(0.23s) = %f s\n", dtime[1]);
	printf("transpose, no openmp \t\t(0.27s) = %f s\n", dtime[2]);
	printf("transpose, openmp \t\t(0.08s) = %f s\n", dtime[3]);

    return 0;

}
