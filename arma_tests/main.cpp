/*
1. cygwin             0.326226s,0.333001s
2. mingw 0.3.9        0.235554s,0.237938s
2. msys2 mingw 0.3.7  0.23445s, 0.235741s

1. Matrix product in python (using numpy):    c = a*b for the size (2000x2000): 0.25625334978103637/0.26406580209732056 s norm =  0.0
2. Matrix product in matlab:                  c = a*b for the size (2000x2000): 0.285133s  norm =  8.4623e-12
   Intel(R) Math Kernel Library Version 11.2.3 Product Build 20150413 for Intel(R) 64 architecture applications, CNR branch AVX2
3. Matrix product in gcc+openblas (cygwin):   C = A*B for the size (2000x2000): 0.3190s (old 0.374518s 2016y)
   32s with blas, lapack
4. Matrix product in fortran (matmul):: C = A * B for the size (2000 x2000 ):  1.1319999694824219s
   gfortran "$(FILE_NAME)" -o "$(NAME_PART)" -O3 -lblas -llapack -pthread

 */
// C = A*B, where A,B are complex matrices with size (2000x2000): 1.34465s

//#define COMPLEX       /*use complex operations*/
//#define LOADMATRICES   /*load produced matrices for test with matalb*/
//#define SAVEMATRICES   /*save produced matrices for test in matalb*/
//#define TESTMATRIXOPS /*use for test of matrix operation in armadillo in not - only matrix product*/
//#define ARMA_DONT_USE_WRAPPER /*directly with BLAS and LAPACK without using the Armadillo run-time library*/
//#define USE_OPENMP



#include <iostream>
#include <armadillo>
#include <sstream>
#include <unordered_map>
#include <cstdio>
#include <string>
#include "openblas_config.h"
#include "cblas.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif // USE_OPENMP


using namespace std;
using namespace arma;

int
main(int argc, char** argv)
{
    {
    unordered_map<unsigned,string> map{{200505,"2.5"},{200805,"3.0"},{201107,"3.1"},{201307,"4.0"},{201511,"4.5"}};

    cout << "Win" << sizeof(void*)*8 << ":" << endl;


#ifdef USE_OPENMP
    cout << "OpenMP " << map.at(_OPENMP).c_str() << ", max number of threads = " << omp_get_max_threads();
    omp_set_dynamic(0);     // Explicitly disable dynamic teams
    omp_set_num_threads(omp_get_max_threads());
    cout << ", uses "<< omp_get_num_threads() << "thread(s)" << endl;
#endif // USE_OPENMP

    cout << "Armadillo version: " << arma_version::as_string() << endl << endl;
    // Armadillo documentation is available at:
    // http://arma.sourceforge.net/docs.html

//  cout << OPENBLAS_VERSION << " with " << OPENBLAS_GEMM_MULTITHREAD_THRESHOLD << " thread(s)" << endl;
    openblas_set_num_threads(openblas_get_num_procs());
    cout << OPENBLAS_VERSION << " with " << openblas_get_num_threads() << " thread(s) (physical="<< openblas_get_num_procs() <<")" << endl;
    cout << "openblas config: " << openblas_get_config() << endl;
    cout << "openblas parallel: " << openblas_get_parallel() << "(0 - SEQUENTIAL; 1 - THREAD; 2 - OPENMP)" << endl;

    cout << "Date : "<< __DATE__  << endl;
    cout << "Time : "<< __TIME__  << endl;
    cout << "====================================================\n\n" << endl;
    }


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

/*
    arma_rng::set_seed_random();
    mat(3, 3, fill::randu).print();
    cx_mat(3, 3, fill::randu).print("random");
*/

#ifndef TESTMATRIXOPS
    wall_clock timer;

    timer.tic();
    mat C;
#ifdef LOADMATRICES
    mat A,B,C2;
    string path="D:\\_install\\_Programming\\parallel_computing\\_matrix_product_tests\\";
    cout << "Loading matrix A..." << endl;
    A.load(path+"A.dat");
    cout << "Loading matrix B..." << endl;
    B.load(path+"B.dat");
    cout << "Loading matrix C2..." << endl;
    C2.load(path+"C.dat");
    cout << "matrices are loaded in ";
#else // generate matrices
    cout << "Create random matrices: ";
#ifndef COMPLEX
    mat A(n, n, fill::randu), B(n, n, fill::randu);
//    mat A = randu(n, n), B = randu(n, n), C;
#else
    cx_mat A = randu<cx_mat>(n, n), B = randu<cx_mat>(n, n), C;
//    cx_mat A = cx_mat(n, n, fill::randu), B = cx_mat(n, n, fill::randu), C;
#endif // COMPLEX
#endif // LOADMATRICES
    cout << timer.toc() << "s" << endl;


//	cout << "Check one of elements: A(1,3) = " << A(1,3)<< endl;
//	cout << "Check one of elements: B(1,3) = " << B(1,3)<< endl;

#ifndef COMPLEX
    cout << "Matrix product in gcc+openblas: C = A*B for the size ("<<n<<"x"<<n<<"): ";
#else
    cout << "C = A*B, where A,B are complex matrices with the size ("<<n<<"x"<<n<<"): ";
#endif // COMPLEX
    int N = 20;
    timer.tic();
    for(int n=1;n<N;n++){
        C = A*B;
    }
    cout << timer.toc()/double(N) << "s" << endl;

#ifdef LOADMATRICES
    timer.tic();
    double CCnorm = norm(C-C2,2);
    cout << "norm2(C-C2)=" << CCnorm << ": " << timer.toc() << "s" << endl;
#endif // LOADMATRICES

#ifdef SAVEMATRICES
    cout << "Storing matrix A..." << endl;
    A.save(path+"A.dat", raw_ascii);
    cout << "Storing matrix B..." << endl;
    B.save(path+"B.dat", raw_ascii);
    cout << "Storing matrix C..." << endl;
    C.save(path+"C.dat", raw_ascii);
#endif // SAVEMATRICES

#else // TESTMATRIXOPS
    cout << "0. Directly specify the matrix with the size (2,3),  elements are uninitialised" << endl;
    mat A(2,3);  // directly specify the matrix size (elements are uninitialised)

    cout << "A.n_rows: " << A.n_rows << endl;  // .n_rows and .n_cols are read only
    cout << "A.n_cols: " << A.n_cols << endl;

    A(1,2) = 456.0;  // direct access an element (indexing starts at 0)
    cout << endl;
    A.print("1. print content of the matrix A (A(1,2) = 456.0 other are randoms):");
    cout << "A(1,2) = " << A(1,2) << endl;

    cout << endl;
    A = 5.0;         // scalars are treated as a 1x1 matrix
    A.print("2. Reinitialize A as a scalar:");

    cout << endl;
    A.set_size(4,5); // change the size (data is not preserved)
    A.fill(5.0);     // set all elements to a particular value
    A.print("3.  Reinitialize A as a matrix (4,5) and fill by 5:");

    cout << endl;
    // endr indicates "end of row"
    A << 0.165300 << 0.454037 << 0.995795 << 0.124098 << 0.047084 << endr
    << 0.688782 << 0.036549 << 0.552848 << 0.937664 << 0.866401 << endr
    << 0.348740 << 0.479388 << 0.506228 << 0.145673 << 0.491547 << endr
    << 0.148678 << 0.682258 << 0.571154 << 0.874724 << 0.444632 << endr
    << 0.245726 << 0.595218 << 0.409327 << 0.367827 << 0.385736 << endr;
    A.print("4. Fill content of A by the operator '<<':");

    cout << endl;
    // determinant
    cout << "5. det(A): " << det(A) << endl;

    cout << endl;
    // inverse
    cout << "6. inv(A): " << endl << inv(A) << endl;

    cout << endl;
    mat A2 = A*inv(A);
    A2.print("7. check: A*invA");

    cout << endl;
    cout << "8. save/load A to a file 'A0.dat'" << endl;
    // save matrix as a text file
    A.save("A0.dat", raw_ascii);

    // load from file
    mat B;
    B.load("A0.dat");

    cout << endl;
    cout << "9. submatrices" << endl;
    // submatrices
    cout << "B( span(0,2), span(3,4) ):" << endl << B( span(0,2), span(3,4) ) << endl;

    cout << "B( 0,3, size(3,2) ):" << endl << B( 0,3, size(3,2) ) << endl;

    cout << "B.row(0): " << endl << B.row(0) << endl;

    cout << "B.col(1): " << endl << B.col(1) << endl;

    // transpose
    cout << "B.t(): " << endl << B.t() << endl;

    // maximum from each column (traverse along rows)
    cout << endl;
    cout << "10. maximum from each column (traverse along rows)" << endl;
    cout << "max(B): " << endl << max(B) << endl;

    // maximum from each row (traverse along columns)
    cout << "max(B,1): " << endl << max(B,1) << endl;

    // maximum value in B
    cout << "max(max(B)) = " << max(max(B)) << endl;

    // sum of each column (traverse along rows)
    cout << "sum(B): " << endl << sum(B) << endl;

    // sum of each row (traverse along columns)
    cout << "sum(B,1) =" << endl << sum(B,1) << endl;

    // sum of all elements
    cout << "accu(B): " << accu(B) << endl;

    // trace = sum along diagonal
    cout << "trace(B): " << trace(B) << endl;

    // generate the identity matrix
    mat C = eye<mat>(4,4);

    // random matrix with values uniformly distributed in the [0,1] interval
    cout << endl;
    mat D = randu<mat>(4,4);
    D.print("11. D (4,4) with randoms :");

    // row vectors are treated like a matrix with one row
    cout << endl;
    rowvec r;
    r << 0.59119 << 0.77321 << 0.60275 << 0.35887 << 0.51683;
    r.print("12 row vector r:");

    // column vectors are treated like a matrix with one column
    cout << endl;
    vec q;
    q << 0.14333 << 0.59478 << 0.14481 << 0.58558 << 0.60809;
    q.print("13. column vector q:");

    // convert matrix to vector; data in matrices is stored column-by-column
    cout << endl;
    vec v = vectorise(A);
    v.print("14 v=A(:):");

    // dot or inner product
    cout << endl;
    cout << "14. as_scalar(r*q): " << as_scalar(r*q) << endl;

    // outer product
    cout << endl;
    cout << "15. q*r: " << endl << q*r << endl;

    // multiply-and-accumulate operation (no temporary matrices are created)
    cout << endl;
    cout << "16. multiply-and-accumulate operation (no temporary matrices are created)" << endl;
    cout << "accu(A % B) = " << accu(A % B) << endl;

    cout << endl;
    // example of a compound operation
    B += 2.0 * A.t();
    B.print("17. example of a compound operation B=2.0 * A.t():");

    // imat specifies an integer matrix
    imat AA;
    imat BB;

    AA << 1 << 2 << 3 << endr << 4 << 5 << 6 << endr << 7 << 8 << 9;
    BB << 3 << 2 << 1 << endr << 6 << 5 << 4 << endr << 9 << 8 << 7;

    // comparison of matrices (element-wise); output of a relational operator is a umat
    cout << endl;
    umat ZZ = (AA >= BB);
    ZZ.print("18. comparison of matrices (element-wise) ZZ = (AA >= BB):");

    // cubes ("3D matrices")
    cout << endl;
    cube Q( B.n_rows, B.n_cols, 2 );

    Q.slice(0) = B;
    Q.slice(1) = 2.0 * B;

    Q.print("19. cubes (3D matrices) Q:");

    // 2D field of matrices; 3D fields are also supported
    cout << endl;
    field<mat> F(4,3);

    for(uword col=0; col < F.n_cols; ++col)
    for(uword row=0; row < F.n_rows; ++row){
        F(row,col) = randu<mat>(2,3);  // each element in field<mat> is a matrix
    }

    F.print("20. each element in F is a matrix:");

#endif // TESTMATRIXOPS
}

