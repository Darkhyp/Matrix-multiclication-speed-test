#include <SDKDDKVer.h>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <sstream>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <omp.h>
#include <random>


using namespace std;

struct Result {
    vector< vector<int> > A;
    vector< vector<int> > B;
};

/*
Result read(string filename) {
    vector< vector<int> > A, B;
    Result ab;
    string line;
    ifstream infile;
    infile.open (filename.c_str());

    int i = 0;
    while (getline(infile, line) &amp;&amp; !line.empty()) {
        istringstream iss(line);
        A.resize(A.size() + 1);
        int a, j = 0;
        while (iss >> a) {
            A[i].push_back(a);
            j++;
        }
        i++;
    }

    i = 0;
    while (getline(infile, line)) {
        istringstream iss(line);
        B.resize(B.size() + 1);
        int a;
        int j = 0;
        while (iss >> a) {
            B[i].push_back(a);
            j++;
        }
        i++;
    }

    infile.close();
    ab.A = A;
    ab.B = B;
    return ab;
}
*/

template<typename T> vector< vector<T> > ijkalgorithm(vector< vector<T> > A, vector< vector<T> > B) {
    int n = A.size();
    T Ctmp;

    // initialise transpose matrix B
    vector< vector<T> > BT(n, vector<T>(n));
    {
        int i,j;
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                BT[i][j] = B[j][i];
            }
        }
    }

    // initialise C with 0s
//    vector<T> tmp(n, 0);
//    vector< vector<T> > C(n, tmp);
    vector< vector<T> > C(n, vector<T>(n));

    {
        int i,j,k;
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                Ctmp = 0;
//                for (k = 0; k < n; k++) Ctmp += A[i][k] * B[k][j];
                for (k = 0; k < n; k++) Ctmp += A[i][k] * BT[j][k];
                C[i][j] = Ctmp;
            }
        }
    }
    return C;
}
template<typename T> void ijkalgorithm_omp(vector< vector<T> > A, vector< vector<T> > B, vector< vector<T> > C) {
    int n = A.size();
    T Ctmp;

    // initialise transpose matrix B
    vector< vector<T> > BT(n, vector<T>(n));
    #pragma omp parallel
    {
        int i,j;
        #pragma omp for private(i,j)
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                BT[i][j] = B[j][i];
            }
        }
    }

    #pragma omp parallel
    {
        int i,j,k;
        #pragma omp for private(i,j,k)
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                Ctmp = 0;
//                for (k = 0; k < n; k++) Ctmp += A[i][k] * B[k][j];
                for (k = 0; k < n; k++) Ctmp += A[i][k] * BT[j][k];
                C[i][j] = Ctmp;
            }
        }
    }
}

template<typename T> vector< vector<T> > ikjalgorithm(vector< vector<T> > A, vector< vector<T> > B) {
    size_t n = A.size();

    // initialise C with 0s
    vector<T> tmp(n, 0);
    vector< vector<T> > C(n, tmp);

    for (size_t i = 0; i < n; i++) {
        for (size_t k = 0; k < n; k++) {
            for (size_t j = 0; j < n; j++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}
template<typename T> vector< vector<T> > ikjalgorithm_omp(vector< vector<T> > A, vector< vector<T> > B) {
    int n = A.size();

    // initialise C with 0s
    vector<T> tmp(n, 0);
    vector< vector<T> > C(n, tmp);

    #pragma omp parallel
    {
        int i,j,k;
        #pragma omp for private(i,j,k)
        for (i = 0; i < n; i++) {
            for (k = 0; k < n; k++) {
                for (j = 0; j < n; j++) {
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }
    }
    return C;
}

void printMatrix(vector< vector<int> > matrix) {
    vector< vector<int> >::iterator it;
    vector<int>::iterator inner;
    for (it=matrix.begin(); it != matrix.end(); it++) {
        for (inner = it->begin(); inner != it->end(); inner++) {
            cout << *inner;
            if(inner+1 != it->end()) {
                cout << "\t";
            }
        }
        cout << endl;
    }
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

	printf("C = A*B, where A,B are matrices with size (%ix%i): \n\n",n,n);

	printf("max number of threads = %i\n",omp_get_max_threads());
    omp_set_num_threads(omp_get_max_threads());

	dtime[0] = omp_get_wtime();
    vector< vector<double> > A,B,C;
    A.resize(n);
    B.resize(n);

    long maxRandom = 10000000000000;
    random_device rd;   // non-deterministic generator
    mt19937 gen(rd());  // to seed mersenne twister. replace the call to rd() with a constant value to get repeatable results.
    uniform_int_distribution<> dist(0,maxRandom); // distribute results between 0 and maxRandom inclusive.

    cout<<"Initializing random matrices (A,B) ...";
    #pragma omp parallel
    {
        int i,j;
        #pragma omp for private(i,j)


        for(i=0; i<n; i++) {
            for(j=0; j<n; j++){
                A[i].push_back((double)dist(gen)/(double)maxRandom);
                B[i].push_back((double)dist(gen)/(double)maxRandom);
            }
        }
    }
	dtime[0] = omp_get_wtime() - dtime[0];
	cout<<dtime[0]<<"s"<<endl;
	printf("Check one of elements: A[3][1] = %g \n", A[3][1]);


    //Result result = read (filename);

    /*
	cout<<"ijk multiplication ... ";
    dtime[0] = omp_get_wtime();
    C = ijkalgorithm(A, B);
	dtime[0] = omp_get_wtime() - dtime[0];
	cout<<dtime[0]<<"s"<<endl;

	cout<<"ikj multiplication ... ";
    dtime[1] = omp_get_wtime();
    C = ikjalgorithm(A, B);
	dtime[1] = omp_get_wtime() - dtime[1];
	cout<<dtime[1]<<"s"<<endl;
    */

	cout<<"ijk omp multiplication ... ";
    dtime[2] = omp_get_wtime();
//    C = ijkalgorithm_omp(A, B);
    ijkalgorithm_omp(A, B, C);
	dtime[2] = omp_get_wtime() - dtime[2];
	cout<<dtime[2]<<"s"<<endl;

	cout<<"ikj omp multiplication ... ";
    dtime[3] = omp_get_wtime();
    C = ikjalgorithm_omp(A, B);
	dtime[3] = omp_get_wtime() - dtime[3];
	cout<<dtime[3]<<"s"<<endl;

    //printMatrix(C);
    return 0;
}
