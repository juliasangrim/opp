#include<mpi.h>
#include<stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#define n 5

void printMatrix(const double* A, int size) {
    printf("A:\n");
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            printf("%.3f ", A[i*size + j]);
        }
        printf("\n");
    }
    printf("\n");
}

void printVector(const double* B, const char* name, int procRank, int procNum, int size) {
    for (int numProc = 0; numProc < procNum; ++numProc) {
        if (procRank == numProc) {
            printf("%s in rank %d:\n", name, procRank);
            for (int i = 0; i < n; ++i) {
                printf("%.3f\n", B[i]);
            }
            printf("\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

void mul(const double* A, const double* B, double* result, int sizePartArray) { //change size and matrix
    for (int i = 0; i < sizePartArray; ++i) {
        result[i] = 0;
        for (int j = 0; j < n; ++j) {
            result[i] += A[i * sizePartArray + j] * B[j];
        }
    }
}

void mulVector(double constant, double* A, int size) {
    for (int i = 0; i < n; ++i) {
        A[i] = constant * A[i];
    }
}

void sub(const double* A, const double* B, double* result, int size) {
    for (int i = 0; i < size ; ++i) {
        result[i] = A[i] - B[i];
    }
}

double scalarMul(const double* A, const double* B) {
    double result = 0;
    for (int i = 0; i < n; ++i) {
        result += A[i] * B[i];
    }
    return result;
}

double absVector(const double* A) {
    double res = 0;
    double sum = 0;
    for (int i = 0; i < n; ++i) {
        sum += pow(A[i], 2);
    }
    res = sqrt(sum);
    return res;
}

void matrixInit(double* A, double* B, double* X,  int procRank) {
    for (int i = 0; i < n; ++i){  //initialisation of matrix A, vector B
        B[i] = (double)(rand() % 2000 - 1000) / 3;
        X[i] = (double)(rand() % 2000 - 1000) / 3;
        for (int j = i; j < n; ++j) {
            double randValue = (double)(rand() % (2000 - 1000) / 3.0);
            if (i == j) {
                randValue += 50;
            }
            if (procRank== 0) {
                A[i * n + j] = randValue;
                A[j * n + i] = randValue;
            }
        }
    }
}

void distribArray(int* shiftIndex, int* numElem, int procNum) {
    int rowNum = n/procNum;
    int restRows = n;
    shiftIndex[0] = 0;
    numElem[0] = n * rowNum;
    for (int i = 1; i < procNum; ++i) {
        restRows -= rowNum;
        rowNum = restRows / (procNum - i);
        numElem[i] = n * rowNum;
        shiftIndex[i] = numElem[i - 1] + shiftIndex[i - 1];
    }

}

void distribVector(int* shiftIndex, int* numElem, int procNum) {
    for (int i = 0; i < procNum; ++i) {
        numElem[i] /= n;
        shiftIndex[i] /= n;
    }
}

//void freeMemory(double* A, double* B, double* X, double* Y, double* tmp, int* sendShiftIndex, int* sendNumElem, int procRank) {
//    if (procRank == 0) {
//        free(A);
//    }
//    free(B);
//    free(X);
//    //free(Y);
//    //free(tmp);
//    free(sendShiftIndex);
//    free(sendNumElem);
//}

int main(int argc, char* argv[]) {
    srand(0);
    double* A;
    double* B = (double*)malloc(n * sizeof(double));
    double* X = (double*)malloc(n * sizeof(double));
    int* shiftIndex = (int*)malloc(n * sizeof(int)); //array contain the shift from beginning of main array for each process
    int* numElem = (int*)malloc(n * sizeof(int)); //array contain amount of elements in each process
    double t = 0;
    MPI_Init(&argc, &argv);
    int procNum;
    int procRank;
    /////////////print and init//////////////////
    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    int rowsNum;
    if (procRank == 0) {
        A = (double*)malloc(n * n * sizeof(double));
    }
    matrixInit(A, B, X, procRank);
    if (procRank == 0) {
        printMatrix(A, n);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    printVector(B, "B", procRank, procNum, n);
    //////////distribute//////////
    distribArray(shiftIndex, numElem, procNum);
    double* partArrayA = (double*)malloc(numElem[procRank] * sizeof(double)); //array for part of array A
    MPI_Scatterv(A, numElem, shiftIndex, MPI_DOUBLE,
                 partArrayA, numElem[procRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    distribVector(shiftIndex, numElem, procNum);
    ///////////calculate/////////////////////
    double* tmp = (double*)malloc(numElem[procRank] * sizeof(double));
    double* Y = (double*)malloc(numElem[procRank] * sizeof(double));
    printVector(X, "X", procRank, procNum, n);
    MPI_Barrier(MPI_COMM_WORLD);
    mul(partArrayA, X, tmp, numElem[procRank]/n);
//    sub(tmp, B, Y); //calculate Y(0)
//    double valueCheck = absVector(Y) / absVector(B); //the value for checking when we should stop calculate
//    double prevValue = 0;
//    double epsilon = 0.00001;
//    while (valueCheck >= epsilon) {
//        prevValue = valueCheck;
//        mul(A, Y, tmp); // multiplex A and Y
//        t = scalarMul(Y, tmp) / scalarMul(tmp, tmp); //(Y, A*Y)/(A*Y,A*Y) calculate the T(n)
//        mulVector(t, Y); //
//        sub(X, Y, X); //calculate the X(n+1)
//        mul(A, X, tmp);
//        sub(tmp, B, Y); // calculate Y(n)
//        valueCheck = absVector(Y) / absVector(B); //calculate new value for checking
//        if (prevValue == valueCheck) {  //in case when matrix have no limits
//            printf("no limits\n");
//            return 0;
//        }
//    }
//    printVector(X, "X");
   // freeMemory(A, B, X, Y, tmp, shiftIndex, numElem, procRank);
    MPI_Finalize();
    return 0;
}
