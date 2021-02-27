#include<mpi.h>
#include<stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#define n 3

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
            result[i] += A[i * n + j] * B[j];
        }
    }
}

void mulVector(double constant, double* A, int size) {
    for (int i = 0; i < size; ++i) {
        A[i] = constant * A[i];
    }
}

void sub(const double* A, const double* B, double* result, int size) {
    for (int i = 0; i < size ; ++i) {
        result[i] = A[i] - B[i];
    }
}

double scalarMul(const double* A, const double* B, int size) {
    double result = 0;
    double sum = 0;
    for (int i = 0; i < size; ++i) {
        sum += A[i] * B[i];
    }
    MPI_Allreduce(&sum, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return result;
}



double absVector(const double* A, int size) {
    double res = 0;
    double sum = 0;
    for (int i = 0; i < size; ++i) {
        sum += pow(A[i], 2);
    }
    MPI_Allreduce(&sum, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    res = sqrt(res);
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

void freeArrays(double* A, double* B, double* X, double* Y, double* tmp, int* shiftIndex, int* numElem, double* partArrayA, double* partArrayB, double* partArrayY, double* partArrayX, int procRank) {
    if (procRank == 0) {
        free(A);
    }
    free(B);
    free(X);
    free(Y);
    free(tmp);
    free(shiftIndex);
    free(numElem);
    free(partArrayA);
    free(partArrayB);
    free(partArrayX);
    free(partArrayY);
}


int main(int argc, char* argv[]) {
    srand(0);
    /////////////init everything//////////////////
    double* A;
    double* B = (double*)malloc(n * sizeof(double));
    double* X = (double*)malloc(n * sizeof(double));
    int* shiftIndex = (int*)malloc(n * sizeof(int)); //array contain the shift from beginning of main array for each process
    int* numElem = (int*)malloc(n * sizeof(int)); //array contain amount of elements in each process
    double* Y = (double*)malloc(n * sizeof(double));
    double t = 0;
    MPI_Init(&argc, &argv);
    int procNum;
    int procRank;
    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    int rowsNum;
    if (procRank == 0) {
        A = (double*)malloc(n * n * sizeof(double));
    }
    matrixInit(A, B, X, procRank);
//    if (procRank == 0) {
//        printMatrix(A, n);
//    }
    MPI_Barrier(MPI_COMM_WORLD);
    //printVector(B, "B", procRank, procNum, n);
    //////////distribute//////////
    double start = MPI_Wtime();
    distribArray(shiftIndex, numElem, procNum);
    double* partArrayA = (double*)malloc(numElem[procRank] * sizeof(double)); //array for part of array A
    MPI_Scatterv(A, numElem, shiftIndex, MPI_DOUBLE,
                 partArrayA, numElem[procRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    distribVector(shiftIndex, numElem, procNum);
    //////making part of arrays///////////////////
    double* partArrayB = (double*)malloc(numElem[procRank] * sizeof(double));
    double* partArrayY = (double*)malloc(numElem[procRank] * sizeof(double));
    double* partArrayX = (double*)malloc(numElem[procRank] * sizeof(double));
    MPI_Scatterv(B, numElem, shiftIndex, MPI_DOUBLE, partArrayB, numElem[procRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(Y, numElem, shiftIndex, MPI_DOUBLE, partArrayY, numElem[procRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    ///////////calculate/////////////////////
    double* tmp = (double*)malloc(numElem[procRank] * sizeof(double));
    mul(partArrayA, X, tmp, numElem[procRank]);
    sub(tmp, partArrayB, partArrayY, numElem[procRank]); //calculate Y(0)
    double valueCheck = absVector(partArrayY, numElem[procRank]) / absVector(partArrayB, numElem[procRank]);//the value for checking when we should stop calculate
    double prevValue = 0;
    double epsilon = 0.00001;
    int counter = 0;
    while (valueCheck >= epsilon) {
        prevValue = valueCheck;
        MPI_Allgatherv(partArrayY, numElem[procRank], MPI_DOUBLE, Y, numElem, shiftIndex, MPI_DOUBLE, MPI_COMM_WORLD);
        mul(partArrayA, Y, tmp, numElem[procRank]); // multiplex A and Y
        t = scalarMul(partArrayY, tmp, numElem[procRank]) / scalarMul(tmp, tmp, numElem[procRank]); //(Y, A*Y)/(A*Y,A*Y) calculate the T(n)
        mulVector(t, partArrayY, numElem[procRank]); //
        MPI_Scatterv(X, numElem, shiftIndex, MPI_DOUBLE, partArrayX, numElem[procRank], MPI_DOUBLE, 0, MPI_COMM_WORLD); //TODO check
        sub(partArrayX, partArrayY, partArrayX, numElem[procRank]); //calculate the X(n+1)
        MPI_Allgatherv(partArrayX, numElem[procRank], MPI_DOUBLE, X, numElem, shiftIndex, MPI_DOUBLE, MPI_COMM_WORLD);
        mul(partArrayA, X, tmp, numElem[procRank]);
        sub(tmp, partArrayB, partArrayY, numElem[procRank]); // calculate Y(n)
        valueCheck = absVector(partArrayY, numElem[procRank]) / absVector(partArrayB, numElem[procRank]); //calculate new value for checking
        if (prevValue <= valueCheck) {  //in case when matrix have no limits
            counter++;
            if (counter >= 6) {
                if (procRank == 0) {
                    printf("No limits");
                    freeArrays(A, B, X, Y, tmp, shiftIndex, numElem, partArrayA, partArrayB, partArrayY, partArrayX, procRank);
                    return 1;
                }
            }
        }
        else {
            counter = 0;
        }
    }
    double end = MPI_Wtime();
    if (procRank == 0){
        printf("Calculating time:%f\n", end - start);
    }
    printVector(X, "X", procRank, procNum, n);
    ///////////////free///////////////
    freeArrays(A, B, X, Y, tmp, shiftIndex, numElem, partArrayA, partArrayB, partArrayY, partArrayX, procRank);
    MPI_Finalize();
    return 0;
}
