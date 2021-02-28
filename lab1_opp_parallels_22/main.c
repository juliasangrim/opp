#include<mpi.h>
#include<stdio.h>
#include <stdlib.h>
#include <math.h>
#define n 3

void printMatrix(const double* A, int size) {
    printf("A:\n");
    for (int i = 0; i < n; ++i) {
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
            for (int i = 0; i < size; ++i) {
                printf("%.3f\n", B[i]);
            }
            printf("\n");
        }
    }
}

void mul(const double* A, const double* B, double* result, int sizePartVector) {
    double* tmp = malloc(n * sizePartVector * sizeof(double));
    double* tmpVector = malloc(n * sizeof(double));
    for (int i = 0; i < n; ++i) {
        tmpVector[i] = 0;
    }
    for (int i = 0; i < sizePartVector; ++i) {
        for (int j = 0; j < n; ++j) {
            tmpVector[j] += A[i * n + j] * B[i];
        }
    }
    MPI_Reduce(tmpVector, result, n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    free(tmp);
    free(tmpVector);
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

void matrixInit(double* A, int procRank) {
    for (int i = 0; i < n; ++i){  //initialisation of matrix A, vector B
        for (int j = i; j < n; ++j) {
            double randValue = (double)(rand() % (200 - 100) / 3.0);
            if (i == j) {
                randValue += 50;
            }
            A[i * n + j] = randValue;
            A[j * n + i] = randValue;
        }
    }
}


void vectorInit(double* X, double* B) {
    for (int i = 0; i < n; ++i) {
        B[i] = (double)(rand() % (2000 - 1000) / 3.0);
        X[i] = (double)(rand() % (2000 - 1000) / 3.0);
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
    double* A = NULL;
    double* B = NULL;
    double* X = NULL;
    double* Y = NULL;
    double* temp = NULL;
    int* shiftIndex = (int*)malloc(n * sizeof(int)); //array contain the shift from beginning of main array for each process
    int* numElem = (int*)malloc(n * sizeof(int)); //array contain amount of elements in each process
    double t = 0;
    MPI_Init(&argc, &argv);
    int procNum;
    int procRank;
    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    if (procRank == 0) {
        A = (double*)malloc(n * n * sizeof(double));
        B = (double*)malloc(n * sizeof(double));
        X = (double*)malloc(n * sizeof(double));
        Y = (double*)malloc(n * sizeof(double));
        temp = (double*)malloc(n * sizeof(double));
        matrixInit(A, procRank);
        printMatrix(A, n);
        vectorInit(X, B);
        printVector(B, "B", procRank, 1, n);
        //printVector(X, "X", procRank, 1, n);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    distribArray(shiftIndex, numElem, procNum);
    double* partArrayA = (double*)malloc(numElem[procRank] * sizeof(double)); //array for part of array A
    MPI_Scatterv(A, numElem, shiftIndex, MPI_DOUBLE, partArrayA, numElem[procRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    distribVector(shiftIndex, numElem, procNum);

    double* partArrayB = (double*)malloc(numElem[procRank] * sizeof(double));
    MPI_Scatterv(B, numElem, shiftIndex, MPI_DOUBLE, partArrayB, numElem[procRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double* partArrayY = (double*)malloc(numElem[procRank] * sizeof(double));
    MPI_Scatterv(Y, numElem, shiftIndex, MPI_DOUBLE, partArrayY, numElem[procRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double* partArrayX = (double*)malloc(numElem[procRank] * sizeof(double));
    MPI_Scatterv(X, numElem, shiftIndex, MPI_DOUBLE, partArrayX, numElem[procRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    double start = MPI_Wtime();
    mul(partArrayA, partArrayX, temp, numElem[procRank]);

    double* partArrayTmp = (double*)malloc(numElem[procRank] * sizeof(double));
    MPI_Scatterv(temp, numElem, shiftIndex, MPI_DOUBLE, partArrayTmp, numElem[procRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    sub(partArrayTmp, partArrayB, partArrayY, numElem[procRank]); //calculate Y(0)
    double valueCheck = absVector(partArrayY, numElem[procRank]) / absVector(partArrayB, numElem[procRank]); //the value for checking when we should stop calculate
    double prevValue = 0;
    double epsilon = 0.00001;
    int count = 0;
    while (valueCheck >= epsilon) {
        prevValue = valueCheck;
        mul(partArrayA, partArrayY, temp, numElem[procRank]); // multiplex A and Y

        MPI_Scatterv(temp, numElem, shiftIndex, MPI_DOUBLE, partArrayTmp, numElem[procRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

        t = scalarMul(partArrayY, partArrayTmp, numElem[procRank]) / scalarMul(partArrayTmp, partArrayTmp, numElem[procRank]); //(Y, A*Y)/(A*Y,A*Y) calculate the T(n)
        mulVector(t, partArrayY, numElem[procRank]); //
        sub(partArrayX, partArrayY, partArrayX, numElem[procRank]); //calculate the X(n+1)
        mul(partArrayA, partArrayX, temp, numElem[procRank]);

        MPI_Scatterv(temp, numElem, shiftIndex, MPI_DOUBLE, partArrayTmp, numElem[procRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
        sub(partArrayTmp, partArrayB, partArrayY, numElem[procRank]); // calculate Y(n)
        valueCheck = absVector(partArrayY, numElem[procRank]) / absVector(partArrayB, numElem[procRank]); //calculate new value for checking
        if (prevValue <= valueCheck) {
            count++;  //in case when matrix have no limits
            if (count >= 6) {
                printf("no limits\n");
                MPI_Abort(MPI_COMM_WORLD, 1);
           }
        }
        else {
            count = 0;
        }
    }
    double end = MPI_Wtime();
    if (procRank == 0){
        printf("Calculating time:%f\n", end - start);
    }
    MPI_Gatherv(partArrayX, numElem[procRank], MPI_DOUBLE, X, numElem, shiftIndex, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (procRank == 0){
        printVector(X, "X", procRank, 1, n);
    }
    ///////////////free///////////////
  //  freeArrays(A, B, X, Y, tmp, shiftIndex, numElem, partArrayA, partArrayB, partArrayY, partArrayX, procRank);
    MPI_Finalize();
    return 0;
}
