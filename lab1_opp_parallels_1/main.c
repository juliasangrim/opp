#include<mpi.h>
#include<stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#define n 5000


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
    }
}

void mul(const double* A, const double* B, double* result, int sizePartArray) { //change size and matrix
//    for (int i = 0; i < n; ++i) {
//        result[i] = 0;
//    }
    for (int i = 0; i < sizePartArray; ++i) {
        result[i] = 0; /// changed
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

void matrixInit(double* A, int procRank) {
    for (int i = 0; i < n; ++i){  //initialisation of matrix A, vector B
        for (int j = i; j < n; ++j) {
            double randValue = (double)(rand() % (600 - 300) / 3.0);
            if (i == j) {
                randValue += 4000;
            }
            A[i * n + j] = randValue;
            A[j * n + i] = randValue;
        }
    }
}


void vectorInit(double* X, int flag) {
    for (int i = 0; i < n; ++i) {
        if (flag == 1) {
            X[i] = (double) (rand() % (6000 - 3000) / 3.0);
        } else {
            X[i] = 0;
        }
    }
}

void distribMatrix(double* A, double** partA, int* shiftIndex, int* numElem, int procNum, int procRank) {
    int rowNum = n / procNum;
    int restRows = n;
    shiftIndex[0] = 0;
    numElem[0] = n * rowNum;
    for (int i = 1; i < procNum; ++i) {
        restRows -= rowNum;
        rowNum = restRows / (procNum - i);
        numElem[i] = n * rowNum;
        shiftIndex[i] = numElem[i - 1] + shiftIndex[i - 1];
    }
    *partA = (double *) malloc(numElem[procRank] * sizeof(double));
    MPI_Scatterv(A, numElem, shiftIndex, MPI_DOUBLE,
                 *partA, numElem[procRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void distribVector(double* B, double* Y, double** partB, double** partY, int* shiftIndex, int* numElem, int procNum, int procRank) {
    for (int i = 0; i < procNum; ++i) {
        numElem[i] /= n;
        shiftIndex[i] /= n;
    }
    *partB = (double*)malloc(numElem[procRank] * sizeof(double));
    *partY = (double*)malloc(numElem[procRank] * sizeof(double));
    MPI_Scatterv(B, numElem, shiftIndex, MPI_DOUBLE, *partB, numElem[procRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(Y, numElem, shiftIndex, MPI_DOUBLE, *partY, numElem[procRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void freeArrays(double* A, double* B, double* X, double* Y, int* shiftIndex, int* numElem, double* partArrayA, double* partArrayB, double* partArrayY, double* partArrayX, int procRank) {
    if (procRank == 0) {
        free(A);
    }
    free(B);
    free(X);
    free(Y);
    free(shiftIndex);
    free(numElem);
    free(partArrayB);
    free(partArrayX);
    free(partArrayY);
}


int calculate(double* partArrayA, double* X, double* partArrayX, double* B, double* partArrayB,
              double* Y, double* partArrayY, int* numElem, int* shiftIndex, int procRank) {
    double t = 0;
    double* partTmp = (double*)malloc(numElem[procRank] * sizeof(double));
    mul(partArrayA, X, partTmp, numElem[procRank]);

    sub(partTmp, partArrayB, partArrayY, numElem[procRank]); //calculate Y(0)
    double valueCheck = absVector(partArrayY, numElem[procRank])
                      / absVector(partArrayB, numElem[procRank]);//the value for checking when we should stop calculate
    double prevValue = 0;
    double epsilon = 0.000000001;
    int counter = 0;
    while (valueCheck >= epsilon) {
        prevValue = valueCheck;
        MPI_Allgatherv(partArrayY, numElem[procRank], MPI_DOUBLE,
                       Y, numElem, shiftIndex, MPI_DOUBLE, MPI_COMM_WORLD);

        mul(partArrayA, Y, partTmp, numElem[procRank]); // multiplex A and Y
        t = scalarMul(partArrayY, partTmp, numElem[procRank])
                / scalarMul(partTmp, partTmp, numElem[procRank]); //(Y, A*Y)/(A*Y,A*Y) calculate the T(n)
        mulVector(t, partArrayY, numElem[procRank]); //
        MPI_Scatterv(X, numElem, shiftIndex, MPI_DOUBLE, /// можно не резать каждый раз, а сделать это за циклом
                     partArrayX, numElem[procRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
        sub(partArrayX, partArrayY, partArrayX, numElem[procRank]); //calculate the X(n+1)

        MPI_Allgatherv(partArrayX, numElem[procRank], MPI_DOUBLE,
                       X, numElem, shiftIndex, MPI_DOUBLE, MPI_COMM_WORLD);
        mul(partArrayA, X, partTmp, numElem[procRank]);
        sub(partTmp, partArrayB, partArrayY, numElem[procRank]); // calculate Y(n)

        valueCheck = absVector(partArrayY, numElem[procRank])
                / absVector(partArrayB, numElem[procRank]); //calculate new value for checking
                printf("%f\n ", valueCheck);
        if (prevValue <= valueCheck) {  //in case when matrix have no limits
            counter++;
            if (counter >= 6) {
                if (procRank == 0) {
                    printf("No limits");
                    return 1;
                }
            }
        }
        else {
            counter = 0;
        }
    }
    free(partTmp);
    return  0;
}



int main(int argc, char* argv[]) {
    srand(4);
    /////////////init everything//////////////////
    double* A;
    double* B = (double*)malloc(n * sizeof(double));
    double* X = (double*)malloc(n * sizeof(double));
    double* Y = (double*)malloc(n * sizeof(double));
    double* partArrayA = NULL;
    double* partArrayB = NULL;
    double* partArrayY = NULL;
    int* shiftIndex = (int*)malloc(n * sizeof(int)); //array contain the shift from beginning of main array for each process
    int* numElem = (int*)malloc(n * sizeof(int)); //array contain amount of elements in each process
    MPI_Init(&argc, &argv);
    int procNum;
    int procRank;
    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    if (procRank == 0) {
        A = (double*)malloc(n * n * sizeof(double));
        matrixInit(A, procRank);
       // printMatrix(A, n);
        vectorInit(B, 1);
       // printVector(B, "B", 0, 1, n);
    }
    vectorInit(X, 0);
    MPI_Bcast(B, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    //////////distribute//////////
    distribMatrix(A, &partArrayA, shiftIndex, numElem, procNum, procRank);
    distribVector(B, Y, &partArrayB, &partArrayY, shiftIndex, numElem, procNum, procRank);

    ///////////calculate/////////////////////
    double* partArrayX = (double*)malloc(numElem[procRank] * sizeof(double));
    double start = MPI_Wtime();
    if (calculate(partArrayA, X, partArrayX, B, partArrayB, Y,
                  partArrayY, numElem, shiftIndex, procRank) == 1) {
        freeArrays(A, B, X, Y, shiftIndex, numElem, partArrayA, partArrayB, partArrayY, partArrayX, procRank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    double end = MPI_Wtime();
    if (procRank == 0) {
       // printVector(X, "X", procRank, 1, n);
        printf("Calculating time:%f\n", end - start);
    }
    ///////////////free///////////////
    freeArrays(A, B, X, Y, shiftIndex, numElem, partArrayA, partArrayB, partArrayY, partArrayX, procRank);
    MPI_Finalize();
    return 0;
}
