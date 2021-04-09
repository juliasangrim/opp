#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#define n 1000

void printMatrix(const double* A) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            printf("%f ", A[i*n + j]);
        }
        printf("\n");
    }
}

void printVector(const double* B, const char* name) {
    printf("%s:\n", name);
    for (int i = 0; i < n; ++i) {
        printf("%f\n", B[i]);
    }
    printf("\n");
}

void matrixInit(double* A, double* B) {
    for (int i = 0; i < n; i++){  //initialisation of matrix A, vector B
        B[i] = (double)(rand() % 2000 - 1000) / 3;
        for (int j = i; j < n; j++) {
            double randValue = (double)(rand() % (200 - 100) / 3.0);
            if (i == j) {
                randValue += 50;
            }
            A[i * n + j] = randValue;
            A[j * n + i] = randValue;
        }
    }
}

void mul(const double* A, const double* B, double* result) {
    for (int i = 0; i < n; ++i) {
        result[i] = 0;
        for (int j = 0; j < n; ++j) {
            result[i] += A[i * n + j] * B[j];
        }
    }
}

void mulVector(double constant, double* A) {
    for (int i = 0; i < n; ++i) {
        A[i] = constant * A[i];
    }
}

void sub(const double* A, const double* B, double* result) {
    for (int i = 0; i < n ; ++i) {
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

int calculate(double* A, double* X, double* B, double* Y, double* temp) {
    double t = 0;
    mul(A, X, temp);
    sub(temp, B, Y); //calculate Y(0)
    double valueCheck = absVector(Y) / absVector(B); //the value for checking when we should stop calculate
    double prevValue = 0;
    double epsilon = 0.0001;
    int count = 0;
    while (valueCheck >= epsilon) {
        prevValue = valueCheck;
        mul(A, Y, temp); // multiplex A and Y
        t = scalarMul(Y, temp) / scalarMul(temp, temp); //(Y, A*Y)/(A*Y,A*Y) calculate the T(n)
        mulVector(t, Y); //
        sub(X, Y, X); //calculate the X(n+1)
        mul(A, X, temp);
        sub(temp, B, Y); // calculate Y(n)
        valueCheck = absVector(Y) / absVector(B); //calculate new value for checking
        if (prevValue <= valueCheck) {
            count++;  //in case when matrix have no limits
            if (count >= 6) {
                printf("no limits\n");
                return 1;
            }
        }
        else {
            count = 0;
        }
    }
    return 0;
}

int main(int argc, char* argv[]) {
    srand(0);
    /////////////init////////////////
    double A[n * n];
    double B[n];
    double X[n] = {0};
    double Y[n] = {0};
    double temp[n] = {0};
    double t = 0;
    /////////////////////////////////
    matrixInit(A, B);
    struct timespec start, end;
    //printMatrix(A);
    //printVector(B, "B");
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    if (calculate(A, X, B, Y, temp) == 1) {
        return 1;
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    //printVector(X, "X");
    printf("Time taken: %lf sec.\n",end.tv_sec-start.tv_sec+ 0.000000001*(end.tv_nsec-start.tv_nsec));
    return 0;
}
