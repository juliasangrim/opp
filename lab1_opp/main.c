#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#define n 2

void printMatrix(const double* A) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
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

void mul(const double* A, const double* B, double* result) {
    for (int i = 0; i < n; i++) {
        result[i] = 0;
        for (int j = 0; j < n; j++) {
            result[i] += A[i * n + j] * B[j];
        }
    }
}

void mulVector(double constant, double* A) {
    for (int i = 0; i < n; i++) {
        A[i] = constant * A[i];
    }
}

void sub(const double* A, const double* B, double* result) {
    for (int i = 0; i < n ; i++) {
        result[i] = A[i] - B[i];
    }
}

double scalarMul(const double* A, const double* B) {
    double result = 0;
    for (int i = 0; i < n; i++) {
        result += A[i] * B[i];
    }
    return result;
}

double absVector(const double* A) {
    double res = 0;
    double sum = 0;
    for (int i = 0; i < n; i++) {
        sum += pow(A[i], 2);
    }
    res = sqrt(sum);
    return res;
}

int main(int argc, char* argv[]) {
    srand(time(NULL));
    double A[n * n];
    double B[n];
    double X[n] = {0};
    double Y[n] = {0};
    double temp[n] = {0};
    double t = 0;
    for (int i = 0; i < n; i++){
        B[i] = (double)(rand() % 2000 - 1000) / 3;
        for (int j = i; j < n; j++) {
            double randValue = (double)(rand() % (2000 - 1000) / 3.0);
            if (i == j) {
                randValue += 50;
            }
            A[i * n + j] = randValue;
            A[j * n + i] = randValue;
        }
    }
    printMatrix(A);
    printVector(B, "B");
    mul(A, X, temp);
    sub(temp, B, Y);
    double valueCheck = absVector(Y) / absVector(B);
    double prevValue = 0;
    double epsilon = 0.00001;
    while (valueCheck >= epsilon) {
        prevValue = valueCheck;
        mul(A, Y, temp);
        t = scalarMul(Y, temp) / scalarMul(temp, temp);
        mulVector(t, Y);
        sub(X, Y, X);
        mul(A, X, temp);
        sub(temp, B, Y);
        valueCheck = absVector(Y) / absVector(B);
        if (prevValue == valueCheck) {
            printf("No limits\n");
            return 0;
        }
    }
    printVector(X, "X");
    return 0;
}
