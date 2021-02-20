#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int* mul(const int* A, const int* B, int n) {
    int C[n];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            C[i] = 0;
            C[i] += A[i * n + j] * B[j];
        }
    }
    return (int *) C[n];
}

int* sub(const int* A, const int* B, int n) {
    int C[n][n];
    for (int i = 0; i < n ; i++) {
        for (int j = 0; j < n; j++) {
            C[i * n + j] = A[i * n + j] - B[i];
        }
    }
    return (int *)C[n];
}

int main(int argc, char* argv[]) {
    int n = atoi( argv[1]);
    srand(time(NULL));
    int A[n][n];
    int B[n], X[n], Y[n], T[n];
    for (int i = 0; i < n; i++){
        B[i] = rand() % 2000 - 1000;
        X[i] = rand() % 2000 - 1000;
        for (int j = i; j < n; j++) {
            A[i][j] = rand() % 200 - 100;
            A[j][i] = A[i][j];
        }
        if (A[i][i] < 0) {
            A[i][i] = abs(A[i][i]);
            A[i][i] += 50;
        }
    }
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++) {
            printf("%3d ", A[i][j]);
        }
        printf("\n");
    }
    for (int i = 0; i < n; i++) {
        printf("%d ", B[i]);

    }

    Y[n] = (int) mul(A[n], (const int *) B[n], n);
    A[n] = (int) sub(A[n], (const int *) B[n], n);

    for (int i = 0; i < n; i++) {
        printf("%d ", Y[i]);

    }
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++) {
            printf("%3d ", T[i][j]);
        }
        printf("\n");
    }
    return 0;
}
