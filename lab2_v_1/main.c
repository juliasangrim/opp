#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#define n 5000

void printMatrix(const double* A)
{
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            printf("%f ", A[i*n + j]);
        }
        printf("\n");
    }
}

void printVector(const double* B, const char* name)
{
    printf("%s:\n", name);
    for (int i = 0; i < n; ++i)
    {
        printf("%f\n", B[i]);
    }
    printf("\n");
}

void arrayInit(double* A, double* B)
{
    for (int i = 0; i < n; i++){  //initialisation of matrix A, vector B
        B[i] = (double)(rand() % 2000 - 1000) / 3;
        for (int j = i; j < n; j++)
        {
            double randValue = (double)(rand() % (200 - 100) / 3.0);
            if (i == j)
            {
                randValue += 3000;
            }
            A[i * n + j] = randValue;
            A[j * n + i] = randValue;
        }
    }
}

void mul(const double* A, const double* B, double* result)
{
    double globalSum;
    int i, j;
    for (i = 0; i < n; ++i)
    {
        globalSum = 0;
        #pragma omp parallel for default(shared) reduction(+:globalSum)
        for (j = 0; j < n; ++j)
        {
            globalSum += A[i * n + j] * B[j];
        }
        result[i] = globalSum;
    }
}

void mulVector(double constant, double* A)
{
     int i;
    #pragma omp parallel for default(shared)
     for (i = 0; i < n; ++i) {
         A[i] = constant * A[i];
     }
}

void sub(const double* A, const double* B, double* result)
{
    int i;
    #pragma omp parallel for default(shared)
    for (i = 0; i < n; ++i) {
        result[i] = A[i] - B[i];
    }
}

double scalarMul(const double* A, const double* B)
{
    double result = 0;
    int i = 0;
    double tmpResult = 0;
    #pragma omp parallel for default(shared) private(i) reduction(+:tmpResult)
        for (i = 0; i < n; ++i) {
            tmpResult += A[i] * B[i];
        }
        result = tmpResult;
    return result;
}

double absVector(const double* A)
{
    double res = 0;
    double tmpSum = 0;
    double sumResult;
    int i;
    #pragma omp parallel for default(shared) private(i) reduction(+:tmpSum)
        for (i = 0; i < n; ++i) {
            tmpSum += pow(A[i], 2);
        }
        sumResult = tmpSum;
        res = sqrt(sumResult);
    return res;
}

int calculate(double* A, double* X, double* B, double* Y, double* temp)
{
    double t = 0;
    double prevValue = 0;
    double epsilon = 0.000000001;
    int count = 0;
    int cnt = 0;
    mul(A, X, temp);
    sub(temp, B, Y); //calculate Y(0)
    double valueCheck = absVector(Y) / absVector(B); //the value for checking when we should stop calculate
    while (valueCheck >= epsilon)
    {
        cnt++;
        prevValue = valueCheck;
        mul(A, Y, temp); // multiplex A and Y
        t = scalarMul(Y, temp) / scalarMul(temp, temp); //(Y, A*Y)/(A*Y,A*Y) calculate the T(n)
        mulVector(t, Y); //
        sub(X, Y, X); //calculate the X(n+1)
        mul(A, X, temp);
        sub(temp, B, Y); // calculate Y(n)
        valueCheck = absVector(Y) / absVector(B); //calculate new value for checking
        if (prevValue <= valueCheck)
        {
            count++;  //in case when matrix have no limits
            if (count >= 6)
            {
                printf("no limits\n");
                return 1;
            }
        }
        else {
            count = 0;
        }
    }
    printf("%d", cnt);
    return 0;
}

int main(int arg, char* argv[])
{
    int threadsNum = atoi(argv[1]);
    omp_set_num_threads(threadsNum);
    srand(4);
    /////////////init////////////////
    double* A = (double*)malloc(sizeof(double) * n * n);
    double* B = (double*)malloc(sizeof(double) * n);
    double* X = (double*)calloc(n, sizeof(double));
    double* Y  = (double*)calloc( n, sizeof(double));
    double* temp = (double*)calloc(n, sizeof(double));
    double t = 0;
    /////////////////////////////////
    arrayInit(A, B);
    struct timespec start, end;
  //  printMatrix(A);
  //  printVector(B, "B");
  
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    if (calculate(A, X, B, Y, temp) == 1)
    {
        return 1;
    }
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);

   // printVector(X, "X");
    printf("Time taken: %lf sec.\n",end.tv_sec-start.tv_sec + 0.000000001 * (end.tv_nsec-start.tv_nsec));
    free(A);
    free(X);
    free(Y);
    free(B);
    free(temp);
    return 0;
}
