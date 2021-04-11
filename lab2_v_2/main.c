#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#include <stdbool.h>
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

void freeArray(double* A, double* B, double* Y, double* X, double* temp)
{
    free(A);
    free(B);
    free(Y);
    free(X);
    free(temp);
}

void calculate(double* A, double* X, double* B, double* Y, double* temp, bool* flag)
{
    double globalSum = 0;
    int i, j;
    double t = 0;
    double absY = 0;
    double absB = 0;
    double valueCheck = 0;
    double prevValue = 0;
    double epsilon = 0.000000001;
    int count = 0;
    int cnt = 0;
    double YTemp = 0;
    double TempTemp = 0;
    bool isItLimit = false;
#pragma omp parallel default(shared) private(cnt)
    {
        ////////////mul(A, X, temp);///////////////
        #pragma omp for private(i, j)
        for (i = 0; i < n; ++i) {
            temp[i] = 0;
            for (j = 0; j < n; ++j) {
                X[j] = 1;
                temp[i] += A[i * n + j] * X[j];
            }
            ///////////sub(temp, B, Y); /////////////calculate Y(0)
            Y[i] = temp[i] - B[i];
        }
       //printVector(temp, "temp");
        //printVector(Y, "Y");
        ////////////////absVector(Y)/////////////////
        #pragma omp single
        globalSum = 0;
        #pragma omp for
        for (i = 0; i < n; ++i) {
            #pragma omp atomic
            globalSum += pow(Y[i], 2);
        }
        #pragma omp single
        {
            absY = sqrt(globalSum);
            printf("%lf ", absY);
            globalSum = 0;
        }
        //////////////////absVector(B)////////////////
        #pragma omp for
        for (i = 0; i < n; ++i) {
            #pragma omp atomic
            globalSum += pow(B[i], 2);
        }

        #pragma omp single
        {
        valueCheck = absY / absB;
        absB = sqrt(globalSum);
        //printf("%lf ", absB);
        }
        /////////////////double valueCheck = absVector(Y) / absVector(B); ///////////////////the value for checking when we should stop calculate
        //printf("%lf ", valueCheck);
        if (valueCheck >= epsilon) {
            isItLimit = true;
        }
        while (isItLimit)
        {
            cnt++;
            #pragma omp single
            prevValue = valueCheck;
            /////////////mul(A, Y, temp); /////////// multiplex A and Y
            #pragma omp for private(i, j)
            for (i = 0; i < n; ++i) {
                temp[i] = 0;
                for (j = 0; j < n; ++j) {
                    temp[i] += A[i * n + j] * Y[j];
                }
            }
            //////////////////t = scalarMul(Y, temp) / scalarMul(temp, temp); ////////////(Y, A*Y)/(A*Y,A*Y) calculate the T(n)
            #pragma omp single
            globalSum = 0;
            #pragma omp for reduction(+:globalSum)
            for (i = 0; i < n; ++i) {
                globalSum += Y[i] * temp[i];
            }
            #pragma omp single
            YTemp = globalSum;
            #pragma omp single
            globalSum = 0;
            #pragma omp for reduction(+:globalSum)
            for (i = 0; i < n; ++i) {
                globalSum += temp[i] * temp[i];
            }
            #pragma omp single
            TempTemp = globalSum;
            #pragma omp single
            t = YTemp / TempTemp;
            #pragma omp for private(i)
            for (i = 0; i < n; ++i) {
                Y[i] = t * Y[i];
                ///////////////sub(X, Y, X); ////////////calculate the X(n+1)
                X[i] = X[i] - Y[i];
            }
            //////////////mul(A, X, temp);//////////////////
            #pragma omp for private(i, j)
            for (i = 0; i < n; ++i) {
                temp[i] = 0;
                for (j = 0; j < n; ++j) {
                    temp[i] += A[i * n + j] * X[j];
                }
                ////////////////sub(temp, B, Y); //////////////// calculate Y(n)
                Y[i] = temp[i] - B[i];
            }
            ///////////////////absVector(Y)//////////////////
            #pragma omp single
            globalSum = 0;
            #pragma omp for private(i) reduction(+:globalSum)
            for (i = 0; i < n; ++i) {
                globalSum += pow(Y[i], 2);
            }
            #pragma omp single
            absY = sqrt(globalSum);
            /////////////////absVector(B)///////////////////////
            #pragma omp single
            globalSum = 0;
            #pragma omp for private(i) reduction(+:globalSum)
            for (i = 0; i < n; ++i) {
                globalSum += pow(B[i], 2);
            }
            #pragma omp single
            absB = sqrt(globalSum);
            ///////////////////double valueCheck = absVector(Y) / absVector(B); ////////////////the value for checking when we should stop calculate
            #pragma omp single
	    {
            valueCheck = absY / absB;
            if (valueCheck >= epsilon) {
                isItLimit = true;
            } else {
                isItLimit = false;
            }
            if (prevValue <= valueCheck) {
                count++;  //in case when matrix have no limits
                if (count >= 6) {
                    printf("no limits\n");
                    *flag = true;
                    isItLimit = false;
                }
            } else {
                count = 0;
            }
	    }

        }
        printf("COUNT %d\n", cnt);
    }
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
    /////////////////////////////////
    arrayInit(A, B);
    double start, end;
    //  printMatrix(A);
      //printVector(B, "B");
    bool flag = 0;
    start = omp_get_wtime();
    calculate(A, X, B, Y, temp, &flag);
    if (flag == true)
    {
        freeArray(A, B, Y, X, temp);
        return 1;
    }
    end = omp_get_wtime();

   //  printVector(X, "X");
    printf("Time taken: %lf sec.\n",end - start);
    freeArray(A, B, Y, X, temp);
    return 0;
}
