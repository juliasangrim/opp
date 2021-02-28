#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define N 8
#define eps 0.00001

int size;
int rank;

void printMatrix(const double *matrix, int rows, int cols)
{
    for (int y = 0; y < rows; ++y)
    {
        for (int x = 0; x < cols; ++x)
        {
            printf("%f ", matrix[y*cols + x]);
        }
        printf("\n");
    }
}

void printVector(const double *vector, const char *name, int border)
{
    printf("----------\n%s:\n", name);
    for (int i = 0; i < border; ++i)
    {
        printf("%16.12f\n", vector[i]);
    }
    printf("----------\n");
}

void printEquation(const double *A, const double *b)
{
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            printf("%.5f ", A[i*N + j]);
        }
        printf("%.5f\n", b[i]);
    }
}

void matrix_vector_mul(const double *matrix, const double *vector, double *result, int rows, int cols)
/**
 * 'result' can't be the same vector as 'vector'
 */
{
    for (int i = 0; i < rows; ++i)
    {
        result[i] = 0;
        for (int j = 0; j < cols; ++j)
        {
            result[i] += matrix[i*cols + j] * vector[j];
        }
    }
}

void vector_vector_sub(const double *first, const double *second, double *result, int border, int offset)
/**
 * 'result' can be the same vector as 'first' or 'second'
 */
{
    for (int i = 0; i < border; ++i)
    {
        result[i] = first[i] - second[i+offset];
    }
}

double scalar_mul(const double *first, const double *second, int border)
{
    double result = 0;
    double gatherResult;
    for (int i = 0; i < border; ++i)
    {
        result += first[i] * second[i];
    }
    MPI_Allreduce(&result, &gatherResult, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return gatherResult;
}

void double_vector_mul(const double f, const double *vector, double *result)
/**
 * 'result' can be the same vector as 'vector'
 */
{
    for (int i = 0; i < N; ++i)
    {
        result[i] = f * vector[i];
    }
}

double vector_module(const double *vector)
{
    double result = 0;
    double gatherResult;
    for (int i = 0; i < N; ++i)
    {
        result += vector[i] * vector[i];
    }
    MPI_Allreduce(&result, &gatherResult, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    gatherResult = sqrt(gatherResult);
    return gatherResult;
}

void dataInitialization(double **A_mat, double **b_vec, double **x_vec, double **y_vec)
{
    *b_vec = (double*) malloc (N * sizeof(double)); // каждый процесс обладает копией вектора b
    *x_vec = (double*) malloc (N * sizeof(double));
    *y_vec = (double *) malloc (N * sizeof(double));
    if (rank == 0)
    {
        long seed = 1614096481;//time(NULL);
        printf("seed: %ld\n\n", seed); // 1614096381(20, 0.00001) - ~16.7sec; 1614096481(30, 0.00001) - ~30sec
        srand(seed);

        *A_mat = (double*) malloc (N * N * sizeof(double)); // матрицей A обладает только нулевой процесс
        for (int y = 0; y < N; ++y)
        {
            for (int x = y; x < N; ++x)
            {
                double val = (rand() % 2000 - 1000) / 13.0;
                if (x == y) (*A_mat)[y*N + x] = val +50;
                else
                {
                    (*A_mat)[y*N + x] = val;
                    (*A_mat)[x*N + y] = val;
                }
            }
            (*b_vec)[y] = (rand() % 2000 - 1000) / 13.0;           // нулевой процесс инициализирует вектор b
            (*x_vec)[y] = 0;
        }
    }
    MPI_Bcast(*b_vec, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);    //...и раздает его всем процессам
    MPI_Bcast(*x_vec, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void dataDistribution(double *A_mat,
                      double **A_mat_forProc,
                      int **elemForProc,
                      int **indForProc)
{
    /// определяем с какого индекса раздавать матрицу A каждому процессу и по сколько элементов
    (*elemForProc) = (int*) malloc (size * sizeof(int));
    (*indForProc) = (int*) malloc (size * sizeof(int));
    int restRows = N;
    int rowNum_forProc = (restRows / size);
    (*elemForProc)[0] = rowNum_forProc * N;
    (*indForProc)[0] = 0;
    for (int i = 1; i < size; ++i)
    {
        restRows -= rowNum_forProc;
        rowNum_forProc = restRows / (size - i);
        (*elemForProc)[i] = rowNum_forProc * N;
        (*indForProc)[i] = (*indForProc)[i - 1] + (*elemForProc)[i - 1];
    }

    *A_mat_forProc = (double*) malloc ((*elemForProc)[rank] * sizeof(double));
    MPI_Scatterv(A_mat, (*elemForProc), (*indForProc), MPI_DOUBLE, *A_mat_forProc, (*elemForProc)[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for (int i = 0; i < size; ++i)
    {
        (*elemForProc)[i] /= N;
        (*indForProc)[i] /= N;
    }
}

void dataDeletion(double *A_mat, double *A_mat_forProc, double *b_vec, double *x_vec, double *y_vec, int *elemForProc, int *indForProc)
{
    if (rank == 0)
    {
        free(A_mat);
    }
    free(A_mat_forProc);
    free(b_vec);
    free(x_vec);
    free(y_vec);
    free(elemForProc);
    free(indForProc);
}

int calculate(double *A_mat_forProc,
              double *b_vec,
              double *x_vec,
              double *y_vec,
              double t,
              int *vecElemForProc,
              int *vecIndForProc)
{
    int rowNum_forProc = vecElemForProc[rank];
    double *x_vec_forProc = (double*) malloc (rowNum_forProc * sizeof(double));
    double *b_vec_forProc = (double*) malloc (rowNum_forProc * sizeof(double));
    double *y_vec_forProc = (double*) malloc (rowNum_forProc * sizeof(double));
    double *tmp_vec_forProc = (double*) malloc (rowNum_forProc * sizeof(double));
    double check;
    double prev_check;
    int divergenceCounter = 0;
    int iterationsCounter = 0;
    MPI_Scatterv(x_vec, vecElemForProc, vecIndForProc, MPI_DOUBLE, x_vec_forProc, vecElemForProc[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatterv(b_vec, vecElemForProc, vecIndForProc, MPI_DOUBLE, b_vec_forProc, vecElemForProc[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    matrix_vector_mul(A_mat_forProc, x_vec, y_vec_forProc, rowNum_forProc, N); // вычисляем y(0)
    vector_vector_sub(y_vec_forProc, b_vec_forProc, y_vec_forProc, rowNum_forProc, 0);
    MPI_Allgatherv(y_vec_forProc, vecElemForProc[rank], MPI_DOUBLE, y_vec, vecElemForProc, vecIndForProc, MPI_DOUBLE, MPI_COMM_WORLD);

    double b_vec_module = vector_module(b_vec); // вычисляем check
    double y_vec_module = vector_module(y_vec);
    check = y_vec_module / b_vec_module;

    while (check >= eps)
    {
        matrix_vector_mul(A_mat_forProc, y_vec, tmp_vec_forProc, rowNum_forProc, N); // вычисляем t(n)
        double scalar_mul_yTmp = scalar_mul(y_vec_forProc, tmp_vec_forProc, rowNum_forProc);
        double scalar_mul_tmpTmp = scalar_mul(tmp_vec_forProc, tmp_vec_forProc, rowNum_forProc);
        t = scalar_mul_yTmp / scalar_mul_tmpTmp;

        printf("PRINT THIS!\n");
        double_vector_mul(t, y_vec_forProc, y_vec_forProc); // вычисляем x(n+1)
        vector_vector_sub(x_vec_forProc, y_vec_forProc, x_vec_forProc, rowNum_forProc, 0);
        MPI_Allgatherv(x_vec_forProc, vecElemForProc[rank], MPI_DOUBLE, x_vec, vecElemForProc, vecIndForProc, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        printf("AND DONT PRINT THIS!\n");

        matrix_vector_mul(A_mat_forProc, x_vec, y_vec_forProc, rowNum_forProc, N); // вычисляем y(0)
        vector_vector_sub(y_vec_forProc, b_vec_forProc, y_vec_forProc, rowNum_forProc, 0);
        MPI_Allgatherv(y_vec_forProc, vecElemForProc[rank], MPI_DOUBLE, y_vec, vecElemForProc, vecIndForProc, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

        prev_check = check;
        y_vec_module = vector_module(y_vec); // перевычисляем число для критерия завершения счета
        check = y_vec_module / b_vec_module;

        if (check >= prev_check)
        {
            divergenceCounter++;
            if (divergenceCounter >= 5)
            {
                return -1;
            }
        } else
        {
            divergenceCounter = 0;
        }
        ++iterationsCounter;
//        if (rank == 0) printf("check = %0.12lf\n", check);
    }
//    free(x_vec_forProc);
//    free(b_vec_forProc);
//    free(y_vec_forProc);
//    free(tmp_vec_forProc);
    return iterationsCounter;
}

int main(int argc, char *argv[]) {
    double start, end;
    double *A_mat = NULL;
    double *A_mat_forProc = NULL;
    double *b_vec = NULL;
    double *x_vec = NULL;
    double *y_vec = NULL;
    double t;
    int *elemForProc = NULL;
    int *indForProc = NULL;

    MPI_Init(&argc, &argv);
    start = MPI_Wtime();
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    dataInitialization(&A_mat, &b_vec, &x_vec, &y_vec);
    dataDistribution(A_mat, &A_mat_forProc, &elemForProc, &indForProc);
//    if (rank == 0) printEquation(A_mat, b_vec);

    int iterationsCounter = calculate(A_mat_forProc, b_vec, x_vec, y_vec, t, elemForProc, indForProc);
    end = MPI_Wtime();
    if (rank == 0)
    {
        printf("Calculations have done in %d iterations.\nTime taken: %f sec\n", iterationsCounter, end-start);
        printf("Result is:\n");
        printVector(x_vec, "x_vec", N);
    }
//    dataDeletion(A_mat, A_mat_forProc, b_vec, x_vec, y_vec, elemForProc, indForProc);
    MPI_Finalize();

    return 0;
}
