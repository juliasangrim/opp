#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define columnsA (4 * 2)
#define rowsA (8 * 2)
#define columnsB (3 * 3)
#define rowsB (4 * 2)

#define procRows 8
#define procColumns 3

void printMatrix(double* matrix, int rows, int columns, char* name)
{
    printf("%s\n", name);
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < columns; ++j)
        {
            printf("%f ", matrix[i * columns + j]);
        }
        printf("\n");
    }
}

void ininMatrix(double* matrix, int rows, int columns)
{
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < columns; ++j)
        {
            matrix[i * columns + j] = (double)(rand() % (200 - 100));
        }
    }
}

void mul(const double* matrixA, const double* matrixB, double* matrixC, int rowsNumA, int rowsNumB, int columnsNumA, int columnsNumB)
{
    for (int i = 0; i < rowsNumA; ++i)
    {
        for (int j = 0; j < columnsNumB; ++j)
        {
            for (int k = 0; k < rowsNumB; ++k)
            {
                matrixC[i * columnsNumB + j] += matrixA[i * columnsNumA + k] * matrixB[k * columnsNumB + j];
            }
        }
    }
}

void destributeMatrix(double* matrixA, double* matrixB, double* vecRowsA, double* vecColumnsB, MPI_Comm commRows, MPI_Comm commColumns, int rankProc) {
    int procRankRows;
    MPI_Comm_rank(commRows, &procRankRows);
    int rowsNumForProc = rowsA / procRows;
    if (procRankRows == 0)
    {
        MPI_Scatter(matrixA, columnsA * rowsNumForProc, MPI_DOUBLE, vecRowsA,
                    columnsA * rowsNumForProc, MPI_DOUBLE, 0, commColumns);
    }
    MPI_Bcast(vecRowsA, columnsA * rowsNumForProc, MPI_DOUBLE, 0, commRows);

    MPI_Datatype column, columnType;
    int columnNumForProc = columnsB / procColumns;
    MPI_Type_vector(rowsB, columnNumForProc, columnsB, MPI_DOUBLE, &column);
    MPI_Type_commit(&column);
    MPI_Type_create_resized(column, 0, columnNumForProc * sizeof(double), &columnType);
    MPI_Type_commit(&columnType);
    int procRankColumns;
    MPI_Comm_rank(commColumns, &procRankColumns);

    if (procRankColumns == 0)
    {
        MPI_Scatter(matrixB, 1, columnType, vecColumnsB,
                    rowsB * columnNumForProc, MPI_DOUBLE, 0, commRows);
    }
    MPI_Bcast(vecColumnsB, rowsB * columnNumForProc, MPI_DOUBLE, 0, commColumns);
}

void gatherMatrix(double* matrixC, double* partC)
{
    int rowsNumForProc = rowsA / procRows;
    int columnNumForProc = columnsB / procColumns;
    int procNum;
    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
       MPI_Datatype block, blockType;
    MPI_Type_vector(rowsNumForProc, columnNumForProc, columnsB, MPI_DOUBLE, &block);
    MPI_Type_commit(&block);
    MPI_Type_create_resized(block, 0, columnNumForProc * sizeof(double), &blockType);
    MPI_Type_commit(&blockType);
    int* recvCount = (int*)malloc(procNum * sizeof(int));
    int* shifts = (int*)malloc(procNum * sizeof(int));
    int count = columnsB;
    int shift = 0;
    int countSkippedLine = 1;
    for (int i = 0; i < procNum; ++i)
    {
        recvCount[i] = 1;
        if (count == 0)
        {
            count = columnsB;
            shift = countSkippedLine * (columnsB / columnNumForProc) * rowsNumForProc;
            countSkippedLine++;
        }
        shifts[i] = shift;
        count -= columnNumForProc;
        shift++;
       // printf("%d ", shifts[i]);
    }
    MPI_Gatherv(partC, rowsNumForProc * columnNumForProc, MPI_DOUBLE, matrixC,
    recvCount, shifts, blockType, 0, MPI_COMM_WORLD);
}

void calculate(double* matrixA, double* matrixB, double* matrixC, int rankProc)
{
    int rowsNumForProc = rowsA / procRows;
    int columnNumForProc = columnsB / procColumns;
    MPI_Comm commRows, commColumns;
    double* vecRowsA = (double*)calloc(rowsNumForProc * columnsA, sizeof(double));
    double* vecColumnsB = (double*)calloc(columnNumForProc * rowsB, sizeof(double));
    int colour = 0;
    colour = rankProc / procColumns;
    MPI_Comm_split(MPI_COMM_WORLD, colour, rankProc, &commRows);
    colour = rankProc % procColumns;
    MPI_Comm_split(MPI_COMM_WORLD, colour, rankProc, &commColumns);

    destributeMatrix(matrixA, matrixB, vecRowsA, vecColumnsB, commRows, commColumns, rankProc);

    double* partMatrixC = (double*)calloc(rowsNumForProc * columnNumForProc, sizeof(double));
    mul(vecRowsA, vecColumnsB, partMatrixC, rowsNumForProc, rowsB, columnsA, columnNumForProc);
    gatherMatrix(matrixC, partMatrixC);
}


int main(int argc, char* argv[])
{
    srand(1);
    MPI_Init(&argc, &argv);
    int procNum, procRank;
    double* matrixA = NULL;
    double* matrixB = NULL;
    double* matrixC = NULL;
    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    if (procRank == 0) {
        matrixA = (double *) malloc(columnsA * rowsA * sizeof(double));
        matrixB = (double *) malloc(columnsB * rowsB * sizeof(double));
        matrixC = (double *) calloc(rowsA * columnsB, sizeof(double));
        ininMatrix(matrixA, rowsA, columnsA);
        ininMatrix(matrixB, rowsB, columnsB);
        printMatrix(matrixA, rowsA, columnsA, "A");
        printMatrix(matrixB, rowsB, columnsB, "B");
    }
    calculate(matrixA, matrixB, matrixC, procRank);
    if (procRank == 0) {
        printMatrix(matrixC, rowsA, columnsB, "C");
        free(matrixA);
        free(matrixB);
        free(matrixC);
    }
    return 0;
}
