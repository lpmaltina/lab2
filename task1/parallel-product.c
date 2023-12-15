#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "utils.h"

void parallelProduct(double* matrix, double* vec, double* result, int matrixRows, int matrixCols, int rank, int size) {

    int blockRows = matrixRows / size;
    int remainderRows = matrixRows % size;
    int localBlockRows = blockRows;
    if (rank < remainderRows) {
        localBlockRows++;
    }

    double* localMatrixBlock = (double*)malloc(sizeof(double) * localBlockRows * matrixCols);
    int* displacements = (int*)malloc(sizeof(int) * size);
    int* counts = (int*)malloc(sizeof(int) * size);

    int offset = 0;
    for (int i = 0; i < size; ++i) {
        displacements[i] = offset * matrixCols;
        counts[i] = localBlockRows * matrixCols;
        if (i < remainderRows) {
            counts[i] += matrixCols;
        }
        offset += counts[i] / matrixCols;
    }

    MPI_Scatterv(matrix, counts, displacements, MPI_DOUBLE, localMatrixBlock, localBlockRows * matrixCols, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double* localVector = (double*)malloc(sizeof(double) * matrixCols);
    MPI_Scatter(vec, matrixCols, MPI_DOUBLE, localVector, matrixCols, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (int i = 0; i < localBlockRows; ++i) {
        result[i] = 0.0;
        for (int j = 0; j < matrixCols; ++j) {
            result[i] += localMatrixBlock[i * matrixCols + j] * localVector[j];
        }
    }

    double* allResults = (double*)malloc(sizeof(double) * matrixRows);
    int* recvCounts = (int*)malloc(sizeof(int) * size);
    int* displacementsGather = (int*)malloc(sizeof(int) * size);

    for (int i = 0; i < size; ++i) {
        recvCounts[i] = blockRows + (i < remainderRows ? 1 : 0);
        displacementsGather[i] = i * blockRows;
    }

    MPI_Gatherv(result, localBlockRows, MPI_DOUBLE, allResults, recvCounts, displacementsGather, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        for (int i = 0; i < matrixRows; ++i) {
            result[i] = allResults[i];
        }
    }

    free(localMatrixBlock);
    free(localVector);
    free(allResults);
    free(displacements);
    free(counts);
    free(recvCounts);
    free(displacementsGather);
}