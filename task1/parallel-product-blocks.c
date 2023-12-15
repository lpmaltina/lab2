#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void parallelProduct(double* matrix, double* vector, double* result, int rows, int cols, int rank, int size) {
    int blockSize = rows / size;
    int remainingRows = rows % size;
    int myRows = rank < remainingRows ? blockSize + 1 : blockSize;

    double* myMatrix = (double*)malloc(myRows * cols * sizeof(double));
    double* myResult = (double*)malloc(myRows * sizeof(double));

    int* sendCounts = (int*)malloc(size * sizeof(int));
    int* displacements = (int*)malloc(size * sizeof(int));

    if (rank == 0) {
        int offset = 0;
        for (int i = 0; i < size; i++) {
            sendCounts[i] = (i < remainingRows ? blockSize + 1 : blockSize) * cols;
            displacements[i] = offset;
            offset += sendCounts[i];
        }
    }

    MPI_Scatterv(matrix, sendCounts, displacements, MPI_DOUBLE, myMatrix, myRows * cols, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Bcast(vector, cols, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (int i = 0; i < myRows; ++i) {
        myResult[i] = 0.0;
        for (int j = 0; j < cols; ++j) {
            myResult[i] += myMatrix[i * cols + j] * vector[j];
        }
    }

    if (rank == 0) {
        for (int i = 0; i < size; i++) {
            sendCounts[i] = i < remainingRows ? blockSize + 1 : blockSize;
            displacements[i] = (i > 0) ? displacements[i - 1] + sendCounts[i - 1] : 0;
        }
    }

    MPI_Gatherv(myResult, myRows, MPI_DOUBLE, result, sendCounts, displacements, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    free(myMatrix);
    free(myResult);
    free(sendCounts);
    free(displacements);
}