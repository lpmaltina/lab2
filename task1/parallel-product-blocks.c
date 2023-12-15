#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void parallelProduct(double* matrix, double* vector, double* result, int rows, int cols, int rank, int size, int blockSize) {
    int i;
    int j;
    int b;
    int myStartBlock, myEndBlock;
    int numBlocks = cols / blockSize;
    if (cols % blockSize != 0) {
        numBlocks++;
    }

    myStartBlock = rank * numBlocks / size;
    myEndBlock = (rank + 1) * numBlocks / size;
    if (rank == size - 1) {
        myEndBlock = numBlocks;
    }

    double* myMatrixBlock = (double*)malloc(rows * blockSize * sizeof(double));
    double* myVectorBlock = (double*)malloc(blockSize * sizeof(double));

    for (int b = myStartBlock; b < myEndBlock; b++) {
        int myBlockOffset = b * blockSize;
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < blockSize; j++) {
                int matrixIndex = i * cols + myBlockOffset + j;
                if (matrixIndex < rows * cols) {
                    myMatrixBlock[i * blockSize + j] = matrix[matrixIndex];
                } else {
                    myMatrixBlock[i * blockSize + j] = 0.0;
                }
            }
        }

        for (int i = 0; i < blockSize; i++) {
            int vectorIndex = myBlockOffset + i;
            if (vectorIndex < cols) {
                myVectorBlock[i] = vector[vectorIndex];
            } else {
                myVectorBlock[i] = 0.0;
            }
        }
    }

    double* myResultBlock = (double*)malloc(rows * sizeof(double));
    for (int i = 0; i < rows; i++) {
        myResultBlock[i] = 0.0;
        for (int j = 0; j < blockSize; j++) {
            myResultBlock[i] += myMatrixBlock[i * blockSize + j] * myVectorBlock[j];
        }
    }

    int* recvcounts = (int*)malloc(size * sizeof(int));
    int* displs = (int*)malloc(size * sizeof(int));
    for (int i = 0; i < size; i++) {
        int startBlock = i * numBlocks / size;
        int endBlock = (i + 1) * numBlocks / size;
        if (i == size - 1) {
            endBlock = numBlocks;
        }
        recvcounts[i] = rows * (endBlock - startBlock);
        displs[i] = rows * startBlock;
    }

    MPI_Gatherv(myResultBlock, rows, MPI_DOUBLE, result, recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    free(myMatrixBlock);
    free(myVectorBlock);
    free(myResultBlock);
    free(recvcounts);
    free(displs);
}