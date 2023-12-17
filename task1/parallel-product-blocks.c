#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void parallelProductBlock(double* matrix, double* vector, double* result, int rows, int cols, int rank, int size, int blockSize) {
    int i, j, b;
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

    int blockMatrixSize = rows * blockSize;

    double* myMatrixBlock = (double*)malloc(blockMatrixSize * sizeof(double));
    double* myVectorBlock = (double*)malloc(blockSize * sizeof(double));
    MPI_Bcast(vector, cols, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (b = myStartBlock; b < myEndBlock; b++) {
        int myBlockOffset = b * blockSize;
        MPI_Scatter(matrix + myBlockOffset, blockSize, MPI_DOUBLE, myMatrixBlock, blockSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        for (i = 0; i < rows; i++) {
            for (j = 0; j < blockSize; j++) {
                int matrixIndex = i * cols + myBlockOffset + j;
                if (matrixIndex < rows * cols) {
                    myMatrixBlock[i * blockSize + j] = matrix[matrixIndex];
                } else {
                    myMatrixBlock[i * blockSize + j] = 0.0;
                }
            }
        }

        for (i = 0; i < blockSize; i++) {
            int vectorIndex = myBlockOffset + i;
            if (vectorIndex < cols) {
                myVectorBlock[i] = vector[vectorIndex];
            } else {
                myVectorBlock[i] = 0.0;
            }
        }

        double* myResultBlock = (double*)malloc(rows * sizeof(double));
        for (i = 0; i < rows; i++) {
            myResultBlock[i] = 0.0;
            for (j = 0; j < blockSize; j++) {
                myResultBlock[i] += myMatrixBlock[i * blockSize + j] * myVectorBlock[j];
            }
        }

        int* recvcounts = (int*)malloc(size * sizeof(int));
        int* displs = (int*)malloc(size * sizeof(int));
        for (i = 0; i < size; i++) {
            int startBlock = i * numBlocks / size;
            int endBlock = (i + 1) * numBlocks / size;
            if (i == size - 1) {
                endBlock = numBlocks;
            }
            recvcounts[i] = rows * (endBlock - startBlock);
            displs[i] = rows * startBlock;
        }

        MPI_Gatherv(myResultBlock, rows, MPI_DOUBLE, result, recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        free(myResultBlock);
        free(recvcounts);
        free(displs);
    }

    free(myMatrixBlock);
    free(myVectorBlock);
}
