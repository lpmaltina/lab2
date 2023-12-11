#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void parallelProduct(double* matrix, double* vector, double* result, int rows, int cols, int rank, int size) {
    int i;
    int j;
    int myRows = rows / size;
    double* myMatrix = (double*) malloc(myRows * cols * sizeof(double));
    double* myResult = (double*) malloc(myRows * sizeof(double));

    MPI_Bcast(vector, cols, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(matrix, myRows * cols, MPI_DOUBLE, myMatrix, myRows * cols, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    for (i = 0; i < myRows; ++i) {
        myResult[i] = 0.0;
        for (j = 0; j < cols; ++j) {
            myResult[i] += myMatrix[i * cols + j] * vector[j];
        }
    }
    
    MPI_Gather(myResult, myRows, MPI_DOUBLE, result, myRows, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    free(myMatrix);
    free(myResult);
}
