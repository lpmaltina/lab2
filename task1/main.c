#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "utils.h"
#include "parallel.h"

int main(int argc, char** argv) {
    double* matrix;
    double* vec;
    double* product;
    int rank;
    int size;
    int matrixRows;
    int matrixCols;
    char fileNameMatrix[30];
    char fileNameVec[30];
    char fileNameResult[40];
    double start = 0;
    double finish = 0;
    double elapsed = 0;
    double time_elapsed = 0;
    int dim = strtol(argv[1], NULL, 10);

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    matrixRows = dim;
    matrixCols = dim;
    matrix = (double*) malloc(matrixRows * matrixCols * sizeof(double));
    vec = (double*) malloc(matrixCols * sizeof(double));
    product = (double*) malloc(matrixCols * sizeof(double));

    if (rank == 0){
        sprintf(fileNameMatrix, "input/matrix-%d-%d.txt", matrixRows, matrixCols);
        readArray(fileNameMatrix, matrix, matrixRows * matrixCols);

        sprintf(fileNameVec, "input/vec-%d.txt", matrixCols);
        readArray(fileNameVec, vec, matrixCols);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();
    parallelProduct(matrix, vec, product, matrixRows, matrixCols, rank, size);
    finish = MPI_Wtime();
    MPI_Barrier(MPI_COMM_WORLD);

    time_elapsed = finish - start;
    MPI_Reduce(&time_elapsed, &elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0){
        sprintf(fileNameResult, "output/parallel-product-%d-%d.txt", matrixRows, matrixCols);
        writeArray(fileNameResult, product, matrixRows);
        printf(
                "Time (dim=%d, n_threads=%d): %lf\n",
                dim,
                size,
                time_elapsed
        );
    }

    free(matrix);
    free(vec);
    free(product);

    MPI_Finalize();
    return 0;
}