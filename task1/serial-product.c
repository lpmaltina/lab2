#include <stdio.h>
#include <stdlib.h>
#include "utils.h"

void serialProduct(
        double* matrix,
        double* vec,
        double* product,
        int matrixRows,
        int matrixCols
    ) {
    int i;
    int j;

    for (i = 0; i < matrixRows; ++i) {
        product[i] = 0.0;
        for (j = 0; j < matrixCols; ++j) {
            product[i] += matrix[i * matrixCols + j] * vec[j];
        }
    }
}

int main() {
    int matrixRows;
    int matrixCols;
    double* matrix;
    double* vec;
    double* product;
    int dim;
    int i;
    int j;
    char fileNameMatrix[30];
    char fileNameVec[30];
    char fileNameResult[40];

    for (dim = 512; dim <= 8192; dim *= 2) {
        matrixRows = dim;
        matrixCols = dim;

        matrix = (double*) malloc(matrixRows * matrixCols * sizeof(double));
        vec = (double*) malloc(matrixCols * sizeof(double));
        product = (double*) malloc(matrixRows * sizeof(double));
        
        sprintf(fileNameMatrix, "input/matrix-%d-%d.txt", matrixRows, matrixCols);
        readArray(fileNameMatrix, matrix, matrixRows * matrixCols);

        sprintf(fileNameVec, "input/vec-%d.txt", matrixCols);
        readArray(fileNameVec, vec, matrixCols);

        sprintf(fileNameResult, "output/serial-product-%d-%d.txt", matrixRows, matrixCols);
        serialProduct(matrix, vec, product, matrixRows, matrixCols);
        writeArray(fileNameResult, product, matrixRows);

        free(matrix);
        free(vec);
        free(product);
    }
    return 0;
}
