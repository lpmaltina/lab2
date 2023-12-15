#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "utils.h"

double generateDouble(int min, int max, int n){
    int denominator = pow(10, n);
    int maxNum = max * denominator;
    int minNum = min * denominator;
    return (minNum + rand() % ((maxNum - minNum) + 1)) / (double) denominator;
}

void randomFill(double* array, int n){
    int i;
    for (i = 0; i < n; ++i){
        array[i] = generateDouble(-100, 100, 1);
    }
}

int main(){
    srand(time(NULL));
    int matrixRows;
    int matrixCols;
    int dim;
    double* matrix;
    double* vec;
    char fileNameMatrix[30];
    char fileNameVec[30];

    for (dim = 512; dim <= 8192; dim *= 2) {
        matrixRows = dim;
        matrixCols = dim;

        matrix = (double*) malloc(matrixRows * matrixCols * sizeof(double*));
        vec = (double*) malloc(matrixCols * sizeof(double));
        
        randomFill(matrix, matrixRows * matrixCols);
        randomFill(vec, matrixCols);

        sprintf(fileNameMatrix, "input/matrix-%d-%d.txt", matrixRows, matrixCols);
        writeArray(fileNameMatrix, matrix, matrixRows * matrixCols);

        sprintf(fileNameVec, "input/vec-%d.txt", matrixCols);
        writeArray(fileNameVec, vec, matrixCols);

        free(matrix);
        free(vec);
    }
    return 0;
}
