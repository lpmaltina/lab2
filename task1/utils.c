#include <stdio.h>
#include <stdlib.h>
#include "utils.h"

int readArray(char* fileName, double* array, int n) {
    int i;
    FILE* fp = fopen(fileName, "r");
    if (fp == NULL) { return 1; }

    for (i = 0; i < n; ++i){
        fscanf(fp, "%lf", &array[i]);
    }
    fclose(fp);
    return 0;
}

int writeArray(char* fileName, double* array, int n) {
    int i;
    FILE *fp = fopen(fileName, "w+");
    if (fp == NULL) { return 1; }

    for (i = 0; i < n; ++i) {
        fprintf(fp, "%.6f\n", array[i]);
    }

    fclose(fp);
    return 0;
}
