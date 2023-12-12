#include <stdio.h>
#include <stdlib.h>
#include "utils.h"

void readArray(char* fileName, double* array, int n) {
    int i;
    FILE* fp = fopen(fileName, "r");
    for (i = 0; i < n; ++i){
        fscanf(fp, "%lf", &array[i]);
    }
    fclose(fp);
}

void writeArray(char* fileName, double* array, int n) {
    FILE *fp = fopen(fileName, "w+");
    
    int i;
    for (i = 0; i < n; ++i) {
        fprintf(fp, "%.6f\n", array[i]);
    }

    fclose(fp);
}
