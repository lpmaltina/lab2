#ifndef PARALLEL
#define PARALLEL

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>


void parallelProduct(
    double* matrix,
    double* vector,
    double* result,
    int rows,
    int cols,
    int rank,
    int size
);

#endif //PARALLEL
