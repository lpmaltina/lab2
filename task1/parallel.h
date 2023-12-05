#ifndef PARALLEL
#define PARALLEL

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
