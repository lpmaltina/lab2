#ifndef PARALLEL
#define PARALLEL

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>


int parallelProduct(
	double* matrix,
	double* vector,
	double* result,
	int ord);

#endif //PARALLEL
