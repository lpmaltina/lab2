#include "parallel.h"


/* MPI_Bcast( */
/* 	,  /\* INOUT buffer *\/ */
/* 	,  /\* IN count *\/ */
/* 	,  /\* IN datatype *\/ */
/* 	,  /\* IN root *\/ */
/* 	); /\* IN comm *\/ */

void printMatSqr(double const *const p_mat, int const ord);
double *transMatSqr(double const *const p_mat, int const ord);
double *incrMatSqr(
	double const *const p_mat, int const ord, int const ordNew);

/* Multiply matrix p_matL and p_matR, write result to p_rslt. */
int parallelProduct(
	double* p_matL, /* root */
	double* p_matR, /* root */
	double* p_rslt, /* root */
	int ord)        /* root */
{
	int i, j;

	int *p_elmsNs = NULL;
	int *p_elmsOfsts = NULL;
	int ordNew = 0;
	double *p_matLNew = NULL;
	double *p_matRNewT = NULL;

	double *p_lcl_matLElms = NULL;
	double *p_lcl_matRElms = NULL;

	int rank;
	int size;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	MPI_Bcast(&ord, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);


	/* Adapt matrices to number of processes. */
	ordNew = ord % size ?
		ord - (ord % size) + size :
		ord;

	if (rank == 0) {
		double *p_matRNew = incrMatSqr(p_matR, ord, ordNew);
		p_matLNew = incrMatSqr(p_matL, ord, ordNew);
		p_matRNewT = transMatSqr(p_matRNew, ordNew);

		printMatSqr(p_matR, ord);
		printMatSqr(p_matRNew, ordNew);
		/* p_matRNewT = transMatSqr(p_matRNew, ordNew); */
		/* free(p_matRNew); */

		/* printMatSqr(p_matL, ord); */
		/* printMatSqr(p_matRNewT, ordNew); */
	}

	/* Offset matrices */
	/* if (rank == 0) { */

	/* 	for (i = 0; i < size; i++) { */
			
	/* 	} */

	/* 	for (i = 0; i < ordNew; i++) { */
	/* 		int tmp =  */
	/* 	} */
	/* } */

	/* if (rank == 0) { */
		/* Offset left matrix by columns */
		/* for (i = 0; i < ord; ++i) { */
			
		/* } */

		/* Transpose right matrix */
		/* Offset right transposed matrix by columns */

	printf(
		"just for test:\nord     is %d\nordNew is %d\n\n",
		ord, ordNew);

	

	return 0;
}


double *transMatSqr(double const *const p_mat, int const ord)
{
	int i;
	int j;

	double *p_matRslt = malloc(ord * sizeof *p_matRslt);

	for (i = 0; i < ord; i++) {
		for (j = 0; j < ord; j++) {
			p_matRslt[(i * ord) + j] = p_mat[(j * ord) + i];
		}
	}
	return p_matRslt;
}

double *incrMatSqr(
	double const *const p_mat, int const ord, int const ordNew)
{
	double *p_matRslt = calloc(ordNew * ordNew, sizeof *p_matRslt);
	for (int i = 0; i < ord; i++) {
		for (int j = 0; j < ord; j++) {
			p_matRslt[i * ord + j] = p_mat[i * ordNew + j];
		}
	}
	return p_matRslt;
}

void printMatSqr(double const *const p_mat, int const ord)
{
	int i;
	int j;

	for (i = 0; i < ord; i++) {
		for (j = 0; j < ord; j++) {
			printf("%.1lf\t", p_mat[i * ord + j]);
		}
		printf("\n");
	}
	printf("\n");
}
