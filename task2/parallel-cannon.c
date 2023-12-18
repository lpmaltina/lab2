#include "parallel.h"


void printMatSq(double const *const p_mat, int const ord);
void mulMatsSq(double *p_matRslt, double *p_matL, double *p_matR, int ord);
void addToMatSq(double *p_matRslt, double *p_matAddend, int ord);
void rearngMatSq(double *p_matRslt, double *p_mat, int ordNew, int ord);
void shftMatLtThr(
	double *p_mat, int rows, int clms, int rowBgg, int rowEnd, int shft);
void shftMatUpThr(
	double *p_mat, int rows, int clms, int clmBgg, int clmEnd, int shft);


/* Multiply matrix p_matL and p_matR, write result to p_rslt. */
int parallelProduct(
	double* p_matL_input,
	double* p_matR_input,
	double* p_rslt,
	int matOrd) /* matrixes' order */
{
	int i = 0;
	int j = 0;
	int k = 0;

	int ord = 0;

	double *p_matL = NULL;
	double *p_matR = NULL;
	double **pp_matRslts = NULL;
	double *p_matRslt = NULL;

	double *p_l_matL = NULL;
	double *p_l_matR = NULL;
	double *p_l_matRslt = NULL;

	int rnk = 0;
	int sz = 0;

	MPI_Comm_rank(MPI_COMM_WORLD, &rnk);
	MPI_Comm_size(MPI_COMM_WORLD, &sz);

	MPI_Bcast(&matOrd, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);

	ord = matOrd % sz ?
		matOrd - (matOrd % sz) + sz :
		matOrd;

	/* Prepare matrices for Cannon: rearrange to make more suitable for */
	/* number of processes, then shift em as Cannon said. */
	if (rnk == 0) {

		p_matL = malloc(ord * ord * sizeof *p_matL);
		rearngMatSq(p_matL, p_matL_input, ord, matOrd);
		for (i = 0; i < sz; i++) {
			shftMatLtThr(
				p_matL,
				ord,
				ord,
				(ord/sz) * i,
				(ord/sz) * (i + 1),
				(ord/sz) * i);
		}

		p_matR = malloc(ord * ord * sizeof *p_matR);
		rearngMatSq(p_matR, p_matR_input, ord, matOrd);
		for (i = 0; i < sz; i++) {
			shftMatUpThr(
				p_matR,
				ord,
				ord,
				(ord/sz) * i,
				(ord/sz) * (i + 1),
				(ord/sz) * i);
		}
	}

	p_l_matL = malloc((ord/sz) * (ord/sz) * sizeof *p_l_matL);
	p_l_matR = malloc((ord/sz) * (ord/sz) * sizeof *p_l_matR);
	p_l_matRslt = malloc((ord/sz) * (ord/sz) * sizeof *p_l_matRslt);

	pp_matRslts = malloc(sz * sizeof *pp_matRslts);
	for (i = 0; i < sz; i++) {
		pp_matRslts[i] = malloc(ord * ord * sizeof *pp_matRslts);
	}

	/* Compute all result matrices (that ones to sum together). */
	for (k = 0; k < sz; k++) {

		/* Compute yet another result matrix. Iterate over rows of */
		/* submatrices in input matrices (not submatrices' rows). */
		for (j = 0; j < sz; j++) {

			/* Scatter yet another two rows of submatrices in */
			/* input matrices (not submatrices' row). Iterate */
			/* over submatrices' rows. */
			for (i = 0; i < ord/sz; i++) {

				/* Scatter submatrices' rows. Mean */
				/* submatrices of left input matrix. */
				MPI_Scatter(
					p_matL + ord*(ord/sz)*j + ord*i, /* IN sendbuf */
					ord/sz,                          /* IN sendcount */
					MPI_DOUBLE,                      /* IN sendtype */
					p_l_matL + (ord/sz) * i,         /* OUT recvbuf */
					ord/sz,                          /* IN recvcount */
					MPI_DOUBLE,                      /* IN recvtype */
					0,                               /* IN root */
					MPI_COMM_WORLD);                 /* IN comm */

				/* Scatter submatrices' rows. Mean */
				/* submatrices of right input matrix. */
				MPI_Scatter(
					p_matR + ord*(ord/sz)*j + ord*i, /* IN sendbuf */
					ord/sz,                          /* IN sendcount */
					MPI_DOUBLE,                      /* IN sendtype */
					p_l_matR + (ord/sz) * i,         /* OUT recvbuf */
					ord/sz,                          /* IN recvcount */
					MPI_DOUBLE,                      /* IN recvtype */
					0,                               /* IN root */
					MPI_COMM_WORLD);                 /* IN comm */
			}

			mulMatsSq(p_l_matRslt, p_l_matL, p_l_matR, ord/sz);

			/* Gather row of result submatrices (not */
			/* submatrices' row). */
			for (i = 0; i < ord/sz; i++) {
				/* Gather result submatrices' rows. */
				MPI_Gather(
					p_l_matRslt + (ord/sz) * i,              /* IN sendbuf */
					ord/sz,                                  /* IN sendcount */
					MPI_DOUBLE,                              /* IN sendtype */
					pp_matRslts[k] + ord*(ord/sz)*j + ord*i, /* OUT recvbuf */
					ord/sz,                                  /* IN recvcount */
					MPI_DOUBLE,                              /* IN recvtype */
					0,                                       /* IN root */
					MPI_COMM_WORLD);                         /* IN comm */
			}
		}
		if (rnk == 0) {
			shftMatLtThr(p_matL, ord, ord, 0, ord, ord/sz);
			shftMatUpThr(p_matR, ord, ord, 0, ord, ord/sz);
		}
	}

	/* Sum result matrices to get final result matrix. */
	if (rnk == 0) {
		p_matRslt = calloc(ord * ord, sizeof *p_matRslt);
		for (i = 0; i < sz; i++) {
			addToMatSq(p_matRslt, pp_matRslts[i], ord);
		}
		memset(p_rslt, 0, matOrd * matOrd * sizeof *p_rslt);
		rearngMatSq(p_rslt, p_matRslt, matOrd, ord);

		/* puts("here am I"); */
		/* printMatSq(p_rslt, matOrd); */
	}

	for (i = 0; i < sz; i++) {
		free(pp_matRslts[i]);
	}

	free(p_matL);
	free(p_matR);
	free(p_matRslt);
	free(p_l_matL);
	free(p_l_matR);
	free(p_l_matRslt);

	return 0;
}


void shftMatLtThr(
	double *p_mat, int rows, int clms, int rowBgg, int rowEnd, int shft)
{
	double* p_matTmp = NULL;
	int shftNew = shft % clms;

	if (shftNew == 0) { return; }

	p_matTmp = malloc(rows * clms * sizeof *p_matTmp);
	memcpy(p_matTmp, p_mat, rows * clms * sizeof *p_matTmp);

	for (int i = rowBgg; i < rowEnd; i++) {
		for (int j = 0; j < clms; j++) {
			int jNew = (j - shftNew + clms) % clms;
			p_mat[i * clms + jNew] = p_matTmp[i * clms + j];
		}
	}

	free(p_matTmp);
}

void shftMatUpThr(
	double *p_mat, int rows, int clms, int clmBgg, int clmEnd, int shft)
{
	double *p_matTmp = NULL;
	int shftNew = -shft % rows;

	if (shftNew == 0) { return; }

	p_matTmp = malloc(rows * clms * sizeof *p_matTmp);
	memcpy(p_matTmp, p_mat, rows * clms * sizeof *p_matTmp);

	for (int i = 0; i < rows; i++) {
		int iNew = (i - shftNew + rows) % rows;
		for (int j = clmBgg; j < clmEnd; j++) {
			p_mat[i * clms + j] = p_matTmp[iNew * clms + j];
		}
	}

	free(p_matTmp);
}

void mulMatsSq(double *p_matRslt, double *p_matL, double *p_matR, int ord)
{
	for (int i = 0; i < ord; i++) {
		for (int j = 0; j < ord; j++) {
			double sum = 0.;
			for (int k = 0; k < ord; k++) {
				sum += p_matL[i * ord + k] * p_matR[k * ord + j];
			}
			p_matRslt[i * ord + j] = sum;
		}
	}
}

void addToMatSq(double *p_matRslt, double *p_matAddend, int ord)
{
	for (int i = 0; i < ord; i++) {
		for (int j = 0; j < ord; j++) {
			p_matRslt[i * ord + j] += p_matAddend[i * ord + j];
		}
	}
}

void rearngMatSq(double *p_matRslt, double *const p_mat, int ordNew, int ord)
{
	memset(p_matRslt, 0, ordNew * ordNew * sizeof *p_matRslt);
	for (int i = 0; i < ord && i < ordNew; i++) {
		for (int j = 0; j < ord && j < ordNew; j++) {
			p_matRslt[i * ordNew + j] = p_mat[i * ord + j];
		}
	}
}

void printMatSq(double const *const p_mat, int const ord)
{
	for (int i = 0; i < ord; i++) {
		for (int j = 0; j < ord; j++) {
			printf("%.0lf \t", p_mat[i * ord + j]);
		}
		printf("\n");
	}
	printf("\n");
}
