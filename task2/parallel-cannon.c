#include "parallel.h"


/* MPI_Bcast( */
/* 	,  /\* INOUT buffer *\/ */
/* 	,  /\* IN count *\/ */
/* 	,  /\* IN datatype *\/ */
/* 	,  /\* IN root *\/ */
/* 	); /\* IN comm *\/ */

void printMatSq(double const *const p_mat, int const ord);
void printMat(double const *const p_mat, int const rows, int const clms);
double *transMatSq(double const *const p_mat, int const ord);
double *transMat(double const *const p_mat, int const rows, int const clms);
double *incrMatSq(
	double const *const p_mat, int const ord, int const ordNew);
double* multMatsSq(double* p_matL, double* p_matR, int ord);


void shftMatLtThr(
	double *const p_mat,
	int rows,
	int clms,
	int rowBgg,
	int rowEnd,
	int shft);

void shftMatUpThr(
	double *const p_mat,
	int rows,
	int clms,
	int clmBgg,
	int clmEnd,
	int shft);

/* Multiply matrix p_matL and p_matR, write result to p_rslt. */
int parallelProduct(
	double* p_matL, /* root */
	double* p_matR, /* root */
	double* p_rslt, /* root */
	int ord)        /* root */
{
	int i = 0;
	int j = 0;

	int ordNew = 0;
	int rowsPerPrcsN = 0;

	int *p_elmsNs = NULL;
	int *p_elmsOfsts = NULL;
	double *p_matLNew = NULL;
	double *p_matRNew = NULL;
	double *p_matRNewT = NULL;

	double *p_l_matLElms = NULL;
	double *p_l_matRElms = NULL;
	double *p_l_matRElmsT = NULL;
	double *p_l_matL = NULL;
	double *p_l_matR = NULL;
	double *p_l_matRslt = NULL;

	int rnk = 0;
	int sz = 0;

	MPI_Comm_rank(MPI_COMM_WORLD, &rnk);
	MPI_Comm_size(MPI_COMM_WORLD, &sz);

	MPI_Bcast(&ord, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);


	/* Adapt matrices to number of processes. */
	ordNew = ord % sz ?
		ord - (ord % sz) + sz :
		ord;

	if (rnk == 0) {

		puts("p_matLNew");
		p_matLNew = incrMatSq(p_matR, ord, ordNew);
		printMatSq(p_matLNew, ordNew);

		for (i = 0; i < sz; i++) {
			shftMatLtThr(
				p_matLNew,
				ordNew,
				ordNew,
				(ordNew / sz) * i,
				(ordNew / sz) * (i + 1),
				(ordNew / sz) * i);
		}

		puts("p_matLNew shifted");
		printMatSq(p_matLNew, ordNew);

		puts("p_matRNew");
		p_matRNew = incrMatSq(p_matR, ord, ordNew);
		printMatSq(p_matRNew, ordNew);

		for (i = 0; i < sz; i++) {
			shftMatUpThr(
				p_matRNew,
				ordNew,
				ordNew,
				(ordNew / sz) * i,
				(ordNew / sz) * (i + 1),
				(ordNew / sz) * i);
		}

		puts("p_matRNew shifted");
		printMatSq(p_matRNew, ordNew);
	}

	p_l_matL = malloc((ordNew / sz) * (ordNew / sz) * sizeof *p_l_matL);
	p_l_matR = malloc((ordNew / sz) * (ordNew / sz) * sizeof *p_l_matR);

	for (i = 0; i < ordNew / sz; i++) {
		MPI_Scatter(
			p_matLNew + ordNew * i,       /* IN sendbuf */
			ordNew / sz,                  /* IN sendcount */
			MPI_DOUBLE,                   /* IN sendtype */
			p_l_matL + (ordNew / sz) * i, /* OUT recvbuf */
			ordNew / sz,                  /* IN recvcount */
			MPI_DOUBLE,                   /* IN recvtype */
			0,                            /* IN root */
			MPI_COMM_WORLD);              /* IN comm */
		MPI_Scatter(
			p_matRNew + ordNew * i,       /* IN sendbuf */
			ordNew / sz,                  /* IN sendcount */
			MPI_DOUBLE,                   /* IN sendtype */
			p_l_matR + (ordNew / sz) * i, /* OUT recvbuf */
			ordNew / sz,                  /* IN recvcount */
			MPI_DOUBLE,                   /* IN recvtype */
			0,                            /* IN root */
			MPI_COMM_WORLD);              /* IN comm */
	}

	p_l_matRslt = multMatsSq(p_l_matL, p_l_matR, ordNew / sz);

	/* if (rnk == 0) { */
	/* 	printMat(p_l_matL, ordNew / sz, ordNew / sz); */
	/* 	printMat(p_l_matR, ordNew / sz, ordNew / sz); */
	/* 	printMat(p_l_matRslt, ordNew / sz, ordNew / sz); */
	/* } */

	/* MPI_Scatter( */
	/* 	p_matLNew,              /\* IN sendbuf *\/ */
	/* 	ordNew * (ordNew / sz), /\* IN sendcount *\/ */
	/* 	MPI_DOUBLE,             /\* IN sendtype *\/ */
	/* 	p_l_matLElms,           /\* OUT recvbuf *\/ */
	/* 	ordNew * (ordNew / sz), /\* IN recvcount *\/ */
	/* 	MPI_DOUBLE,             /\* IN recvtype *\/ */
	/* 	0,                      /\* IN root *\/ */
	/* 	MPI_COMM_WORLD);        /\* IN comm *\/ */

	/* MPI_Scatter( */
	/* 	p_matRNewT,             /\* IN sendbuf *\/ */
	/* 	ordNew * (ordNew / sz), /\* IN sendcount *\/ */
	/* 	MPI_DOUBLE,             /\* IN sendtype *\/ */
	/* 	p_l_matRElmsT,           /\* OUT recvbuf *\/ */
	/* 	ordNew * (ordNew / sz), /\* IN recvcount *\/ */
	/* 	MPI_DOUBLE,             /\* IN recvtype *\/ */
	/* 	0,                      /\* IN root *\/ */
	/* 	MPI_COMM_WORLD);        /\* IN comm *\/ */

	/* Every process got its macro row and macro column, but macro */
	/* column shall be transposed to act in matrix multiplying. */

	/* p_l_matRElms = transMat(p_l_matRElmsT, ordNew / sz, ordNew); */

	/* if (rnk == 0) { */
	/* 	puts("p_l_matLElms"); */
	/* 	printMat(p_l_matLElms, ordNew / sz, ordNew); */

	/* 	puts("p_l_matRElms"); */
	/* 	printMat(p_l_matRElms, ordNew, ordNew / sz); */
	/* } */

	/* for (i = 0; i < (rnk + 1) * 2000000000; i++); */

	/* printf("p_l_matRElmsT %d\n", rnk); */
	/* printMat(p_l_matRElmsT, ordNew / sz, ordNew); */

	/* printf("p_l_matRElmsT %d shfted\n", rnk); */
	/* shftMatLtThr(p_l_matRElmsT, ordNew / sz, ordNew, (ordNew / sz) * rnk); */
	/* printMat(p_l_matRElmsT, ordNew / sz, ordNew); */

	/* puts("p_l_matRElms"); */
	/* printMat(p_l_matRElms, ordNew, ordNew / sz); */

	/* shftMatLtThr(p_l_matRElms, ordNew, ordNew / sz, rnk); */

	/* if (rnk == 0) { */
	/* 	for (i = 0; i < 2000000000; i++); */
	/* 	puts("p_l_matRElms shfted"); */
	/* 	printMat(p_l_matRElms, ordNew, ordNew / sz); */
	/* } */

	return 0;
}


void shftMatLtThr(
	double *const p_mat,
	int rows,
	int clms,
	int rowBgg,
	int rowEnd,
	int shft)
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
	return;
}

void shftMatUpThr(
	double *const p_mat,
	int rows,
	int clms,
	int clmBgg,
	int clmEnd,
	int shft)
{
	double* p_matTmp = NULL;
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
	return;
}

double* multMatsSq(double* p_matL, double* p_matR, int ord)
{
	double *p_matRslt = malloc(ord * ord * sizeof *p_matRslt);

	for (int i = 0; i < ord; i++) {
		for (int j = 0; j < ord; j++) {
			double sum = 0.;
			for (int k = 0; k < ord; k++) {
				sum += p_matL[i * ord + k] * p_matR[k * ord + j];
			}
			p_matRslt[i * ord + j] = sum;
		}
	}

	return p_matRslt;
}

double *transMat(double const *const p_mat, int const rows, int const clms)
{
	double *const p_matRslt = malloc(rows * clms * sizeof *p_matRslt);
	for (int i = 0; i < clms; i++) {
		for (int j = 0; j < rows; j++) {
			p_matRslt[i * rows + j] = p_mat[j * clms + i];
		}
	}
	return p_matRslt;
}

double *transMatSq(double const *const p_mat, int const ord)
{
	return transMat(p_mat, ord, ord);
}

double *incrMatSq(
	double const *const p_mat, int const ord, int const ordNew)
{
	double *p_matRslt = calloc(ordNew * ordNew, sizeof *p_matRslt);
	for (int i = 0; i < ord; i++) {
		for (int j = 0; j < ord; j++) {
			p_matRslt[i * ordNew + j] = p_mat[i * ord + j];
		}
	}
	return p_matRslt;
}

void printMat(double const *const p_mat, int const rows, int const clms)
{
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < clms; j++) {
			printf("%.1lf\t", p_mat[i * clms + j]);
		}
		printf("\n");
	}
	printf("\n");

	return;
}

void printMatSq(double const *const p_mat, int const ord)
{
	printMat(p_mat, ord, ord);
	return;
}
