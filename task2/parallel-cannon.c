#include "parallel.h"


void printMatSq(double const *const p_mat, int const ord);
double *mulMatsSq(double *p_matL, double *p_matR, int ord);
void addToMatSq(double *p_matRslt, double *p_matAddend, int ord);
double *rearngMatSq(double *const p_mat, int ord, int ordNew);

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
	int k = 0;

	int ordNew = 0;
	int rowsPerPrcsN = 0;

	int *p_elmsNs = NULL;
	int *p_elmsOfsts = NULL;
	double *p_matLNew = NULL;
	double *p_matRNew = NULL;
	double *p_matRNewT = NULL;
	double **pp_matRslts = NULL;
	double *p_matRslt = NULL;

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
		p_matLNew = rearngMatSq(p_matR, ord, ordNew);
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

		/* puts("p_matLNew shifted"); */
		/* printMatSq(p_matLNew, ordNew); */

		puts("p_matRNew");
		p_matRNew = rearngMatSq(p_matR, ord, ordNew);
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

		/* puts("p_matRNew shifted"); */
		/* printMatSq(p_matRNew, ordNew); */
	}

	p_l_matL = malloc((ordNew / sz) * (ordNew / sz) * sizeof *p_l_matL);
	p_l_matR = malloc((ordNew / sz) * (ordNew / sz) * sizeof *p_l_matR);

	pp_matRslts = calloc(sz, sizeof *pp_matRslts);
	for (k = 0; k < sz; k++) {
		pp_matRslts[k] = calloc(ordNew * ordNew, sizeof *pp_matRslts);
		for (j = 0; j < sz; j++) {
			for (i = 0; i < ordNew / sz; i++) {
				MPI_Scatter(
					p_matLNew + ordNew*(ordNew / sz)*j + ordNew*i, /* IN sendbuf */
					ordNew / sz,                  /* IN sendcount */
					MPI_DOUBLE,                   /* IN sendtype */
					p_l_matL + (ordNew / sz) * i, /* OUT recvbuf */
					ordNew / sz,                  /* IN recvcount */
					MPI_DOUBLE,                   /* IN recvtype */
					0,                            /* IN root */
					MPI_COMM_WORLD);              /* IN comm */
				MPI_Scatter(
					p_matRNew + ordNew*(ordNew / sz)*j + ordNew*i, /* IN sendbuf */
					ordNew / sz,                  /* IN sendcount */
					MPI_DOUBLE,                   /* IN sendtype */
					p_l_matR + (ordNew / sz) * i, /* OUT recvbuf */
					ordNew / sz,                  /* IN recvcount */
					MPI_DOUBLE,                   /* IN recvtype */
					0,                            /* IN root */
					MPI_COMM_WORLD);              /* IN comm */
			}

			p_l_matRslt = mulMatsSq(p_l_matL, p_l_matR, ordNew / sz);

			for (i = 0; i < ordNew / sz; i++) {
				MPI_Gather(
					p_l_matRslt + (ordNew / sz) * i, /* IN sendbuf */
					ordNew / sz,                     /* IN sendcount */
					MPI_DOUBLE,                      /* IN sendtype */
					pp_matRslts[k] + ordNew*(ordNew / sz)*j + ordNew*i, /* OUT recvbuf */
					ordNew / sz,                     /* IN recvcount */
					MPI_DOUBLE,                      /* IN recvtype */
					0,                               /* IN root */
					MPI_COMM_WORLD);                 /* IN comm */
			}

			if (rnk == 0) {
				puts("p_l_matL, p_l_matR, p_l_matRslt");
				printMatSq(p_l_matL, ordNew / sz);
				printMatSq(p_l_matR, ordNew / sz);
				printMatSq(p_l_matRslt, ordNew / sz);
			}
		}
		if (rnk == 0) {
			printf("pp_matRslts[%d]\n", k);
			printMatSq(pp_matRslts[k], ordNew);
			shftMatLtThr(p_matLNew, ordNew, ordNew, 0, ordNew, ordNew / sz);
			shftMatUpThr(p_matRNew, ordNew, ordNew, 0, ordNew, ordNew / sz);
		}
	}

	if (rnk == 0) {
		p_matRslt = calloc(ordNew * ordNew, sizeof *p_matRslt);
		for (i = 0; i < sz; i++) {
			addToMatSq(p_matRslt, pp_matRslts[i], ordNew);
		}
	
		puts("here am I");
		printMatSq(
			rearngMatSq(p_matRslt, ordNew, ord),
			ord);
	}

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
	return;
}

double *mulMatsSq(double *p_matL, double *p_matR, int ord)
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

void addToMatSq(double *p_matRslt, double *p_matAddend, int ord)
{
	for (int i = 0; i < ord; i++) {
		for (int j = 0; j < ord; j++) {
			p_matRslt[i * ord + j] += p_matAddend[i * ord + j];
		}
	}
}

double *rearngMatSq(double *const p_mat, int ord, int ordNew)
{
	double *p_matRslt = calloc(ordNew * ordNew, sizeof *p_matRslt);
	for (int i = 0; i < ord && i < ordNew; i++) {
		for (int j = 0; j < ord && j < ordNew; j++) {
			p_matRslt[i * ordNew + j] = p_mat[i * ord + j];
		}
	}
	return p_matRslt;
}

void printMatSq(double const *const p_mat, int const ord)
{
	for (int i = 0; i < ord; i++) {
		for (int j = 0; j < ord; j++) {
			printf("%.2lf    \t", p_mat[i * ord + j]);
		}
		printf("\n");
	}
	printf("\n");

	return;
}
