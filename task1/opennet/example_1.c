#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpi.h"


#define MAT_IN_ROWS_N 4
#define MAT_IN_CLMS_N 4
#define MAT_IN ((double []) { \
	-37.1, 4.9, -28.9, -54.2, \
	74.1, 32.8, 85.4, -97.8, \
	32.1, -76.6, -75.9, 99.2, \
	64.5, -17.7, 37.8, -12.1 })

#define VEC_IN_LEN 4
#define VEC_IN ((double []) { \
	22.1, -43.3, 57.3, 54.8 })

#define VEC_OUT_LEN 4


int main(int argc, char **argv)
{
	/* Necessary variables. */
	int i = 0;
	int j = 0;
	int rank = 0;
	int size = 0;

	/* Root's variables. */
	int *p_elms_ns = NULL;
	int *p_elms_ofsts = NULL;
	double *p_mat_in = NULL;
	double *p_vec_in = NULL;
	double *p_vec_out = NULL;

	/* Any process's (local) variables. */
	int lcl_elms_n = 0;
	int lcl_elms_ofst = 0;
	double *p_lcl_mat_elms = NULL;
	double *p_lcl_vec_elms = NULL;
	double lcl_rslt_elm = 0.;

	/* Initialization de MPI. */
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	if (size == 1) { puts("You want parallel stuff or what? Anyway."); }

	/* Initialize matrix and vector. */
	if (rank == 0) {
		p_mat_in = calloc(
			MAT_IN_ROWS_N * MAT_IN_CLMS_N, sizeof(*p_mat_in));
		memcpy(
			p_mat_in,
			MAT_IN,
			MAT_IN_ROWS_N * MAT_IN_CLMS_N * sizeof(*p_mat_in));
		p_vec_in = calloc(VEC_IN_LEN, sizeof(*p_vec_in));
		memcpy(p_vec_in, VEC_IN, VEC_IN_LEN * sizeof(*p_vec_in));
		p_vec_out = calloc(MAT_IN_ROWS_N, sizeof(*p_vec_out));
	}

	/* A process defines how much matrix elements it should work with */
	/* and offset to these elements in a matrix row. */
	lcl_elms_n = MAT_IN_CLMS_N / size;
	lcl_elms_ofst = rank * lcl_elms_n;
	if (size > 1 && rank == size - 1) {
		lcl_elms_n = MAT_IN_CLMS_N - lcl_elms_n * (size - 1);
	}

	/* Gathering numbers of elements processes will work with */
	/* and offsets to these elements in a matrix row. */
	if (rank == 0) {
		p_elms_ns = calloc(size, sizeof(*p_elms_ns));
		p_elms_ofsts = calloc(size, sizeof(*p_elms_ofsts));
	}
	MPI_Gather(
		&lcl_elms_n,     /* IN sendbuf */
		1,               /* IN sendcount */
		MPI_INTEGER,     /* IN sendtype */
		p_elms_ns,       /* OUT recvbuf */
		1,               /* IN recvcount */
		MPI_INTEGER,     /* IN recvtype */
		0,               /* IN root */
		MPI_COMM_WORLD); /* IN comm */
	MPI_Gather(
		&lcl_elms_ofst,  /* IN sendbuf */
		1,               /* IN sendcount */
		MPI_INTEGER,     /* IN sendtype */
		p_elms_ofsts,    /* OUT recvbuf */
		1,               /* IN recvcount */
		MPI_INTEGER,     /* IN recvtype */
		0,               /* IN root */
		MPI_COMM_WORLD); /* IN comm */

	p_lcl_mat_elms = calloc(lcl_elms_n, sizeof(*p_lcl_mat_elms));
	p_lcl_vec_elms = calloc(lcl_elms_n, sizeof(*p_lcl_vec_elms));

	/* Scattering vector elements among processes. */
	MPI_Scatterv(
		p_vec_in,        /* IN sendbuf */
		p_elms_ns,       /* IN sendcounts */
		p_elms_ofsts,    /* IN displs */
		MPI_DOUBLE,      /* IN sendtype */
		p_lcl_vec_elms,  /* OUT recvbuf */
		lcl_elms_n,      /* IN recvcount */
		MPI_DOUBLE,      /* IN recvtype */
		0,               /* IN root */
		MPI_COMM_WORLD); /* IN comm */

	for (i = 0; i < MAT_IN_ROWS_N; i++) {

		/* Scattering matrix elements among processes. */
		MPI_Scatterv(
			p_mat_in + MAT_IN_CLMS_N * i, /* IN sendbuf */
			p_elms_ns,                    /* IN sendcounts */
			p_elms_ofsts,                 /* IN displs */
			MPI_DOUBLE,                   /* IN sendtype */
			p_lcl_mat_elms,               /* OUT recvbuf */
			lcl_elms_n,                   /* IN recvcount */
			MPI_DOUBLE,                   /* IN recvtype */
			0,                            /* IN root */
			MPI_COMM_WORLD);              /* IN comm */

		/* Calculate local products, sum and send em to root. */
		lcl_rslt_elm = 0.;
		for (j = 0; j < lcl_elms_n; j++) {
			lcl_rslt_elm += p_lcl_mat_elms[j] * p_lcl_vec_elms[j];
		}

		MPI_Reduce(
			&lcl_rslt_elm,   /* IN sendbuf */
			p_vec_out + i,   /* OUT recvbuf */
			1,               /* IN count */
			MPI_DOUBLE,      /* IN datatype */
			MPI_SUM,         /* IN operation */
			0,               /* IN root */
			MPI_COMM_WORLD); /* IN comm */

		if (rank == 0) {
			printf("here it is: %0.3lf\n", p_vec_out[i]);
		}
	}

	/* Clean up and finish. */
	free(p_elms_ns);
	free(p_elms_ofsts);
	free(p_mat_in);
	free(p_vec_in);
	free(p_vec_out);

	free(p_lcl_mat_elms);
	free(p_lcl_vec_elms);;

	MPI_Finalize();
	return 0;
}
