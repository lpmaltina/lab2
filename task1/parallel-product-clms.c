#include "parallel.h"


/* TODO Consider using const in function definition like this: */
/* void parallelProduct( */
/* 	double const *const p_mat, */
/* 	double const *const p_vec, */
/* 	double *const p_rslt, */
/* 	int const mat_rows_n, */
/* 	int const mat_clms_n, */
/* 	int const rank, */
/* 	int const size) */

/* Multiply matrix p_mat on vector p_vec, write result to p_rslt */
void parallelProduct(
	double *p_mat,  /* root */
	double *p_vec,  /* root */
	double *p_rslt, /* root */
	int mat_rows_n, /* local, but always the same */
	int mat_clms_n, /* local, but always the same */
	int rank,       /* local */
	int size)       /* local, but always the same */
{
	/* Temporary variables. */
	int i = 0;
	int j = 0;

	/* Root's variables. */
	int *p_elms_ns = NULL;
	int *p_elms_ofsts = NULL;

	/* Local variables. */
	int lcl_elms_n = 0;
	int lcl_elms_ofst = 0;
	double *p_lcl_mat_elms = NULL;
	double *p_lcl_vec_elms = NULL;
	double lcl_rslt_elm = 0.;

	if (size == 1) { puts("You want parallel stuff or what? Anyway."); }

	/* A process defines how much matrix elements it should work with */
	/* and offset to these elements in a matrix row. */
	lcl_elms_n = mat_clms_n / size;
	lcl_elms_ofst = rank * lcl_elms_n;
	if (size > 1 && rank == size - 1) {
		lcl_elms_n = mat_clms_n - lcl_elms_n * (size - 1);
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
		p_vec,           /* IN sendbuf */
		p_elms_ns,       /* IN sendcounts */
		p_elms_ofsts,    /* IN displs */
		MPI_DOUBLE,      /* IN sendtype */
		p_lcl_vec_elms,  /* OUT recvbuf */
		lcl_elms_n,      /* IN recvcount */
		MPI_DOUBLE,      /* IN recvtype */
		0,               /* IN root */
		MPI_COMM_WORLD); /* IN comm */

	for (i = 0; i < mat_rows_n; i++) {

		/* Scattering matrix elements among processes. */
		MPI_Scatterv(
			p_mat + mat_clms_n * i, /* IN sendbuf */
			p_elms_ns,              /* IN sendcounts */
			p_elms_ofsts,           /* IN displs */
			MPI_DOUBLE,             /* IN sendtype */
			p_lcl_mat_elms,         /* OUT recvbuf */
			lcl_elms_n,             /* IN recvcount */
			MPI_DOUBLE,             /* IN recvtype */
			0,                      /* IN root */
			MPI_COMM_WORLD);        /* IN comm */

		/* Calculate local products, sum and send em to the root. */
		lcl_rslt_elm = 0.;
		for (j = 0; j < lcl_elms_n; j++) {
			lcl_rslt_elm += p_lcl_mat_elms[j] * p_lcl_vec_elms[j];
		}

		MPI_Reduce(
			&lcl_rslt_elm,   /* IN sendbuf */
			p_rslt + i,      /* OUT recvbuf */
			1,               /* IN count */
			MPI_DOUBLE,      /* IN datatype */
			MPI_SUM,         /* IN operation */
			0,               /* IN root */
			MPI_COMM_WORLD); /* IN comm */
	}

	/* Clean up and return */
	free(p_elms_ns);
	free(p_elms_ofsts);
	free(p_lcl_mat_elms);
	free(p_lcl_vec_elms);

	return;
}
