#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "utils.h"

#include "parallel.h"


#define FILE_NAME_LEN 128


int main(int argc, char** argv)
{
	double* p_matA = NULL;
	double* p_matB = NULL;
	double* p_rslt = NULL;

	double timeStrt = 0.;
	double timeFnsh = 0.;
	double timeElapsed = 0.;

	int dim = strtol(argv[1], NULL, 10);
	int matRowsN = dim;
	int matClmsN = dim;

	int rank;
	int size;

	MPI_Init(&argc, &argv); 
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if (rank == 0) {
		char fnameMatA[FILE_NAME_LEN] = {0};
		char fnameMatB[FILE_NAME_LEN] = {0};

		p_matA = malloc(matRowsN * matClmsN * sizeof(*p_matA));
		p_matB = malloc(matRowsN * matClmsN * sizeof(*p_matB));
		p_rslt = malloc(matRowsN * matClmsN * sizeof(*p_rslt));

		sprintf(fnameMatA, "input/matA-%d-%d.txt", matRowsN, matClmsN);
		if (readArray(fnameMatA, p_matA, matRowsN * matClmsN)) {
			puts("KU! INPUT MATRIX A FILENAME ERROR!");
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		sprintf(fnameMatB, "input/matB-%d-%d.txt", matRowsN, matClmsN);
		if (readArray(fnameMatB, p_matB, matRowsN * matClmsN)) {
			puts("KU! INPUT MATRIX B FILENAME ERROR!");
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		if (matRowsN == matClmsN) { printMatSq(p_matA, matRowsN); }
		if (matRowsN == matClmsN) { printMatSq(p_matB, matRowsN); }
	}

	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0) { timeStrt = MPI_Wtime(); }
	parallelProduct(p_matA, p_matB, p_rslt, matRowsN);
	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0) {
		timeFnsh = MPI_Wtime();
		timeElapsed = timeFnsh - timeStrt;
	}

	if (rank == 0) {
		char fnameRslt[FILE_NAME_LEN] = {0};

		if (matRowsN == matClmsN) { printMatSq(p_rslt, matRowsN); }

		sprintf(
			fnameRslt,
			"output/parallel-cannon-%d-%d.txt",
			matRowsN,
			matClmsN);

		if (writeArray(fnameRslt, p_rslt, matRowsN)) {
			puts("KU! OUTPUT FILENAME ERROR!");
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		printf(
			"Time (rows_n=%d, clms_n=%d, threads=%d): %lf\n",
			matRowsN, matClmsN, size, timeElapsed);

		free(p_matA);
		free(p_matB);
		free(p_rslt);
	}

	MPI_Finalize();
	return 0;
}
