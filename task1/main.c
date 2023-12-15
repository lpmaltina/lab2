#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "utils.h"
#include "parallel.h"


#define FILE_NAME_LEN 128


int main(int argc, char** argv)
{
	int lclTmp;

	double* p_mat = NULL;
	double* p_vec = NULL;
	double* p_rslt = NULL;

	double timeStrt = 0.;
	double timeFnsh = 0.;
	double timeElapsed = 0.;

	int dim = strtol(argv[1], NULL, 10);
	int matRowsN = dim;
	int matClmsN = dim;

	int rank;
	int size;

    FILE* fp;
    char* mode = argv[2];

	MPI_Init(&argc, &argv); 
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	p_vec = malloc(matClmsN * sizeof(*p_vec));

	if (rank == 0) {
		char fnameMat[FILE_NAME_LEN] = {0};
		char fnameVec[FILE_NAME_LEN] = {0};

		p_mat = malloc(matRowsN * matClmsN * sizeof(*p_mat));
		p_rslt = malloc(matClmsN * sizeof(*p_rslt));

		sprintf(
			fnameMat,
			"input/mat-%d-%d.txt",
			matRowsN,
			matClmsN);

		if (readArray(fnameMat, p_mat, matRowsN * matClmsN)) {
			puts("KU! INPUT MATRIX FILENAME ERROR!");
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		sprintf(fnameVec, "input/vec-%d.txt", matClmsN);
		if (readArray(fnameVec, p_vec, matClmsN)) {
			puts("KU! INPUT VECTOR FILENAME ERROR!");
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0) { timeStrt = MPI_Wtime(); }
	parallelProduct(p_mat, p_vec, p_rslt, matRowsN, matClmsN, rank, size);
	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0) {
		timeFnsh = MPI_Wtime();
		timeElapsed = timeFnsh - timeStrt;
	}

	if (rank == 0) {
		char fnameRslt[FILE_NAME_LEN] = {0};
        char fnameTimings[FILE_NAME_LEN] = {0};

		sprintf(
			fnameRslt,
			"output/parallel-product-%d-%d.txt",
			matRowsN,
			matClmsN);

		if (writeArray(fnameRslt, p_rslt, matRowsN)) {
			puts("KU! OUTPUT FILENAME ERROR!");
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
        
        sprintf(fnameTimings, "timings/timings-%s.txt", mode);
        fp = fopen(fnameTimings, "a");
        fprintf(
            fp,
            "%d, %d, %d, %lf\n",
            matRowsN,
            matClmsN,
            size,
            timeElapsed
        );
        fclose(fp);

		free(p_mat);
		free(p_rslt);
	}

	free(p_vec);

	MPI_Finalize();
	return 0;
}
