#include <stdio.h>
#include <stdlib.h>
#include <string.h>


void shiftMatL(double* p_mat, int rows, int clms, int shift)
{
	double* p_matTmp = NULL;
	int shiftNew = shift % clms;

	if (shiftNew == 0) { return; }

	p_matTmp = malloc(rows * clms * sizeof *p_matTmp);
	memcpy(p_matTmp, p_mat, rows * clms * sizeof *p_matTmp);

	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < clms; j++) {
			int jNew = (j - shiftNew + clms) % clms;
			p_mat[i * clms + jNew] = p_matTmp[i * clms + j];
		}
	}

	free(p_matTmp);
	return;
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

int main() {

	double mat[] = {0., 1., 2., 3., 4., 5.}; /* rows 3 x clms 2 */
	int rows = 3;
	int clms = 2;

	printMat((double *) &mat, 3, 2);
	double *p_matT = transMat((double *) &mat, 3, 2); /* rows 2 x clms 3 */
	printMat((double *) p_matT, 2, 3);
	shiftMatL(p_matT, 2, 3, 3);
	printMat((double *) p_matT, 2, 3);

	return 0;
}
