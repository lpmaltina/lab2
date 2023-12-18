#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "parallel.h"

#define ERROR_STRING_LENGTH 150

void handleMPIError(int mpiError){
    if (mpiError != MPI_SUCCESS){
        printf("MPI message : %d\n", mpiError);
        char* errorMessage = malloc(ERROR_STRING_LENGTH * sizeof(char));
        int errorMessageLength;
        MPI_Error_string(mpiError, errorMessage, &errorMessageLength);
        printf("Message : %s, Length: %d\n", errorMessage, errorMessageLength);
        free(errorMessage);
    }
}
void parallelProduct(double* matrix, double* vector, double* result, int numRows, int numCols, int myRank, int commSize) {
    int error;
    int commSizeRows = sqrt(commSize);
    int commSizeCols = commSize/commSizeRows;
    int localNumRows = numRows/commSizeRows;
    int localNumCols = numCols/commSizeCols;
    double* localMatrix = (double*)malloc(localNumRows * localNumCols * sizeof(double));
    double* localVector = (double*)malloc(localNumCols * sizeof(double));
    MPI_Datatype mpiVec;
    MPI_Type_vector(localNumRows, localNumCols, numCols, MPI_DOUBLE, &mpiVec);
    MPI_Type_commit(&mpiVec);

    if (myRank == MAIN_PROCESS){
        int position;
        int packSize;
        MPI_Pack_size(1, mpiVec, MPI_COMM_WORLD, &packSize);
        for (int i = 0; i < commSizeRows; ++i){
            for (int j = 0; j < commSizeCols; ++j){
                if (i == 0 && j == 0)
                    continue;
                position = 0;
                error = MPI_Pack(&matrix[i * localNumRows * numCols + j * localNumCols], 1, mpiVec, localMatrix, packSize, &position, MPI_COMM_WORLD);
                handleMPIError(error);
                error = MPI_Send(localMatrix, localNumRows * localNumCols, MPI_DOUBLE, i * commSizeCols + j, SUBMATR_TAG, MPI_COMM_WORLD);
                handleMPIError(error);
                error = MPI_Send(&vector[j * localNumCols], localNumCols, MPI_DOUBLE, i * commSizeCols + j, SUBVEC_TAG, MPI_COMM_WORLD);
                handleMPIError(error);
            }
        }
        position = 0;
        error = MPI_Pack(&matrix[0], 1, mpiVec, localMatrix, packSize, &position, MPI_COMM_WORLD);
        handleMPIError(error);
        memcpy(localVector, &vector[0], localNumCols * sizeof(double));
    } else {
        error = MPI_Recv(localMatrix, localNumRows * localNumCols, MPI_DOUBLE, MAIN_PROCESS, SUBMATR_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        handleMPIError(error);
        error = MPI_Recv(localVector, localNumCols, MPI_DOUBLE, MAIN_PROCESS, SUBVEC_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        handleMPIError(error);
    }
    error = MPI_Type_free(&mpiVec);
    handleMPIError(error);

    for (long int j = 0; j < localNumCols; ++j){
        for (long int i = 0; i < numRows; ++i){
            localMatrix[i * localNumCols + j] *= vector[j];
        }
    }
    double* columns = (double*) malloc(numCols * sizeof(double));
    for (long int i = 0; i < numRows; ++i){
        double sum = 0.0;
        for (long int j = 0; j < localNumCols; ++j){
            sum += localMatrix[i * localNumCols + j];
        }
        columns[i] = sum;
    }

    MPI_Reduce(columns, result, numRows, MPI_DOUBLE, MPI_SUM, MAIN_PROCESS, MPI_COMM_WORLD);
    if (myRank == MAIN_PROCESS){
        for (long int i = 0; i < numRows; ++i){
            result[i] = 0.0;
        }
    }
    if (myRank != MAIN_PROCESS){
        error = MPI_Send(localVector, localNumRows, MPI_DOUBLE, MAIN_PROCESS, SUBVEC_TAG, MPI_COMM_WORLD);
        handleMPIError(error);
    } else {
        MPI_Status status;
        double* recvBuffer = (double*) malloc(localNumRows * sizeof(double));
        for (int i = 0; i < commSize; ++i){
            int src;
            if (i != 0){
                error = MPI_Recv(recvBuffer, localNumRows, MPI_DOUBLE, MPI_ANY_SOURCE, SUBVEC_TAG, MPI_COMM_WORLD, &status);
                handleMPIError(error);
                src = status.MPI_SOURCE;
            } else {
                memcpy(recvBuffer, localVector, localNumRows * sizeof(double));
                src = 0;
            }
            for (long int j = 0; j < localNumRows; ++j){
                result[(src/commSizeCols) * localNumRows + j] += recvBuffer[j];
            }
        }
    }
    free(localMatrix);
    free(localVector);
    free(columns);
}
