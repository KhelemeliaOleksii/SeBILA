
#ifndef COMMUNICATOR_H
#define COMMUNICATOR_H

int communicator_mpi(int rank, int size, int series, int dimention, int (*function)(int, float*, float*), float boundaryX_min[], float boundaryX_max[], float parameters[], float *result, int *counter); 
#endif