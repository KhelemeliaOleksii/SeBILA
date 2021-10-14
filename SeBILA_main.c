// A Program calculates secondary Bohr term of 
// erenergy losses of the projectile particle 
// in electron gas with anisotropic temperature 
// within quantum field apprximation
// Mathematic calulation one can find in 
// SeBILA projects files: SeBILA_work direct

//Author: Khelemelia O.V.
//Institute of applied physics of NAS of Ukraine, Sumy
//Start: 12 april 2020

// Algorithm include 6-dimentional Simpson's quadrature
// for integral with infinity boundaries: 
// int d^6 p f(p) \limits_{-\infinity}^{\infinity}
// Test function: f(p) = exp(-p^2),
// unswer: Pi^3

//version 1.01
// finish on:	1) integrand --- OK
//				2) simpson integrator --- not OK

// version 1.02 
// finish on: stack overflow

// version 1.03
// finish on: 1D Simson quadrature work
	// interval - half
// version 1.04
// finish on: 1D Simson quadrature work
		// interval - quater
// version 1.05
	// add function files
	// add integrator files
// version 1.08 is bad
// version 1.09
	// precalculated result take in
// version 1.10
// version 1.11 is bad
// version 1.12 2D integration
// version 1.13 3D integration
	// precalculated result is taken in
	// do not work fo N-Dimetional integrand
// version 1.14 3D integration
	// verify boundary nodes in integrator procedure
	// for N-Dimention
// version 1.15 ND integration
	// do ND integration
	// do not do task_creator for N-Dimention (only 1D, 2D and 3D)
// version 1.16 ND integration
	// do task_creator for 1-Dimention and 2-Dimention
	// add task_creator fo 3 Dimention
	// don't work for older dimention (4, 5, 6)
		// need to use recursion algorithm
// version 1.17 ND intefrator works
	// one can choose n-dimention function (to 6D)
	// write n
	// test-time for 6d integrand exp(-\vec{R}^2) is above 6 hours
		// R = (x1,x2,...xn);
// version 1.18 function is commented
	// add parameters of function set in main file
	// add boundaries files
		// boundaries function is transfered to boundaries.c file
	// add task files
		// tasks function is transfered to task.c file
	// Qeuestion: how will global variables be excluded (DIM, NUMPOINTS, ...)

// version 2... MPI add
// version 2.05 
	// add datas (tasks) sender function
	// add a series number of the job 
#include<stdio.h>
#include<conio.h>
#include<math.h>
#include <stdlib.h> //malloc
#include "Function.h"
#include "Communicator.h"
#include "Boundaries.h"
#include "mpi.h"
//#include <stdlib.h> //malloc

#define PARAMETERSNUMBER 1


void main (int argc, char **argv) {
	float result = 0.;
	float boundaryX_max[DIM];
	float boundaryX_min[DIM];
	float parameters[PARAMETERSNUMBER];
	int flagMPI_Init;
	int size, rank;
	int counter = 0;
	int series = 1; //a serial number of the work
 	MPI_Init(&argc, &argv); // MPI initialization
	MPI_Initialized (&flagMPI_Init); // test MPI initialization

	MPI_Comm_size(MPI_COMM_WORLD, &size); 
	//definition of indexes of work processors
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); // index number of processor
	
	if (rank == 0) {
		if(size <= 1) {
			fprintf(stdout, "Warning: Start with a 1 (one) processor. \n");fflush(stdout);
//			MPI_Abort(MPI_COMM_WORLD, MPIErrorCode);
		}
	}

	//testfunc
	//bous: -10..10
	//ans3D: 5,5683279968317078452848179821188 Pi^3/2
	//ans2D: 3,1415926535897932384626433832795 PI
	//testfunc1
	//bous: -10..10
	//ans3D: 76847.35652

	// define boundaries 
	// if values of boundaroes is correct call communicator
	if (boundary_Create(DIM, boundaryX_min, boundaryX_max)) {
		communicator_mpi(rank, size, series, DIM, integrand_function_test, boundaryX_min, boundaryX_max, parameters, &result, &counter);
	}
	fprintf(stdout, "proc %d do %d iterations\n", rank, counter); fflush(stdout);
	if (rank == 0) {
		fprintf(stdout,"result of integration is %.9e\n", result);fflush(stdout);
	}
	MPI_Finalize();
	//_getch();
}

