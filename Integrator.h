// Function.h содержит определение прототипа подинтегральной функции
#ifndef INTEGRATOR_H
#define INTEGRATOR_H
#include "Interval_sruct.h"

#define NUMPOINTS1D 5

// DEFINITION: Function calculates summ of the task interval within Simpson's quadrature
// USE: void insideintegrator();
// IN:	int (*function)(int, float*, float*) - an integrand
	//	interval task - a task for integration
// OUT: 
// if precition is not sufficient
// return 1
	// float coord_set[][NUMPOINTS1D] - set of boundaries
	//	float funcvalue_set[] - set of function values in nodes position
//else
// return 0
	//	float *result - result of summation
int Simpson_ND_integrator(int (*function)(int, float*, float*), interval task, float coord_set[][NUMPOINTS1D], float funcvalue_set[], float *result);

#endif