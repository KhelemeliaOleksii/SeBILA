#ifndef TASK_H
#define TASK_H

#include "Interval_sruct.h"
#include "Integrator.h"

// DEFINITION: Function transforms a set of boundaries
	// to the task format
// IN:	float pointmin_set[] - a set low values of boundaries
	//	float pointmax_set[] - a set up values of boundaries
// OUT: int *count_task - counter of tasks
	//	interval list_task[] - a set  of task
void firsttask_Creator(int dimention, float pointmin_set[], float pointmax_set[], int *count_task, interval list_task[]);

int firsttask_CreatorND_mpi(int rank, int size, int *dimention_key, int dim_iter[], float boundaryX_min[], float boundaryX_max[], int *count_task, interval list_task[]);

// DEFINITION: Function transforms sets of boundaries and function values
	// to the task format
// USE: reccurtion
// IN:	float coord_set[][NUMPOINTS1D]	- a set of boundaries
	//	float funcvalue_set[] - a set of function values	
// OUT: int *count_task - counter of tasks
	//	interval list_task[] - a set  of task
// WITHIN: Function uses one_TaskcreatorND() to write function values 
void task_CreatorND(float coord_set[][NUMPOINTS1D], float funcvalue_set[], int *dimention_key, int *step, 
			int dim_iter[], int *count_task, interval list_task[]);

// DEFINITION: Function add inparameter value to outparameter
// IN:	float subsumm - value of partial summ
// OUT: float *sum - value of summ 
void summator(float subsumm, float *sum);

float totalsummator_mpi(int rank, int size, float subSumm);

void task_sender(int rank, int size, int series, int *count_task, interval list_task[]);

int task_asker(int rank, int size, int series, int *count_task, interval list_task[], int *currentsize, int IncreaseProcId[], int INcountcompare[]);
//int task_asker(int rank, int size, int series, int *count_task, interval list_task[]);
#endif