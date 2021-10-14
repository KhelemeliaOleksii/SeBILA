#include<stdio.h>
#include<conio.h>
#include<math.h>
#include <stdlib.h> //malloc
#include "Integrator.h"
#include "Interval_sruct.h"
#include "Task.h"
#include <mpi.h>
#include <string.h>	//memcpy
// procedure creates list of tasks* (Simpson points)
// procedure call Simpson quadrature function
// result from Simpson quadrature function summ in summation function

int communicator_mpi(int myid, int size, int series, int dimention, int (*function)(int, float*, float*), float boundaryX_min[], float boundaryX_max[], float parameters[], float *result, int *counter) {
	int KEY_exit;  // key to exit from conversation

	// tags for communication
	// tag10 - for count of recieve (send) additional tasks
	// tag11 - for recieve (send) additional tasks
	// tag12 - ask for additional tasks
	// tag13 - no task
	// tag14 - stop the work
	MPI_Status *status13r;
	MPI_Status status12r;
	MPI_Status status12s;
	MPI_Status status11r;
	MPI_Status status11s;
	MPI_Status status10r;
	MPI_Status status10s;
	MPI_Request *req13r;
	MPI_Request req12s = MPI_REQUEST_NULL;
	MPI_Request req12r;
	MPI_Request req11r = MPI_REQUEST_NULL;
	MPI_Request req11s = MPI_REQUEST_NULL;
	MPI_Request req10r = MPI_REQUEST_NULL;
	MPI_Request req10s = MPI_REQUEST_NULL;
	int *flag13r;
	int flag12s = 1;
	int flag12r = 1; // I'm ready to receive a request for additional tasks
	int flag11r = 1;
	int flag11s = 1;
	int flag10r = 0;
	int flag10s;

	int countsend;		// count tasks to send
	int k, i, j;		// iterators
	int *flagexitIn;	//flag of exit increasing list
	int *IncreaseProcId;//list of coprocessors Increasing
	int *LeaveProcId;	//list of coprocessors which leave communication
	int LeaveN = 0;		
	int LeaveJ = 0;
	int bye = 1;
	int InN;		
	int InNTemp;
	int Inj;
	int Inprev;
	int flagExitProc = 0;// my flag of exit
	int signexitIncom;	
	int *signallIncom;
	int *signallAsk;
	int *INcountcompare;	// количество задач в возрастающем потоке
	int seriestmp12 = 0;

	int dim_key = 0;
	int count_task=0;
	float result_tmp, result1 = 0;
	float coord_set[DIM][NUMPOINTS1D];
	int numpoints;
	int step;
	int *dim_iter;
	interval *list_task; 
	float *funcvalue_set;
	int currentsize = size;


	numpoints = (int)pow((float)NUMPOINTS1D, dimention);

	funcvalue_set = (float*)malloc(numpoints*sizeof(float));
	dim_iter = (int*)malloc(DIM*sizeof(int));
	list_task = (interval*)malloc(10000*sizeof(interval));

	// list of coprocessors
	IncreaseProcId = (int*)malloc(size*sizeof(int));
	// list of tasks in processors that we know
	// at the start NULL tasks. We know nothing about.
	INcountcompare = (int*)malloc(size*sizeof(int));


	// create a list of processors for a conversation
	// myid processor is root one
	// example: size 6, myid 4
	// IncreaseProcId[0] = 4, IncreaseProcId[1] = 5..0..1..2..3  
	for(i = myid; i < size; i++) { IncreaseProcId[i-myid] = i;}
	for(i = 0; i < myid; i++) { IncreaseProcId[size-myid+i] = i;}

	// create list of tasks in processors 
	// we assume at the start there are NULL tasks.  
	for (i = 0; i < size; i++) {INcountcompare[i] = 0;}

	
	if (!firsttask_CreatorND_mpi(myid, size, &dim_key, dim_iter, boundaryX_min, boundaryX_max, &count_task, list_task)){
		fprintf(stdout, "ERROR: procedure \"firsttask_CreatorND_mpi\" didn't created tasks\n"); fflush(stdout);
		return 0;
	}
	fprintf(stdout, "I`m %d proc has %d task at the start", myid, count_task);
	KEY_exit = 0;		
	//KEY_exit = 1;		
	while (!KEY_exit){
		while (count_task) {
			*counter = *counter+1;
			if (Simpson_ND_integrator (function, list_task[count_task-1], coord_set,  funcvalue_set, &result_tmp))
			{
				step = 0;
				dim_key = 0;
				count_task--;
				task_CreatorND(coord_set, funcvalue_set, &dim_key, &step, dim_iter, &count_task, list_task);
			} else {
				count_task--;
				summator(result_tmp, result);
			}
			//KEY_exit = 1;
/*program tested ot this point*/			

/*next time i have test send-receive procedure*/
			 //procedure to recieve additional tasks
			if (count_task >= 5) {
				task_sender(myid, size, series, &count_task, list_task);
			}
		}
		// procedure to ask and recive addional tasks
		KEY_exit = task_asker(myid, size, series, &count_task, list_task, &currentsize, IncreaseProcId, INcountcompare);
		
		//fprintf(stdout, "I`m %d, current size is %d\n", myid, currentsize);fflush(stdout);

	}
	*result = totalsummator_mpi(myid, size, *result);

	free(dim_iter);
	free(list_task); 
	free(funcvalue_set);
	return 1;	
}

