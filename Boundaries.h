#ifndef BOUNDARIES_H
#define BOUNDARIES_H
// DEFINITION: Function create sets of max and min values 
	//of boundaries for integration procedure
// IN:	int dimention - dimention of task
// OUT: float boundaryX_max[] - a set of max values of boundaries 
	//	float boundaryX_min[] - a set of min values of boundaries 
// WITHIN: Function tests values
// return 0 - uncorrect values;
// return 1	- correct values;
int boundary_Create (int dimention, float *boumin_set, float*boumax_set);
// DEFINITION: Function tests values of boundaries for integration procedure
// IN:	int dimention - dimention of task
// OUT: float boundaryX_max[] - a set of max values of boundaries 
	//	float boundaryX_min[] - a set of min values of boundaries 
// return 0 - uncorrect values;
// return 1	- correct values;
int boundary_Test(int dimention, float *boumin_set, float*boumax_set);


#endif