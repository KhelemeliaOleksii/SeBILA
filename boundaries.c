#include<stdio.h>

// DEFINITION: Function tests values of boundaries for integration procedure
// IN:	int dimention - dimention of task
// OUT: float boundaryX_max[] - a set of max values of boundaries 
	//	float boundaryX_min[] - a set of min values of boundaries 
// return 0 - uncorrect values;
// return 1	- correct values;
int boundary_Test(int dimention, float *boumin_set, float*boumax_set) {
	int i;
	for (i=0; i < dimention; i++) {
		if (boumax_set[i] < boumin_set[i]) {
			fprintf(stdout, "Error: boundary_Test.function->Uncorrect values of boundaries. max < min\n");
			return 0;
		} else if (boumax_set[i] == boumin_set[i]) {
			fprintf(stdout, "Warninig: boundary_Test.function->Boundaries of integaration process are equal. Result is 0\n" );
			return 0;
		} 
	}
	return 1;
}

// DEFINITION: Function create sets of max and min values 
	//of boundaries for integration procedure
// IN:	int dimention - dimention of task
// OUT: float boundaryX_max[] - a set of max values of boundaries 
	//	float boundaryX_min[] - a set of min values of boundaries 
// WITHIN: Function tests values
// return 0 - uncorrect values;
// return 1	- correct values;
int boundary_Create (int dimention, float boundaryX_min[], float boundaryX_max[]) {
	int i;
	for(i = 0; i < dimention; i++) {
		boundaryX_min[i] = -10;
		boundaryX_max[i] = 10;
	}
	if (boundary_Test(dimention, boundaryX_min, boundaryX_max)) {
		return 1;
	} else {
		fprintf(stdout, "Error: boundary_Create.function->Values of boundaries are Uncorrect. See boundary_Test.function\n");
		return 0;
	}

}
