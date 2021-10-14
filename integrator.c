#include<stdio.h>
#include<math.h>
#include<conio.h>
#include "Interval_sruct.h"
#include "Integrator.h"
#include <stdlib.h> // malloc

// DEFINITION: Function calculate Simpson coefficient
// IN: int parameter			--- iterator
// OUT: as income parameter
	//int* parameter			---	coefficient1
	//int* parameter			---	coefficient2
void Simpson_coefficient (int iterator, int* coefficient1, int* coefficient2) {
		if ((iterator == 0)||(iterator == 4)) {
			*coefficient1 = 1;
			*coefficient2 = 1;
		}
		if ((iterator == 2)) {
			*coefficient1=2;
			*coefficient2=4;
		}
		if ((iterator == 1)||(iterator == 3)) {
			*coefficient1=4;
			*coefficient2=0;
		}
}

// DEFINITION: Function defines where are boundary nodes
// USE: reccurtion
// IN: int dimention			
// OUT: int bounodes[] - a set of boundary nodes
void fbounodes(int dimention, int bounodes[]) {
	int boupoints;
	int i;
	int numpoints, numpoints_minus1D;
	boupoints = (int)pow((float)2,dimention);
	numpoints = (int)pow((float)NUMPOINTS1D,dimention);
	numpoints_minus1D = (int)pow((float)NUMPOINTS1D,dimention-1);
	// the first half of nodes in previous dimention
	if (dimention == 1) {
		bounodes[0] = 0;
		bounodes[1] = 4;
	} else {
		fbounodes(dimention-1, bounodes); 
		for (i = boupoints/2; i < boupoints; i++) {
			bounodes[i] = bounodes[i-boupoints/2] + numpoints-numpoints_minus1D;  
		}
	}
}

// DEFINITION: Function test Is a node in boundary nodes set
// IN:	int previouscounter - start from previous position
	//	int counter 
	//	int bounodes[] - set of a boundary nodes
// OUT: int *bounodescounter - position of the boundary node
int boufunctest(int previouscounter, int counter, int bounodes[], int *bounodescounter) {
	int i;
	int numpoints;
	numpoints = (int)pow(2.,DIM);

	for(i=previouscounter; i < numpoints; i++) {
		if (counter == bounodes[i]) {
			*bounodescounter = i;
			return 1;
		}
	}
	return 0;
}

// DEFINITION: Function calculates summ and subsumm
	// of the task interval within Simpson's quadrature
// USE: reccurtion
// IN:	int (*function)(int, float*, float*) - an integrand
	//	interval task - a task for integration
// OUT: float coord_set[][NUMPOINTS1D] - set of boundaries
	//	float funcvalue_set[] - set of function values in nodes position
	//	float *summ_subBubble - half step summation; (5 nodes per dimention)
	//	float *summ_Bubble - full step summation; (3 nodes per dimention)
void insideintegrator (int (*function)(int, float*, float*), int *node_iter, int dimkey, float stepOcoef_product, float stepPcoef_product, node *nodeone,  
	interval task, float coord_set[][NUMPOINTS1D], float funcvalue_set[], float *summ_subBubble, float *summ_bubble) {
	int i;
	int o_coeffSim, p_coeffSim;
	int j_prev = 0;
	int bounodenum; // number of boundary nodes
	int boucounter = 0;
	float step;
	float term_subBubble  = 0;
	float term_bubble  = 0;
	float tmp_stepOcoef_product;
	float tmp_stepPcoef_product;
	int *bounodes; // list of boundary nodes
	
	if  ((DIM-dimkey) == 1){
		step = (float)((task.bou[dimkey][1] - task.bou[dimkey][0])/4.0);
		bounodenum = (int)pow((float)2,DIM);
		bounodes = (int*)malloc(bounodenum*sizeof(node));
		for (i = 0; i < 5; i++) {
			coord_set[dimkey][i] = task.bou[dimkey][0]+(float)i*step;
			Simpson_coefficient(i, &o_coeffSim, &p_coeffSim);
			tmp_stepOcoef_product = (float)(stepOcoef_product*step*o_coeffSim/3.);
			tmp_stepPcoef_product = (float)(stepPcoef_product*step*p_coeffSim*2/3.);
			// if a function value precalculated
			if (task.key == 1) {
				// define boundary nodes
				fbounodes(DIM, bounodes);
				// test is nodes boundary one
				if (boufunctest(j_prev, *node_iter, bounodes, &boucounter)) {
					funcvalue_set[*node_iter] = task.fvalue[boucounter];
					j_prev = *node_iter;
				} else {
					nodeone->point[dimkey] = coord_set[dimkey][i];
					function(DIM, nodeone->point, &funcvalue_set[*node_iter] );
				}
			// if a function value didn`t precalculated
			} else {
				nodeone->point[dimkey] = coord_set[dimkey][i];
				function(DIM, nodeone->point, &funcvalue_set[*node_iter] );
			} 

				term_subBubble = funcvalue_set[*node_iter]*(float)(tmp_stepOcoef_product );
				term_bubble = funcvalue_set[*node_iter]*(float)(tmp_stepPcoef_product);
				*summ_subBubble += term_subBubble;
				*summ_bubble += term_bubble;
			
				*node_iter = *node_iter+1;
		}
		free(bounodes);
	} else {
		step = (float)((task.bou[dimkey][1] - task.bou[dimkey][0])/4.0);
		for (i = 0; i < 5; i++) {
			coord_set[dimkey][i] = task.bou[dimkey][0]+(float)i*step;
			Simpson_coefficient(i, &o_coeffSim, &p_coeffSim);
			tmp_stepOcoef_product = (float)(stepOcoef_product*step*o_coeffSim/3.);
			tmp_stepPcoef_product = (float)(stepPcoef_product*step*p_coeffSim*2/3.);
			nodeone->point[dimkey] = coord_set[dimkey][i];
			insideintegrator (function, node_iter, dimkey+1, tmp_stepOcoef_product, tmp_stepPcoef_product, nodeone,  task, coord_set, funcvalue_set, summ_subBubble, summ_bubble);
		}
	}
}

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
int Simpson_ND_integrator(int (*function)(int, float*, float*), interval task, float coord_set[][NUMPOINTS1D], float funcvalue_set[], float *result) {
	float summ_subBubble, summ_bubble;
	float precision;
	float tmpsum = 0;
	float numpoints;
	float stepOcoef_product = 1;
	float stepPcoef_product = 1;
	int dimkey = 0;
	int node_iter = 0;
	node nodeone;
	*result = 0;
	precision = (float)0.00001;
	numpoints = (int)pow((float)NUMPOINTS1D, DIM);
	summ_subBubble = 0;
	summ_bubble = 0;

	// to calulate summ_subBubble, summ_bubble and coord_set, funcvalue_set
	insideintegrator(function, &node_iter, dimkey, stepOcoef_product, stepPcoef_product, &nodeone, task, coord_set, funcvalue_set, &summ_subBubble, &summ_bubble);

	if (fabs((summ_bubble-summ_subBubble)) > precision ) {
		return 1;
	} else {
		*result += summ_subBubble;
	}
	return 0;
}
