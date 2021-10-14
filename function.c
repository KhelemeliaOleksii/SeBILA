#include<math.h>
#include "Interval_sruct.h"

int integrand_function_test (int argue_count, float* agrue_values, float *result){
// boundaries -infinity..infinity
// answer pow(sqrt(PI), argue_count)
	if (argue_count == 1) {
		*result = (float)exp(-(
			agrue_values[0]*agrue_values[0]));
			//+agrue_values[1]*agrue_values[1]+
			//agrue_values[2]*agrue_values[2]+agrue_values[3]*agrue_values[3]+
			//agrue_values[4]*agrue_values[4]+agrue_values[5]*agrue_values[5]));
		return 0;
	} else if (argue_count == 2) {
		*result = (float)exp(-(
			agrue_values[0]*agrue_values[0]
			+agrue_values[1]*agrue_values[1]));
			//agrue_values[2]*agrue_values[2]+agrue_values[3]*agrue_values[3]+
			//agrue_values[4]*agrue_values[4]+agrue_values[5]*agrue_values[5]));
		return 0;
	} else if (argue_count == 3) {
		*result = (float)exp(-(
			agrue_values[0]*agrue_values[0]
			+agrue_values[1]*agrue_values[1]+
			agrue_values[2]*agrue_values[2]));
			//+agrue_values[3]*agrue_values[3]+
			//agrue_values[4]*agrue_values[4]+agrue_values[5]*agrue_values[5]));
		return 0;
	} else if (argue_count == 4) {
		*result = (float)exp(-(
			agrue_values[0]*agrue_values[0]+
			agrue_values[1]*agrue_values[1]+
			agrue_values[2]*agrue_values[2]+
			agrue_values[3]*agrue_values[3]));
			//+agrue_values[3]*agrue_values[3]+
			//agrue_values[4]*agrue_values[4]+agrue_values[5]*agrue_values[5]));
		return 0;
	} else if (argue_count == 5) {
		*result = (float)exp(-(
			agrue_values[0]*agrue_values[0]+
			agrue_values[1]*agrue_values[1]+
			agrue_values[2]*agrue_values[2]+
			agrue_values[3]*agrue_values[3]+
			agrue_values[4]*agrue_values[4]));
			//+agrue_values[3]*agrue_values[3]+
			//agrue_values[4]*agrue_values[4]+agrue_values[5]*agrue_values[5]));
		return 0;
	} else if (argue_count == 6) {
		*result = (float)exp(-(
			agrue_values[0]*agrue_values[0]+
			agrue_values[1]*agrue_values[1]+
			agrue_values[2]*agrue_values[2]+
			agrue_values[3]*agrue_values[3]+
			agrue_values[4]*agrue_values[4]+
			agrue_values[5]*agrue_values[5]));
			//+agrue_values[3]*agrue_values[3]+
			//agrue_values[4]*agrue_values[4]+agrue_values[5]*agrue_values[5]));
		return 0;
	} else {
		return 1;
	}
}
int integrand_function_test1 (int argue_count, float* agrue_values, float *result){
// boundaries -10..10
// answer3D 76847.35652
	float t;
	if (argue_count == 1) {
		t = pow((agrue_values[0]*agrue_values[0]),(float)0.5);
		*result = t;
		return 0;
	} else if (argue_count == 2) {
		t = pow((agrue_values[0]*agrue_values[0]+
			agrue_values[1]*agrue_values[1]),(float)0.5);
		*result = t;
		return 0;
	} else if (argue_count == 3) {
		t = pow((agrue_values[0]*agrue_values[0]+
			agrue_values[1]*agrue_values[1]+
			agrue_values[2]*agrue_values[2]),(float)0.5);
		*result = t;
		//*result = t/pow((1-t*t),(float)0.5);
		return 0;
	} else {
		return 1;
	}
}
int integrand_function_test2 (int argue_count, float* agrue_values, float *result){
		if (argue_count <= DIM  ) {
			if (agrue_values[0] < 0){
				*result = 0;
			}
			*result = 1/((float)(1+agrue_values[0]*agrue_values[0]));
			return 0;
		}

}
