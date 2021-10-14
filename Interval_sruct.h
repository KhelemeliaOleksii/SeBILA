
#ifndef INTERVAL_STRUCT_H
#define INTERVAL_STRUCT_H
#include "Function.h"

#if (DIM == 1)
typedef struct {
	int key;	// 0 - do calculation of function value
				// 1 - calculation have been done
	float bou[DIM][2];
	float fvalue[2];
} interval;
#elif (DIM == 2)
typedef struct {
	int key;	// 0 - do calculation of function value
				// 1 - calculation have been done
	float bou[DIM][2];
	float fvalue[4]; // 2^DIM
} interval;
#elif (DIM == 3)
typedef struct {
	int key;	// 0 - do calculation of function value
				// 1 - calculation have been done
	float bou[DIM][2];
	float fvalue[8]; // 2^DIM
} interval;
#elif (DIM == 4)
typedef struct {
	int key;	// 0 - do calculation of function value
				// 1 - calculation have been done
	float bou[DIM][2];
	float fvalue[16]; // 2^DIM
} interval;
#elif (DIM == 5)
typedef struct {
	int key;	// 0 - do calculation of function value
				// 1 - calculation have been done
	float bou[DIM][2];
	float fvalue[32]; // 2^DIM
} interval;
#else (DIM == 6)
typedef struct {
	int key;	// 0 - do calculation of function value
				// 1 - calculation have been done
	float bou[DIM][2];
	float fvalue[64]; // 2^DIM
} interval;
#endif

typedef struct {
	float point[DIM];
} node;

#endif