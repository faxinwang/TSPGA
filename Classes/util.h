#pragma once

#include <cstdlib>
#define eps (1e-6)

inline int randInt(int low,int high) {
	if (low >= high) return low;
	return rand() % (high - low + 1) + low;
	
}

inline float randFloat() {
	return rand() / (RAND_MAX + 1.0);
}

inline int dcmp(double x) {
	if (fabs(x) <= eps) return 0;
	if (x > 0) return 1;
	return -1;
}