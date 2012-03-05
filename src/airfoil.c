#include <math.h>
#include <unistd.h>
#include "airfoil.h"

void x_straight(struct airfoil * a){
	double * xPoints = (*a).x_c;
	int numPoints = (*a).numPoints;
	double dx = 1.0 / (numPoints-1);
	int i;

	for(i = 0; i < numPoints; i++){
		*xPoints = dx * i;
		xPoints++;
	}
}
