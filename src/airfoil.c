#include <math.h>
#include <unistd.h>
#include "airfoil.h"

void x_straight(double * xPoints){
	double dx = 1.0 / (numPoints-1);
	int i;

	for(i = 0; i < numPoints; i++){
		*xPoints = dx * i;
		xPoints++;
	}
}
