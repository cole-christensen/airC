/*
 * airc.h
 *
 *  Created on: Apr 5, 2011
 *      Author: chaos
 */

#ifndef AIRC_H_
#define AIRC_H_

void help();
void usage();
struct airfoil * NACA4(int);

struct airfoil {
	char * desc; //short text description
	double * x_c;  //mean camber line
	double * y_c;
	double * x_U;  //upper surface
	double * y_U;
	double * x_L;  // lower surface
	double * y_L;
};


#endif /* AIRC_H_ */
