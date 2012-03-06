/*
 ============================================================================
 Name        : airC.c
 Author      : Cole Christensen 
 Version     : 0.1
 Copyright   : 2012
 Description : 
 ============================================================================
 */

//TODO: PLplot http://plplot.sourceforge.net/

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "airc.h"
#include "airfoil.h"
#include "NACA4.h"


void help();
void usage();

void usage(){
	printf("Airfoil generator example usage:\n\n");
	printf("   airC NACA 2412\n");
}

void x_curved   (struct airfoil *);

// x_c, y_c, x_U, y_U, x_L, y_L
void airfoil_gen(double *,double *,double *, double *,double *,double *);
//void airfoil_gen(double *);

int main(int argc, char *argv[]) {
	for(int i=0; i<argc;i++){
		printf("%s ",argv[i]);
	} printf("\n");

	if(argc > 1){
		if (strcmp(argv[1],"NACA") == 0){
			//NACA4();
			NACA4(atoi(argv[2]),11);
		} else {
			usage();
			return EXIT_FAILURE;
		}
	} else {
		usage();
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

/*
struct airfoil NACA4a(int num){
	struct airfoil a;

	double m,p,t;

	t = num%100;
	num = num - t;
	p = num%1000;
	num=num-p;
	m = num;

	//allocate_airfoil(a);

	x_curved(&a);

	meanLineNACA4(m, p, a.x_c, a.y_c);

	airfoil_gen(a.x_c, a.y_c, a.x_U, a.y_U, a.x_L, a.y_L);


	return a;
}
*/

/**
 * concentrate points towards the leading and trailing edges
 */
void x_curved(struct airfoil * a){
	double dx = 1.0 / (numPoints-1);
	double * xPoints = a->x_c;
	//double * xPoints_ = xPoints;

	int i;

	for(i = 0; i < numPoints; i++){
		*xPoints = dx * i;
		*xPoints = 0.5*cos(*xPoints * M_PI + M_PI) + 0.5;
		xPoints++;
	}
	//return xPoints;
}

