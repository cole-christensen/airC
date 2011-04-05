/*
 ============================================================================
 Name        : airC.c
 Author      : 
 Version     :
 Copyright   : Your copyright notice
 Description : Hello World in C, Ansi-style
 ============================================================================
 */

//TODO: PLplot http://plplot.sourceforge.net/

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void help();
void usage();
struct airfoil NACA4(int);

struct airfoil {
	char * desc; //short text description
	long * x_c;  //mean camber line
	long * y_c;
	long * x_U;  //upper surface
	long * y_U;
	long * x_L;  // lower surface
	long * y_L;
};

size_t numPoints;
void x_straight (double *);
void x_curved   (double *);

// x_c, y_c, x_U, y_U, x_L, y_L
void airfoil_gen(double *,double *,double *, double *,double *,double *);
//void airfoil_gen(double *);

double * dx(double *);
void meanLineNACA4(double, double, double *, double *);

int main(void) {
	int i;
	puts("airC airfoil generator");

	numPoints = 30;
	double m = 0.5;
	double p = 0.6;
	double * x_c = malloc(numPoints*sizeof(double));
	double * y_c = malloc(numPoints*sizeof(double));
	double * x_L = malloc(numPoints*sizeof(double));
	double * x_U = malloc(numPoints*sizeof(double));
	double * y_L = malloc(numPoints*sizeof(double));
	double * y_U = malloc(numPoints*sizeof(double));

	x_curved(x_c);
	meanLineNACA4(m, p, x_c, y_c);
	//airfoil_gen (double * x_c, double * y_c, double * x_U, double * y_U, double * x_L, double * y_L);

	airfoil_gen(x_c, y_c, x_U, y_U, x_L, y_L);
	printf("i \t x_c \t y_c \t x_U \t y_U \t x_L \t y_L\n");
	printf("------\t-----\t------\t------\t------\t------\t------\n");
	for (i=0; i<numPoints;i++){
		printf("%d\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\n",i, *x_c, *y_c, *x_U, *y_U, *x_L, *y_L);
		x_c++;
		y_c++;
		x_U++;
		y_U++;
		x_L++;
		y_L++;
	}

	return EXIT_SUCCESS;
}

double * diff(double * x){
	int i;
	double * dx = malloc(numPoints*sizeof(double));
	dx[0] = x[1] - x[0];
	for(i=1; i<numPoints-1;i++){
		dx[i] = (x[i+1] - x[i-1])/2;
	}
	dx[numPoints-1] = x[numPoints-1] - x[numPoints-2];
	return dx;
}

void meanLineNACA4(double m, double p, double * x, double * y){
	int i;

	for(i=0; i < numPoints; i++){
		if(*x<=p){
			*y = m * *x / pow(p,2) * (2*p - *x);
		} else {
			*y = m * (1 - *x) / pow((1 - p),2) * (1 + *x - 2*p);
		}
		x++;
		y++;
	}
}

/**
 * evenly distributed x points
 */
void x_straight(double * xPoints){
	double dx = 1.0 / (numPoints-1);
	int i;

	for(i = 0; i < numPoints; i++){
		*xPoints = dx * i;
		xPoints++;
	}
}

/**
 * concentrate points towards the leading and trailing edges
 */
void x_curved(double* xPoints){
	double dx = 1.0 / (numPoints-1);
	//double * xPoints;// = malloc(numPoints*sizeof(double));
	//double * xPoints_ = xPoints;

	int i;

	for(i = 0; i < numPoints; i++){
		*xPoints = dx * i;
		*xPoints = 0.5*cos(*xPoints * M_PI + M_PI) + 0.5;
		xPoints++;
	}
	//return xPoints;
}

// x_c, y_c, x_U, y_U, x_L, y_L
/**
 *
 */
void airfoil_gen (double * x_c,double * y_c, double * x_U, double * y_U, double * x_L, double * y_L){
	double y;
	double * dy_c = diff(y_c);
	double * dx   = diff(x_c);
	double theta;
	double t = 0.12;
	int i;
	for(i = 0; i < numPoints; i++){
		y = t/0.2*(0.2969 * sqrt(x_c[i]) - 0.1260 * x_c[i] - 0.3516 * x_c[i] * x_c[i] + 0.2843 * x_c[i] * x_c[i] * x_c[i] - 0.1015 * x_c[i] * x_c[i] * x_c[i] * x_c[i]);
		theta = atan2(dy_c[i],dx[i]);
		x_L[i] = x_c[i] + y * sin(theta);
		y_L[i] = y_c[i] - y * cos(theta);
		x_U[i] = x_c[i] - y * sin(theta);
		y_U[i] = y_c[i] + y * cos(theta);
	}
}

