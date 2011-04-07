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
//#include <argp.h>
#include "airc.h"


void help();
void usage();
struct airfoil NACA4(int);

struct airfoil {
	char * desc; //short text description
	double * x_c;  //mean camber line
	double * y_c;
	double * x_U;  //upper surface
	double * y_U;
	double * x_L;  // lower surface
	double * y_L;
};


size_t numPoints;
void x_straight (double *);
void x_curved   (struct airfoil *);
void allocate_airfoil(struct airfoil *);
void print_airfoil(struct airfoil * a);

// x_c, y_c, x_U, y_U, x_L, y_L
void airfoil_gen(double *,double *,double *, double *,double *,double *);
//void airfoil_gen(double *);

void meanLineNACA4(double, double, double *, double *);

int main(int argc, char *argv[]) {
	int i;
	int num = atoi(argv[2]);
	numPoints = 30;

	puts("airC airfoil generator");
	printf("%s %s %s\n",argv[0],argv[1],argv[2]);

	struct airfoil a = NACA4(num);

	print_airfoil(&a);
	return EXIT_SUCCESS;
}

struct airfoil NACA4(int num){
	struct airfoil a;

	a.x_c = malloc(numPoints*sizeof(double));
	a.y_c = malloc(numPoints*sizeof(double));
	a.x_L = malloc(numPoints*sizeof(double));
	a.x_U = malloc(numPoints*sizeof(double));
	a.y_L = malloc(numPoints*sizeof(double));
	a.y_U = malloc(numPoints*sizeof(double));

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

void allocate_airfoil(struct airfoil * a){
	double hi;
	a->x_c =  &hi;//malloc(100*sizeof(double));
	(*a).y_c = malloc(numPoints*sizeof(double));
	(*a).x_L = malloc(numPoints*sizeof(double));
	(*a).x_U = malloc(numPoints*sizeof(double));
	(*a).y_L = malloc(numPoints*sizeof(double));
	(*a).y_U = malloc(numPoints*sizeof(double));
}

void print_airfoil(struct airfoil * a){
	printf("i \t x_c \t y_c \t x_U \t y_U \t x_L \t y_L\n");
	printf("------\t-----\t------\t------\t------\t------\t------\n");
	int i;
	for (i=0; i<numPoints;i++){
		printf("%d\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\t%6.3f\n",i, *a->x_c, *a->y_c, *a->x_U, *a->y_U, *a->x_L, *a->y_L);
		(*a).x_c++;
		(*a).y_c++;
		(*a).x_U++;
		(*a).y_U++;
		(*a).x_L++;
		(*a).y_L++;
	}
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

