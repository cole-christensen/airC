#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "airfoil.h"

void x_curved   (struct airfoil *);
void meanLineNACA4(double, double, struct airfoil *);
void airfoil_gen(struct airfoil *, double);


double * diff(double * x, int numPoints){
	int i;
	double * dx = malloc(numPoints*sizeof(double));
	dx[0] = x[1] - x[0];
	for(i=1; i<numPoints-1;i++){
		dx[i] = (x[i+1] - x[i-1])/2;
	}
	dx[numPoints-1] = x[numPoints-1] - x[numPoints-2];
	return dx;
}

void allocate_airfoil(struct airfoil * a, int numPoints){
	(*a).x_c = malloc(numPoints*sizeof(double));
	(*a).y_c = malloc(numPoints*sizeof(double));
	(*a).x_L = malloc(numPoints*sizeof(double));
	(*a).x_U = malloc(numPoints*sizeof(double));
	(*a).y_L = malloc(numPoints*sizeof(double));
	(*a).y_U = malloc(numPoints*sizeof(double));
	(*a).numPoints = numPoints;
}

void free_airfoil(struct airfoil * a){
	free((*a).x_c);
	free((*a).y_c);
	free((*a).x_L);
	free((*a).y_L);
	free((*a).x_U);
	free((*a).y_U);
	free(a);
}

void airfoil_gen (struct airfoil * a, double t){
	double * x_L = (*a).x_L;
	double * y_L = (*a).y_L;
	double * x_U = (*a).x_U;
	double * y_U = (*a).y_U;
	double * x_c = (*a).x_c;
	double * y_c = (*a).y_c;
	int numPoints = (*a).numPoints;

	double y;
	double * dy_c = diff(y_c,numPoints);
	double * dx   = diff(x_c,numPoints);
	double theta;

	for(int i = 0; i < numPoints; i++){
		y = t/0.2*(0.2969 * sqrt(x_c[i]) - 0.1260 * x_c[i] - 0.3516 * x_c[i] * x_c[i] + 0.2843 * x_c[i] * x_c[i] * x_c[i] - 0.1015 * x_c[i] * x_c[i] * x_c[i] * x_c[i]);
		theta = atan2(dy_c[i],dx[i]);
		x_L[i] = x_c[i] + y * sin(theta);
		y_L[i] = y_c[i] - y * cos(theta);
		x_U[i] = x_c[i] - y * sin(theta);
		y_U[i] = y_c[i] + y * cos(theta);
	}
}

void print_airfoil(struct airfoil * a){
	printf("i \t x_c \t y_c \t x_U \t y_U \t x_L \t y_L\n");
	printf("------\t-----\t------\t------\t------\t------\t------\n");
	int i;
	for (i=0; i<(*a).numPoints;i++){
		printf("%d\t%6.2f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%7.4f\n",i, *a->x_c, *a->y_c, *a->x_U, *a->y_U, *a->x_L, *a->y_L);
		(*a).x_c++;
		(*a).y_c++;
		(*a).x_U++;
		(*a).y_U++;
		(*a).x_L++;
		(*a).y_L++;
	}
}

void meanLineNACA4(double m, double p, struct airfoil * a){
	double * x = (*a).x_c;
	double * y = (*a).y_c;

	for(int i=0; i < (*a).numPoints; i++){
		if(*x<p){
			*y = m * *x / (p*p) * (2*p - *x);
		} else {
			*y = m * (1 - *x) / ((1 - p)*(1 - p)) * (1 + (*x) - 2*p);
		}
		x++;
		y++;
	}
}

int NACA4parse(int num, double * m, double * p, double * t){
	*t  = num%100;
	num = num - *t;
	*t /= 100;
	*p  = num%1000;
	num = num-*p;
	*p /= 1000;
	*m  = num;
	*m /= 100000;
	return 0;
}

int NACA4(int num, int numPoints){
	double m, p, t;
	NACA4parse(num,&m,&p,&t);
	printf("%f %f %f\n",m,p,t);
	struct airfoil a;
	allocate_airfoil(&a, numPoints);
	x_straight(&a);
	meanLineNACA4(m,p,&a);
	airfoil_gen (&a, t); //TODO: insert t
	print_airfoil(&a);
	return 0;
}
