#include <math.h>
#include <stdlib.h>
#include <stdio.h>


int numPoints = 41;

struct airfoil {
	char * desc; //short text description
	double * x_c;  //mean camber line
	double * y_c;
	double * x_U;  //upper surface
	double * y_U;
	double * x_L;  // lower surface
	double * y_L;
};

void x_straight (double *);
void x_curved   (struct airfoil *);
void meanLineNACA4(double, double, double *, double *);
void airfoil_gen(double *,double *,double *, double *,double *,double *);


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

void allocate_airfoil(struct airfoil * a){
	(*a).x_c = malloc(numPoints*sizeof(double));
	(*a).y_c = malloc(numPoints*sizeof(double));
	(*a).x_L = malloc(numPoints*sizeof(double));
	(*a).x_U = malloc(numPoints*sizeof(double));
	(*a).y_L = malloc(numPoints*sizeof(double));
	(*a).y_U = malloc(numPoints*sizeof(double));
}

void x_straight(double * xPoints){
	double dx = 1.0 / (numPoints-1);
	int i;

	for(i = 0; i < numPoints; i++){
		*xPoints = dx * i;
		xPoints++;
	}
}

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

void print_airfoil(struct airfoil * a){
	printf("i \t x_c \t y_c \t x_U \t y_U \t x_L \t y_L\n");
	printf("------\t-----\t------\t------\t------\t------\t------\n");
	int i;
	for (i=0; i<numPoints;i++){
		printf("%d\t%6.2f\t%6.4f\t%6.4f\t%6.4f\t%6.4f\t%7.4f\n",i, *a->x_c, *a->y_c, *a->x_U, *a->y_U, *a->x_L, *a->y_L);
		(*a).x_c++;
		(*a).y_c++;
		(*a).x_U++;
		(*a).y_U++;
		(*a).x_L++;
		(*a).y_L++;
	}
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

int main(){
	struct airfoil a;
	allocate_airfoil(&a);
	x_straight(a.x_c);
	meanLineNACA4(0.02,0.4,a.x_c,a.y_c);
	airfoil_gen (a.x_c,a.y_c,a.x_U,a.y_U,a.x_L,a.y_L);
	print_airfoil(&a);
	return 0;
}
