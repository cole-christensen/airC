struct airfoil {
        char * desc;   //short text description
        double * x_c;  //mean camber line
        double * y_c;
        double * x_U;  //upper surface
        double * y_U;
        double * x_L;  // lower surface
        double * y_L;
};

void x_straight (double *); /* evenly spaced points */
void x_curved   (struct airfoil *); /* points concentrated at L.E. and T.E. */
void allocate_airfoil(struct airfoil *);
