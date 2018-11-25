#include "stdafx.h"
#include "ImmersedBoundaryOperators.h"


void IB_Operators::InterpolationOperator(double *u, double *v, double *uInterp, double *vInterp, Body *B, double dx, double dy, int nx, int ny) {
	//Interpolates the velocity field into the Lagrangian coordinate system where the body sits
	//A 3-grid-cell delta function is used

	int Bi_min, Bi_max, Bj_min, Bj_max;
	double curr_x, curr_y;
	double tempVal;

	//====Evaluates uInterp=====
	for (int k = 0; k < B->ns; k++) { //Runs through all points in the body
		//Finds the bounds for interpolation
		Bi_min = (int) round(((B->x[k] / dx) - 1.) - 1.5);
		Bi_max = (int) round(((B->x[k] / dx) - 1.) + 1.5);
		Bj_min = (int) round(((B->y[k] / dy) - 0.5) - 1.5);
		Bj_max = (int) round(((B->y[k] / dy) - 0.5) + 1.5);

		tempVal = 0;
		for (int j = Bj_min; j <= Bj_max; j++) {
			for (int i = Bi_min; i <= Bi_max; i++) {
				curr_x = dx * (i + 1.); //Current x,y coordinates in u vector coordinates
				curr_y = dy * (j + 0.5);

				tempVal += IB_Operators::DiscreteDeltaFunction2D(curr_x - B->x[k], curr_y - B->y[k], dx, dy) * u[j*(nx - 1) + i];
			}
		}

		uInterp[k] = tempVal;
	}


	//====Evaluates vInterp=====
	for (int k = 0; k < B->ns; k++) { //Runs through all points in the body
								   //Finds the bounds for interpolation
		Bi_min = (int)round(((B->x[k] / dx) - 0.5) - 1.5);
		Bi_max = (int)round(((B->x[k] / dx) - 0.5) + 1.5);
		Bj_min = (int)round(((B->y[k] / dy) - 1.) - 1.5);
		Bj_max = (int)round(((B->y[k] / dy) - 1.) + 1.5);

		tempVal = 0;
		for (int j = Bj_min; j <= Bj_max; j++) {
			for (int i = Bi_min; i <= Bi_max; i++) {
				curr_x = dx * (i + 0.5); //Current x,y coordinates in u vector coordinates
				curr_y = dy * (j + 1);

				tempVal += IB_Operators::DiscreteDeltaFunction2D(curr_x - B->x[k], curr_y - B->y[k], dx, dy) * v[j*nx + i];
			}
		}

		vInterp[k] = tempVal;
	}
}


void IB_Operators::RegularizationOperator(double *FuRegularized, double *FvRegularized, double *Fu, double *Fv, Body *B, double dx, double dy, int nx, int ny) {
	//Interpolates the body forces into the Eulerian coordinate system where the fluid sits
	//A 3-grid-cell delta function is used


	int Bi_min, Bi_max, Bj_min, Bj_max;
	double curr_x, curr_y;

	//====Evaluates FuRegularized=====
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i <= (nx - 1); i++) { //Zeros all forces in fluid domain
			FuRegularized[j*(nx - 1) + i] = 0;
		}
	}

	for (int k = 0; k < B->ns; k++) { //Runs through all points in the body
		//Finds the bounds where regularization is needed
		Bi_min = (int) round(((B->x[k] / dx) - 1.) - 1.5);
		Bi_max = (int) round(((B->x[k] / dx) - 1.) + 1.5);
		Bj_min = (int) round(((B->y[k] / dy) - 0.5) - 1.5);
		Bj_max = (int) round(((B->y[k] / dy) - 0.5) + 1.5);

		for (int j = Bj_min; j <= Bj_max; j++) {
			for (int i = Bi_min; i <= Bi_max; i++) {
				curr_x = dx * (i + 1.); //Current x,y coordinates in u vector coordinates
				curr_y = dy * (j + 0.5);

				FuRegularized[j*(nx - 1) + i] += IB_Operators::DiscreteDeltaFunction2D(curr_x - B->x[k], curr_y - B->y[k], dx, dy) * Fu[k]; //Updates regularized fluid forces
			}
		}
	}

	//====Evaluates FvRegularized=====
	for (int j = 0; j < (ny - 1); j++) {
		for (int i = 0; i <= nx; i++) { //Zeros all forces in fluid domain
			FvRegularized[j*nx + i] = 0;
		}
	}

	for (int k = 0; k < B->ns; k++) { //Runs through all points in the body
								   //Finds the bounds where regularization is needed
		Bi_min = (int) round(((B->x[k] / dx) - 0.5) - 1.5);
		Bi_max = (int) round(((B->x[k] / dx) - 0.5) + 1.5);
		Bj_min = (int) round(((B->y[k] / dy) - 1.) - 1.5);
		Bj_max = (int) round(((B->y[k] / dy) - 1.) + 1.5);

		for (int j = Bj_min; j <= Bj_max; j++) {
			for (int i = Bi_min; i <= Bi_max; i++) {
				curr_x = dx * (i + 0.5); //Current x,y coordinates in v vector coordinates
				curr_y = dy * (j + 1.);

				FvRegularized[j*nx + i] += IB_Operators::DiscreteDeltaFunction2D(curr_x - B->x[k], curr_y - B->y[k], dx, dy) * Fv[k]; //Updates regularized fluid forces
			}
		}
	}
}



double IB_Operators::DiscreteDeltaFunction2D(double x, double y, double dx, double dy) {
	//Implements the 3-grid-cell discrete delta function for a point (x,y). (x,y)=(0,0) is the delta function peak
	double x_dimless = fabs(x / dx);
	double y_dimless = fabs(y / dy);
	double Delta_x, Delta_y;

	//Delta function in x direction
	if (x_dimless <= 0.5) {
		Delta_x = (1 + sqrt(1 - 3 * x_dimless * x_dimless)) / 3;
	}
	else if (x_dimless <= 1.5) {
		Delta_x = (5 - 3 * x_dimless - sqrt(1 - 3 * (1 - x_dimless) * (1 - x_dimless))) / 6;
	}
	else {
		Delta_x = 0;
	}

	//Delta function in y direction
	if (y_dimless <= 0.5) {
		Delta_y = (1 + sqrt(1 - 3 * y_dimless * y_dimless)) / 3;
	}
	else if (y_dimless <= 1.5) {
		Delta_y = (5 - 3 * y_dimless - sqrt(1 - 3 * (1 - y_dimless) * (1 - y_dimless))) / 6;
	}
	else {
		Delta_y = 0;
	}

	//Returns 2D delta function
	return Delta_x * Delta_y;
}