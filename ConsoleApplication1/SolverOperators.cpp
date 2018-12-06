#include "stdafx.h"
#include "SolverOperators.h"

#define Inverse_Order 1 //2 uses inverse definition with L². 1 truncates to L.

//Global variables are useful here because these functions get called many times. Therefore reallocation is not necessary because the variables are overwritten.

//Global variables for SolveIntermVelocity2D
double* ru;
double* rv;
double* HLu; //HL is a dummy variable that stores H and L values to compute r
double* HLv;
bool SolveInterm_alloc = 0;

//Global variables for ROperator
double* Lu;
double* Lv;
double* Lu2;
double* Lv2;
double dtOver2Re;
bool R_Op_alloc = 0;

//Global variables for ProjectionStep
double *Gup;
double *Gvp;
double *Rup_inv;
double *Rvp_inv;
double *Fxp;
double *Fyp;
bool Projection_alloc = 0;

//Global variables for ConvectiveOperator
double du2dx, dv2dy, duvdx, duvdy;
double uAvgLeft, uAvgRight, vAvgTop, vAvgBottom;
double uAvgTop, uAvgBottom, vAvgLeft, vAvgRight;

//Global variables for LaplaceOperator
double d2udx2, d2udy2, d2vdx2, d2vdy2;

//Global variables for DivergenceOperator
double dudx, dvdy;


void SolverOperators::SolveIntermVelocity2D(double *uPrev, double *vPrev, double *uCurrent, double *vCurrent, BC *B, double *uStar, double *vStar, double Re, double dx, double dy, double dt, int nx, int ny, double reltol) {
	//Solves for the interm. velocity Av*=r. Uses the Conjugate Gradient method.
	//Uses the previous (2 timesteps before) velocity field vPrev and the current (1 timestep before) velocity field vCurrent to
	//infer vStar. Re is the Reynolds number of the flow and dx,dy,dt are the grid and time spacings. nx and ny are used to interpret the velocity field matrix.
	//reltol and set the relative and absolute tolerances.
	//B contains the bounday conditions object

	//Preallocates
	if (!SolveInterm_alloc) {
		ru = new double[(nx - 1)*ny]; //ru is the same size as u
		rv = new double[nx*(ny - 1)]; //rv is the same size as v
		HLu = new double[(nx - 1)*ny];
		HLv = new double[nx*(ny - 1)];
		SolveInterm_alloc = 1;
	}

	//===========================================================================================
	//Computes the r vector

	//Computes first H_prev to save memory.
	ConvectiveOperator2D(uPrev, vPrev, HLu, HLv, B, dx, dy, nx, ny);

	//Remembering expression for r, given Adams-Bashfort2:
	//r=velCurrent/dt + (1/2Re)*L(velCurrent) - (1.5*H_Current - 0.5*H_Prev)
	//Start storing r:
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < (nx - 1); i++) {
			ru[j * (nx - 1) + i] = 0.5 * HLu[j * (nx - 1) + i] * dt;
		}
	}
	for (int j = 0; j < (ny - 1); j++) {
		for (int i = 0; i < nx; i++) {
			rv[j * nx + i] = 0.5 * HLv[j * nx + i] * dt;
		}
	}


	//Performs the same for the H_Current
	ConvectiveOperator2D(uCurrent, vCurrent, HLu, HLv, B, dx, dy, nx, ny);

	//Stores the current H in r:
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < (nx - 1); i++) {
			ru[j * (nx - 1) + i] -= 1.5 * HLu[j * (nx - 1) + i] * dt;
		}
	}
	for (int j = 0; j < (ny - 1); j++) {
		for (int i = 0; i < nx; i++) {
			rv[j * nx + i] -= 1.5 * HLv[j * nx + i] * dt;
		}
	}

	//Performs the same for L_Current
	LaplaceOperator2D(uCurrent, vCurrent, HLu, HLv, dx, dy, nx, ny);

	//Stores the current L*dt/2Re in r:
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < (nx - 1); i++) {
			ru[j * (nx - 1) + i] += HLu[j * (nx - 1) + i] * dt / (2 * Re) + uCurrent[j * (nx - 1) + i];
		}
	}
	for (int j = 0; j < (ny - 1); j++) {
		for (int i = 0; i < nx; i++) {
			rv[j * nx + i] += HLv[j * nx + i] * dt / (2 * Re) + vCurrent[j * nx + i];
		}
	}

	//Gets the boundary conditions for Laplace at current timestep
	B->Laplacian(HLu, HLv);
	//Stores the current L*dt/2Re in r:
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < (nx - 1); i++) {
			ru[j * (nx - 1) + i] += HLu[j * (nx - 1) + i] * dt / (2 * Re);
		}
	}
	for (int j = 0; j < (ny - 1); j++) {
		for (int i = 0; i < nx; i++) {
			rv[j * nx + i] += HLv[j * nx + i] * dt / (2 * Re);
		}
	}


	//Advances boundary conditions 1 time step and computes the laplacian again
	B->ComputeConvectiveBCs();
	B->Laplacian(HLu, HLv);
	//Stores the current L*dt/2Re in r:
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < (nx - 1); i++) {
			ru[j * (nx - 1) + i] += HLu[j * (nx - 1) + i] * dt / (2 * Re);
		}
	}
	for (int j = 0; j < (ny - 1); j++) {
		for (int i = 0; i < nx; i++) {
			rv[j * nx + i] += HLv[j * nx + i] * dt / (2 * Re);
		}
	}

	//Solves using Conjugate Gradient
	CGMomentum2D(&R_Operator, uCurrent, vCurrent, ru, rv, uStar, vStar, Re, dx, dy, dt, nx, ny, reltol);

	
}

void SolverOperators::R_Operator(double *u, double *v, double *res_u, double *res_v, double Re, double dx, double dy, double dt, int nx, int ny) {
	//R operator is the matrix that multiplies v* (fractional velocity field).
	//Results output in res_u, res_v. Form of u and v

	//Allocates result of Laplace operator
	if (!R_Op_alloc) {
		Lu = new double[(nx - 1) * ny];
		Lv = new double[nx * (ny - 1)];
		Lu2 = new double[(nx - 1) * ny];
		Lv2 = new double[nx * (ny - 1)];
		R_Op_alloc = 1;
	}

	//Gets result from Laplace
	LaplaceOperator2D(u, v, Lu, Lv, dx, dy, nx, ny);

	//Forms the R matrix:
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < (nx - 1); i++) {
			res_u[j * (nx - 1) + i] = u[j * (nx - 1) + i] - (dt / (2 * Re)) * Lu[j * (nx - 1) + i];
		}
	}
	for (int j = 0; j < (ny - 1); j++) {
		for (int i = 0; i < nx; i++) {
			res_v[j * nx + i] = v[j * nx + i] - (dt / (2 * Re)) * Lv[j * nx + i];
		}
	}
	
}

void SolverOperators::R_Operator_Inv(double *u, double *v, double *res_u, double *res_v, double Re, double dx, double dy, double dt, int nx, int ny) {
	//Outputs the inverted R operator
	

	//Gets result from Laplace
	LaplaceOperator2D(u, v, Lu, Lv, dx, dy, nx, ny);
#if Inverse_Order == 2
	LaplaceOperator2D(Lu, Lv, Lu2, Lv2, dx, dy, nx, ny);
#endif


	//Forms the Rinv matrix:
	dtOver2Re = (dt / (2 * Re));
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < (nx - 1); i++) {
#if Inverse_Order == 1
			res_u[j * (nx - 1) + i] = u[j * (nx - 1) + i] + dtOver2Re * Lu[j * (nx - 1) + i];
#elif Inverse_Order == 2
			res_u[j * (nx - 1) + i] = u[j * (nx - 1) + i] + dtOver2Re * Lu[j * (nx - 1) + i] + dtOver2Re * dtOver2Re * Lu2[j * (nx - 1) + i];
#endif
		}
	}
	for (int j = 0; j < (ny - 1); j++) {
		for (int i = 0; i < nx; i++) {
#if Inverse_Order == 1
			res_v[j * nx + i] = v[j * nx + i] + dtOver2Re * Lv[j * nx + i];
#elif Inverse_Order == 2
			res_v[j * nx + i] = v[j * nx + i] + dtOver2Re * Lv[j * nx + i] + dtOver2Re * dtOver2Re * Lv2[j * nx + i];
#endif
		}
	}

}

void SolverOperators::ProjectionStep(double *uStar, double *vStar, double *uNext, double *vNext, double *P, Body *B, double Re, double dx, double dy, double dt, int nx, int ny) {
	//Performs the projection step uNext=uStar-dt Rinv G P
	
	//Preallocates
	if (!Projection_alloc) {
		Gup = new double[(nx - 1)*ny];
		Gvp = new double[nx*(ny - 1)];
		Rup_inv = new double[(nx - 1)*ny];
		Rvp_inv = new double[nx*(ny - 1)];
		Fxp = new double[(nx - 1)*ny];
		Fyp = new double[nx*(ny - 1)];
		Projection_alloc = 1;
	}


	SolverOperators::GradientOperator2D(P, Gup, Gvp, dx, dy, nx, ny);
	IB_Operators::RegularizationOperator(Fxp, Fyp, B->Fx, B->Fy, B, dx, dy, nx, ny);

	//Updates u component of Gu to include immersed boundary forces
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < (nx - 1); i++) {
			Gup[j * (nx - 1) + i] = (Gup[j * (nx - 1) + i] - Fxp[j * (nx - 1) + i]) ;
		}
	}

	//Updates v component of Gv to include immersed boundary forces
	for (int j = 0; j < (ny - 1); j++) {
		for (int i = 0; i < nx; i++) {
			//---v component----
			Gvp[j * nx + i] = (Gvp[j * nx + i] - Fyp[j * nx + i]); //
		}
	}


	SolverOperators::R_Operator_Inv(Gup, Gvp, Rup_inv, Rvp_inv, Re, dx, dy, dt, nx, ny); //Applies R_inv as always

	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < (nx - 1); i++) {
			uNext[j * (nx - 1) + i] = uStar[j * (nx - 1) + i] - Rup_inv[j * (nx - 1) + i];
		}
	}
	for (int j = 0; j < (ny - 1); j++) {
		for (int i = 0; i < nx; i++) {
			vNext[j * nx + i] = vStar[j * nx + i] - Rvp_inv[j * nx + i];
		}
	}

}


void SolverOperators::CGMomentum2D(void (*A_Fn) (double *u, double *v, double *res_u, double *res_v, double Re, double dx, double dy, double dt, int nx, int ny), double *uGuess, double *vGuess, double *ru, double *rv,
															double *uResult, double *vResult, double Re, double dx, double dy, double dt, int nx, int ny, double reltol) {
	//===========================================================================================
	//===========================================================================================
	//===============================CONJUGATE GRADIENT METHOD===================================
	//===========================================================================================
	//===========================================================================================
	//Solves the Av*=r equation using the beautiful Conjugate Gradient Method
	//A_Fn is a pointer to a function that is evaluated (matrix A)
	//uGuess and vGuess are the initial guesses to start the algorithm
	//r is a pointer to the right-hand-side vector
	//B is a boundary condition object
	//uResult and vResult are the results output
	//Re, dx, dy, nx, ny, reltol are self-evident.


	//Preallocates
	double *d_u, *d_v;
	double *r_u, *r_v;

	double *uTemp, *vTemp;
	d_u = new double[(nx - 1)*ny];
	d_v = new double[nx*(ny - 1)];
	r_u = new double[(nx - 1)*ny];
	r_v = new double[nx*(ny - 1)];

	uTemp = new double[(nx - 1)*ny];
	vTemp = new double[nx*(ny - 1)];

	int iter = 0;
	double alpha, alphaNum, alphaDen;
	double beta, betaNum, betaDen;
	double delta, delta0, deltaOld;

	//Define initial residue
	delta = 0;
	A_Fn(uGuess, vGuess, uTemp, vTemp, Re, dx, dy, dt, nx, ny);
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < (nx - 1); i++) {
			r_u[j * (nx - 1) + i] = ru[j * (nx - 1) + i] - uTemp[j * (nx - 1) + i];
			d_u[j * (nx - 1) + i] = r_u[j * (nx - 1) + i];
			delta += r_u[j * (nx - 1) + i] * r_u[j * (nx - 1) + i];
			uResult[j * (nx - 1) + i] = uGuess[j * (nx - 1) + i];
		}
	}
	for (int j = 0; j < (ny - 1); j++) {
		for (int i = 0; i < nx; i++) {
			r_v[j * nx + i] = rv[j * nx + i] - vTemp[j * nx + i];
			d_v[j * nx + i] = r_v[j * nx + i];
			delta += r_v[j * nx + i] * r_v[j * nx + i];
			vResult[j * nx + i] = vGuess[j * nx + i];
		}
	}

	delta0 = delta;

	while (delta > (reltol * reltol * delta0)) { //Stop condition
		iter++;

		//Updates alpha
		alphaDen = 0;
		A_Fn(d_u, d_v, uTemp, vTemp, Re, dx, dy, dt, nx, ny);
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < (nx - 1); i++) {
				alphaDen += d_u[j * (nx - 1) + i] * uTemp[j * (nx - 1) + i];
			}
		}
		for (int j = 0; j < (ny - 1); j++) {
			for (int i = 0; i < nx; i++) {
				alphaDen += d_v[j * nx + i] * vTemp[j * nx + i];
			}
		}
		alpha = delta / alphaDen;

		//Updates uStar, vStar
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < (nx - 1); i++) {
				uResult[j * (nx - 1) + i] = uResult[j * (nx - 1) + i] + alpha * d_u[j * (nx - 1) + i];
			}
		}
		for (int j = 0; j < (ny - 1); j++) {
			for (int i = 0; i < nx; i++) {
				vResult[j * nx + i] = vResult[j * nx + i] + alpha * d_v[j * nx + i];
			}
		}

		
		//Quick formula to update the residuals
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < (nx - 1); i++) {
				r_u[j * (nx - 1) + i] = r_u[j * (nx - 1) + i] - alpha * uTemp[j * (nx - 1) + i];
			}
		}
		for (int j = 0; j < (ny - 1); j++) {
			for (int i = 0; i < nx; i++) {
				r_v[j * nx + i] = r_v[j * nx + i] - alpha * vTemp[j * nx + i];
			}
		}

		//Updates delta
		deltaOld = delta;
		delta = 0;
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < (nx - 1); i++) {
				delta += r_u[j * (nx - 1) + i] * r_u[j * (nx - 1) + i];
			}
		}
		for (int j = 0; j < (ny - 1); j++) {
			for (int i = 0; i < nx; i++) {
				delta += r_v[j * nx + i] * r_v[j * nx + i];
			}
		}


		//Updates Beta
		beta = delta / deltaOld;


		//Updates d
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < (nx - 1); i++) {
				d_u[j * (nx - 1) + i] = r_u[j * (nx - 1) + i] + beta * d_u[j * (nx - 1) + i];
			}
		}
		for (int j = 0; j < (ny - 1); j++) {
			for (int i = 0; i < nx; i++) {
				d_v[j * nx + i] = r_v[j * nx + i] + beta * d_v[j * nx + i];
			}
		}
	}
	
	delete[] d_u;//deletes to save memory
	delete[] d_v;
	delete[] r_u;//deletes to save memory
	delete[] r_v;

	delete[] uTemp;//deletes to save memory
	delete[] vTemp;
}

void SolverOperators::ConvectiveOperator2D(double *u, double *v, double *Hu, double *Hv, BC *B, double dx, double dy, int nx, int ny) {
	//Outputs in H the 2D convective operator Hu = ddx(uu) + ddy(uv),  Hv = ddx(uv) + ddy(vv)
	//Uses the central difference from staggered grid for all derivatives. Boundaries have to be dealt with differntly for u and v. Terms from boundaries come from B.


	//Inits Hu and Hv (Copy this snippet where the call is made)
	//Hu = new double[(nx - 1)*ny];
	//Hv = new double[nx*(ny - 1)];


	//Evaluates Hu
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < (nx - 1); i++) {
			//u derivatives
			if (i == 0) { //Uses left boundary u
				uAvgLeft = (u[j * (nx - 1) + i] + B->GetVal(BC::VEC_U_LEFT, j)) / 2;
				uAvgRight = (u[j * (nx - 1) + (i + 1)] + u[j * (nx - 1) + i]) / 2;
			}
			else if (i == (nx - 2)) { //Uses right boundary u
				uAvgLeft = (u[j * (nx - 1) + i] + u[j * (nx - 1) + (i - 1)]) / 2;
				uAvgRight = (B->GetVal(BC::VEC_U_RIGHT, j) + u[j * (nx - 1) + i]) / 2;
			}
			else {
				uAvgLeft = (u[j * (nx - 1) + i] + u[j * (nx - 1) + (i - 1)]) / 2;
				uAvgRight = (u[j * (nx - 1) + (i + 1)] + u[j * (nx - 1) + i]) / 2;
			}
			
			if (j == 0) { //Uses bottom boundary v
				vAvgTop = (v[j * nx + i] + v[j * nx + (i + 1)]) / 2;
				vAvgBottom = (B->GetVal(BC::VEC_V_BOTTOM, i) + B->GetVal(BC::VEC_V_BOTTOM, i + 1)) / 2;

				uAvgTop = (u[j * (nx - 1) + i] + u[(j + 1) * (nx - 1) + i]) / 2;
				uAvgBottom = (u[j * (nx - 1) + i] + B->GetVal(BC::VEC_U_BOTTOM, i)) / 2;
			}
			else if (j == (ny - 1)) { //Uses top boundary v
				vAvgTop = (B->GetVal(BC::VEC_V_TOP, i) + B->GetVal(BC::VEC_V_TOP, i + 1)) / 2;
				vAvgBottom = (v[(j - 1) * nx + i] + v[(j - 1) * nx + (i + 1)]) / 2;

				uAvgTop = (u[j * (nx - 1) + i] + B->GetVal(BC::VEC_U_TOP, i)) / 2;
				uAvgBottom = (u[j * (nx - 1) + i] + u[(j - 1) * (nx - 1) + i]) / 2;
			}
			else {
				vAvgTop = (v[j * nx + i] + v[j * nx + (i + 1)]) / 2;
				vAvgBottom = (v[(j - 1) * nx + i] + v[(j - 1) * nx + (i + 1)]) / 2;

				uAvgTop = (u[j * (nx - 1) + i] + u[(j + 1) * (nx - 1) + i]) / 2;
				uAvgBottom = (u[j * (nx - 1) + i] + u[(j - 1) * (nx - 1) + i]) / 2;
			}
			
			//Evaluates Hu
			du2dx = (uAvgRight * uAvgRight - uAvgLeft * uAvgLeft) / dx;
			duvdy = (uAvgTop * vAvgTop - uAvgBottom * vAvgBottom) / dy;
			Hu[j * (nx - 1) + i] = du2dx + duvdy;
		}
	}

	//Evaluates Hv
	for (int j = 0; j < (ny - 1); j++) {
		for (int i = 0; i < nx; i++) {
			//v derivatives
			if (i == 0) {
				vAvgLeft = (B->GetVal(BC::VEC_V_LEFT, j) + v[j * nx + i]) / 2;
				vAvgRight = (v[j * nx + (i + 1)] + v[j * nx + i]) / 2;

				uAvgLeft = (B->GetVal(BC::VEC_U_LEFT, j) + B->GetVal(BC::VEC_U_LEFT, j + 1)) / 2;
				uAvgRight = (u[j * (nx - 1) + i] + u[(j + 1) * (nx - 1) + i]) / 2;
			}
			else if (i == (nx - 1)) {
				vAvgLeft = (v[j * nx + (i - 1)] + v[j * nx + i]) / 2;
				vAvgRight = (B->GetVal(BC::VEC_V_RIGHT, j) + v[j * nx + i]) / 2;

				uAvgLeft = (u[j * (nx - 1) + (i - 1)] + u[(j + 1) * (nx - 1) + (i - 1)]) / 2;
				uAvgRight = (B->GetVal(BC::VEC_U_RIGHT, j) + B->GetVal(BC::VEC_U_RIGHT, j + 1)) / 2;
			}
			else {
				vAvgLeft = (v[j * nx + (i - 1)] + v[j * nx + i]) / 2;
				vAvgRight = (v[j * nx + (i + 1)] + v[j * nx + i]) / 2;

				uAvgLeft = (u[j * (nx - 1) + (i - 1)] + u[(j + 1) * (nx - 1) + (i - 1)]) / 2;
				uAvgRight = (u[j * (nx - 1) + i] + u[(j + 1) * (nx - 1) + i]) / 2;
			}

			if (j == 0) {
				vAvgTop = (v[j * nx + i] + v[(j + 1) * nx + i]) / 2;
				vAvgBottom = (B->GetVal(BC::VEC_V_BOTTOM, i) + v[j * nx + i]) / 2;
			}
			else if (j == (ny - 2)) {
				vAvgTop = (v[j * nx + i] + B->GetVal(BC::VEC_V_TOP, i)) / 2;
				vAvgBottom = (v[(j - 1) * nx + i] + v[j * nx + i]) / 2;
			}
			else {
				vAvgTop = (v[j * nx + i] + v[(j + 1) * nx + i]) / 2;
				vAvgBottom = (v[(j - 1) * nx + i] + v[j * nx + i]) / 2;
			}
						
			//Evaluates Hv
			dv2dy = (vAvgTop * vAvgTop - vAvgBottom * vAvgBottom) / dy;
			duvdx = (uAvgRight * vAvgRight - uAvgLeft * vAvgLeft) / dx;
			Hv[j * nx + i] = duvdx + dv2dy;
			}
		}
	
	}	

void SolverOperators::LaplaceOperator2D(double *u, double *v, double *Lu, double *Lv, double dx, double dy, int nx, int ny) {
	//Outputs in L the 2D Laplace operator Lu = d2udx2 + d2udy2,  Lv = d2vdx2 + d2vdy2
	//Uses the central difference for all derivatives. Leaves zero at boundary nodes. For BCs use the LaplaceBC2D function


	//Inits Lu and Lv (Copy this snippet where the call is made)
	//Lu = new double[(nx - 1) * ny];
	//Lv = new double[nx * (ny - 1)];

	//Evaluates Lu
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < (nx - 1); i++) {
			//---u derivatives----
			//   d2udx2
			if (i == 0) {
				d2udx2 = (-2 * u[j * (nx - 1) + i] + u[j * (nx - 1) + (i + 1)]) / (dx * dx); //Needs boundary node			
			}
			else if (i == (nx - 2)) {
				d2udx2 = (u[j * (nx - 1) + (i - 1)] - 2 * u[j * (nx - 1) + i]) / (dx * dx); //Needs boundary node
			}
			else {
				d2udx2 = (u[j * (nx - 1) + (i - 1)] - 2 * u[j * (nx - 1) + i] + u[j * (nx - 1) + (i + 1)]) / (dx * dx); //Midpoint two-sided
			}
			//   d2udy2
			if (j == 0) {
				d2udy2 = (- 2 * u[j * (nx - 1) + i] + u[(j + 1) * (nx - 1) + i]) / (dy * dy); //Needs boundary node
			}
			else if (j == (ny - 1)) {
				d2udy2 = (u[(j - 1) * (nx - 1) + i] - 2 * u[j * (nx - 1) + i]) / (dy * dy); //Needs boundary node
			}
			else {
				d2udy2 = (u[(j - 1) * (nx - 1) + i] - 2 * u[j * (nx - 1) + i] + u[(j + 1) * (nx - 1) + i]) / (dy * dy); //Midpoint two-sided
			}

			//Evaluates Lu
			Lu[j * (nx - 1) + i] = d2udx2 + d2udy2; //Laplacian is calculated at the edges (as vel. vectors)
		}
	}

	//Evaluate Lv
	for (int j = 0; j < (ny - 1); j++) {
		for (int i = 0; i < nx; i++) {
			//---v derivatives----
			//   d2vdx2
			if (i == 0) {
				d2vdx2 = (- 2 * v[j * nx + i] + v[j * nx + (i + 1)]) / (dx * dx); //Needs boundary node
			}
			else if (i == (nx - 1)) {
				d2vdx2 = (v[j * nx + (i - 1)] - 2 * v[j * nx + i]) / (dx * dx); //Needs boundary node
			}
			else { 
				d2vdx2 = (v[j * nx + (i - 1)] - 2 * v[j * nx + i] + v[j * nx + (i + 1)]) / (dx * dx); //Midpoint two-sided
			}
			//   d2vdy2
			if (j == 0) {
				d2vdy2 = (- 2 * v[j * nx + i] + v[(j + 1) * nx + i]) / (dy * dy); //Needs boundary node
			}
			else if (j == (ny - 2)) {
				d2vdy2 = (v[(j - 1) * nx + i] - 2 * v[j * nx + i]) / (dy * dy); //Needs boundary node
			}
			else {
				d2vdy2 = (v[(j - 1) * nx + i] - 2 * v[j * nx + i] + v[(j + 1) * nx + i]) / (dy * dy); //Midpoint two-sided
			}

			//Evaluates Lv
			Lv[j * nx + i] = d2vdx2 + d2vdy2; //Laplacian is calculated at the edges (as vel. vectors)
		}
	}
}

void SolverOperators::GradientOperator2D(double *P, double *Gu, double *Gv, double dx, double dy, int nx, int ny) {
	//Performs the discrete 2D gradient operator and returns the vectors in Gu and Gv. No bcs needed
	int i;
	int j;

	//P[ny * nx - 1] = 0; //Enforces pinned point

	for (j = 0; j < ny; j++) {
		for (i = 0; i < (nx - 1); i++) {
			//---u component----
			Gu[j * (nx - 1) + i] = (P[j * nx + (i + 1)] - P[j * nx + i]) / dx; //Remember grid indices are different than in class!
		}
	}	
	
	for (j = 0; j < (ny - 1); j++) {
		for (i = 0; i < nx; i++) {
			//---v component----
			Gv[j * nx + i] = (P[(j + 1) * nx + i] - P[j * nx + i]) / dy; //
		}
	}

}

void SolverOperators::DivergenceOperator2D(double *u, double *v, double *D, double dx, double dy, int nx, int ny) {
	//Evaluates the discrete divergence of a vector field. Result in a P-centered grid, output in D.
	//BC vector needs to be added 

	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {
			//---u component----
			if (i == 0) {
				dudx = (u[j * (nx - 1) + i]) / dx; //BC
			}
			else if (i == (nx - 1)) {
				dudx = (- u[j * (nx - 1) + (i - 1)]) / dx; //BC
			}
			else {
				dudx = (u[j * (nx - 1) + i] - u[j * (nx - 1) + (i - 1)]) / dx; //Interior point
			}

			//---v component----
			if (j == 0) {
				dvdy = (v[j * nx + i]) / dy; //BC
			}
			else if (j == (ny - 1)) {
				dvdy = (- v[(j - 1) * nx + i]) / dy; //BC
			}
			else {
				dvdy = (v[j * nx + i] - v[(j - 1) * nx + i]) / dy; //Interior point
			}

			//Divergence
			D[j * nx + i] = dudx + dvdy;
		}
	}

	D[ny * nx - 1] = 0; //Pinned point at top-right corner. Won't be used, but if it happens to be used will generate a division by zero! 
}
