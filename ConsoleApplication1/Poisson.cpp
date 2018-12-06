#include "stdafx.h"
#include "Poisson.h"
#include "SolverOperators.h"

//Global variables for CG
double *d_p;
double *r_p;
double *s_p;
double *pStarOld;
double *pTemp;
double *fxTemp;
double *fyTemp;
double *d_fx;
double *d_fy;
double *r_fx;
double *r_fy;
double *s_fx;
double *s_fy;
double *fxOld;
double *fyOld;
bool CG_alloc = 0;

//Global variables for Laplacian
double *Gu;
double *Gv;
double *Fx;
double *Fy;
double *Ru_inv;
double *Rv_inv;
bool Lap_alloc = 0;

void Poisson::CGModifiedPressurePoisson2D(double *pGuess, double *rp, double *pResult, Body* B, double *uInterp, double *vInterp, double dx, double dy, double dt, double Re, int nx, int ny, double reltol) {
	//===========================================================================================
	//===========================================================================================
	//===============================CONJUGATE GRADIENT METHOD===================================
	//===========================================================================================
	//===========================================================================================
	//Solves the [QT Rinv Q Lambda=QT u* - r2] equation using the amazing Conjugate Gradient Method
	//matrix A that is the evaluation function is the Laplacian_Q function below
	//pGuess is the initial guess to start the algorithm in the Pressure. Guesses for Fx and Fy are put inside the Body B
	//rp is a pointer to the right-hand-side vector in pressure
	//uInterp and vInterp are pointers to the right-hand-side vector of interpolated velocities
	//B is a Body object that contains the set of lagrangian forces
	//uResult and vResult are the results output, as well as the body forces in B->Fx and B->Fy
	//Re, dx, dy, nx, ny, reltol are self-evident.


	//Preallocates once
	if (!CG_alloc) {
		d_p = new double[nx*ny];
		pStarOld = new double[nx*ny];
		r_p = new double[nx*ny];
		s_p = new double[nx*ny];
		pTemp = new double[nx*ny];
		fxTemp = new double[B->ns];
		fyTemp = new double[B->ns];
		d_fx = new double[B->ns];
		d_fy = new double[B->ns];
		r_fx = new double[B->ns];
		r_fy = new double[B->ns];
		s_fx = new double[B->ns];
		s_fy = new double[B->ns];
		fxOld = new double[B->ns];
		fyOld = new double[B->ns];
		CG_alloc = 1;
	}
	

	int iter = 0;
	double alpha, alphaDen;
	double beta;
	double delta, delta0, deltaOld;
	
	//Define initial residue
	delta = 0;
	Poisson::Laplacian_Q(pGuess, B->Fx, B->Fy, pTemp, fxTemp, fyTemp, B, dx, dy, dt, Re, nx, ny);
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {
			r_p[j * nx + i] = rp[j * nx + i] - pTemp[j * nx + i];
			//d_p[j * nx + i] = r_p[j * nx + i];
			pResult[j * nx + i] = pGuess[j * nx + i]; //Inits result
		}
	}
	for (int k = 0; k < B->ns; k++) {
		r_fx[k] = (uInterp[k] - B->vel_x[k]) - fxTemp[k]; //RHS of force equation is [Eu* - u_b]
		r_fy[k] = (vInterp[k] - B->vel_y[k]) - fyTemp[k];
		//d_fx[k] = r_fx[k];
		//d_fy[k] = r_fy[k];
		//Doesn't init B->Fx and Fy as it is not necessary
	}
	Preconditioner(r_p, r_fx, r_fy, d_p, d_fx, d_fy, dx, dy, nx, ny, B); //Applies the preconditioner
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {
			delta += r_p[j * nx + i] * d_p[j * nx + i];
		}
	}
	for (int k = 0; k < B->ns; k++) {
		delta += r_fx[k] * d_fx[k] + r_fy[k] * d_fy[k];
	}

	delta0 = delta;

	while (delta > (reltol * reltol * delta0)) {
		iter++;

		//Updates alpha
		alphaDen = 0;
		Poisson::Laplacian_Q(d_p, d_fx, d_fy, pTemp, fxTemp, fyTemp, B, dx, dy, dt, Re, nx, ny);
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				alphaDen += d_p[j * nx + i] * pTemp[j * nx + i];
			}
		}
		for (int k = 0; k < B->ns; k++) {
			alphaDen += d_fx[k] * fxTemp[k] + d_fy[k] * fyTemp[k];
		}
		alpha = delta / alphaDen;

		//Updates Result
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				pResult[j * nx + i] += alpha * d_p[j * nx + i];
			}
		}
		for (int k = 0; k < B->ns; k++) {
			B->Fx[k] += alpha * d_fx[k];
			B->Fy[k] += alpha * d_fy[k];
		}

		//Quick formula to update the residuals
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				r_p[j * nx + i] = r_p[j * nx + i] - alpha * pTemp[j * nx + i];
			}
		}
		for (int k = 0; k < B->ns; k++) {
			r_fx[k] = r_fx[k] - alpha * fxTemp[k];
			r_fy[k] = r_fy[k] - alpha * fyTemp[k];
		}

		//Updates s with preconditioner
		Preconditioner(r_p, r_fx, r_fy, s_p, s_fx, s_fy, dx, dy, nx, ny, B); //Applies the preconditioner
		
		//Updates delta
		deltaOld = delta;
		delta = 0;
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				delta += r_p[j * nx + i] * s_p[j * nx + i];
			}
		}
		for (int k = 0; k < B->ns; k++) {
			delta += r_fx[k] * s_fx[k] + r_fy[k] * s_fy[k];
		}


		//Updates Beta
		beta = delta / deltaOld;
		//Updates d
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				d_p[j * nx + i] = s_p[j * nx + i] + beta * d_p[j * nx + i];
			}
		}
		for (int k = 0; k < B->ns; k++) {
			d_fx[k] = s_fx[k] + beta * d_fx[k];
			d_fy[k] = s_fy[k] + beta * d_fy[k];
		}
	}	

}

void Poisson::Laplacian_Q(double *p, double *BFx, double *BFy, double *res_p, double *res_uInterp, double *res_vInterp, Body* B, double dx, double dy, double dt, double Re, int nx, int ny) {
	//Generates the Laplace operator from the pressure field divergence

	//preallocates
	if (!Lap_alloc) {
		Gu = new double[(nx - 1)*ny];
		Gv = new double[nx*(ny - 1)];
		Fx = new double[(nx - 1)*ny];
		Fy = new double[nx*(ny - 1)];
		Ru_inv = new double[(nx - 1)*ny];
		Rv_inv = new double[nx*(ny - 1)];
		Lap_alloc = 1;
	}
	
	SolverOperators::GradientOperator2D(p, Gu, Gv, dx, dy, nx, ny);
	IB_Operators::RegularizationOperator(Fx, Fy, BFx, BFy, B, dx, dy, nx, ny);

	//Updates u component of Gu to include immersed boundary forces
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < (nx - 1); i++) {
			Gu[j * (nx - 1) + i] = (Gu[j * (nx - 1) + i] - Fx[j * (nx - 1) + i]);
		}
	}

	//Updates v component of Gv to include immersed boundary forces
	for (int j = 0; j < (ny - 1); j++) {
		for (int i = 0; i < nx; i++) {
			//---v component----
			Gv[j * nx + i] = (Gv[j * nx + i] - Fy[j * nx + i]); //
		}
	}


	SolverOperators::R_Operator_Inv(Gu, Gv, Ru_inv, Rv_inv, Re, dx, dy, dt, nx, ny); //Applies R_inv as always

	SolverOperators::DivergenceOperator2D(Ru_inv, Rv_inv, res_p, dx, dy, nx, ny);//Divergence of G*Rinv
	IB_Operators::InterpolationOperator(Ru_inv, Rv_inv, res_uInterp, res_vInterp, B, dx, dy, nx, ny); //Interpolates velocity in Lagrangian grid
}

void Poisson::Preconditioner(double* rp, double* rfx, double* rfy, double* sp, double* sfx, double* sfy, double dx, double ds, int nx, int ny, Body* B) {
	//Implements the CG preconditioner to accelerate convergence.
	//Simple scaling r_p by dx^2 and r_f by ds. Returns in s_p and s_f

	double dx2 = dx * dx;

	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {
			sp[j * nx + i] = rp[j * nx + i] * dx2;
		}
	}
	for (int k = 0; k < B->ns; k++) {
		sfx[k] = rfx[k];
		sfy[k] = rfy[k];
	}

}