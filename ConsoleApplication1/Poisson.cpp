#include "stdafx.h"
#include "Poisson.h"
#include "SolverOperators.h"

//Global variables for CG
double *d_p;
double *r_p;
double *pStarOld;
double *pTemp;
bool CG_alloc = 0;

//Global variables for Laplacian
double *Gu;
double *Gv;
double *Ru_inv;
double *Rv_inv;
bool Lap_alloc = 0;

void Poisson::CGPressurePoisson2D(void(*A_Fn) (double *p, double *res_p, double dx, double dy, double dt, double Re, int nx, int ny), double *pGuess, double *rp,
	double *pResult, double dx, double dy, double dt, double Re, int nx, int ny, double reltol) {
	//===========================================================================================
	//===========================================================================================
	//===============================CONJUGATE GRADIENT METHOD===================================
	//===========================================================================================
	//===========================================================================================
	//Solves the LP=r equation using the amazing Conjugate Gradient Method
	//A_Fn is a pointer to a function that is evaluated (matrix A)
	//uGuess and vGuess are the initial guesses to start the algorithm
	//r is a pointer to the right-hand-side vector
	//B is a boundary condition object
	//uResult and vResult are the results output
	//Re, dx, dy, nx, ny, reltol are self-evident.

	//Further comments: This function could be the same for P and uStar (intermediate vel gradient). I separated the functions to consider the array sizing requirements. 

	//Preallocates once
	if (!CG_alloc) {
		d_p = new double[nx*ny];
		pStarOld = new double[nx*ny];
		r_p = new double[nx*ny];
		pTemp = new double[nx*ny];
		CG_alloc = 1;
	}
	

	int iter = 0;
	double alpha, alphaDen;
	double beta;
	double delta, delta0, deltaOld;
	
	//Define initial residue
	delta = 0;
	A_Fn(pGuess, pTemp, dx, dy, dt, Re, nx, ny);
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {
			r_p[j * nx + i] = rp[j * nx + i] - pTemp[j * nx + i];
			d_p[j * nx + i] = r_p[j * nx + i];
			delta += r_p[j * nx + i] * r_p[j * nx + i];
			pResult[j * nx + i] = pGuess[j * nx + i]; //Inits result
		}
	}


	delta0 = delta;

	while (delta > (reltol * reltol * delta0)) {
		iter++;

		//Updates alpha
		alphaDen = 0;
		A_Fn(d_p, pTemp, dx, dy, dt, Re, nx, ny);
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				alphaDen += d_p[j * nx + i] * pTemp[j * nx + i];
			}
		}
		alpha = delta / alphaDen;

		//Updates pStar
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				pResult[j * nx + i] += alpha * d_p[j * nx + i];
			}
		}
		//pResult[nx * ny - 1] = 0; //Enforces the pinned point (not needed because D(nx*ny+1)=0 but just making sure)

		//Quick formula to update the residuals
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				r_p[j * nx + i] = r_p[j * nx + i] - alpha * pTemp[j * nx + i];
			}
		}

		//Updates delta
		deltaOld = delta;
		delta = 0;
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				delta += r_p[j * nx + i] * r_p[j * nx + i];
			}
		}


		//Updates Beta
		beta = delta / deltaOld;
		//Updates d
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				d_p[j * nx + i] = r_p[j * nx + i] + beta * d_p[j * nx + i];
			}
		}
	}	

}

void Poisson::LaplacianPressure(double *p, double *res_p, double dx, double dy, double dt, double Re, int nx, int ny) {
	//Generates the Laplace operator from the pressure field divergence

	//preallocates
	if (!Lap_alloc) {
		Gu = new double[(nx - 1)*ny];
		Gv = new double[nx*(ny - 1)];
		Ru_inv = new double[(nx - 1)*ny];
		Rv_inv = new double[nx*(ny - 1)];
		Lap_alloc = 1;
	}
	
	SolverOperators::GradientOperator2D(p, Gu, Gv, dx, dy, nx, ny);
	SolverOperators::R_Operator_Inv(Gu, Gv, Ru_inv, Rv_inv, Re, dx, dy, dt, nx, ny);

	SolverOperators::DivergenceOperator2D(Ru_inv, Rv_inv, res_p, dx, dy, nx, ny);//Divergence of G*Rinv
	//SolverOperators::DivergenceOperator2D(Gu, Gv, res_p, dx, dy, nx, ny);//Just the Laplacian (for debug purposes)



}

void Poisson::SORPressurePoisson2D(double *pGuess, double *rp, double *pResult, double dx, double dy, double dt, double Re, int nx, int ny, double reltol) {
	//Solves the Poisson equations with a successive over-relaxation method.
	//O(150x) slower than CG so it's not being used!
	//Implementation does not include Rinv.
	
	//Even though it is not active, it was still very useful during debug.

	double dx2 = dx * dx;
	double dy2 = dy * dy;
	double pLeft, pRight, pUp, pDn, pMid;
	double Beta = 1.6; //SOR coefficient 
	double changeMag;
	double absMag;
	double relChange = 1e6;
	double pOld;

	int i, j;

	for (int i = 0; i < nx*ny; i++) { pResult[i] = pGuess[i]; }//Inits pResult
	

	int iter;

	iter = 0;
	while (relChange > reltol) { //Stops if eiher condition is satisfied
		iter++;

		changeMag = 0; //Zeros relative change
		absMag = 0; //Zeros abs magnitude

		//Uses Gauss-Seidel + SOR
		for (j = 1; j < (ny - 1); j++) {
			for (i = 1; i < (nx - 1); i++) {
				//Interior points
						
				pLeft = pResult[j * nx + (i - 1)];
				pRight = pResult[j * nx + (i + 1)];
				pUp = pResult[(j + 1) * nx + i];
				pDn = pResult[(j - 1) * nx + i];

				pOld = pResult[j * nx + i];
				pResult[j * nx + i] = Beta * ((((pLeft + pRight) * dy2) + ((pUp + pDn) * dx2) - (dx2 * dy2 * rp[j * nx + i])) / (2 * (dx2 + dy2))) + (1 - Beta) * pOld;
				
				changeMag += (pResult[j * nx + i] - pOld) * (pResult[j * nx + i] - pOld);
				absMag += pResult[j * nx + i] * pResult[j * nx + i];
			}
		}


		i = 0;
		for (j = 1; j < (ny - 1); j++) {
			//Left boundary

			pLeft = 0;
			pRight = pResult[j * nx + (i + 1)];
			pUp = pResult[(j + 1) * nx + i];
			pDn = pResult[(j - 1) * nx + i];

			pOld = pResult[j * nx + i];
			pResult[j * nx + i] = Beta * ((((pLeft + pRight) * dy2) + ((pUp + pDn) * dx2) - (dx2 * dy2 * rp[j * nx + i])) / (dx2 + 2 * dy2)) + (1 - Beta) * pOld; //Note denominator is not 2*(dx2 + dy2), that's because of BC

			changeMag += (pResult[j * nx + i] - pOld) * (pResult[j * nx + i] - pOld);
			absMag += pResult[j * nx + i] * pResult[j * nx + i];
		}
		i = nx - 1;
		for (j = 1; j < (ny - 1); j++) {
			//Right boundary

			pLeft = pResult[j * nx + (i - 1)];
			pRight = 0;
			pUp = pResult[(j + 1) * nx + i];
			pDn = pResult[(j - 1) * nx + i];

			pOld = pResult[j * nx + i];
			pResult[j * nx + i] = Beta * ((((pLeft + pRight) * dy2) + ((pUp + pDn) * dx2) - (dx2 * dy2 * rp[j * nx + i])) / (dx2 + 2 * dy2)) + (1 - Beta) * pOld;

			changeMag += (pResult[j * nx + i] - pOld) * (pResult[j * nx + i] - pOld);
			absMag += pResult[j * nx + i] * pResult[j * nx + i];
		}
		j = 0;
		for (i = 1; i < (nx - 1); i++) {
			//Bottom boundary

			pLeft = pResult[j * nx + (i - 1)];
			pRight = pResult[j * nx + (i + 1)];
			pUp = pResult[(j + 1) * nx + i];
			pDn = 0;

			pOld = pResult[j * nx + i];
			pResult[j * nx + i] = Beta * ((((pLeft + pRight) * dy2) + ((pUp + pDn) * dx2) - (dx2 * dy2 * rp[j * nx + i])) / (2 * dx2 + dy2)) + (1 - Beta) * pOld;

			changeMag += (pResult[j * nx + i] - pOld) * (pResult[j * nx + i] - pOld);
			absMag += pResult[j * nx + i] * pResult[j * nx + i];
		}
		j = ny - 1;
		for (i = 1; i < (nx - 1); i++) {
			//Top boundary

			pLeft = pResult[j * nx + (i - 1)];
			pRight = pResult[j * nx + (i + 1)];
			pUp = 0;
			pDn = pResult[(j - 1) * nx + i];

			pOld = pResult[j * nx + i];
			pResult[j * nx + i] = Beta * ((((pLeft + pRight) * dy2) + ((pUp + pDn) * dx2) - (dx2 * dy2 * rp[j * nx + i])) / (2 * dx2 + dy2)) + (1 - Beta) * pOld;

			changeMag += (pResult[j * nx + i] - pOld) * (pResult[j * nx + i] - pOld);
			absMag += pResult[j * nx + i] * pResult[j * nx + i];
		}

		//Bottom-Left boundary corner
		i = 0; j = 0;
		pLeft = 0;
		pRight = pResult[j * nx + (i + 1)];
		pUp = pResult[(j + 1) * nx + i];
		pDn = 0;

		pOld = pResult[j * nx + i];
		pResult[j * nx + i] = Beta * ((((pLeft + pRight) * dy2) + ((pUp + pDn) * dx2) - (dx2 * dy2 * rp[j * nx + i])) / (dx2 + dy2)) + (1 - Beta) * pOld;

		changeMag += (pResult[j * nx + i] - pOld) * (pResult[j * nx + i] - pOld);
		absMag += pResult[j * nx + i] * pResult[j * nx + i];


		//Bottom-Right boundary corner
		i = nx - 1; j = 0;
		pLeft = pResult[j * nx + (i - 1)];
		pRight = 0;
		pUp = pResult[(j + 1) * nx + i];
		pDn = 0;

		pOld = pResult[j * nx + i];
		pResult[j * nx + i] = Beta * ((((pLeft + pRight) * dy2) + ((pUp + pDn) * dx2) - (dx2 * dy2 * rp[j * nx + i])) / (dx2 + dy2)) + (1 - Beta) * pOld;

		changeMag += (pResult[j * nx + i] - pOld) * (pResult[j * nx + i] - pOld);
		absMag += pResult[j * nx + i] * pResult[j * nx + i];

		//Top-Left boundary corner
		i = 0; j = ny - 1;
		pLeft = 0;
		pRight = pResult[j * nx + (i + 1)];
		pUp = 0;
		pDn = pResult[(j - 1) * nx + i];

		pOld = pResult[j * nx + i];
		pResult[j * nx + i] = Beta * ((((pLeft + pRight) * dy2) + ((pUp + pDn) * dx2) - (dx2 * dy2 * rp[j * nx + i])) / (dx2 + dy2)) + (1 - Beta) * pOld;

		changeMag += (pResult[j * nx + i] - pOld) * (pResult[j * nx + i] - pOld);
		absMag += pResult[j * nx + i] * pResult[j * nx + i];

		//Top-Right boundary corner
		i = nx - 1; j = ny - 1;
		/*pLeft = pResult[j * nx + (i - 1)];
		pRight = 0;
		pUp = 0;
		pDn = pResult[(j - 1) * nx + i];

		pOld = pResult[j * nx + i];
		pResult[j * nx + i] = Beta * ((((pLeft + pRight) * dy2) + ((pUp + pDn) * dx2) - (dx2 * dy2 * rp[j * nx + i])) / (dx2 + dy2)) + (1 - Beta) * pOld;
		
		changeMag += (pResult[j * nx + i] - pOld) * (pResult[j * nx + i] - pOld);
		absMag += pResult[j * nx + i] * pResult[j * nx + i];*/

		//Pinned point
		pResult[j * nx + i] = 0; //Pinned point has zero pressure


		//relative change
		relChange = sqrt(changeMag / absMag);

	}

	//Outputs the number of iterations
	std::cout << "Poisson iter: " << iter << std::endl;
}
