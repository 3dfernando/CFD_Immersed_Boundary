#include "stdafx.h"
#include "Poisson.h"
#include "Plotting.h"
#include "SolverOperators.h"
#include "ImmersedBoundaryOperators.h"
#include "BC_2D.h"
#include "Body.h"
#include "Setup_Interface.h"
#include <iostream>
#include <time.h>
#include <tchar.h>
#include <math.h>
#include <thread>
#include <ctime>
#include <cmath>
#include <direct.h> //Gets current path
#include <stdio.h>
#define GetCurrentDir _getcwd

using namespace std;

void sleep(int milliseconds);

#pragma region Variables
	double dx;
	double dy;
	const double pi = 3.14159;
#pragma endregion

BOOL ctrl_handler(DWORD event)
{
	if (event == CTRL_CLOSE_EVENT) {
		//Had some file saving routines here but they're not necessary anymore.
		return TRUE;
	}
	return FALSE;
}

int main()
{
	//Main routine does both initialization of the parameters for this simulation AND the time step integration.

	//======================================================================
	//=====================INITIALIZATION OF PARAMETERS=====================
	//======================================================================
	double* x; double*y;
	double dx; //x spacing of grid cells
	double dy; //y spacing of grid cells
	double dt; //time step increment
	double Re;  //Simulation Reynolds number

	double reltol; //Relative tolerance to stop CG algorithm

	int nx; //nx, ny are the number of grid points
	int ny; 
	int last_t; //t is the time step (not the true simulation time).
	int t_decimation; //Saves a field file only every t_decimation time steps

	//-------Boundary Conditions Variables--------
	double* uTop;	double* uBottom;
	double* uLeft;	double* uRight;
	double* vTop;	double* vBottom;
	double* vRight;	double* vLeft;

	int uLK, uRK, uBK, uTK;
	int vLK, vRK, vBK, vTK;

	//-----------Domain variables----------
	double* uOld; //Stores previous time step u, v, P
	double* vOld;
	double* POld;

	double* u; //Stores current time step u, v, P
	double* v;
	double* P;

	//================================================
	//=====Opens the simulation parameters file=======
	//================================================
	char cCurrentPath[FILENAME_MAX];
	GetCurrentDir(cCurrentPath, sizeof(cCurrentPath));
	std::string f(cCurrentPath);
	f = f + "\\Initialization.mat";
	Setup_Interface::GetSimulationParameters(f, &dx, &dy, &dt, &nx, &ny, &Re, &reltol, &last_t, &t_decimation, &uLeft, &uRight, &uTop, &uBottom, &vLeft, &vRight, &vTop, &vBottom, &u, &v, &P,
		&uLK, &uRK, &uTK, &uBK, &vLK, &vRK, &vTK, &vBK);

	//Initializes the old variables with the current values
	uOld = new double[(nx - 1) * ny];
	vOld = new double[nx * (ny - 1)];
	POld = new double[nx*ny];
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < (nx - 1); i++) {
			uOld[j*(nx - 1) + i] = u[j*(nx - 1) + i];
		}
	}
	for (int j = 0; j < (ny - 1); j++) {
		for (int i = 0; i < nx; i++) {
			vOld[j*nx + i] = v[j*nx + i];
		}
	}
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {
			POld[j*nx + i] = P[j*nx + i];
		}
	}

	//===================================================================
	//=====================Inits Boundary conditions=====================
	//===================================================================
	BC* Boundary;
	Boundary = new BC(nx, ny, dx, dy, dt, u, v, uOld, vOld);

	Boundary->SetVector(BC::VEC_U_TOP, uTop);
	Boundary->SetVector(BC::VEC_U_BOTTOM, uBottom);
	Boundary->SetVector(BC::VEC_U_RIGHT, uRight);
	Boundary->SetVector(BC::VEC_U_LEFT, uLeft);

	Boundary->SetKind(BC::VEC_U_LEFT, uLK);
	Boundary->SetKind(BC::VEC_U_RIGHT, uRK);
	Boundary->SetKind(BC::VEC_U_TOP, uTK);
	Boundary->SetKind(BC::VEC_U_BOTTOM, uBK);

	Boundary->SetVector(BC::VEC_V_TOP, vTop);
	Boundary->SetVector(BC::VEC_V_BOTTOM, vBottom);
	Boundary->SetVector(BC::VEC_V_RIGHT, vRight);
	Boundary->SetVector(BC::VEC_V_LEFT, vLeft);

	Boundary->SetKind(BC::VEC_V_TOP, vTK);
	Boundary->SetKind(BC::VEC_V_BOTTOM, vBK);
	Boundary->SetKind(BC::VEC_V_RIGHT, vRK);
	Boundary->SetKind(BC::VEC_V_LEFT, vLK);


	//===============================================
	//========Initializes the body in the IBM========
	//===============================================
	int ns;
	double ds;
	double* xBody; 
	double* yBody;
	double* vxBody;
	double* vyBody;
	double* xForceBody;
	double* yForceBody;
	double xCenter; //Center of rotation
	double yCenter;
	double* xTraj;
	double* yTraj;
	double* thetaTraj;

	std::string f2(cCurrentPath);
	f2 = f2 + "\\Body.mat";
	Setup_Interface::LoadBodyAndTrajectory(f2, &ns, &ds, &xBody, &yBody, &xCenter, &yCenter, &xTraj, &yTraj, &thetaTraj);

	//Inits forces and velocities to zero in the first time step
	vxBody = new double[ns];
	vyBody = new double[ns];
	xForceBody = new double[ns];
	yForceBody = new double[ns];
	for (int k = 0; k < ns; k++) {
		vxBody[k] = 0;
		vyBody[k] = 0;
		xForceBody[k] = 0;
		yForceBody[k] = 0;
	}

	Body* ImmersedBody = new Body(xBody, yBody, vxBody, vyBody, xForceBody, yForceBody, ds, ns, xTraj, yTraj, thetaTraj, xCenter, yCenter);


	//==========Defines x, y coordinates of the grid =================
	//================================================================
	x = new double[nx*ny]; //x coordinate of the cell centers
	y = new double[nx*ny]; //y coordinate of the cell centers
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < (nx + 1); i++) {
			x[j*nx + i] = ((double)i)*dx + (dx / 2); //Half cell added to refer to cell center 
			y[j*nx + i] = ((double)j)*dy + (dy / 2);
		}
	}

	

	//
	//=====================Other variables for the integrator=====================
	//
	HBITMAP B;
	HBITMAP Bu;
	HBITMAP Bv;
	HBITMAP Bp;
	double* uStar;
	double* vStar;
	double* Gu;
	double* Gv;
	double* D;
	double* Dbc;
	double* uInterp;
	double* vInterp;
	uStar = new double[(nx - 1)*ny];
	vStar = new double[nx*(ny - 1)];
	Gu = new double[nx*ny];
	Gv = new double[nx*ny];
	D = new double[nx*ny];
	Dbc = new double[nx*ny];
	uInterp = new double[ImmersedBody->ns];
	vInterp = new double[ImmersedBody->ns];


	//========================================================================
	//========================================================================
	//========================INTEGRATES IN TIME==============================
	//========================================================================
	//========================================================================
	for (int t = 0; t < last_t; t++) {
				
		Plotting::Print_CFL_Re(Re, dx, dy, dt, t, u, v, nx, ny); //Calculates maximum Courant and Reynolds and displays in Console window
		
		//Outputs the tecplot file
		if ((t % t_decimation == 0) && (t > 0)) {
			Plotting::OutputTecplot(x, y, u, v, P, Boundary, ImmersedBody, dx, dy, nx, ny, dt, (double)t * dt, t, false);
		}

		//Updates the u, v, P plots in console window
		Bu = Plotting::ArrayToBitmap(u, nx - 1, ny);
		Plotting::PlotBMP(450, 0, Bu);

		Bv = Plotting::ArrayToBitmap(v, nx, ny - 1);
		Plotting::PlotBMP(450, ny + 5, Bv);

		Bp = Plotting::ArrayToBitmap(P, nx, ny);
		Plotting::PlotBMP(450, ny * 2 + 10, Bp);
		
		//Solves Intermediary Velocity v*
		SolverOperators::SolveIntermVelocity2D(uOld, vOld, u, v, Boundary, uStar, vStar, Re, dx, dy, dt, nx, ny, reltol);

		//Calculates divergence of v*
		SolverOperators::DivergenceOperator2D(uStar, vStar, D, dx, dy, nx, ny);

		//Adds BC for divergence (Don't divide anymore by delta-t as that's done in the Q operator)
		Boundary->Divergence(Dbc);
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				D[j*nx + i] = D[j*nx + i] + Dbc[j*nx + i];
			}
		}

		//Interpolates uStar and vStar to the Lagrangian grid
		IB_Operators::InterpolationOperator(uStar, vStar, uInterp, vInterp, ImmersedBody, dx, dy, nx, ny);

		//Solves the Poisson equation
		//clock_t start = clock(); //Times the routine
		Poisson::CGModifiedPressurePoisson2D(POld, D, P, ImmersedBody, uInterp, vInterp, dx, dy, dt, Re, nx, ny, reltol);
		//clock_t end = clock();
		//double time = (double)(end - start) / CLOCKS_PER_SEC * 1000.0;
						
		//Copies P to Pold
		for (int i = 0; i < nx*ny; i++) { POld[i] = P[i]; } //POld could be removed and save time but potential savings are small

		//Calculates the gradient of P
		SolverOperators::ProjectionStep(uStar, vStar, u, v, P, ImmersedBody, Re, dx, dy, dt, nx, ny);
					
		//Moves u and v to uOld and vOld
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < (nx - 1); i++) {
				uOld[j*(nx - 1) + i] = u[j*(nx - 1) + i];
			}
		}
		for (int j = 0; j < (ny - 1); j++) {
			for (int i = 0; i < nx; i++) {
				vOld[j*nx + i] = v[j*nx + i];
			}
		}

		//Updates the body trajectory
		ImmersedBody->UpdateTrajectory(t, dt);

		//Clears the clutter in the console window.
		if (t % 10 == 0) {
			system("CLS");
		}
	}

	//========================================================================
	//========================================================================
	//========================================================================
	//========================================================================
	

	//After finishing all time steps.
    return 0;

}


