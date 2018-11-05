#include "stdafx.h"
#include "Poisson.h"
#include "Plotting.h"
#include "SolverOperators.h"
#include "BC_2D.h"
#include <iostream>
#include <time.h>
#include <tchar.h>
#include <math.h>
#include <thread>
#include <ctime>
#include <cmath>

using namespace std;

void sleep(int milliseconds);

#pragma region Variables
	double dx;
	double dy;
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
	double dx = 0.05; //x spacing of grid cells
	double dy = 0.05; //y spacing of grid cells
	double dt = 0.01; //time step increment
	double Re = 200;  //Simulation Reynolds number

	double reltol = 1e-10; //Relative tolerance to stop CG algorithm

	int nx = 100; //nx, ny are the number of grid points
	int ny = 99; //nx different than ny to spot swapping mistakes
	int last_t = 200000; //t is the time step (not the true simulation time).
	int t_decimation = 10; //Saves a field file only every t_decimation time steps

	SetConsoleCtrlHandler((PHANDLER_ROUTINE)(ctrl_handler), TRUE); //Kills console window gracefully

	x = new double[nx*ny]; //x coordinate of the cell centers
	y = new double[nx*ny]; //y coordinate of the cell centers
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < (nx + 1); i++) {
			x[j*nx + i] = ((double)i)*dx + (dx / 2); //Half cell added to refer to cell center 
			y[j*nx + i] = ((double)j)*dy + (dy / 2);
		}
	}

	double* uOld; //Stores previous time step u, v, P
	double* vOld;
	double* POld;

	double* u; //Stores current time step u, v, P
	double* v;
	double* P;

	u = new double[(nx - 1) * ny];
	v = new double[nx * (ny - 1)];
	uOld = new double[(nx - 1) * ny];
	vOld = new double[nx * (ny - 1)];
	P = new double[nx*ny];
	POld = new double[nx*ny];
	
	//=====================Inits u and v =====================
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < (nx - 1); i++) {
			u[j*(nx - 1) + i] = 0;
			uOld[j*(nx - 1) + i] = u[j*(nx - 1) + i];
		}
	}

	for (int j = 0; j < (ny - 1); j++) {
		for (int i = 0; i < nx; i++) {
			v[j*nx + i] = 0;
			vOld[j*nx + i] = v[j*nx + i];
		}
	}
	//=====================Inits P =====================
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {
			P[j*nx + i] = 0;
			POld[i] = P[j*nx + i];
		}
	}

	//===================================================================
	//=====================Inits Boundary conditions=====================
	//===================================================================
	BC* Boundary;
	Boundary = new BC(nx, ny, dx, dy, dt, u, v, uOld, vOld);

	double* uTop;	double* uBottom;
	uTop = new double[nx - 1]; uBottom = new double[nx - 1];
	for (int i = 0; i < (nx - 1); i++) { uTop[i] = 0; uBottom[i] = 1;}
	Boundary->SetVector(BC::VEC_U_TOP, uTop);
	Boundary->SetVector(BC::VEC_U_BOTTOM, uBottom);

	double* uLeft;	double* uRight;
	uLeft = new double[ny]; uRight = new double[ny];
	for (int i = 0; i < ny; i++) {
		if (i < (ny / 3)) {
			uLeft[i] = 0;
			uRight[i] = 0;
		}
		else if (i < (2 * ny / 3)) {
			uLeft[i] = 0;
			uRight[i] = 0;
		}
		else {
			uLeft[i] = 0;
			uRight[i] = 0;
		}
	}
	Boundary->SetVector(BC::VEC_U_RIGHT, uRight);
	Boundary->SetVector(BC::VEC_U_LEFT, uLeft);
	Boundary->SetKind(BC::VEC_U_LEFT,BC::BC_Dirichlet);
	Boundary->SetKind(BC::VEC_U_RIGHT, BC::BC_Dirichlet);
	
	double* vTop;	double* vBottom;
	vTop = new double[nx]; vBottom = new double[nx];
	for (int i = 0; i < nx; i++) {
		if (i < (nx / 3)) {
			vTop[i] = 0;
			vBottom[i] = 0;
		}
		else if (i < (2 * nx / 3)) {
			vTop[i] = 0;
			vBottom[i] = 0;
		}
		else {
			vTop[i] = 0;
			vBottom[i] = 0;
		}
	}
	Boundary->SetVector(BC::VEC_V_TOP, vTop);
	Boundary->SetVector(BC::VEC_V_BOTTOM, vBottom);
	Boundary->SetKind(BC::VEC_V_TOP, BC::BC_Dirichlet);
	Boundary->SetKind(BC::VEC_V_BOTTOM, BC::BC_Dirichlet);

	double* vRight;	double* vLeft;
	vRight = new double[ny - 1]; vLeft = new double[ny - 1];
	for (int i = 0; i < (ny - 1); i++) { vRight[i] = 0; vLeft[i] = 0; }
	Boundary->SetVector(BC::VEC_V_RIGHT, vRight);
	Boundary->SetVector(BC::VEC_V_LEFT, vLeft);

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
	uStar = new double[(nx - 1)*ny];
	vStar = new double[nx*(ny - 1)];
	Gu = new double[nx*ny];
	Gv = new double[nx*ny];
	D = new double[nx*ny];
	Dbc = new double[nx*ny];


	//========================================================================
	//========================================================================
	//========================INTEGRATES IN TIME==============================
	//========================================================================
	//========================================================================
	for (int t = 0; t < last_t; t++) {
				
		Plotting::Print_CFL_Re(Re, dx, dy, dt, t, u, v, nx, ny); //Calculates maximum Courant and Reynolds and displays in Console window
		
		//Outputs the tecplot file
		if (t % t_decimation == 0) {
			Plotting::OutputTecplot(x, y, u, v, P, Boundary, nx, ny, (double)t * dt, t, false);
		}

		//Updates the u, v, P plots in console window
		Bu = Plotting::ArrayToBitmap(u, nx - 1, ny);
		Plotting::PlotBMP(400, 0, Bu);

		Bv = Plotting::ArrayToBitmap(v, nx, ny - 1);
		Plotting::PlotBMP(400, ny + 5, Bv);

		Bp = Plotting::ArrayToBitmap(P, nx, ny);
		Plotting::PlotBMP(400, ny * 2 + 10, Bp);
		
		//Solves Intermediary Velocity v*
		SolverOperators::SolveIntermVelocity2D(uOld, vOld, u, v, Boundary, uStar, vStar, Re, dx, dy, dt, nx, ny, reltol);

		//Calculates divergence of v*
		SolverOperators::DivergenceOperator2D(uStar, vStar, D, dx, dy, nx, ny);

		//Adds BC for divergence and divides by delta t
		Boundary->Divergence(Dbc);
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				D[j*nx + i] = (D[j*nx + i] + Dbc[j*nx + i]) / dt;
			}
		}

		//Solves the Poisson equation
		//clock_t start = clock(); //Times the routine
		Poisson::CGPressurePoisson2D(Poisson::LaplacianPressure, POld, D, P, dx, dy, dt, Re, nx, ny, reltol);
		//clock_t end = clock();
		//double time = (double)(end - start) / CLOCKS_PER_SEC * 1000.0;
						
		//Copies P to Pold
		for (int i = 0; i < nx*ny; i++) { POld[i] = P[i]; } //POld could be removed and save time but potential savings are small

		//Calculates the gradient of P
		SolverOperators::ProjectionStep(uStar, vStar, u, v, P, Re, dx, dy, dt, nx, ny);
					

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



#pragma region RandomRoutines
	void sleep(int milliseconds)
	{
		clock_t time_end;
		time_end = clock() + milliseconds * CLOCKS_PER_SEC / 1000;
		while (clock() < time_end)
		{
		}
	}
#pragma endregion

