#pragma once
#include "Plotting.h"
#include "BC_2D.h"
#include "ImmersedBoundaryOperators.h"
#include "Body.h"
#include <time.h>

class SolverOperators
{
public:

	static void SolveIntermVelocity2D(double *uPrev, double *vPrev, double *uCurrent, double *vCurrent, BC *B, double *uStar, double *vStar, double Re, double dx, double dy, double dt, int nx, int ny, double reltol);
	static void ConvectiveOperator2D(double *u, double *v, double *Hu, double *Hv, BC *B, double dx, double dy, int nx, int ny);
	static void LaplaceOperator2D(double *u, double *v, double *Lu, double *Lv, double dx, double dy, int nx, int ny);
	static void GradientOperator2D(double *P, double *Gu, double *Gv, double dx, double dy, int nx, int ny);
	static void DivergenceOperator2D(double *u, double *v, double *D, double dx, double dy, int nx, int ny); 
	static void CGMomentum2D(void(*A_Fn) (double *u, double *v, double *res_u, double *res_v, double Re, double dx, double dy, double dt, int nx, int ny), double *uGuess, double *vGuess, double *ru, double *rv,
		double *uResult, double *vResult, double Re, double dx, double dy, double dt, int nx, int ny, double reltol);
	static void R_Operator(double *u, double *v, double *res_u, double *res_v, double Re, double dx, double dy, double dt, int nx, int ny);
	static void R_Operator_Inv(double *u, double *v, double *res_u, double *res_v, double Re, double dx, double dy, double dt, int nx, int ny);
	static void ProjectionStep(double *uStar, double *vStar, double *uNext, double *vNext, double *P, Body *B, double Re, double dx, double dy, double dt, int nx, int ny);
};

