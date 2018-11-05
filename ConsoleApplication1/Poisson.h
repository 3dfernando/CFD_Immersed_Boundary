#pragma once
#include <thread>
#include <vector>
#include <algorithm>

class Poisson
{
public:
	static void CGPressurePoisson2D(void(*A_Fn) (double *p, double *res_p, double dx, double dy, double dt, double Re, int nx, int ny), double *pGuess, double *rp,
		double *pResult, double dx, double dy, double dt, double Re, int nx, int ny, double reltol);
	static void LaplacianPressure(double *p, double *res_p, double dx, double dy, double dt, double Re, int nx, int ny);
	static void SORPressurePoisson2D(double *pGuess, double *rp, double *pResult, double dx, double dy, double dt, double Re, int nx, int ny, double reltol);

};

