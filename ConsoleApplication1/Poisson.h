#pragma once
#include <thread>
#include <vector>
#include <algorithm>
#include "Body.h"
#include "ImmersedBoundaryOperators.h"

class Poisson
{
public:
	static void CGModifiedPressurePoisson2D(double *pGuess, double *rp, double *pResult, Body* B, double *uInterp, double *vInterp, double dx, double dy, double dt, double Re, int nx, int ny, double reltol);
	static void Laplacian_Q(double *p, double *BFx, double *BFy, double *res_p, double *res_uInterp, double *res_vInterp, Body* B, double dx, double dy, double dt, double Re, int nx, int ny);
};

