#pragma once
#include <time.h>
#include <math.h> 
#include "Body.h"


class IB_Operators
{
public:
	static void InterpolationOperator(double *u, double *v, double *uInterp, double *vInterp, Body *B, double dx, double dy, int nx, int ny);
	static void RegularizationOperator(double *FuRegularized, double *FvRegularized, double *Fu, double *Fv, Body *B, double dx, double dy, int nx, int ny);
	static double DiscreteDeltaFunction2D(double x, double y, double dx, double dy);
};

