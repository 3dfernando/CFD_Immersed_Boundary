#include "stdafx.h"
#include "Body.h"



Body::Body(double* xPos, double* yPos, double* xVel, double* yVel, double* xForce, double* yForce, double delta_s, double n_s)
{
	x = xPos;
	y = yPos;
	vel_x = xVel;
	vel_y = yVel;
	Fx = xForce;
	Fy = yForce;
	ds = delta_s;
	ns = n_s;
}


Body::~Body()
{
	delete[] x;
	delete[] y;
	delete[] Fx;
	delete[] Fy;
	delete[] vel_x;
	delete[] vel_y;
}
