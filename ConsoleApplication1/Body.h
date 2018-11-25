#pragma once
class Body
{
public:
	Body(double* xPos, double* yPos, double* xVel, double* yVel, double* xForce, double* yForce, double delta_s, double n_s);
	~Body();

	//Stores body coordinates and velocities
	double* x;
	double* y;
	double* vel_x;
	double* vel_y;
	double* Fx;
	double* Fy;
	double ds;
	int ns;
};

