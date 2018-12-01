#pragma once
#include <math.h>

class Body
{
public:
	Body(double* xPos, double* yPos, double* xVel, double* yVel, double* xForce, double* yForce, double delta_s, double n_s, double* Tx, double* Ty, double* Tt, double xC, double yC);
	void UpdateTrajectory(int t, double dt);
	~Body();

	//Stores body coordinates and velocities
	double* x;
	double* y;
	double* vel_x;
	double* vel_y;
	double* Fx;
	double* Fy;
	double ds;
	double* Trajectory_x;
	double* Trajectory_y;
	double* Trajectory_theta;
	double xCenter;
	double yCenter;
	int ns;
};

