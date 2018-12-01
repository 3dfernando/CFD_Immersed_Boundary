#include "stdafx.h"
#include "Body.h"



Body::Body(double* xPos, double* yPos, double* xVel, double* yVel, double* xForce, double* yForce, double delta_s, double n_s, double* Tx, double* Ty, double* Tt, double xC, double yC)
{
	x = xPos;
	y = yPos;
	vel_x = xVel;
	vel_y = yVel;
	Fx = xForce;
	Fy = yForce;
	ds = delta_s;
	ns = n_s;
	Trajectory_x = Tx;
	Trajectory_y = Ty;
	Trajectory_theta = Tt;
	xCenter = xC;
	yCenter = yC;
}

void Body::UpdateTrajectory(int t, double dt) {
	if (t > 0) {
		//Updates the position and velocity of the points in the body
		double dx = Trajectory_x[t] - Trajectory_x[t - 1];
		double dy = Trajectory_y[t] - Trajectory_y[t - 1];
		double dtheta = Trajectory_theta[t] - Trajectory_theta[t - 1];
		double cost = cos(dtheta);
		double sint = sin(dtheta);

		double xOld;
		double yOld;

		for (int i = 0; i < ns; i++) {
			xOld = x[i];
			yOld = y[i];

			//Applies rotation
			x[i] = xCenter + cost * (x[i] - xCenter) - sint * (y[i] - yCenter);
			y[i] = yCenter + sint * (x[i] - xCenter) + cost * (y[i] - yCenter);

			//Applies displacement
			x[i] += dx;
			y[i] += dy;

			//Gets point velocity
			vel_x[i] = (x[i] - xOld) / dt;
			vel_y[i] = (y[i] - yOld) / dt;
		}

		//Displaces center, too!
		xCenter += dx;
		yCenter += dy;
	}
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
