#pragma once
#include "mat.h"
#include <string>
#include <iostream>     // std::ios, std::istream, std::cout
#include <fstream>      // std::filebuf
#include <ctime>
#include <cstdlib> 

class Setup_Interface
{
public:
	static void Setup_Interface::GetSimulationParameters(std::string file, double* delta_x, double* delta_y, double* delta_t, int* n_x, int* n_y, double* Reynolds, double* rel_tol, int* lastTime, int* decimation, double** uLeft, double** uRight, double** uTop, double** uBottom, double** vLeft, double** vRight, double** vTop, double** vBottom, double** u, double** v, double** p, int* uLK, int* uRK, int* uTK, int* uBK, int* vLK, int* vRK, int* vTK, int* vBK);
	static void Setup_Interface::sleep(int milliseconds);
	static void Setup_Interface::LoadBodyAndTrajectory(std::string file, int* ns, double* ds, double** xBody, double** yBody, double* xCenter, double* yCenter, double** xTraj, double** yTraj, double** thetaTraj);
};

