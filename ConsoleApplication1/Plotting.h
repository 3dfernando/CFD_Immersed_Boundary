#pragma once
#include <windows.h>
#include <iostream>
#include <algorithm>
#include <string>
#include "TECIO.h"
#include "MASTER.h" 
#include "BC_2D.h"

class Plotting
{
public:
	static void PlotBMP(int X, int Y, HBITMAP B);
	static HBITMAP ArrayToBitmap(double *A, int W, int H);
	static void OutputTecplot(double *x, double *y, double *u, double *v, double *p, BC *B, int nx, int ny, double t, int Strand, bool EF);
	static void Print_CFL_Re(double Re, double dx, double dy, double dt, double t, double *u, double *v, int nx, int ny);
};

