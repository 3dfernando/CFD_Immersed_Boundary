#include "stdafx.h"
#include "Plotting.h"
//Global variables 
BITMAP bm;

void Plotting::PlotBMP(int X, int Y, HBITMAP B)
{
	//Plots a bitmap in the screen at X,Y coordinates from the command console

	//Gets the bitmap from the handle
	GetObject(B, sizeof(BITMAP), &bm);

	HDC hdc = CreateCompatibleDC(NULL);

	SelectObject(hdc, B);

	HDC hdc_x = GetDC(GetConsoleWindow());

	BitBlt(hdc_x, X, Y, bm.bmWidth, bm.bmHeight, hdc, 0, 0, SRCCOPY);
	ReleaseDC(HWND_DESKTOP, hdc_x);

}

HBITMAP Plotting::ArrayToBitmap(double *A, int W, int H)
{
	//Transfers the array to a HBITMAP Windows API object such that it can be drawn in the screen.
	double* maxA;
	double* minA;
	maxA = std::max_element(A, A + W*H);
	minA = std::min_element(A, A + W*H);

	double Avg;
	double Diff;
	double Level;

	Avg = (*maxA + *minA) / 2;
	Diff = (*maxA - *minA) / 2;

	//Makes array width be a multiple of 4(otherwise bitmap plots improperly for some weird Windows reason)
	int newW = ceil((double)W / 4) * 4;

	BYTE* RGB;
	int nC = 3; //Number of channels
	int idx, idx2;
	RGB = new BYTE[newW * H * nC];

	for (int j = 0; j < H; j++) {
		for (int i = 0; i < newW; i++) {
			idx = j * newW + i;
			idx2 = j * W + i;
			if (i < W) {				
				if (A[idx2] > Avg) {
					Level = 255 * (1 - ((A[idx2] - Avg) / Diff));
					RGB[nC * idx] = (BYTE)Level;
					RGB[nC * idx + 1] = (BYTE)Level;
					RGB[nC * idx + 2] = 255;
				}
				else {
					Level = 255 * (1 - ((Avg - A[idx2]) / Diff));
					RGB[nC * idx] = 255;
					RGB[nC * idx + 1] = (BYTE)Level;
					RGB[nC * idx + 2] = (BYTE)Level;
				}
			}
			else {
				RGB[nC * idx] = 0;
				RGB[nC * idx + 1] = 0;
				RGB[nC * idx + 2] = 0;
			}
		}
	}

	
	BITMAPINFOHEADER bmih;
	bmih.biSize = sizeof(BITMAPINFOHEADER);
	bmih.biWidth = newW;
	bmih.biHeight = -H;
	bmih.biPlanes = 1;
	bmih.biBitCount = 24;
	bmih.biCompression = BI_RGB;
	bmih.biSizeImage = 0;
	bmih.biXPelsPerMeter = 10;
	bmih.biYPelsPerMeter = 10;
	bmih.biClrUsed = 0;
	bmih.biClrImportant = 0;

	BITMAPINFO dbmi;
	ZeroMemory(&dbmi, sizeof(dbmi));
	dbmi.bmiHeader = bmih;
	dbmi.bmiColors->rgbBlue = 0;
	dbmi.bmiColors->rgbGreen = 0;
	dbmi.bmiColors->rgbRed = 0;
	dbmi.bmiColors->rgbReserved = 0;

	HDC hdc = ::GetDC(NULL);

	HBITMAP hbmp = CreateDIBitmap(hdc, &bmih, CBM_INIT, RGB, &dbmi, DIB_RGB_COLORS);
	//Cleanup
	delete[] RGB;

	return hbmp;
}

void Plotting::OutputTecplot(double *x, double *y, double *u, double *v, double *p, BC *B, Body* BD, double dx, double dy, int nx, int ny, double dt, double t, int Strand, bool EF) {
	//Outputs the result to tecplot as a *.plt file
	INTEGER4 Debug = 1;
	INTEGER4 VIsDouble = 1;
	INTEGER4 FileType = 0;
	INTEGER4 FileFormat = 1; // 0 == PLT, 1 == SZPLT
	INTEGER4 I = 0; /* Used to track return codes */

	char fName_c[30];
	sprintf(fName_c, "R_%d.szplt", Strand);
	
	I = TECINI142((char*)"CFD Results", (char*)"X Y U V P Fx Fy", fName_c, (char*)".", &FileFormat, &FileType, &Debug, &VIsDouble);

	//Inits all the variables to write the proper format
	INTEGER4 ICellMax = 0;
	INTEGER4 JCellMax = 0;
	INTEGER4 KCellMax = 0;
	INTEGER4 DIsDouble = 0;
	double SolTime = t;
	INTEGER4 StrandID = (INTEGER4)Strand;
	INTEGER4 ParentZn = 0;
	INTEGER4 IsBlock = 1; /* Block */
	INTEGER4 NFConns = 0;
	INTEGER4 FNMode = 0;
	INTEGER4 TotalNumFaceNodes = 1;
	INTEGER4 TotalNumBndryFaces = 1;
	INTEGER4 TotalNumBndryConnections = 1;
	INTEGER4 ShrConn = 0;
	/*Ordered Zone Parameters*/
	INTEGER4 IMax = nx;
	INTEGER4 JMax = ny;
	INTEGER4 KMax = 1;

	//Computes the integral of drag forces:
	double TotalFx, TotalFy;
	TotalFx = 0; TotalFy = 0;
	for (int k = 0; k < BD->ns; k++) {
		TotalFx += BD->Fx[k];
		TotalFy += BD->Fy[k];
	}
	TotalFx = -TotalFx * dx * dy / dt;
	TotalFy = -TotalFy * dx * dy / dt;
	std::cout << "Fx=" << TotalFx;
	std::cout << "Fy=" << TotalFy << std::endl;

	//Saves Tecplot File
	float* X1;
	float* Y1;
	float* U1;
	float* V1;
	float* P1;
	float* Fx1;
	float* Fy1;
	double* Fx2;
	double* Fy2;
	X1 = new float[nx*ny];
	Y1 = new float[nx*ny];
	U1 = new float[nx*ny];
	V1 = new float[nx*ny];
	P1 = new float[nx*ny];
	Fx1 = new float[nx*ny];
	Fy1 = new float[nx*ny];
	Fx2 = new double[nx * ny];
	Fy2 = new double[nx * ny];

	//Calls the regularization operator and produces the body force projected in the fluid
	IB_Operators::RegularizationOperator(Fx2, Fy2, BD->Fx, BD->Fy, BD, dx, dy, nx, ny);



	for (int i = 0; i < nx*ny; i++) {
		X1[i] = (float)x[i]; //Typecasts variables
		Y1[i] = (float)y[i];
	}

	//Performs interpolation of the u, v fields to the center of the cell
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {
			P1[j*nx + i] = (float)p[j*nx + i];
		}
	}
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {
			if (i == 0) {
				U1[j*nx + i] = (float)((B->GetVal(BC::VEC_U_LEFT, j) + u[j*(nx - 1) + i]) / 2); //BC
				Fx1[j*nx + i] = (float)(Fx2[j*(nx - 1) + i] / 2); //Force at boundary is 0
			}
			else if (i == (nx - 1)) {
				U1[j*nx + i] = (float)((u[j*(nx - 1) + (i - 1)] + B->GetVal(BC::VEC_U_RIGHT, j)) / 2); //BC
				Fx1[j*nx + i] = (float)(Fx2[j*(nx - 1) + (i - 1)] / 2); //Force at boundary is 0
			}
			else {
				U1[j*nx + i] = (float)((u[j*(nx - 1) + (i - 1)] + u[j*(nx - 1) + i]) / 2); //Interpolates to get U
				Fx1[j*nx + i] = (float)((Fx2[j*(nx - 1) + (i - 1)] + Fx2[j*(nx - 1) + i]) / 2); //Interpolates to get Fx
			}
			Fx1[j*nx + i] = -Fx1[j*nx + i] * dx * dy / dt; //Corrects scale factors beta/alpha (H=E^T) and -dt (from lambda)
		}
	}
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {
			if (j == 0) {
				V1[j*nx + i] = (float)((B->GetVal(BC::VEC_V_BOTTOM, i) + v[j*nx + i]) / 2); //BC
				Fy1[j*nx + i] = (float)(Fy2[j*nx + i] / 2); //Force at boundary is 0
			}
			else if (j == (ny - 1)) {
				V1[j*nx + i] = (float)((v[(j - 1)*nx + i] + B->GetVal(BC::VEC_V_TOP, i)) / 2); //BC
				Fy1[j*nx + i] = (float)(Fy2[(j - 1)*nx + i] / 2); //Force at boundary is 0
			}
			else {
				V1[j*nx + i] = (float)((v[(j - 1)*nx + i] + v[j*nx + i]) / 2); //Interpolates to get V
				Fy1[j*nx + i] = (float)((Fy2[(j - 1)*nx + i] + Fy2[j*nx + i]) / 2); //Interpolates to get Fy
			}
			Fy1[j*nx + i] = -Fy1[j*nx + i] * dx * dy / dt; //Corrects scale factors beta/alpha (H=E^T) and -dt (from lambda)
		}
	}
	

	/* Ordered Zone */
	INTEGER4 ZoneType = 0;
	I = TECZNE142((char*) "Flow Field",
		&ZoneType,
		&IMax,
		&JMax,
		&KMax,
		&ICellMax,
		&JCellMax,
		&KCellMax,
		&SolTime,
		&StrandID,
		&ParentZn,
		&IsBlock,
		&NFConns,
		&FNMode,
		&TotalNumFaceNodes,
		&TotalNumBndryFaces,
		&TotalNumBndryConnections,
		NULL,
		NULL,
		NULL,
		&ShrConn);
	INTEGER4 III = IMax * JMax * KMax;
	I = TECDAT142(&III, X1, &DIsDouble);
	I = TECDAT142(&III, Y1, &DIsDouble);
	I = TECDAT142(&III, U1, &DIsDouble);
	I = TECDAT142(&III, V1, &DIsDouble);
	I = TECDAT142(&III, P1, &DIsDouble);
	I = TECDAT142(&III, Fx1, &DIsDouble);
	I = TECDAT142(&III, Fy1, &DIsDouble);


	/* FE Zone for body */



	I = TECEND142(); //Ends file

	delete[] X1;
	delete[] Y1;
	delete[] U1;
	delete[] V1;
	delete[] P1;
	delete[] Fx1;
	delete[] Fy1;
	delete[] Fx2;
	delete[] Fy2;
}

void Plotting::Print_CFL_Re(double Re, double dx, double dy, double dt, double t, double *u, double *v, int nx, int ny) {
	//Just prints the CFL and cell Re in the screen.
	double CFL = 0;
	double CellRe = 0;

	CFL = 0; CellRe = 0;
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < (nx - 1); i++) {
			if (abs(u[j*(nx - 1) + i] * dt / dx) > CFL) { CFL = abs(u[j*(nx - 1) + i] * dt / dx); }
			if (abs(u[j*(nx - 1) + i] * dx * Re) > CellRe) { CellRe = abs(u[j*(nx - 1) + i] * dx * Re); }
		}
	}
	for (int j = 0; j < (ny - 1); j++) {
		for (int i = 0; i < nx; i++) {
			if (abs(v[j*nx + i] * dt / dy) > CFL) { CFL = abs(v[j*nx + i] * dt / dy); }
			if (abs(v[j*nx + i] * dy * Re) > CellRe) { CellRe = abs(v[j*nx + i] * dy * Re); }
		}
	}
	std::cout << "CFL=" << CFL << "; CellRe=" << CellRe;
	std::cout << "iter " << t << std::endl;
}