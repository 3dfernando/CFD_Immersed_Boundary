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

	BYTE* RGB;
	int nC = 3;
	RGB = new BYTE[W * H * nC];

	for (int i = 0; i<(W*H); i++) {
		if (A[i]>Avg) {
			Level = 255 * (1 - ((A[i] - Avg) / Diff));
			RGB[nC * i] = (BYTE)Level;
			RGB[nC * i + 1] = (BYTE)Level;
			RGB[nC * i + 2] = 255;
		}
		else {
			Level = 255 * (1 - ((Avg - A[i]) / Diff));
			RGB[nC * i] = 255;
			RGB[nC * i + 1] = (BYTE)Level;
			RGB[nC * i + 2] = (BYTE)Level;
		}
	}

	
	BITMAPINFOHEADER bmih;
	bmih.biSize = sizeof(BITMAPINFOHEADER);
	bmih.biWidth = W;
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

void Plotting::OutputTecplot(double *x, double *y, double *u, double *v, double *p, BC *B, int nx, int ny, double t, int Strand, bool EF) {
	//Outputs the result to tecplot as a *.plt file
	INTEGER4 Debug = 1;
	INTEGER4 VIsDouble = 1;
	INTEGER4 FileType = 0;
	INTEGER4 FileFormat = 1; // 0 == PLT, 1 == SZPLT
	INTEGER4 I = 0; /* Used to track return codes */

	char fName_c[30];
	sprintf(fName_c, "R_%d.szplt", Strand);
	
	I = TECINI142((char*)"CFD Results", (char*)"X Y U V P", fName_c, (char*)".", &FileFormat, &FileType, &Debug, &VIsDouble);

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


	float* X1;
	float* Y1;
	float* U1;
	float* V1;
	float* P1;
	X1 = new float[nx*ny];
	Y1 = new float[nx*ny];
	U1 = new float[nx*ny];
	V1 = new float[nx*ny];
	P1 = new float[nx*ny];
	for (int i = 0; i < nx*ny; i++) {
		X1[i] = (float)x[i]; //Typecasts variables
		Y1[i] = (float)y[i];
	}

	//Performs interpolation of the u, v fields to the center of the cell
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {
			P1[j*nx + i] = p[j*nx + i];
		}
	}
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {
			if (i == 0) {
				U1[j*nx + i] = (B->GetVal(BC::VEC_U_LEFT, j) + u[j*(nx - 1) + i]) / 2; //BC
			}
			else if (i == (nx - 1)) {
				U1[j*nx + i] = (u[j*(nx - 1) + (i - 1)] + B->GetVal(BC::VEC_U_RIGHT, j)) / 2; //BC
			}
			else {
				U1[j*nx + i] = (u[j*(nx - 1) + (i - 1)] + u[j*(nx - 1) + i]) / 2; //Interpolates to get U
			}
		}
	}
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {
			if (j == 0) {
				V1[j*nx + i] = (B->GetVal(BC::VEC_V_BOTTOM, i) + v[j*nx + i]) / 2; //BC
			}
			else if (j == (ny - 1)) {
				V1[j*nx + i] = (v[(j - 1)*nx + i] + B->GetVal(BC::VEC_V_TOP, i)) / 2; //BC
			}
			else {
				V1[j*nx + i] = (v[(j - 1)*nx + i] + v[j*nx + i]) / 2; //Interpolates to get V
			}
		}
	}
	

	/* Ordered Zone */
	INTEGER4 ZoneType = 0;
	I = TECZNE142((char*) "Ordered Zone",
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

	I = TECEND142(); //Ends file

	delete[] X1;
	delete[] Y1;
	delete[] P1;
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