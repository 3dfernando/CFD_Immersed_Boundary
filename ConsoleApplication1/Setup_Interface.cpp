#include "stdafx.h"
#include "Setup_Interface.h"

using namespace std;

void Setup_Interface::GetSimulationParameters(std::string file, double* delta_x, double* delta_y, double* delta_t, int* n_x, int* n_y, double* Reynolds, 
			double* rel_tol, int* lastTime, int* decimation, double** uLeft, double** uRight, double** uTop, double** uBottom,
	double** vLeft, double** vRight, double** vTop, double** vBottom, double** u, double** v, double** p, int* uLK, int* uRK, int* uTK, int* uBK, int* vLK, int* vRK, int* vTK, int* vBK) {
	//Opens the file that has the simulation parameters (in csv format) and parses the variables out

	char* f;
	int n = file.length();
	f = new char[n+1];
	strcpy(f, file.c_str());

	MATFile *pmat;
	mxArray *scalar;

	pmat = matOpen(f, "r");
	if (pmat == NULL) {
		printf("Error reopening file %s\n", file);
		Setup_Interface::sleep(10000);
		exit(0);
	}

	//Copies the initialization data to the output variables


	//================GRID DEFINITION================
	double temp;
	int intTemp;

	scalar = matGetVariable(pmat, "dx");
	if (scalar == NULL) {
		printf("Error reading dx\n");
		Setup_Interface::sleep(10000);
		exit(0);
	}
	temp = mxGetScalar(scalar);
	memcpy(delta_x, &temp, sizeof(double));

	scalar = matGetVariable(pmat, "dy");
	if (scalar == NULL) {
		printf("Error reading dy\n");
		Setup_Interface::sleep(10000);
		exit(0);
	}
	temp = mxGetScalar(scalar);
	memcpy(delta_y, &temp, sizeof(double));
	
	scalar = matGetVariable(pmat, "dt");
	if (scalar == NULL) {
		printf("Error reading dt\n");
		Setup_Interface::sleep(10000);
		exit(0);
	}
	temp = mxGetScalar(scalar);
	memcpy(delta_t, &temp, sizeof(double));

	scalar = matGetVariable(pmat, "nx");
	if (scalar == NULL) {
		printf("Error reading nx\n");
		Setup_Interface::sleep(10000);
		exit(0);
	}
	intTemp = (int) mxGetScalar(scalar);
	memcpy(n_x, &intTemp, sizeof(int));

	scalar = matGetVariable(pmat, "ny");
	if (scalar == NULL) {
		printf("Error reading ny\n");
		Setup_Interface::sleep(10000);
		exit(0);
	}
	intTemp = (int)mxGetScalar(scalar);
	memcpy(n_y, &intTemp, sizeof(int));

	scalar = matGetVariable(pmat, "Re");
	if (scalar == NULL) {
		printf("Error reading Re\n");
		Setup_Interface::sleep(10000);
		exit(0);
	}
	temp = mxGetScalar(scalar);
	memcpy(Reynolds, &temp, sizeof(double));

	scalar = matGetVariable(pmat, "reltol");
	if (scalar == NULL) {
		printf("Error reading reltol\n");
		Setup_Interface::sleep(10000);
		exit(0);
	}
	temp = mxGetScalar(scalar);
	memcpy(rel_tol, &temp, sizeof(double));

	scalar = matGetVariable(pmat, "last_t");
	if (scalar == NULL) {
		printf("Error reading last_t\n");
		Setup_Interface::sleep(10000);
		exit(0);
	}
	intTemp = (int)mxGetScalar(scalar);
	memcpy(lastTime, &intTemp, sizeof(int));

	scalar = matGetVariable(pmat, "t_decimation");
	if (scalar == NULL) {
		printf("Error reading t_decimation\n");
		Setup_Interface::sleep(10000);
		exit(0);
	}
	intTemp = (int)mxGetScalar(scalar);
	memcpy(decimation, &intTemp, sizeof(int));


	//================BOUNDARY CONDITION INITIALIZATION================
	mxArray *matArray;
	//size_t nEl;

	double* tempArray1;
	matArray = matGetVariable(pmat, "uLeft");
	if (matArray == NULL) {
		printf("Error reading uLeft\n");
		Setup_Interface::sleep(10000);
		exit(0);
	}
	tempArray1 = mxGetPr(matArray);
	*uLeft = tempArray1;

	double* tempArray2;
	matArray = matGetVariable(pmat, "uRight");
	if (matArray == NULL) {
		printf("Error reading uRight\n");
		Setup_Interface::sleep(10000);
		exit(0);
	}
	tempArray2 = mxGetPr(matArray);
	*uRight = tempArray2;

	double* tempArray3;
	matArray = matGetVariable(pmat, "uTop");
	if (matArray == NULL) {
		printf("Error reading uTop\n");
		Setup_Interface::sleep(10000);
		exit(0);
	}
	tempArray3 = mxGetPr(matArray);
	*uTop = tempArray3;

	double* tempArray4;
	matArray = matGetVariable(pmat, "uBottom");
	if (matArray == NULL) {
		printf("Error reading uBottom\n");
		Setup_Interface::sleep(10000);
		exit(0);
	}
	tempArray4 = mxGetPr(matArray);
	*uBottom = tempArray4;

	double* tempArray5;
	matArray = matGetVariable(pmat, "vLeft");
	if (matArray == NULL) {
		printf("Error reading vLeft\n");
		Setup_Interface::sleep(10000);
		exit(0);
	}
	tempArray5 = mxGetPr(matArray);
	*vLeft = tempArray5;

	double* tempArray6;
	matArray = matGetVariable(pmat, "vRight");
	if (matArray == NULL) {
		printf("Error reading vRight\n");
		Setup_Interface::sleep(10000);
		exit(0);
	}
	tempArray6 = mxGetPr(matArray);
	*vRight = tempArray6;

	double* tempArray7;
	matArray = matGetVariable(pmat, "vTop");
	if (matArray == NULL) {
		printf("Error reading vTop\n");
		Setup_Interface::sleep(10000);
		exit(0);
	}
	tempArray7 = mxGetPr(matArray);
	*vTop = tempArray7;

	double* tempArray8;
	matArray = matGetVariable(pmat, "vBottom");
	if (matArray == NULL) {
		printf("Error reading vBottom\n");
		Setup_Interface::sleep(10000);
		exit(0);
	}
	tempArray8 = mxGetPr(matArray);
	*vBottom = tempArray8;

	//================BOUNDARY CONDITION KIND INITIALIZATION================
	scalar = matGetVariable(pmat, "uLeftKind");
	if (scalar == NULL) {
		printf("Error reading uLeftKind\n");
		Setup_Interface::sleep(10000);
		exit(0);
	}
	intTemp = (int)mxGetScalar(scalar);
	memcpy(uLK, &intTemp, sizeof(int));

	scalar = matGetVariable(pmat, "uRightKind");
	if (scalar == NULL) {
		printf("Error reading uRightKind\n");
		Setup_Interface::sleep(10000);
		exit(0);
	}
	intTemp = (int)mxGetScalar(scalar);
	memcpy(uRK, &intTemp, sizeof(int));

	scalar = matGetVariable(pmat, "uTopKind");
	if (scalar == NULL) {
		printf("Error reading uTopKind\n");
		Setup_Interface::sleep(10000);
		exit(0);
	}
	intTemp = (int)mxGetScalar(scalar);
	memcpy(uTK, &intTemp, sizeof(int));

	scalar = matGetVariable(pmat, "uBottomKind");
	if (scalar == NULL) {
		printf("Error reading uBottomKind\n");
		Setup_Interface::sleep(10000);
		exit(0);
	}
	intTemp = (int)mxGetScalar(scalar);
	memcpy(uBK, &intTemp, sizeof(int));

	scalar = matGetVariable(pmat, "vLeftKind");
	if (scalar == NULL) {
		printf("Error reading vLeftKind\n");
		Setup_Interface::sleep(10000);
		exit(0);
	}
	intTemp = (int)mxGetScalar(scalar);
	memcpy(vLK, &intTemp, sizeof(int));

	scalar = matGetVariable(pmat, "vRightKind");
	if (scalar == NULL) {
		printf("Error reading vRightKind\n");
		Setup_Interface::sleep(10000);
		exit(0);
	}
	intTemp = (int)mxGetScalar(scalar);
	memcpy(vRK, &intTemp, sizeof(int));

	scalar = matGetVariable(pmat, "vTopKind");
	if (scalar == NULL) {
		printf("Error reading vTopKind\n");
		Setup_Interface::sleep(10000);
		exit(0);
	}
	intTemp = (int)mxGetScalar(scalar);
	memcpy(vTK, &intTemp, sizeof(int));

	scalar = matGetVariable(pmat, "vBottomKind");
	if (scalar == NULL) {
		printf("Error reading vBottomKind\n");
		Setup_Interface::sleep(10000);
		exit(0);
	}
	intTemp = (int)mxGetScalar(scalar);
	memcpy(vBK, &intTemp, sizeof(int));



	//================DOMAIN INITIALIZATION================
	double* tempArray10;
	matArray = matGetVariable(pmat, "u");
	if (matArray == NULL) {
		printf("Error reading u\n");
		Setup_Interface::sleep(10000);
		exit(0);
	}
	tempArray10 = mxGetPr(matArray);
	*u = tempArray10;

	double* tempArray11;
	matArray = matGetVariable(pmat, "v");
	if (matArray == NULL) {
		printf("Error reading v\n");
		Setup_Interface::sleep(10000);
		exit(0);
	}
	tempArray11 = mxGetPr(matArray);
	*v = tempArray11;

	double* tempArray12;
	matArray = matGetVariable(pmat, "p");
	if (matArray == NULL) {
		printf("Error reading p\n");
		Setup_Interface::sleep(10000);
		exit(0);
	}
	tempArray12 = mxGetPr(matArray);
	*p = tempArray12;

}

void Setup_Interface::LoadBodyAndTrajectory(std::string file, int* ns, double* ds, double** xBody, double** yBody, double* xCenter, double* yCenter, double** xTraj, double** yTraj, double** thetaTraj) {
	//Loads the body points and the trajectory as a function of time

	char* f;
	int n = file.length();
	f = new char[n + 1];
	strcpy(f, file.c_str());

	MATFile *pmat;
	mxArray *scalar;
	mxArray *matArray;

	pmat = matOpen(f, "r");
	if (pmat == NULL) {
		printf("Error reopening file %s\n", file);
		Setup_Interface::sleep(10000);
		exit(0);
	}



	//Loads variables
	double temp;

	scalar = matGetVariable(pmat, "xC");
	if (scalar == NULL) {
		printf("Error reading xC\n");
		Setup_Interface::sleep(10000);
		exit(0);
	}
	temp = mxGetScalar(scalar);
	memcpy(xCenter, &temp, sizeof(double));

	scalar = matGetVariable(pmat, "yC");
	if (scalar == NULL) {
		printf("Error reading yC\n");
		Setup_Interface::sleep(10000);
		exit(0);
	}
	temp = mxGetScalar(scalar);
	memcpy(yCenter, &temp, sizeof(double));

	scalar = matGetVariable(pmat, "ds");
	if (scalar == NULL) {
		printf("Error reading ds\n");
		Setup_Interface::sleep(10000);
		exit(0);
	}
	temp = mxGetScalar(scalar);
	memcpy(ds, &temp, sizeof(double));

	//Loads matrices

	double* tempArray1;
	matArray = matGetVariable(pmat, "xP");
	if (matArray == NULL) {
		printf("Error reading xP\n");
		Setup_Interface::sleep(10000);
		exit(0);
	}
	tempArray1 = mxGetPr(matArray);
	*xBody = tempArray1;

	double* tempArray2;
	matArray = matGetVariable(pmat, "yP");
	if (matArray == NULL) {
		printf("Error reading yP\n");
		Setup_Interface::sleep(10000);
		exit(0);
	}
	tempArray2 = mxGetPr(matArray);
	*yBody = tempArray2;

	double* tempArray3;
	matArray = matGetVariable(pmat, "xTraj");
	if (matArray == NULL) {
		printf("Error reading xTraj\n");
		Setup_Interface::sleep(10000);
		exit(0);
	}
	tempArray3 = mxGetPr(matArray);
	*xTraj = tempArray3;

	double* tempArray4;
	matArray = matGetVariable(pmat, "yTraj");
	if (matArray == NULL) {
		printf("Error reading yTraj\n");
		Setup_Interface::sleep(10000);
		exit(0);
	}
	tempArray4 = mxGetPr(matArray);
	*yTraj = tempArray4;

	double* tempArray5;
	matArray = matGetVariable(pmat, "thetaTraj");
	if (matArray == NULL) {
		printf("Error reading thetaTraj\n");
		Setup_Interface::sleep(10000);
		exit(0);
	}
	tempArray5 = mxGetPr(matArray);
	*thetaTraj = tempArray5;

	int intTemp;
	scalar = matGetVariable(pmat, "ns");
	if (scalar == NULL) {
		printf("Error reading ns\n");
		Setup_Interface::sleep(10000);
		exit(0);
	}
	intTemp = (int)mxGetScalar(scalar);
	memcpy(ns, &intTemp, sizeof(int));

}

void Setup_Interface::sleep(int milliseconds)
{
	clock_t time_end;
	time_end = clock() + milliseconds * CLOCKS_PER_SEC / 1000;
	while (clock() < time_end)
	{
	}
}