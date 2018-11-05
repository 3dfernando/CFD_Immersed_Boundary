#include "stdafx.h"
#include "BC_2D.h"

//Stores boundary conditions
double* uTop;
double* vTop;
double* uBottom;
double* vBottom;
double* uLeft;
double* uRight;
double* vRight;
double* vLeft;

//Old boundary conditions for the convective terms
double* vTopOld;
double* vBottomOld;
double* uLeftOld;
double* uRightOld;
int OldBCInitialized = 0;

//Stores where the velocity vectors are and the grid spacings
double* u;
double* v;
double* uOld;
double* vOld;
double d_x;
double d_y;
double dt;
int nx; 
int ny;
//Stores which kind of BC and default values
int uTopKind = BC::BC_Dirichlet;
int vTopKind = BC::BC_Dirichlet;
int uBottomKind = BC::BC_Dirichlet;
int vBottomKind = BC::BC_Dirichlet;
int uLeftKind = BC::BC_Dirichlet;
int vLeftKind = BC::BC_Dirichlet;
int uRightKind = BC::BC_Outflow;
int vRightKind = BC::BC_Dirichlet;



//Class defines boundary conditions of the CFD code
BC::BC(int n_el_x, int n_el_y, double Delta_x, double Delta_y, double Delta_t, double* uField, double* vField, double* uFieldOld, double* vFieldOld)
{
	//Defines a constructor
	//n_el_x and n_el_y are the grid size
	//Delta_x, y, t are used for derivative computation
	//uField and vField are pointers to the velocity fields. This way convective BCs are easy to compute
	d_x = Delta_x;
	d_y = Delta_y;
	dt = Delta_t;

	u = uField;
	v = vField;

	uOld = uFieldOld;
	vOld = vFieldOld;

	nx = n_el_x;
	ny = n_el_y;
}

void BC::SetVector(int WhichVector, double* VectorValues)
{
	//Sets a BC vector.
	//Remembering: size(uTop)=nx-1, size(uBottom)=nx-1, size(uLeft)=ny, size(uRight)=ny,
	//size(vTop)=nx, size(vBottom)=nx, size(vLeft)=ny-1, size(vRight)=ny-1,
	//VectorValues is a pointer to the initialized value. The pointer is just copied to the class memory.

	if (WhichVector == VEC_U_TOP) {
		uTop = VectorValues;
	}
	else if (WhichVector == VEC_U_BOTTOM) {
		uBottom = VectorValues;
	}
	else if (WhichVector == VEC_U_LEFT) {
		uLeft = VectorValues;
		if (OldBCInitialized < 4) {
			//Inits old value for convective BC only once (this sizes the vector properly)
			uLeftOld = new double[ny];
			for (int i = 0; i < ny; i++) {
				uLeftOld[i] = VectorValues[i];
			}
			OldBCInitialized++;
		}
	}
	else if (WhichVector == VEC_U_RIGHT) {
		uRight = VectorValues;
		if (OldBCInitialized < 4) {
			//Inits old value for convective BC only once (this sizes the vector properly)
			uRightOld = new double[ny];
			for (int i = 0; i < ny; i++) {
				uRightOld[i] = VectorValues[i];
			}
			OldBCInitialized++;
		}
	}
	else if (WhichVector == VEC_V_TOP) {
		vTop = VectorValues;
		if (OldBCInitialized < 4) {
			//Inits old value for convective BC only once (this sizes the vector properly)
			vTopOld = new double[nx];
			for (int i = 0; i < nx; i++) {
				vTopOld[i] = VectorValues[i];
			}
			OldBCInitialized++;
		}
	}
	else if (WhichVector == VEC_V_BOTTOM) {
		vBottom = VectorValues;
		if (OldBCInitialized < 4) {
			//Inits old value for convective BC only once (this sizes the vector properly)
			vBottomOld = new double[nx];
			for (int i = 0; i < nx; i++) {
				vBottomOld[i] = VectorValues[i];
			}
			OldBCInitialized++;
		}
	}
	else if (WhichVector == VEC_V_LEFT) {
		vLeft = VectorValues;
	}
	else if (WhichVector == VEC_V_RIGHT) {
		vRight = VectorValues;
	}
}

void BC::SetKind(int WhichVector, int Kind)
{
	//Sets a BC vector kind.
	//Kind defines the boundary condition type, i.e. Dirichlet, convective, etc.
	if (WhichVector == VEC_U_TOP) {
		uTopKind = Kind;
	}
	else if (WhichVector == VEC_U_BOTTOM) {
		uBottomKind = Kind;
	}
	else if (WhichVector == VEC_U_LEFT) {
		uLeftKind = Kind;
	}
	else if (WhichVector == VEC_U_RIGHT) {
		uRightKind = Kind;
	}
	else if (WhichVector == VEC_V_TOP) {
		vTopKind = Kind;
	}
	else if (WhichVector == VEC_V_BOTTOM) {
		vBottomKind = Kind;
	}
	else if (WhichVector == VEC_V_LEFT) {
		vLeftKind = Kind;
	}
	else if (WhichVector == VEC_V_RIGHT) {
		vRightKind = Kind;
	}
}

double BC::GetVal(int WhichVector, int id_x)
{
	//Gets a BC vector location.
	//Remembering: size(uTop)=nx-1, size(uBottom)=nx-1, size(uLeft)=ny, size(uRight)=ny,
	//size(vTop)=nx, size(vBottom)=nx, size(vLeft)=ny-1, size(vRight)=ny-1,
	//Returns a double (only one number) for the index id_x
	if (WhichVector == VEC_U_TOP) {
		return uTop[id_x];
	}
	else if (WhichVector == VEC_U_BOTTOM) {
		return uBottom[id_x];
	}
	else if (WhichVector == VEC_U_LEFT) {
		return uLeft[id_x];
	}
	else if (WhichVector == VEC_U_RIGHT) {
		return uRight[id_x];
	}
	else if (WhichVector == VEC_V_TOP) {
		return vTop[id_x];
	}
	else if (WhichVector == VEC_V_BOTTOM) {
		return vBottom[id_x];
	}
	else if (WhichVector == VEC_V_LEFT) {
		return vLeft[id_x];
	}
	else if (WhichVector == VEC_V_RIGHT) {
		return vRight[id_x];
	}
	else {
		throw("You need to provide a vector");
	}
}

void BC::Laplacian(double* Lu, double* Lv) {
	//Returns the laplacian of the boundary conditions. Obviously it's only populated at the domain borders!

	double d2ud_x2, d2ud_y2, d2vd_x2, d2vd_y2;

	//Inits Lu and Lv (Copy this snippet where the call is made)
	//Lu = new double[(nx - 1) * ny];
	//Lv = new double[nx * (ny - 1)];
	
	//Evaluates Lu
	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < (nx - 1); i++) {
			//---u derivatives----
			//   d2ud_x2
			if (i == 0) {//uLeft
				d2ud_x2 = (uLeft[j]) / (d_x * d_x); //Boundary node
			}
			else if (i == (nx - 2)) {//uRight
				d2ud_x2 = (uRight[j]) / (d_x * d_x); //Boundary node
			}
			else {
				d2ud_x2 = 0; //Inside domain
			}
			//   d2ud_y2
			if (j == 0) {//uBottom
				d2ud_y2 = (uBottom[i]) / (d_y * d_y); //Boundary node
			}
			else if (j == (ny - 1)) {//uTop
				d2ud_y2 = (uTop[i]) / (d_y * d_y); //Boundary node
			}
			else {
				d2ud_y2 = 0; //Inside domain
			}

			//Evaluates Lu 
			Lu[j * (nx - 1) + i] = d2ud_x2 + d2ud_y2; 
		}
	}

	//Evaluate Lv
	for (int j = 0; j < (ny - 1); j++) {
		for (int i = 0; i < nx; i++) {
			//---v derivatives----
			//   d2vd_x2
			if (i == 0) {//vLeft
				d2vd_x2 = (vLeft[j]) / (d_x * d_x);
			}
			else if (i == (nx - 1)) {//vRight
				d2vd_x2 = (vRight[j]) / (d_x * d_x);
			}
			else {
				d2vd_x2 = 0;//Inside domain
			}
			//   d2vd_y2
			if (j == 0) {//vBottom
				d2vd_y2 = (vBottom[i]) / (d_y * d_y);
			}
			else if (j == (ny - 2)) {//vTop
				d2vd_y2 = (vTop[i]) / (d_y * d_y);
			}
			else {
				d2vd_y2 = 0; //Inside Domain
			}

			//Evaluates Lv
			Lv[j * nx + i] = d2vd_x2 + d2vd_y2; 
		}
	}
}

void BC::Divergence(double* D) {
	//Returns the divergence of the boundary conditions. Obviously it's only populated at the domain borders!

	double dudx, dvdy;

	for (int j = 0; j < ny; j++) {
		for (int i = 0; i < nx; i++) {
			//---u component----
			if (i == 0) {
				dudx = -uLeft[j] / d_x; //Missing BC
			}
			else if (i == (nx - 1)) {
				dudx = uRight[j] / d_x; //Missing BC
			}
			else {
				dudx = 0; //Interior point
			}

			//---v component----
			if (j == 0) {
				dvdy = -vBottom[i] / d_y; //Interior point
			}
			else if (j == (ny - 1)) {
				dvdy = vTop[i] / d_y; //Interior point
			}
			else {
				dvdy = 0; //Interior point
			}

			//Divergence
			D[j * nx + i] = dudx + dvdy;
		}
	}
	D[ny * nx - 1] = 0; //Pinned point at top-right corner. Won't be used, but sets to zero as a good measure! 
}

void BC::ComputeConvectiveBCs(){
	//Computes convective bcs and stores in the respective vectors.
	//dudt + Uavg * dud_x = 0;
	//u[n+1] =u[n] + dt * (1.5*g[n] - 0.5 * g[n-1]); where g = -uavg * dudx; This uses AB2 discretization.

	double AvgVel;
	double AvgVelOld;
	double ddx;
	double g;
	double gOld;


	if (uLeftKind == BC_Outflow) {
		//Computes average velocity based on old velocity field
		AvgVel = 0;
		AvgVelOld = 0;
		for (int j = 0; j < ny; j++) {
			AvgVel += u[j * (nx - 1) + 0];
			AvgVelOld += uOld[j * (nx - 1) + 0];
		}
		AvgVel = AvgVel / ny;
		AvgVelOld = AvgVelOld / ny;

		//Computes new boundary vector
		for (int j = 0; j < ny; j++) {
			g = -AvgVel * (u[j * (nx - 1) + 0] - uLeft[j]) / d_x;
			gOld = -AvgVelOld * (uOld[j * (nx - 1) + 0] - uLeftOld[j]) / d_x;

			uLeftOld[j] = uLeft[j]; //Updates old value
			uLeft[j] += dt * (1.5 * g - 0.5 * gOld); //Updates new value using AB2
		}
	}
	if (uRightKind == BC_Outflow) {
		AvgVel = 0;
		AvgVelOld = 0;
		for (int j = 0; j < ny; j++) {
			AvgVel += u[j * (nx - 1) + (nx - 2)];
			AvgVelOld += uOld[j * (nx - 1) + (nx - 2)];
		}
		AvgVel = AvgVel / ny;
		AvgVelOld = AvgVelOld / ny;

		//Computes new boundary vector
		for (int j = 0; j < ny; j++) {
			g = -AvgVel * (-u[j * (nx - 1) + (nx - 2)] + uRight[j]) / d_x;
			gOld = -AvgVelOld * (-uOld[j * (nx - 1) + (nx - 2)] + uRightOld[j]) / d_x;

			uRightOld[j] = uRight[j]; //Updates old value
			uRight[j] += dt * (1.5 * g - 0.5 * gOld); //Updates new value using AB2
		}
	}
	if (vTopKind == BC_Outflow) {
		AvgVel = 0;
		AvgVelOld = 0;
		for (int i = 0; i < nx; i++) {
			AvgVel += v[(ny - 2) * nx + i];
			AvgVelOld += vOld[(ny - 2) * nx + i];
		}
		AvgVel = AvgVel / nx;
		AvgVelOld = AvgVelOld / nx;

		//Computes new boundary vector
		for (int i = 0; i < nx; i++) {
			g = -AvgVel * (-v[(ny - 2) * nx + i] + vTop[i]) / d_y;
			gOld = -AvgVelOld * (-vOld[(ny - 2) * nx + i] + vTopOld[i]) / d_y;

			vTopOld[i] = vTop[i]; //Updates old value
			vTop[i] += dt * (1.5 * g - 0.5 * gOld); //Updates new value using AB2
		}
	}
	if (vBottomKind == BC_Outflow) {
		AvgVel = 0;
		AvgVelOld = 0;
		for (int i = 0; i < nx; i++) {
			AvgVel += v[0 * nx + i];
			AvgVelOld += vOld[0 * nx + i];
		}
		AvgVel = AvgVel / nx;
		AvgVelOld = AvgVelOld / nx;

		//Computes new boundary vector
		for (int i = 0; i < nx; i++) {
			g = -AvgVel * (-v[0 * nx + i] + vBottom[i]) / d_y;
			gOld = -AvgVelOld * (-vOld[0 * nx + i] + vBottomOld[i]) / d_y;

			vBottomOld[i] = vBottom[i]; //Updates old value
			vBottom[i] += dt * (1.5 * g - 0.5 * gOld); //Updates new value using AB2
		}
	}


	//The remaining 4 combinations aren't implemented.
	if (uTopKind == BC_Outflow) {
		throw "This boundary condition needs to be set as a Dirichlet zero.";
	}
	if (uBottomKind == BC_Outflow) {
		throw "This boundary condition needs to be set as a Dirichlet zero.";
	}
	if (vLeftKind == BC_Outflow) {
		throw "This boundary condition needs to be set as a Dirichlet zero.";
	}
	if (vRightKind == BC_Outflow) {
		throw "This boundary condition needs to be set as a Dirichlet zero.";
	}
}

BC::~BC()
{	
	//Destructor.
	delete[] uTop;
	delete[] vTop;
	delete[] uBottom;
	delete[] vBottom;
	delete[] uLeft;
	delete[] uRight;
	delete[] vRight;
	delete[] vLeft;
}
