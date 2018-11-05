#pragma once
class BC
{
public:
	BC(int n_el_x, int n_el_y, double Delta_x, double Delta_y, double Delta_t, double* uField, double* vField, double* uFieldOld, double* vFieldOld);
	void SetVector(int WhichVector, double* VectorValues);
	double GetVal(int WhichVector, int idx);
	void SetKind(int WhichVector, int Kind);
	void Laplacian(double* Lu, double* Lv);
	void Divergence(double* D);
	void ComputeConvectiveBCs();
	~BC();

	//Constants for BC identification
	static const int BC_Dirichlet = 0;
	static const int BC_Outflow = 1;
	static const int VEC_U_TOP = 0;
	static const int VEC_U_BOTTOM = 1;
	static const int VEC_U_LEFT = 2;
	static const int VEC_U_RIGHT = 3;
	static const int VEC_V_TOP = 4;
	static const int VEC_V_BOTTOM = 5;
	static const int VEC_V_LEFT = 6;
	static const int VEC_V_RIGHT = 7;
};

