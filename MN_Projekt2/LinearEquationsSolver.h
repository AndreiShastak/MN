#pragma once
#include "Matrix.h"
#include "IMatrix.h"
#include <math.h>

static class LinearEquationsSolver
{
public:
	LinearEquationsSolver();
	~LinearEquationsSolver();

	static Matrix Jacobi(Matrix, Matrix, double);
	static Matrix Gauss_Seidel(Matrix, Matrix, double);
	static Matrix LU_Factorization(Matrix, Matrix);
	static void testMethods();
	

	static double Norm(Matrix);

private:
	static void testHelppingMethod(int N);
};

