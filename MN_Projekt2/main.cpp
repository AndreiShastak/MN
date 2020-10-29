#include "Matrix.h"
#include "IMatrix.h"
#include "LinearEquationsSolver.h"
#include <iostream>
#include <math.h>

#define N_size 962



int main() {
	//Zadanie A
	//Matrix exA_A(10, -1, -1, N_size);
	//Matrix exA_b(5.0, N_size);


	//Zadanie B
	//Matrix resultJacobi = LinearEquationsSolver::Jacobi(exA_A, exA_b, pow(10, -9));
	//Matrix resultGaussSeidel = LinearEquationsSolver::Gauss_Seidel(exA_A, exA_b, pow(10, -9));


	//Zadanie C
	//Matrix exC_A(3, -1, -1, N_size);
	//Matrix resultJacobi2 = LinearEquationsSolver::Jacobi(exC_A, exA_b, pow(10, -9));
	//Matrix resultGaussSeidel2 = LinearEquationsSolver::Gauss_Seidel(exC_A, exA_b, pow(10, -9));


	//Zadanie D
	//Matrix exD_A(3, -1, -1, N_size);
	//Matrix exD_b(5.0, N_size);
	//Matrix resultLU = LinearEquationsSolver::LU_Factorization(exD_A, exD_b);
	//std::cout << LinearEquationsSolver::Norm(exD_A*resultLU - exD_b) << std::endl;


	//Zadanie E
	//LinearEquationsSolver::testMethods();


	return 0;
}