#include "IMatrix.h"



IMatrix::IMatrix(int N) : Matrix(N,N)
{
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (i == j) {
				matrix[i][j] = 1;
			}
			else {
				matrix[i][j] = 0;
			}
		}
	}
}


IMatrix::~IMatrix()
{
}
