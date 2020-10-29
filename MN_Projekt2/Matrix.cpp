#include "Matrix.h"


Matrix::Matrix(int N, int M)
{
	this->N = N;
	this->M = M;
	matrix = new double*[N];
	for (int i = 0; i < N; i++) {
		matrix[i] = new double[M];
	}
}

Matrix::Matrix(double a1, double a2, double a3, int size)
	:Matrix(size,size)
{
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			if (i == j)matrix[i][j] = a1;
			else if (j == i + 1 || j == i - 1) matrix[i][j] = a2;
			else if (j == i + 2 || j == i - 2) matrix[i][j] = a3;
			else matrix[i][j] = 0;
		}
	}
}

Matrix::Matrix(double f, int size)
	:Matrix(size, 1)
{
	for (int i = 0; i < size; i++) {
		matrix[i][0] = sin(i * (f + 1));
	}
}

Matrix::Matrix(const Matrix & copy)
{

	this->N = copy.N;
	this->M = copy.M;
	matrix = new double*[N];
	for (int i = 0; i < N; i++) {
		matrix[i] = new double[M];
	}

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			matrix[i][j] = copy[i][j];
			//matrix[i][j] = copy.getElement(i, j);
		}
	}

}

Matrix::~Matrix()
{
	for (int i = 0; i < N; i++) {
		delete[] matrix[i];
	}
	delete[]matrix;
	matrix = NULL;
}


void Matrix::printMatrix()
{
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			std::cout << matrix[i][j] << " ";
		}
		std::cout << std::endl;
	}
}

void Matrix::populate(int with)
{
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			matrix[i][j] = with;
		}
	}
}

int Matrix::getN()
{
	return N;
}

int Matrix::getM()
{
	return M;
}

Matrix & Matrix::operator=(const Matrix & right)
{
	if (this == &right)
		return *this; //self assignment
	if (matrix != NULL)
		this->~Matrix(); //clean up already allocated memory

	this->N = right.N;
	this->M = right.M;

	matrix = new double*[N];
	for (int i = 0; i < N; i++) {
		matrix[i] = new double[M];
	}

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			matrix[i][j] = right[i][j];
		}
	}

	return *this;
}

double* const  Matrix::operator[](int i)
{
	return matrix[i];
}

double * const Matrix::operator[](int i) const
{
	return matrix[i];
}

const Matrix Matrix::operator*(const Matrix & right) const
{
	Matrix c(this->N,right.M);
	double sum_elems;
	if (this->M == right.N) {
		for (int i = 0; i < this->N; ++i)
		{
			for (int j = 0; j < right.M; ++j)
			{
				sum_elems = 0;
				for (int k = 0; k < right.N; ++k)
				{
					sum_elems += (*this)[i][k] * right[k][j];
				}
				c[i][j] = sum_elems;
			}
		}
	}
	return c;
}

const Matrix Matrix::operator*(const double scalar) const
{
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			matrix[i][j] *= scalar;
		}
	}
	return *this;
}

const Matrix Matrix::operator+(const Matrix & right) const
{
	Matrix c = right;
	if (this->N == right.N && this->M == right.M) {
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < M; j++) {
				c[i][j] = (*this)[i][j] + right[i][j];
			}
		}
	}
	return c;
}

const Matrix Matrix::operator-(const Matrix & right) const
{
	Matrix c = right;
	c = c * (-1);
	return c + (*this);
}
