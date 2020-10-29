#pragma once
#include <stdio.h>
#include <iostream>
#include <math.h>

class Matrix
{
public:
	Matrix(int N, int M);
	//Coeficient matrix
	Matrix(double a1, double a2, double a3, int size);
	//result vector
	Matrix(double f, int size);
	//Copy constructor
	Matrix(const Matrix&);
	~Matrix();


	void printMatrix();
	void populate(int with);
	int getN();
	int getM();

	//Operators overloading 
	Matrix& operator = (const Matrix&);
	double* const operator [](int i);
	double* const operator [](int i) const;
	const Matrix operator * (const Matrix&) const;
	const Matrix operator * (const double) const;
	const Matrix operator + (const Matrix&) const;
	const Matrix operator - (const Matrix&) const;

protected:
	double** matrix;
	int N, M;
};

