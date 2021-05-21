/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : Hyeonkyoung Lee
Created          : 11-05-2021
Modified         : 18-03-2021
Language/ver     : C++ in MSVS2019

Description      : myNM.cpp
----------------------------------------------------------------*/

#include "myNM.h"
// Return the dy/dx results for the input data. (truncation error: O(h^2))
// Move this function to myNM.h and myNM.cpp
Matrix	gradient(Matrix _x, Matrix _y) {
	int n = _x.rows;
	int ny = _y.rows;
	
	if (n != ny) {
		printf("ERROR: length of x and y must be equal");
		return { 0 };
	}

	Matrix Out = zeros(_x.rows, 1);

	//Assuming constant h
	double h = _x.at[1][0] - _x.at[0][0];


	for (int k = 0; k < _x.rows ; k++) {
		if (_x.rows > 2) {
			if (k == 0) {
				Out.at[k][0] = (-3 * _y.at[0][0] + 4 * _y.at[1][0] - _y.at[2][0]) / (2 * h);
			}
			else if (k == _x.rows-1) {
				Out.at[k][0] = (_y.at[k - 2][0] - 4 * _y.at[k - 1][0] + 3 * _y.at[k][0]) / (2 * h);
			}
			else {
				Out.at[k][0] = (_y.at[k + 1][0] - _y.at[k - 1][0]) / (2 * h);
			}
		}
		else if(_x.rows==2){
			Out.at[0][0] = (_y.at[1][0] - _y.at[0][0]) / h;
			Out.at[1][0] = (_y.at[1][0] - _y.at[0][0]) / h;
		}
	}
	return Out;
}

Matrix	arr2Mat(double* _1Darray, int _rows, int _cols)
{
	Matrix Output = createMat(_rows, _cols);

	for (int i = 0; i < _rows; i++)
		for (int j = 0; j < _cols; j++)
			Output.at[i][j] = _1Darray[i * _cols + j];

	return Output;
}

void gradient_1Darray(double _x[], double _y[],  double _dydx[], int m) {
	
	double h = _x[1] - _x[0];

	for (int k = 0; k < m; k++) {
		if (m > 2) {
			if (k == 0) {
				_dydx[k] = (-3 * _y[0] + 4 * _y[1] - _y[2]) / (2 * h);
			}
			else if (k == m - 1) {
				_dydx[k] = (_y[k - 2] - 4 * _y[k - 1] + 3 * _y[k]) / (2 * h);
			}
			else {
				_dydx[k] = (_y[k + 1] - _y[k - 1]) / (2 * h);
			}
		}
		else if (m == 2) {
			_dydx[0] = (_y[1] - _y[0]) / h;
			_dydx[1] = (_y[1] - _y[0]) / h;
		}
	}
	printf("dydx =\n");
	for (int i = 0; i < m; i++) {
		printf("%15.6f\n", _dydx[i]);
	}
}

// Define a function that defines the target equation.
double myFunc(const double x) {
	return 3 * x * x;
}
double mydFunc(const double x) {
	return 6 * x;
}

// Return the dy/dx results for the target equation. (truncation error: O(h^2))
// Move this function to myNM.h and myNM.cpp
Matrix	gradientFunc(double func(const double x), Matrix xin) {
	int n = xin.rows;
	Matrix y = zeros(n, 1);

	for (int i = 0; i < n; i++) {
		y.at[i][0] = func(xin.at[i][0]);
	}
	return gradient(xin,y);
}

double newtonRaphsonFunc(double func(const double x), double dfunc(const double x), float x0, float tol) {
	int k = 0;
	int Nmax = 1000;
	float h;
	float xn = x0;
	float ep = fabs(func(xn));

	do {
		printf("  Iteration:%d \t", k);
		printf("X(n): %f \t", xn);
		printf("Tolerance: %.10f\n", ep);
		if (dfunc(xn) == 0)
		{
			EXIT_FAILURE;
		}
		h = -func(xn) / dfunc(xn);
		xn += h;
		ep = fabs(func(xn));
		k++;
	} while (k<Nmax && ep>tol);

	printf("  Iteration:%d \t", k);
	printf("X(n): %f \t", xn);
	printf("Tolerance: %.10f\n", ep);

	return xn;
}