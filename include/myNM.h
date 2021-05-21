/*----------------------------------------------------------------\
@ Numerical Methods by Young-Keun Kim - Handong Global University

Author           : Hyeonkyoung Lee
Created          : 02-04-2021
Modified         : 18-03-2021
Language/ver     : C++ in MSVS2019

Description      : myNM.h
----------------------------------------------------------------*/

#ifndef		_MY_NM_H		// use either (#pragma once) or  (#ifndef ...#endif)
#define		_MY_NM_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "myMatrix.h"

void tempFunc(int m);


// Return the dy/dx results for the input data. (truncation error: O(h^2))
Matrix	gradient(Matrix _x, Matrix _y);

void gradient_1Darray(double _x[], double _y[], double _dydx[], int m);

// Define a function that defines the target equation.
double myFunc(const double x);

double mydFunc(const double x);

Matrix	arr2Mat(double* _1Darray, int _rows, int _cols);

// Return the dy/dx results for the target equation. (truncation error: O(h^2))
Matrix	gradientFunc(double func(const double x), Matrix xin);

double newtonRaphsonFunc(double func(const double x), double dfunc(const double x), float x0, float tol);
#endif