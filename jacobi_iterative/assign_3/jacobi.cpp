/**
 * @file    jacobi.cpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements matrix vector multiplication and Jacobi's method.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */
#include "jacobi.h"

/*
 * TODO: Implement your solutions here
 */

// my implementation:
#include <iostream>
#include <math.h>

// Calculates y = A*x for a square n-by-n matrix A, and n-dimensional vectors x
// and y
void matrix_vector_mult(const int n, const double* A, const double* x, double* y)
{
    // TODO
    double sum = 0;
    for (int i = 0; i < n; i++)
    {
    	sum = 0;
    	for (int j = 0; j < n; j++)
    	{
    		sum += A[i*n + j] * x[j];
    	}
    	y[i] = sum;
    }
}

// Calculates y = A*x for a n-by-m matrix A, a m-dimensional vector x
// and a n-dimensional vector y
void matrix_vector_mult(const int n, const int m, const double* A, const double* x, double* y)
{
    // TODO
    double sum = 0;
        for (int i = 0; i < n; i++)
    {
    	sum = 0;
    	for (int j = 0; j < m; j++)
    	{
    		sum += A[i*m + j] * x[j];
    	}
    	y[i] = sum;
    }
}

// implements the sequential jacobi method
void jacobi(const int n, double* A, double* b, double* x, int max_iter, double l2_termination)
{
    // initialize x with a 0 vector;
    // find D;
    
    for (int i = 0; i < n; i ++)
    {
    	x[i] = 0;
    }

	// find D = diag(A) and R = A - D;
	double D[n];
    for (int i = 0; i < n; i++)
    {
    	D[i] = A[i*n + i];
    }

    int counter = 0;
	double norm_diff;
	do
	{
		// calculate ||Ax - b||_2;
		norm_diff = 0;
		double y[n];
		matrix_vector_mult(n, A, x, &y[0]);

		for (int i = 0; i < n; i++)
		{
			norm_diff += pow(y[i] - b[i], 2);
		}


		// update x;

		for (int i = 0; i < n; i++)
		{
			x[i] = (b[i] - y[i])/D[i] + x[i];
		}

		counter += 1;

	}
    while (sqrt(norm_diff) > l2_termination && counter < max_iter);

    
}
