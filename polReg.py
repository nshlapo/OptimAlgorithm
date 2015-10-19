from __future__ import division
import numpy as np
import sympy as sp
import string

points = [(1, 1), (1, 0), (2, 2), (12, 5432)]

coeffVals = [1, 1, 1, 1]

def calcError(points, coeffVals):
	'''
	Calculates error by taking sum of error between each 
	point and the polynomial function

	Inputs:
		points (list of tuples): Points to fit
		coeffVals (list): Coefficients of polynomial function
	Returns:
		totalError (int): Sum of errors squared
	'''
	totalError = 0
	for point in points:
		xVal = point[0]
		yVal = point[1]
		funcVal = sum([coeff*(xVal**index) for index, coeff in enumerate(coeffVals)])
		totalError += (yVal - funcVal)**2

	return totalError

def calcDeriv(points, coeffVals, indexCoeff):
	'''
	Calculates derivative of function with respect to given 
	coefficient index at each given error function location

	Inputs: 
		points (list of tuples): Points to fit
		coeffVals (list): Coefficients of polynomial function
		indexCoeff (int): Term of function for which derivative 
		is taken respect to
	Returns: 
		finalDeriv (int): Value of derivative of function with
		respect to the given index
	'''
	for point in points:
		xVal = point[0]
		yVal = point[1]
		xDepDeriv = sum([coeff*(xVal**(index+indexCoeff)) for index, coeff in enumerate(coeffVals)])
		deriv = yVal*(xVal**(indexCoeff)) - xDepDeriv

	finalDeriv = -2*deriv
	return finalDeriv

# In final function, iterate over degrees until satisfactory fit is reached
degree = 5

domain = np.meshgrid(*[np.linspace(-100, 100, 50) for i in range(degree)])
print domain

# We're running into a memory error! 
# Solutions: make it sparse (helps with memory), use smaller types (8bit ints)


# print calcDeriv(points, [1, 1], 0)
# print calcError(points, [1, 1])
