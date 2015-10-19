from __future__ import division
import numpy as np
import sympy as sp
import string
import itertools

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

def calcTotalError(points, coeffMatrices):
	'''
	Calculates total error over entire domain for each set
	of coefficient values
	
	Inputs:
		coeffMatrices (list of np matrices): List of matrices
		representing all possible combinations of coeff vals

	Returns: 
		errorMatrix (np matrix): Error at each given pointl		 	

	'''
	errorMatrix = np.zeros([coeffMatrices[0].shape[i] for i in range(len(coeffMatrices.shape[0]))])

	ranges = [range(coeffMatrix.shape[i]) for i in range(len(coeffMatrix.shape)) for coeffMatrix in coeffMatrices]

	domain = itertools.product(*ranges)

	for location in domain:
		coeffVals = [coeffMatrix[location] for coeffMatrix in coeffMatrices]
		errorAtLoc = calcError(points, coeffVals)
		errorMatrix[location] = errorAtLoc

	return errorMatrix


def indexByTuple(matrix, index):
	'''
	Indexes matrix by location represented as tuples
	i.e. indexByTuple(matrix, (1,2,3)) = matrix[1][2][3]

	Inputs:
		matrix (np array): Matrix to be indexed
		index (tuple): Represents location in matrix

	Returns
		val (int): Value at location 'index'

	'''

	val = matrix
	for i in index:
		val = val[i]

	return val


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
