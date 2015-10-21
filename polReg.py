from __future__ import division
import numpy as np
import itertools


points = [(1, 1), (1, 0), (2, 2), (12, 5432)]

coeffVals = np.meshgrid(*[np.linspace(-100, 100, 2000) for i in range(degree)], sparse=True)
# coeffVals = [1, 1, 1, 1]

def calc_point_error(points, coeffVals):
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

def calc_all_error(points, coeffMatrices):
	'''
	Calculates total error over entire domain for each set
	of coefficient values

	Inputs:
		points (list of tuples): 'data set' of points to fit
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


def calc_gradient(points, coeffVals, indexCoeff):
	'''
	Calculates derivative of function with respect to given
	coefficient index at each given error function location

	Inputs:
		points (list of tuples): 'data set' of points to fit
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

def calc_all_gradient(points, coeffMatrices):
	'''
	Calculates gradient for every point in domain, where each 
	point includes a set of coefficient values describing the 
	function

	Inputs:
		points (list of tuples): 'data set' of points to fit
		coeffMatrices (list of np matrices): List of matrices
		representing all possible combinations of coeff vals

	Returns:
		errorMatrix (np matrix): Error at each given point
	'''
	gradientMatrix = np.zeros([coeffMatrices[0].shape[i] for i in range(len(coeffMatrices.shape[0]))])

	ranges = [range(coeffMatrix.shape[i]) for i in range(len(coeffMatrix.shape)) for coeffMatrix in coeffMatrices]

	domain = itertools.product(*ranges)

	for location in domain:
		coeffVals = [coeffMatrix[location] for coeffMatrix in coeffMatrices]
		gradientAtLoc = []

		for i in range(len(coeffVals)):
			gradientAtLoc[i] = calc_gradient(points, coeffVals, i)

		gradientMatrix[location] = gradientAtLoc

	return gradientMatrix


# In final function, iterate over degrees until satisfactory fit is reached
degree = 22

def gradientDescent(coeffs, degree):
    ''' Performs a gradient descent using coeffs as the domain and returns the
        location of the minimum (in coefficient space), and the number of
        iterations to reach the minumum.

        coeffs: a list of meshgrided matrices representing
        degree: the length of the coefficient arrays in coeffs
    '''

    # initialize variables and arrays we'll need for gr dsc
    grad = calcDeriv(points, coeffVals, indexCoeff)
    iCoeffs = numpy.zeros(degree)
    iz = func(coeffs)
    lambdas = np.logspace(-8, 1, 50)
    it = 0
    grad = numpy.ones(degree)

    # perform the gr dsc
    while (np.linalg.norm(grad) > .00001):
        it  = it + 1
        grad = gradient(coeffs)
        optiLambda = optimalLambda(iCoeffs, grad, lamdbas)
        fCoeffs = iCoeffs - optiLambda*grad;
        fz = func(fCoeffs);
        iCoeffs = fCoeffs;
        iz = fz;

    return it, fCoeffs

def optimalLambda(iCoeffs, grad, lambdas):
    ''' Returns the optimal lambda for a given gradient descent step.

        iCoeffs: starting location (coefficients) for current step
        grad: the gradient vector for the current location
        lambdas: list of potential lambda values
    '''

    pCoeffs = iCoeffs - lambdas*grad
    pZ = func(pCoeffs)
    ind = Zp.argmax(axis=0)
    optiLambda = lambdas[ind];

    return optiLamdba