# Created by Patrick Huston and Nur Shlapobersky on 10/10/15

from __future__ import division
import numpy as np
import itertools
import matplotlib.pyplot as plt
import random

def calc_error(points, coeffVals):
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


def calc_all_gradient(points, coeffVals):
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

    gradientAtLoc = []

    for i in range(len(coeffVals)):
        gradientAtLoc.append(calc_gradient(points, coeffVals, i))

    return gradientAtLoc


def gradientDescent(points, degree):
    '''
    Performs a gradient descent using coeffs as the domain and returns the
    location of the minimum (in coefficient space), and the number of
    iterations to reach the minumum.

    Inputs:
        points (list of tuples): points in input data set
        degree (int): the length of the coefficient arrays in coeffs

    Returns:
        it (int): Number of iterations calculation took
        fCoeffs (list): fitted coefficients
    '''

    # initialize variables and arrays we'll need for gr dsc
    iCoeffs = np.zeros(degree)
    iError = calc_error(points, iCoeffs)
    lambdas = np.logspace(-10, 2, 100)
    it = 0
    grad = 100*np.ones(degree)

    # perform the gr dsc
    while (np.linalg.norm(grad) > (10**(degree - 2))):
        plt.ion()
        it += 1
        grad = calc_all_gradient(points, iCoeffs)
        # grad1 = calc_all_gradient_test(points, iCoeffs)
        if it%200 is 0:
            print "Gradient: ", np.linalg.norm(grad)
            print "Error: ", iError
            plot_results(points, (it, fCoeffs))

        optiLambda = optimalLambda(points, iCoeffs, grad, lambdas)
        fCoeffs = iCoeffs - np.multiply(optiLambda,grad);
        fError = calc_error(points, fCoeffs);
        iCoeffs = fCoeffs;
        iError = fError;

    return it, fCoeffs


def optimalLambda(points, iCoeffs, grad, lambdas):
    '''
    Returns the optimal lambda for a given gradient descent step.

    Inputs:
        points (list of tuples): Points in input data set
        iCoeffs (tuple): starting location (coefficients) for current step
        grad (list): the gradient vector for the current location
        lambdas (list): list of potential lambda values
    '''

    new_coeffs = []

    for lambdo in lambdas:
        new_coeffs.append(iCoeffs - np.multiply(lambdo, grad))

    potential_errors = [calc_error(points, new_coeff_set) for new_coeff_set in new_coeffs]

    return lambdas[potential_errors.index(min(potential_errors))]


def eval_func(points, optimVals):
    '''
    Evaluates the function described by the calculated
    coefficients, puts results into a list which can
    be plotted more conveniently.

    Inputs:
        points (list of tuples): Points in input data set
        optimVals (list): Calculated coefficients of polynomial
    '''

    # Calculate min and max of data set to determine domain
    min_x = min(points, key=lambda point:point[0])[0]
    max_x = max(points, key=lambda point:point[0])[0]

    domain = np.linspace(min_x, max_x, 100)

    func_res = []

    # Calculate function result at each point in domain
    for x in domain:
        res = 0
        for index, coeff in enumerate(optimVals[1]):
            res += coeff*(x**index)

        func_res.append(res)
    return domain, func_res


def plot_results(points, optimVals):
    '''
    Plots the input data set against the polynomial
    function described by the coefficients

    Inputs:
        points (list of tuples): Points in data set
        optimVals (list): Computed coefficients of fitted polynomial

    '''
    domain, function_res = eval_func(points, optimVals)

    xs = [point[0] for point in points]
    ys = [point[1] for point in points]

    plt.clf()
    plt.scatter(xs, ys)
    plt.plot(domain, function_res)
    plt.pause(0.0001)

    plt.show()


def pol_reg(points, TSS, sampleSize):
    '''
    Starts of pipeline of polynomial regression

    Inputs:
        points (list of tuples): Points in data set to be fitted
        TSS (int): Calculated total sum of squares
        sampleSize ()

    Returns:
        optimCoeffs (list): Optimal coefficients for best fit

    '''
    degree = 2
    prevRbar2 = -1
    while degree < 6:
        optimCoeffs = gradientDescent(points, degree)
        RSS = calc_error(points, optimCoeffs[1])

        R2 = 1 - (RSS/TSS)
        Rbar2 = R2 - (1-R2)*degree/(sampleSize - degree - 2)
        print "Degree: ", degree
        print "R2: ", R2
        print "Rbar2: ", Rbar2

        if prevRbar2 > Rbar2:
            print 'Break'
            return prevOptimCoeffs

        prevRbar2 = Rbar2
        degree += 1
        prevOptimCoeffs = optimCoeffs

    return optimCoeffs


if __name__ == '__main__':

    # Generate random data set loosely representing a parabola
    points = [(i,i**2+random.randint(0,100)) for i in range(-10,10)]
    yp = [point[1] for point in points]
    yAvg = sum(yp)/len(yp)

    sampleSize = len(yp)
    TSS = sum([(y-yAvg)**2 for y in yp])

    results = pol_reg(points, TSS, sampleSize)
    print "Iterations, Coeff: ", results

    plt.ioff()
    plot_results(points, results)