"""
Name:
	polyeval

Purpose:
	Takes an array of a set of polynomial coefficients and
	evaluates a given set of xvalues using them

Category:

Calling Example:
	result = polyeval(coeffs, xvals)

Inputs:
	coeffs: an array containing a set of polynomial coefficients
	xvals: locations to evaluate. Its type is the return type

Outputs:
	Returns an array, dimension len(xvals) containing xvals
	evaluated with the polynomial coefficients

History:
	Written by:    Dara Zeehandelaar
	1 Sept 2003:   John Dermody, Setting return type to xvals
	Apr 2021:      Adapted for python by Joshua Fraustro

Created on 4/17/2021$
"""

import numpy as np

def polyeval(coeffs, xvals):

	degree = len(coeffs)
	xvals = np.array(xvals)

	polys = np.zeros(len(xvals))
	polys[0] = coeffs[0]

	if degree > 0:
		for d in range(1, degree):
			polys = polys + coeffs[d] * xvals ** d

	return polys


