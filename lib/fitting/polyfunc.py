"""
Name:

Purpose:

Category:

Calling Example:

Inputs:

Outputs:

History:

Created on 9/19/2021$
"""

import numpy as np
from numpy.polynomial import polynomial as poly
from lib.excep import *


def polyeval(coeffs, xvals):
	# TODO DOCS: polyeval
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

	return poly.polyval(xvals, coeffs)


def polyfunc(xvals, datav, varv, specv, eval, coeffv, deg):
	# TODO DOCS: polyfunc
	"""
	Name:
		polyfunc

	Purpose:
		Function for fitting the data using a polynomial fit

	Category:
		Optimal Spectrum Extraction Package
			- Vector fitting functions

	Calling Example:
		result = polyfunc(xvals, datav, varv, specv, eval, coeffv, deg)

	Inputs:
		xvals:  indices of the datav
		datav:  the data vector to fit, must be the same size as xvals
		varv:   the variance vector
		specv:  the spectrum vector
		deg:    the degree of fit

	Outputs:
		The estimated spectrum

	Optional Outputs:
		coeffv: Will contain the coefficients of the fit

	Procedure:
		If eval then evaluate the coefficients passed in. Otherwise, fit the
		data / spectrum using the variances as weights using a polynomial function.
		If deg is 0 then use mean, if deg is 1 then use linfit.

	History:

	Created on 4/17/2021$
	"""

	# Check Inputs
	nx = len(xvals)

	if not deg:
		deg = 2

	if nx != len(datav):
		raise VectorLengthException(nx, len(datav))
	if nx != len(varv):
		raise VectorLengthException(nx, len(varv))
	if nx != len(specv):
		raise VectorLengthException(nx, len(specv))
	if deg < 0:
		raise ParameterException("Degree cannot be < 0.")
	if nx <= deg:
		raise ParameterException("Number of xvals must be greater than degree.")

	# Evaluate Coefficients
	if eval:  # evalate given coefficients
		fiteval = polyeval(coeffs=coeffv, xvals=xvals)
		zl = np.where(varv == 0.)  # use actual data at varv 0 locations
		if len(zl) > 0:
			fiteval[zl] = (datav[zl] / specv[zl])

		return fiteval, coeffv

	# Fit Data
	nz = np.where(varv != 0.)  # locations where variance is zero
	est = datav / specv  # initial estimate
	if len(nz) > 0:
		if deg == 0:  # use an average over the column
			mn = np.average(datav / specv)
			est[nz] = mn
			coeffv = [mn, 0]  # correct for right length
		else:
			merrors = varv / specv ** 2 > 1E-8  # use all errors for estimation
			if deg == 1:
				coeffv = poly.polyfit(xvals, datav / specv, 1)
				estz = poly.polyval(xvals, coeffv)
			else:
				coeffv = poly.polyfit(xvals, datav / specv, deg)
				estz = poly.polyval(xvals, coeffv)
			est[nz] = estz[nz]

	return est, coeffv
