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

from lib.done.excep import VectorLengthException, ParameterException
from lib.done.polyeval import polyeval
import numpy as np

def polyfunc(xvals, datav, varv, specv, eval, coeffv, deg):

	# Check Inputs
	nx = len(xvals)

	if len(deg) == 0:
		deg = 2

	if nx != len(datav):
		raise VectorLengthException(nx, datav)
	if nx != len(varv):
		raise VectorLengthException(nx, varv)
	if nx != len(specv):
		raise VectorLengthException(nx, specv)
	if deg < 0:
		raise ParameterException("Degree cannot be < 0.")
	if nx <= deg:
		raise ParameterException("Number of xvals must be greater than degree.")

	# Evaluate Coefficients
	if eval: # evalate given coefficients
		fiteval = polyeval(coeffv, xvals)
		zl = np.where(varv == 0.)   # use actual data at varv 0 locations
		if len(zl) > 0:
			fiteval[zl] = (datav[zl] / specv[zl])

		return fiteval

	# Fit Data
	nz = np.where(varv != 0.)   # locations where variance is zero
	est = datav / specv         # initial estimate
	if len(nz > 0):
		if deg == 0:    # use an average over the column
			mn = np.average(datav / specv)
			est[nz] = mn
			coeffv = [mn, 0]    # correct for right length
		else:
			merrors = varv / specv**2 > 1E-8    # use all errors for estimation
			if (deg == 1):
				coeffv = np.polyfit(xvals, datav/specv, 1)
			else:
				coeffv = np.polyfit(xvals, datav/specv, deg)

	# TODO Where does estz in the following from the original come from? It's never declared...
	# coeffv = linfit(xvals, datav/specv, $
	#                       measure_errors = merrors, yfit = estz)
	# coeffv = poly_fit(xvals, datav/specv, deg, /double, $
	#                         yfit = estz, yband = yband, measure_errors = merrors)
	return est