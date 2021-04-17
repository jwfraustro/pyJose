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

def polyfunc(xvals, datav, varv, specv, eval, coeffv, deg):
	#TODO

	return