#TODO gaussfunc docs
"""
Name:
	gaussfunc
Purpose:
	Function for estimating a Gaussian curve to the data.
Category:
	Optimal Spectrum Extraction Package
		- Vector Fitting functions
Calling Example:

Inputs:

Outputs:

History:

Created on 4/17/2021$
"""

from lib.done.excep import VectorLengthException, ParameterException
from math import exp
import numpy as np
from scipy.ndimage.filters import uniform_filter1d as smooth

def gausseval(x,a, f=None, p=None):

	z = (x - a[1]) / a[2]

	f = a[0]*exp(-z**2/2)

	if not p:
		return f
	else:

		p = np.zeros((len(x), 3))

		p[:, 0] = f/a[0]
		p[:, 1] = f*z/a[2]
		p[:, 2] = p[:, 1] * z

	return f, p



def gaussfunc(xvals, datav, varv, specv, eval, coeff, reest=None):

	# Check Inputs
	nx = len(xvals)

	if nx != len(datav):
		raise VectorLengthException(nx, datav)
	if nx != len(varv):
		raise VectorLengthException(nx, varv)
	if nx != len(specv) and len(specv) != 1:
		raise VectorLengthException(nx, specv)
	if nx <= 4:
		raise ParameterException("datav has less than 4 elements.")

	# Evaluate Coefficients
	if eval:
		f = gausseval(xvals, coeff)
		return f

	# Prepare Fitting
	# TODO INT_TABULATED: 5-Point Newton-Cotes implementation
	#if reest:
		specv = int_tabulated(xvals, datav) #reestimate spectrum

	hbw = 4 < nx / 2-1  #half box width for estimation
	sdata = smooth(datav / specv, 2*hbw+1)  #estimate gauss
	shi = max(sdata, chi)
	slo = min(sdata, cli)
	if abs(shi) > abs(slo):
		ci = chi
	else:
		ci = cli
	hi = (datav/specv > sdata)[chi] # estimated gaussian height
	center = xvals[ci]              # estimated gaussian center
	top = np.where(abs(datav/specv) > abs(hi/exp(1)))
	fwhm = max(xvals[top]) - min(xvals[top]) + 1    #estimated gaussian fullwidthhalfmax
	base = np.median(sdata) #estimated base level

	coeff = [hi, center, fwhm/2]    #coefficients for estimate

	nz = np.where(varv != 0)    #locations to fit

	est = datav / specv
	origcoeff = coeff

	# TODO FIT DATA
# 	; FIT
# 	DATA
#
# 	;est[nz] = gaussfit(xvals[nz], datav[nz] / specv, coeff, $
# 	;                     measure_errors = sqrt(varv[nz]), nterms = 4, $
# 	;                     estimates = origcoeff, chisq = chisq)
#
# 	merrors = specv ^ 2 / varv > 1e-8
# 	estz = curvefit(xvals, datav / specv, merrors, $
# 	coeff, iter = fititer, function_name = 'gausseval', $
# 	tol = 1e-6, / double, status = status, $
# 	itmax = 10, chisq = chisq)
# 	est[nz] = estz[nz]
#
# 	if (status ne 0) then begin; estimation wrong so use median instead
# 	mdata = median(datav / specv, 2 * hbw + 1)
# 	estz = curvefit(xvals, mdata, merrors, $
# 	origcoeff, iter = fititer, function_name = 'gausseval', $
# 	/ double, status = status)
# 	est[nz] = estz[nz]
#
#
# endif
#
# return, est
# end


	return