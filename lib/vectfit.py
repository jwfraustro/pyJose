# TODO vectorfit docs
"""
Name:

Purpose:

Category:

Calling Example:

Inputs:

Outputs:

History:

Created on 5/25/2021$
"""

from math import exp
from lib.excep import VectorLengthException, ParameterException
import numpy as np
from scipy.signal import medfilt as median
from scipy.ndimage.filters import uniform_filter1d as smooth
from scipy.interpolate import interp1d as interpol

def gausseval(x, a, f=None, p=None):
	# TODO gaussfunc docs
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

	z = (x - a[1]) / a[2]

	f = a[0] * exp(-z ** 2 / 2)

	if not p:
		return f
	else:

		p = np.zeros((len(x), 3))

		p[:, 0] = f / a[0]
		p[:, 1] = f * z / a[2]
		p[:, 2] = p[:, 1] * z

	return f, p


def gaussfunc(xvals, datav, varv, specv, eval, coeff, reest=None):
	# Check Inputs
	nx = len(xvals)

	if nx != len(datav):
		raise VectorLengthException("x_vals", "data_v")
	if nx != len(varv):
		raise VectorLengthException("x_vals", "var_v")
	if nx != len(specv) and len(specv) != 1:
		raise VectorLengthException("x_vals", "spec_v")
	if nx <= 4:
		raise ParameterException("data_v has less than 4 elements.")

	# Evaluate Coefficients
	if eval:
		f = gausseval(xvals, coeff)
		return f

		# Prepare Fitting
		# TODO INT_TABULATED: 5-Point Newton-Cotes implementation
		# if reest:
		specv = int_tabulated(xvals, datav)  # reestimate spectrum

	hbw = 4 < nx / 2 - 1  # half box width for estimation
	sdata = smooth(datav / specv, 2 * hbw + 1)  # estimate gauss
	shi = max(sdata, chi)
	slo = min(sdata, cli)
	if abs(shi) > abs(slo):
		ci = chi
	else:
		ci = cli
	hi = (datav / specv > sdata)[chi]  # estimated gaussian height
	center = xvals[ci]  # estimated gaussian center
	top = np.where(abs(datav / specv) > abs(hi / exp(1)))
	fwhm = max(xvals[top]) - min(xvals[top]) + 1  # estimated gaussian fullwidthhalfmax
	base = np.median(sdata)  # estimated base level

	coeff = [hi, center, fwhm / 2]  # coefficients for estimate

	nz = np.where(varv != 0)  # locations to fit

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

def boxcarfunc(x_vals, data_v, var_v, spec_v, eval, coeffv, boxcar_hw):
	# TODO boxcarfunc docs
	"""
	Function: boxcarfunc

	Description:
		Function for fitting the data with a boxcar average or median fit.

	Category:
		Optimal Spectrum Extraction Package
			- Vector fitting functions

	Parameters:
		x_vals:     The x values for the data
		data_v:     The vector with cosmic rays to be fitted
		var_v:      The variance vector
		spec_v:     The spectrum vector. Bad pixels have a spec_v of 0
		eval:       Set to True to evaluate the inputs
		coeffv:     A filler in this routine
		boxcar_hw:  The halfwidth of the median or smooth filter

	Outputs:
		Returns an estimation of the data by returning data_v ran
		through a median filter. If eval is set the routine smooths
		over the good pixels, then interpolates the bad pixels.

	Restrictions:
		All vectors must be the same length. Boxcar_hw must be an int >= 1

	Procedure:
		If eval == False, then do a median filter over the data_v / spec_v. Else,
		smooth over the good pixels and interpolate over the bad pixels.

	Example:
		TODO

	Modification History:
		TODO


	"""
	# Check inputs
	nx = len(x_vals)

	if len(data_v) != nx:
		raise VectorLengthException("data_v", "x_vals")
	if len(var_v) != nx:
		raise VectorLengthException("var_v", "x_vals")
	if len(spec_v) != nx:
		raise VectorLengthException("spec_v", "x_vals")

	coeffv = [np.average(data_v / spec_v)]

	# Fit Data
	if not eval:
		return median(data_v / spec_v, boxcar_hw * 2 + 1)

	# Evaluate
	if np.count_nonzero == 0:
		return np.zeros(nx)

	gv = spec_v.nonzero()

	goodx = x_vals[gv]

	estg = smooth(data_v[gv] / spec_v[gv], boxcar_hw * 2 + 1)
	fiteval = np.zeros(nx)
	bv = np.where(spec_v == 0)[0]
	if len(bv) != 0:
		badx = x_vals[bv]
		estb_func = interpol(goodx, estg, kind='linear')
		estb = estb_func(badx)
		fiteval[badx] = estb

	fiteval[goodx] = estg

	return fiteval


def centermass(x_vals, data_v, var_v, spec_v):
	# TODO centermass docs
	"""
	Function: centermass

	Description:
		Function which returns the center of mass of datav

	Category:
		Optimal Spectrum Extraction Package
			- Vector fitting functions

	Parameters:
		x_vals: The x values for the data
		data_v: the sky subtracted data
		var_v:  the variance
		spec_v: the spectrum

	Outputs:
		Returns the center of mass of the data vector / spec vector

	Restrictions:
		All vectors must be the same length.

	Modification History:
		TODO


	"""
	nx = len(x_vals)

	if len(data_v) != nx:
		raise VectorLengthException("data_v", "x_vals")

	multv = sum(data_v / spec_v)
	return sum(x_vals * data_v / spec_v / multv)


def extractfunc(xvals, datav, varv, profv, eval, coeffv, opvar):
	# TODO extractfunc docs
	"""
	Name:
		extractfunc
	Purpose:
		A function which extracts the optimal spectrum from a vector.
	Category:
		Optimal Spectrum Extraction Package
			- Vector Fitting functions
	Calling Example:

	Inputs:

	Outputs:

	History:

	Created on 4/17/2021$
	"""
	# Check Inputs

	nx = len(datav)

	if nx != len(profv):
		raise VectorLengthException("datav", "profv")
	if nx != len(varv):
		raise VectorLengthException("datav", "varv")

	# Always Extract

	gl = np.where(profv != 0)
	if (len(gl) == 0):
		print("No good pixels in datav.")
		return 0

	denom = np.sum((profv[gl] * profv[gl]) / varv[gl])  # avoid recalc
	opt = np.sum((profv[gl] * datav[gl]) / varv[gl]) / denom
	opvar = np.sum(profv[gl]) / denom

	return opt


def polyeval(coeffs, xvals):
	# TODO polyeval docs
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

	return np.polyval(coeffs, xvals)

def polyfunc(xvals, datav, varv, specv, eval, coeffv, deg):
	# TODO polyfunc docs
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
	if eval:  # evalate given coefficients
		fiteval = polyeval(coeffv, xvals)
		zl = np.where(varv == 0.)  # use actual data at varv 0 locations
		if len(zl) > 0:
			fiteval[zl] = (datav[zl] / specv[zl])

		return fiteval

	# Fit Data
	nz = np.where(varv != 0.)  # locations where variance is zero
	est = datav / specv  # initial estimate
	if len(nz > 0):
		if deg == 0:  # use an average over the column
			mn = np.average(datav / specv)
			est[nz] = mn
			coeffv = [mn, 0]  # correct for right length
		else:
			merrors = varv / specv ** 2 > 1E-8  # use all errors for estimation
			if (deg == 1):
				coeffv = np.polyfit(xvals, datav / specv, 1)
			else:
				coeffv = np.polyfit(xvals, datav / specv, deg)

	# TODO Where does estz in the following from the original come from? It's never declared...
	# coeffv = linfit(xvals, datav/specv, $
	#                       measure_errors = merrors, yfit = estz)
	# coeffv = poly_fit(xvals, datav/specv, deg, /double, $
	#                         yfit = estz, yband = yband, measure_errors = merrors)
	return est
