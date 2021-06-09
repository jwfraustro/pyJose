# TODO DOCS: vectfit
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
from scipy.interpolate import CubicSpline

def newton_cotes(xvals, datav):

	def booles_rule(y, h):

		bl = (7*y[0]+
		      32*y[1]+
		      12*y[2]+
		      32*y[3]+
		      7*y[4])*(2*h)/45

		return bl

	# check if arrays are sorted
	# if not, do so
	if not all(xvals[i] <= xvals[i + 1] for i in range(len(xvals) - 1)):

		zipped_lists = zip(xvals, datav)
		sort_pairs = sorted(zipped_lists)
		tups = zip(*sort_pairs)
		xvals, datav = [list(tup) for tup in tups]

	# Divide the xvals into 5 equally spaced points

	a = min(xvals)
	b = max(xvals)
	h = (b - a) / 4

	xvals_new = [a, a+h, a+2*h, a+3*h, b]
	print(xvals_new)

	# Use cubic spline interpolation to calculate datav points at those values
	f = CubicSpline(xvals, datav)
	datav_new = f(xvals_new)
	print(datav_new)

	# Perform Boole's Rule
	bl = booles_rule(datav_new, h)

	return bl

def gausseval(x, a, f=None, p=None):
	# TODO DOCS: gausseval
	"""
	Name:
		gausseval
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


def gaussfunc(RunConfig):
	# Check Inputs
	nx = len(RunConfig.xvals)

	if nx != len(RunConfig.datav):
		raise VectorLengthException("x_vals", "data_v")
	if nx != len(RunConfig.varv):
		raise VectorLengthException("x_vals", "var_v")
	if nx != len(RunConfig.specv) and len(RunConfig.specv) != 1:
		raise VectorLengthException("x_vals", "spec_v")
	if nx <= 4:
		raise ParameterException("data_v has less than 4 elements.")

	# Evaluate Coefficients
	if eval:
		f = gausseval(RunConfig.xvals, RunConfig.coeff)
		return f

	# Prepare Fitting
	 if RunConfig.reest:
		RunConfig.specv = newton_cotes(RunConfig.xvals, RunConfig.datav)  # reestimate spectrum

	RunConfig.hbw = 4 < nx / 2 - 1  # half box width for estimation
	RunConfig.sdata = smooth(RunConfig.datav / RunConfig.specv, 2 * RunConfig.hbw + 1)  # estimate gauss
	RunConfig.shi = max(RunConfig.sdata, RunConfig.chi)
	RunConfig.slo = min(RunConfig.sdata, RunConfig.cli)
	if abs(RunConfig.shi) > abs(RunConfig.slo):
		RunConfig.ci = RunConfig.chi
	else:
		RunConfig.ci = RunConfig.cli
	RunConfig.hi = (RunConfig.datav / RunConfig.specv > RunConfig.sdata)[RunConfig.chi]  # estimated gaussian height
	RunConfig.center = RunConfig.xvals[RunConfig.ci]  # estimated gaussian center
	RunConfig.top = np.where(abs(RunConfig.datav / RunConfig.specv) > abs(RunConfig.hi / exp(1)))
	RunConfig.fwhm = max(RunConfig.xvals[RunConfig.top]) - min(RunConfig.xvals[RunConfig.top]) + 1  # estimated gaussian fullwidthhalfmax
	RunConfig.base = np.median(RunConfig.sdata)  # estimated base level

	RunConfig.coeff = [RunConfig.hi, RunConfig.center, RunConfig.fwhm / 2]  # coefficients for estimate

	RunConfig.nz = np.where(RunConfig.varv != 0)  # locations to fit

	RunConfig.est = RunConfig.datav / RunConfig.specv
	RunConfig.origcoeff = RunConfig.coeff

	# FIXME gaussfit fit data
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

def boxcarfunc(RunConfig):
	# TODO DOCS: boxcarfunc
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


	Modification History:



	"""
	# Check inputs
	nx = len(RunConfig.x_vals)

	if len(RunConfig.datav) != nx:
		raise VectorLengthException("data_v", "x_vals")
	if len(RunConfig.varv) != nx:
		raise VectorLengthException("var_v", "x_vals")
	if len(RunConfig.specv) != nx:
		raise VectorLengthException("spec_v", "x_vals")

	RunConfig.coeffv = [np.average(RunConfig.datav / RunConfig.specv)]

	# Fit Data
	if not eval:
		return median(RunConfig.datav / RunConfig.specv, RunConfig.boxcarhw * 2 + 1)

	# Evaluate
	if np.count_nonzero == 0:
		return np.zeros(nx)

	RunConfig.gv = RunConfig.specv.nonzero()

	RunConfig.goodx = RunConfig.xvals[RunConfig.gv]

	RunConfig.estg = smooth(RunConfig.datav[RunConfig.gv] / RunConfig.spec_v[RunConfig.gv], RunConfig.boxcarhw * 2 + 1)
	RunConfig.fiteval = np.zeros(nx)
	RunConfig.bv = np.where(RunConfig.specv == 0)[0]
	if len(RunConfig.bv) != 0:
		RunConfig.badx = RunConfig.xvals[RunConfig.bv]
		RunConfig.estb_func = interpol(RunConfig.goodx, RunConfig.estg, kind='linear')
		RunConfig.estb = RunConfig.estb_func(RunConfig.badx)
		RunConfig.fiteval[RunConfig.badx] = RunConfig.estb

	RunConfig.fiteval[RunConfig.goodx] = RunConfig.estg

	return RunConfig.fiteval


def centermass(RunConfig):
	# TODO DOCS: centermass
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
	nx = len(RunConfig.xvals)

	if len(RunConfig.datav) != nx:
		raise VectorLengthException("data_v", "x_vals")

	RunConfig.multv = sum(RunConfig.datav / RunConfig.specv)
	return sum(RunConfig.xvals * RunConfig.datav / RunConfig.specv / RunConfig.multv)


def extractfunc(RunConfig):
	# TODO DOCS: extractfunc
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

	nx = len(RunConfig.datav)

	if nx != len(RunConfig.profv):
		raise VectorLengthException("datav", "profv")
	if nx != len(RunConfig.varv):
		raise VectorLengthException("datav", "varv")

	# Always Extract

	RunConfig.gl = np.where(RunConfig.profv != 0)
	if (len(RunConfig.gl) == 0):
		print("No good pixels in datav.")
		return 0

	RunConfig.denom = np.sum((RunConfig.profv[RunConfig.gl] * RunConfig.profv[RunConfig.gl]) / RunConfig.varv[RunConfig.gl])  # avoid recalc
	RunConfig.opt = np.sum((RunConfig.profv[RunConfig.gl] * RunConfig.datav[RunConfig.gl]) / RunConfig.varv[RunConfig.gl]) / RunConfig.denom
	RunConfig.opvar = np.sum(RunConfig.profv[RunConfig.gl]) / RunConfig.denom

	return RunConfig.opt


def polyeval(RunConfig):
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

	return np.polyval(RunConfig.coeffs, RunConfig.xvals)

def polyfunc(RunConfig):
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
	nx = len(RunConfig.xvals)

	if len(RunConfig.deg) == 0:
		deg = 2

	if nx != len(RunConfig.datav):
		raise VectorLengthException(nx, RunConfig.datav)
	if nx != len(RunConfig.varv):
		raise VectorLengthException(nx, RunConfig.varv)
	if nx != len(RunConfig.specv):
		raise VectorLengthException(nx, RunConfig.specv)
	if RunConfig.deg < 0:
		raise ParameterException("Degree cannot be < 0.")
	if nx <= RunConfig.deg:
		raise ParameterException("Number of xvals must be greater than degree.")

	# Evaluate Coefficients
	if eval:  # evalate given coefficients
		RunConfig.fiteval = polyeval(RunConfig)
		RunConfig.zl = np.where(RunConfig.varv == 0.)  # use actual data at varv 0 locations
		if len(RunConfig.zl) > 0:
			RunConfig.fiteval[RunConfig.zl] = (RunConfig.datav[RunConfig.zl] / RunConfig.specv[RunConfig.zl])

		return RunConfig.fiteval

	# Fit Data
	RunConfig.nz = np.where(RunConfig.varv != 0.)  # locations where variance is zero
	RunConfig.est = RunConfig.datav / RunConfig.specv  # initial estimate
	if len(RunConfig.nz > 0):
		if RunConfig.deg == 0:  # use an average over the column
			RunConfig.mn = np.average(RunConfig.datav / RunConfig.specv)
			RunConfig.est[RunConfig.nz] = RunConfig.mn
			RunConfig.coeffv = [RunConfig.mn, 0]  # correct for right length
		else:
			RunConfig.merrors = RunConfig.varv / RunConfig.specv ** 2 > 1E-8  # use all errors for estimation
			if (RunConfig.deg == 1):
				RunConfig.coeffv = np.polyfit(RunConfig.xvals, RunConfig.datav / RunConfig.specv, 1)
			else:
				RunConfig.coeffv = np.polyfit(RunConfig.xvals, RunConfig.datav / RunConfig.specv, RunConfig.deg)

	return RunConfig.est

if __name__ == '__main__':
	x = [0.0, .12, .22, .32, .36, .40, .44, .54, .64, .70, .80]
	f = [0.200000, 1.30973, 1.30524, 1.74339, 2.07490, 2.45600, 2.84299, 3.50730, 3.18194, 2.36302, 0.231964]
	print(newton_cotes(x, f))