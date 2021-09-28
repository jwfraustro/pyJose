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

		bl = ( 7 * y[0] +
		      32 * y[1] +
		      12 * y[2] +
		      32 * y[3] +
		      7 * y[4]) * (2 * h) / 45

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

	xvals_new = [a, a + h, a + 2 * h, a + 3 * h, b]
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
		# 
		p[:, 0] = f / a[0]
		p[:, 1] = f * z / a[2]
		p[:, 2] = p[:, 1] * z

	return f, p


def gaussfunc(rc):
	# Check Inputs
	nx = len(rc.xvals)

	if nx != len(rc.datav):
		raise VectorLengthException("x_vals", "data_v")
	if nx != len(rc.varv):
		raise VectorLengthException("x_vals", "var_v")
	if nx != len(rc.specv) and len(rc.specv) != 1:
		raise VectorLengthException("x_vals", "spec_v")
	if nx <= 4:
		raise ParameterException("data_v has less than 4 elements.")

	# Evaluate Coefficients
	if eval:
		f = gausseval(rc.xvals, rc.coeff)
		return f

	# Prepare Fitting
	if rc.reest:
		rc.specv = newton_cotes(rc.xvals, rc.datav)  # reestimate spectrum

	rc.hbw = 4 < nx / 2 - 1  # half box width for estimation
	rc.sdata = smooth(rc.datav / rc.specv, 2 * rc.hbw + 1)  # estimate gauss
	rc.shi = max(rc.sdata, rc.chi)
	rc.slo = min(rc.sdata, rc.cli)
	if abs(rc.shi) > abs(rc.slo):
		rc.ci = rc.chi
	else:
		rc.ci = rc.cli
	rc.hi = (rc.datav / rc.specv > rc.sdata)[rc.chi]  # estimated gaussian height
	rc.center = rc.xvals[rc.ci]  # estimated gaussian center
	rc.top = np.where(abs(rc.datav / rc.specv) > abs(rc.hi / exp(1)))
	rc.fwhm = max(rc.xvals[rc.top]) - min(rc.xvals[rc.top]) + 1  # estimated gaussian fullwidthhalfmax
	rc.base = np.median(rc.sdata)  # estimated base level

	rc.coeff = [rc.hi, rc.center, rc.fwhm / 2]  # coefficients for estimate

	rc.nz = np.where(rc.varv != 0)  # locations to fit

	rc.est = rc.datav / rc.specv
	rc.origcoeff = rc.coeff

	# FIXME gaussfit fit data
	#
	# 	;est[nz] = gaussfit(xvals[nz], datav[nz] / specv, coeff, $
	# 	;                     measure_errors = sqrt(varv[nz]), nterms = 4, $
	# 	;                     estimates = origcoeff, chisq = chisq)
	#
	# 	merrors = specv ** 2 / varv > 1e-8
	# 	estz = curvefit(xvals, datav / specv, merrors, $
	# 	coeff, iter = fititer, function_name = 'gausseval', $
	# 	tol = 1e-6, / double, status = status, $
	# 	itmax = 10, chisq = chisq)
	# 	est[nz] = estz[nz]
	#
	# 	if (status != 0): then begin; estimation wrong so use median instead
	# 	mdata = median(datav / specv, 2 * hbw + 1)
	# 	estz = curvefit(xvals, mdata, merrors, $
	# 	origcoeff, iter = fititer, function_name = 'gausseval', $
	# 	/ double, status = status)
	# 	est[nz] = estz[nz]
	#
	#

	return # est

def boxcarfunc(rc):
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
	nx = len(rc.xvals)

	if len(rc.datav) != nx:
		raise VectorLengthException("data_v", "x_vals")
	if len(rc.varv) != nx:
		raise VectorLengthException("var_v", "x_vals")
	if len(rc.specv) != nx:
		raise VectorLengthException("spec_v", "x_vals")

	rc.coeffv = [np.average(rc.datav / rc.specv)]

	# Fit Data
	if not eval:
		return median(rc.datav / rc.specv, rc.boxcarhw * 2 + 1)

	# Evaluate
	if np.count_nonzero == 0:
		return np.zeros(nx)

	gv = rc.specv.nonzero()

	goodx = rc.xvals[gv]

	estg = smooth(rc.datav[gv] / rc.specv[gv], rc.boxcarhw * 2 + 1)
	rc.fiteval = np.zeros(nx)
	bv = np.where(rc.specv == 0)[0]
	if len(bv) != 0:
		badx = rc.xvals[bv]
		estb_func = interpol(goodx, estg, kind='linear')
		estb = estb_func(badx)
		rc.fiteval[badx] = estb

	rc.fiteval[goodx] = estg

	return rc.fiteval


def centermass(rc):
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

	"""
	nx = len(rc.xvals)

	if len(rc.datav) != nx:
		raise VectorLengthException("data_v", "x_vals")

	rc.multv = sum(rc.datav / rc.specv)
	return sum(rc.xvals * rc.datav / rc.specv / rc.multv)


def extractfunc(rc):
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

	nx = len(rc.datav)

	if nx != len(rc.profv):
		raise VectorLengthException("datav", "profv")
	if nx != len(rc.varv):
		raise VectorLengthException("datav", "varv")

	# Always Extract

	gl = np.where(rc.profv != 0)
	if (len(gl) == 0):
		print("No good pixels in datav.")
		return 0

	denom = np.sum((rc.profv[gl] * rc.profv[gl]) / rc.varv[gl])  # avoid recalc
	rc.opt = np.sum((rc.profv[gl] * rc.datav[gl]) / rc.varv[gl]) / denom
	rc.opvar = np.sum(rc.profv[gl]) / denom

	return rc.opt


if __name__ == '__main__':
	x = [0.0, .12, .22, .32, .36, .40, .44, .54, .64, .70, .80]
	f = [0.200000, 1.30973, 1.30524, 1.74339, 2.07490, 2.45600, 2.84299, 3.50730, 3.18194, 2.36302, 0.231964]
	print(newton_cotes(x, f))