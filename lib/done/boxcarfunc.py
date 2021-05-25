import numpy as np
from scipy.signal import medfilt as median
from scipy.ndimage.filters import uniform_filter1d as smooth
from scipy.interpolate import interp1d as interpol
from .excep import VectorLengthException

#TODO boxcarfunc docs
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

def boxcarfunc(x_vals, data_v, var_v, spec_v, eval, coeffv, boxcar_hw):
	# Check inputs
	nx = len(x_vals)

	for vector in [data_v, var_v, spec_v]:
		if len(vector) != nx:
			raise VectorLengthException(vector, nx)

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

