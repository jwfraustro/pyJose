from lib.excep import *
import numpy as np
from scipy.ndimage.filters import uniform_filter1d as smooth
from scipy.interpolate import interp1d as interpol

def boxcarfunc(xvals, datav, varv, specv, coeffv, boxcarhw):
	"""
		Function: boxcarfunc

		Description:
			Function for fitting the data with a boxcar average or median fit.

		Category:
			Optimal Spectrum Extraction Package
				- Vector fitting functions

		Parameters:
			xvals:     The x values for the data
			datav:     The vector with cosmic rays to be fitted
			varv:      The variance vector
			specv:     The spectrum vector. Bad pixels have a specv of 0
			eval:       Set to True to evaluate the inputs
			coeffv:     A filler in this routine
			boxcarhw:  The halfwidth of the median or smooth filter

		Outputs:
			Returns an estimation of the data by returning datav ran
			through a median filter. If eval is set the routine smooths
			over the good pixels, then interpolates the bad pixels.

		Restrictions:
			All vectors must be the same length. boxcarhw must be an int >= 1

		Procedure:
			If eval == False, then do a median filter over the datav / specv. Else,
			smooth over the good pixels and interpolate over the bad pixels.

		Example:


		Modification History:



		"""
	# Check inputs
	nx = len(xvals)

	if len(datav) != nx:
		raise VectorLengthException("datav", "xvals")
	if len(varv) != nx:
		raise VectorLengthException("varv", "xvals")
	if len(specv) != nx:
		raise VectorLengthException("specv", "xvals")

	coeffv = [np.average(datav / specv)]

	# Fit Data
	if not eval:
		return np.median(datav / specv, boxcarhw * 2 + 1)

	# Evaluate
	if np.count_nonzero == 0:
		return np.zeros(nx)

	gv = specv.nonzero()

	goodx = xvals[gv]

	estg = smooth(datav[gv] / specv[gv], boxcarhw * 2 + 1)
	fiteval = np.zeros(nx)
	bv = np.where(specv == 0)[0]
	if len(bv) != 0:
		badx = xvals[bv]
		estb_func = interpol(goodx, estg, kind='linear')
		estb = estb_func(badx)
		fiteval[badx] = estb

	fiteval[goodx] = estg

	return fiteval