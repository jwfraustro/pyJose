# TODO DOCS: horne
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
from lib.excep import ParameterException
import numpy as np
from scipy import interpolate

def stdextr(rc):
	# TODO DOCS: stdextr
	"""
	Standard box extraction of spectrum, returns spectrum and variance image, and estimates over bad pixels.

	Calling Example:

	Inputs:
		dataim: Sky-subtracted, processed image to extract spectrum from. Horizontal is pixel position,
				vertical is wavelength, assuming no curvature.
		varim:  Variance image from processed image.
		x1, x2: Boundaries in x (inclusive) over which spectrum will be summed.

	Optional Keywords:
		inmask:     The cosmic ray mask for the data image.
		stdvar:     Variance of extracted spectra, length n where n is the number of wavelengths
					taken of the image (vertical direction).
		adjspec:    Set to linearly interpolate the data, then extract the standard spectrum.
					Only useful if input mask has bad pixels.
	Outputs:
		Returns stdspec, an array of length n containing extracted spectra,
		where n is the number of wavelengths of the image (vertical direction).

	Restrictions:
		dataim, varim, and inmask must be the same size. 0<=x1<=x2<=nx-1
	History:

	Created on 4/17/2021$
	"""

	# Check Inputs
	ny = np.shape(rc.dataim)[0]
	nx = np.shape(rc.dataim)[1]

	if not np.any(rc.inmask):
		rc.inmask = np.ones((nx, ny))

	if rc.x1 < 0 or rc.x1 > rc.x2:
		raise ParameterException("x1 must be greater than 0 and less than x2.")
	if rc.x2 > nx-1 or rc.x2 < rc.x1:
		raise ParameterException("x2 must be less than nx-1 and greater than x1.")
	if np.shape(rc.varim)[1] != nx or np.shape(rc.varim)[0] != ny:
		raise ParameterException("varim must be the same shape as dataim.")
	if np.shape(rc.inmask)[1] != nx or np.shape(rc.inmask)[0] != ny:
		raise ParameterException("inmask must be the same shape as dataim.")

	# Interpolate over bad pixels
	if rc.adjspec:
		rc.adjspec = np.array(ny, np.double)
		for i in range(ny):
			bv = np.where(rc.inmask[i, rc.x1:rc.x2] == 0)
			gv = np.where(rc.inmask[i, rc.x1:rc.x2] == 1)
			datav = rc.dataim[i, rc.x1:rc.x2]
			if len(bv) > 0:
				interpfunc = interpolate.interp1d(gv, datav[gv],  kind='linear')
				datav[bv] = interpfunc(bv)
			rc.adjspec[i] = np.sum(datav)

	rc.stdspec = np.sum((rc.dataim * rc.inmask)[:, rc.x1:rc.x2], 1)
	rc.stdvar = np.sum((rc.varim * rc.inmask)[:, rc.x1:rc.x2], 1)

	return rc.stdspec
