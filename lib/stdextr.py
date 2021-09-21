"""
Name: stdextr

Purpose: Standard box extraction of spectrum, returns spectrum and variance image
	and estimates over bad pixels.

Category:
	Optimal Spectrum Extraction Package
		- Horne Step Functions

Calling Example:

Inputs:
	dataim: sky-subtracted, processed image to extract spectrum from. Horizontal is pixel position, vertical is 
			wavelength, assumes no curavature.
	varim: variance image from processed image
	x1,x2: boundaries in x over which spectrum will be summed
	
	Optional Inputs:
	inmask: the cosmic ray mask for the data image
	adjspec: set to linearly interpolate the data, then extract the standard spectrum. Only useful if input mask has
			bad pixels.
	

Outputs:

	Returns stdspec, an array of length n containing extracted spectra, where n is the number of wavelengths of the 
	image (vertical direction).
	stdvar: variance of extracted spectra, length n where n is the number of wavelengths taken of the image (vertical direction)

	
Restrictions:
		dataim, varim, and inmask must be the same size. 0<=x1<=x2<=nx-1
		
History:

Created on 9/19/2021$
"""

from lib.excep import ParameterException
import numpy as np
from scipy import interpolate

def stdextr(dataim, varim, x1, x2, **kwargs):

	# Check Inputs
	ny = np.shape(dataim)[0]
	nx = np.shape(dataim)[1]

	inmask = kwargs.get("inmask", np.ones((nx, ny)))
	varim = kwargs.get("varim", np.ones((nx, ny)))


	if x1 < 0 or x1 > x2:
		raise ParameterException("x1 must be greater than 0 and less than x2.")
	if x2 > nx-1 or x2 < x1:
		raise ParameterException("x2 must be less than nx-1 and greater than x1.")
	if np.shape(varim)[1] != nx or np.shape(varim)[0] != ny:
		raise ParameterException("varim must be the same shape as dataim.")
	if np.shape(inmask)[1] != nx or np.shape(inmask)[0] != ny:
		raise ParameterException("inmask must be the same shape as dataim.")

	# Interpolate over bad pixels
	adjspec = kwargs.get("adjspec", False)
	if adjspec:
		adjspec = np.array(ny, np.double)
		for i in range(ny):
			bv = np.where(inmask[i, x1:x2] == 0)
			gv = np.where(inmask[i, x1:x2] == 1)
			datav = dataim[i, x1:x2]
			if len(bv) > 0:
				interpfunc = interpolate.interp1d(gv, datav[gv],  kind='linear')
				datav[bv] = interpfunc(bv)
			adjspec[i] = np.sum(datav)

	stdspec = np.sum((dataim * inmask)[:, x1:x2], 1)
	stdvar = np.sum((varim * inmask)[:, x1:x2], 1)

	return stdspec, stdvar, adjspec