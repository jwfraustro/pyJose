"""
Name: fitbg.py

Purpose: Creates a sky background image from a data image, interpolated across the spectrum, to create an image of 
sky lines only.

Category:

Calling Example:

Inputs:
	dataim: The data image, with spectrum
	x1, x2: boundaries in x which contain spectrum (spectrum is inclusive of x1,x2)

Optional Inputs:
	nobgfit: set True to not fit a background and assume the average
	
	Cut into rows and passed to procvect.py:
	inmask: the main mask used by all functions
	varim: the array that contains the variance of the locations
	skyvar: the variance image for the sky subtraction
	
	Passed through to procvect.py:
	bgdeg: degree of polynomial interpolation (default is 1)
	bthresh: the threshold for sigma rejection
	q: the electron per adu of the data
	v0: Read-noise squared of the data
	bpct: the percentage of bad pixels that kicks out of an interation
	
	Made from procvect.py:
	bgres: the residuals for cosmic ray rejection
	bgmask: the output mask of the cosmic rays found
	
	Debugging:
	verbose: level of printed output
	plottype: level of plot to show (0-4)
	gotovect: the row at which to stop the loop
	errvect: will contain the rows that exited with bad pixels

Outputs:
	an array of size dataim in which, for each wavelength, the spectrum
	from x1 to x2 has been removed, interpolated over (with coefficients coeffs) 
	and polynomial has been evaluated at all x values

History:

Created on 9/19/2021$
"""

import numpy as np
from lib.excep import *
from lib.misc import *
from lib.procvect import procvect
from lib.utils import save_fits


def fitbg(dataim, x1, x2, **kwargs):
	dims = np.shape(dataim)

	nx = dims[1]
	ny = dims[0]

	if "verbose" in kwargs:
		if kwargs['verbose'] > 1:
			print("Starting background fitting")

	bgdeg = kwargs.get("bgdeg", 1)
	bthresh = kwargs.get("bthresh", 5)
	verbose = kwargs.get("verbose", 0)
	plottype = kwargs.get("plottype", 0)
	gotovect = kwargs.get("gotovect", -1)
	inmask = kwargs.get("inmask", np.ones((nx, ny)))
	varim = kwargs.get("varim", np.ones((nx, ny)))
	skyvar = kwargs.get("skyvar", np.zeros((nx, ny)))

	if (x2 > nx - 1) or (x2 < x1):
		raise ParameterException("x2 cannot be greater than nx-1 or greater than x1.")
	if (np.shape(varim)[0] != nx) or (np.shape(varim)[1] != ny):
		raise ParameterException("Dimensions of varim do not match dataim.")
	if (np.shape(inmask)[0] != nx) or (np.shape(inmask)[1] != ny):
		raise ParameterException("Dimensions of inmask do not match dataim.")
	if (np.shape(skyvar)[0] != nx) or (np.shape(skyvar)[1] != ny):
		raise ParameterException("Dimensions of skyvar do not match dataim.")

	xvals1 = np.arange(x1)
	xvals2 = np.arange(nx - 1 - x2) + x2 + 1
	xvals = np.array([*xvals1, *xvals2])

	# Subtract bias
	nobgfit = kwargs.get("nobgfit", False)

	if nobgfit:
		bgim = (np.ones((nx, ny))) * np.median(dataim[:, xvals])
		return bgim

	# Prepare for row by row fitting
	yvals = np.arange(ny)
	bgim = np.copy(dataim)
	bgres = np.zeros((nx, ny), np.single)
	bgmask = np.ones((nx, ny), np.byte)
	errvect = np.ones(ny)
	bgim[:, x1:x2] = 0
	allx = np.arange(nx)

	# FIT BY ROW
	# Cut up the data into rows and pass each to procvect. Tell procvect
	# to use polyfunc to estimate the background. Procvect will handle
	# bad pixel rejection, and return the polynomial over all x

	for i in range(ny):
		print("Fitting row ", i)
		datav = dataim[i, :]
		maskv = inmask[i, :]
		varv = varim[i, :]
		bcrv = bgres[i, :]
		skyvarv = skyvar[i, :]

		if (sum(maskv[0:x1]) < 2) or (sum(maskv[x2 + 1:nx]) < 2):
			parm = 0
		else:
			parm = bgdeg

		if i == gotovect:
			verbose = 5

		# plot_fitbg(datav, maskv, varv, skyvarv, kwargs['output_dir'])

		bgim[i, :], errflag = procvect(datav, varv=varv, maskv=maskv, crv=bcrv, thresh=bthresh, xvals=xvals,
		                               skyvarv=skyvarv, vectnum=i, parm=parm, **kwargs)

		if errflag:
			errvect[i] = 0  # There was a problem fitting this row

	return bgim, varim
