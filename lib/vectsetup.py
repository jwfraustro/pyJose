# TODO DOCS: vectsetup
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

import numpy as np
from lib.excep import ParameterException, VectorLengthException
from lib.procvect import procvect
from lib.misc import plot_procvect
from lib.misc import plot_fitbg

def extrspec(rc):
	# TODO DOCS: extrspec
	"""

	Optimally extracts spectra using weighted profiles.

	Calling Example:

	Inputs:
		dataim: Sky-subtracted, processed image to extract spectrum from.
				Horizontal is pixel position, vertical is wavelength, assuming no curvature.
		profim: image of spatial profiles
		varim:  variance image from processed image
		v_0:    root(v_0) is squared readout noise in DN
		q:      effective number of photons per DN
		x1, x2: boundaries in x which contain spectrum (inclusive)
	Optional Keywords:

	Outputs:
		Returns fopt, an array of optimally extracted spectra, and their variances.
	History:

	Created on 4/17/2021$
	"""
	# TODO FUNC: extrspec

	ny = np.shape(rc.dataim)[0]
	nx = np.shape(rc.dataim)[1]

	if rc.verbose > 1:
		print("Starting Optimal Extraction")

	#TODO fix exceptions here
	if rc.x1 < 0 or rc.x1 > rc.x2:
		raise ParameterException("Exception here.")
	if rc.x2 > nx -1 or rc.x2 < rc.x2:
		raise ParameterException("Exception here.")
	#str += (size(varim ))[1] ne nx || (size(varim ))[2] ne ny ? "varim "  : ""
	#str += (size(profim))[1] ne nx || (size(profim))[2] ne ny ? "profim " : ""
	#str += (size(inmask))[1] ne nx || (size(inmask))[2] ne ny ? "inmask " : ""
	#str += (size(bgim  ))[1] ne nx || (size(bgim  ))[2] ne ny ? "bgim "   : ""
	#str += (size(skyvar))[1] ne nx || (size(skyvar))[2] ne ny ? "skyvar " : ""
	#str += 4 ne n_elements(plottype)                          ? "plottype" : ""

	exres = np.zeros((nx, ny))
	exmask = np.ones((nx, ny))
	crv = np.zeros((rc.x2-rc.x1+1))
	optspec = np.zeros((ny))
	opvar = np.zeros((ny))
	errvect = np.ones((ny))
	func = "extractfunc"

	#LOOP THRROUGH ROWS
	#Cut up the images into rows and pass each to procvect.  Tells
	#procvect to use extractfunc.pro which returns the optimal
	#extraction.  It also return in parm the variance of the extraction.
	#Procvect will handle the sigma rejection.

	for i in range(ny):
		varv = rc.varim[i, rc.x1:rc.x2]
		maskv = rc.inmask[i, rc.x1:rc.x2]
		datav = rc.dataim[i, rc.x1:rc.x2]
		multv = rc.profim[i, rc.x1:rc.x2]
		bgv = rc.bgim[i, rc.x1:rc.x2]
		skyvarv = rc.skyvar[i, rc.x1:rc.x2]

	#TODO Plotting
	#if (i eq gotovect) then verbose = 5
    #if plottype[1] or verbose eq 5 then begin
    #device, window_state = ws   ; get the present window state
    #if not ws[12] then window, 12 else wset, 12 ; open window 12
    #!p.multi = [0, 2, 3, 1, 1]  ; fit all six inside
    #!p.charsize = 2.0
    #!x.margin = [8, 2]
    #!y.margin = [2, 2]
    #plot, datav,   title = 'Data Vector for Extraction', /ystyle
    #plot, multv,   title = 'Profile', /ystyle
    #plot, maskv,   title = 'Input Mask', yrange = [-0.1, 1.1]
    #plot, varv,    title = 'Variance', /ystyle
    #plot, skyvarv, title = 'Sky variance', /ystyle
    #plot, bgv,     title = 'Background', /ystyle
    #!p.multi = 0
    #!p.charsize = 1.0
    #!x.margin = [10, 3]
    #!y.margin = [4, 2]
    #wait, 0.01                  ; pause to allow user to view
    #endif

		rc.optspec[i] = procvect(rc)

		if (rc.errflag):
			errvect[i] = 0
		opvar[i] = parm
		exmask[i, rc.x1:rc.x2] = maskv
		varim[i, rc.x1:rc.x2] = varv
		exres[i, rc.x1:rc.x2] = crv*maskv


	return


def fitprof(rc):
	"""
	Creates an image of normalized spatial profiles, enforces positivity, normalization.
		If needed, can also shift and expand the data image prior to fitting. By default fits a
		polynomial in the spacial direction.

	Inputs:
		dataim: Sky-subtracted processed image.
				Horizontal is pixel, vertical is wavelength.
		spec:   Array holding extracted spectra, same dimension as vertical of dataim.
		x1, x2: Boundaries in x of spectrum (inclusive).
	Optional Keywords:
		General:
			varim:  The variance image for weighting of the polyfit.
					default: array of 1s
			skyvar: The corresponding sky variance image.
					default: array of 0s
			bgim:   The background image to be subtracted off the data.
					default: array of 0s
			inmask: The mask used for all functions.
					default: array of 1s
			q:      The gain of the image.
					default: 1
			v0:     The squared read noise of the data.
					default: 0
			bpct:   The percentage of allowable bad pixels before halting.
					default: 0.5
		Geometry Modification:
			adjfunc:    Function to use to conduct any geometry manipulation.
			adjparms:   The parameters of the geometry for this dataset.
			adjoptions: Any options that can be set for the adjfunc.

		Profile Smoothing Options:
			noproffit:  If set, will not smooth profile image.
			fitgauss:   Set to fit profile in spacial direction with a gaussian.
			fitboxcar:  Set to fit the profile with a boxcar filter.
			boxcarhw:   The width of the median in boxcar average.
						default: 3
			pthresh:    The threshold for sigma rejection.
						default: 3
			profdeg:    Degree of smoothing for profile using a polynomial fit.
						default: 2

		Output from Smooth Profile Fitting:
			profmask:   The mask modified by the fit vector routine
			profres:    The final residuals from the fit
			errvect:    Contains recovered vectors.
			difpmask:   The pixels rejected during the profile fitting routine.

		Debugging:
			verbose:    Level of output to screen.
						default: 0
			plottype:   The type of plot to output while running.
						default: 0
			gotovect:   Stop when this vector is reached.
						default: None
	Outputs:
		Returns an image with the same dimensions as dataim containing the spatial profile,
		optionally smoothed, positivity enforced and normalized.
	"""
	# Set defaults

	dims = np.shape(rc.dataim)
	nx = dims[1]
	ny = dims[0]

	if not np.any(rc.varim):
		rc.varim = np.ones((ny, nx), np.double)
	if not np.any(rc.inmask):
		rc.inmask = np.ones((ny, nx), np.byte)
	if not rc.pthresh:
		rc.pthresh = 3
	if not np.any(rc.bgim):
		rc.bgim = np.ones((ny, nx), np.double)
	if not np.any(rc.skyvar):
		rc.skyvar = np.ones((ny, nx), np.double)
	if not rc.profdeg:
		rc.profdeg = 3
	if not rc.boxcarhw:
		rc.boxcarhw = 3
	if not rc.verbose:
		rc.verbose = 0
	if not rc.plottype:
		rc.plottype = 0
	if not rc.pgotovect:
		rc.pgotovect = -1

	# Check inputs
	if len(rc.spec) != ny:
		raise VectorLengthException("spec", "ny")
	if rc.x1 < 0 or rc.x1 > rc.x2:
		raise ParameterException("x1 must be greater than 0 and less than x2.")
	if rc.x2 > nx - 1 or rc.x2 < rc.x1:
		raise ParameterException("x2 must be less than nx-1 and greater than x1")
	if np.shape(rc.varim)[1] != nx or np.shape(rc.varim)[0] != ny:
		raise ParameterException("varim must be the same dimensions as dataim.")
	if np.shape(rc.inmask)[1] != nx or np.shape(rc.inmask)[0] != ny:
		raise ParameterException("inmask must be the same dimensions as dataim.")
	if np.shape(rc.bgim)[1] != nx or np.shape(rc.bgim)[0] != ny:
		raise ParameterException("bgim must be the same dimensions as dataim.")
	if np.shape(rc.skyvar)[1] != nx or np.shape(rc.skyvar)[0] != ny:
		raise ParameterException("skyvar must be the same dimensions as dataim.")

	if rc.verbose > 1:
		print("Beginning profile fitting.")

	# No Profile Smoothing
	# If profile fitting is not chosen, the data image can be divided by the spectrum at that row.
	# This is not recommended as it is sensitive to the standard extraction's bad pixels.
	# The image is then normalized and made greater than zero everywhere.

	if rc.noproffit:
		rc.specim = np.matmul(np.ones(nx), rc.spec)
		rc.profim = rc.dataim / rc.specim
		rc.profim = rc.profim * rc.inmask > 0
		rc.t = np.sum(rc.profim[:, rc.x1:rc.x2], 1)
		rc.t = np.matmul(np.ones(nx), rc.t)
		rc.profim = rc.profim / rc.t
		rc.profmask = rc.inmask
		rc.difpmask = rc.profmask - rc.profmask
		rc.Perrvect = np.ones(ny, np.byte)
		return rc.profim

	if rc.adjfunc:
		inarray = np.zeros((ny, rc.x2-rc.x2+1, 5))
		inarray[:,:,0] = rc.dataim[:, rc.x1:rc.x2]
		inarray[:, :, 1] = rc.varim[:, rc.x1:rc.x2]
		inarray[:, :, 2] = rc.inmask[:, rc.x1:rc.x2]
		inarray[:, :, 3] = rc.skyvar[:, rc.x1:rc.x2]
		inarray[:, :, 4] = rc.bgim[:, rc.x1:rc.x2]

		# TODO setup adjust function calling
		# outarray = call_function(AdjFunc, InArray, OptInfo, AdjParms = AdjParms, $
		#                            AdjOptions = AdjOptions)
		outarray = None

		rc.pdataim = outarray[:,:,0]
		rc.pvarim = outarray[:, :, 1]
		rc.pmask = outarray[:, :, 2]
		rc.pskyvar = outarray[:, :, 3]
		rc.pbgim = outarray[:, :, 4]

		# Estimate the influence of a bad pixel on its neighbors
		# TODO check this
		rc.pmask = rc.pmask > 0.99 and rc.pmask < 1.01

	#INITIALIZE
	#Because the images may be a different size now, new initilization is
	#needed.  If a Gaussian profile is assumed, transverse by rows.  If
	#polynomial fitting or running average fitting is used loop by columns.

	pnx = np.shape(rc.pdataim)[1]
	pny = np.shape(rc.pdataim)[0]

	if rc.fitgauss == True:
		f1 = 0
		f2 = pny-1
		nvect = pnx
		ndata = pny
		func = "gaussfunc"
		parm = 1
		# TODO
		# this is where the fitting function will need to be moved for the main loop
	else:
		f1 = 0
		f2 = pnx - 1
		nvect = pny
		ndata = pnx
		if rc.fitboxcar == True:
			func = "boxcarfunc"
			parm = rc.boxcarhw
		else:
			func = "polyfunc"
			parm = rc.profdeg

	rc.profmask = rc.pmask
	rc.profres = np.zeros((pnx, pny))
	rc.pprofim = np.zeros((pnx, pny))
	rc.pspecim = np.matmul(np.ones((pnx)), rc.spec)
	rc.datav = np.zeros(nvect)
	rc.maskv = np.zeros(nvect)
	rc.varv = np.zeros(nvect)
	rc.multv = np.zeros(nvect)
	rc.skyvarv = np.zeros(nvect)
	rc.bgv = np.zeros(nvect)
	rc.errvect = np.ones(ndata)
	index = np.arange(ndata+1)
	rc.pprofim = rc.pdataim / rc.pspecim
	yvals = np.arange(pny+1)
	xvals = np.arange(pnx+1)








	#TODO FUNC: rest of fitprof

	return


def fitbg(rc):
	# TODO DOCS: fitbg
	"""
	Creates a sky background image from a data image, interpolated across the spectrum, to create an image of sky lines only.

	Calling Example:

	Inputs:

	Outputs:

	History:

	Created on 4/17/2021$
	"""
	# Fill inputs and verify
	dims = np.shape(rc.data)

	nx = dims[1]
	ny = dims[0]

	if rc.verbose > 1:
		print("Starting background fitting")

	if not rc.inmask:
		rc.inmask = np.ones((nx, ny), np.byte)
	if not rc.varim:
		rc.varim = np.ones((nx, ny), np.single)
	if not rc.skyvar:
		rc.skyvar = np.zeros((nx, ny), np.single)

	if (rc.x2 > nx - 1) or (rc.x2 < rc.x1):
		raise ParameterException("x2 cannot be greater than nx-1 or greater than x1.")
	if (np.shape(rc.varim)[0] != nx) or (np.shape(rc.varim)[1] != ny):
		raise ParameterException("Dimensions of varim do not match dataim.")
	if (np.shape(rc.inmask)[0] != nx) or (np.shape(rc.inmask)[1] != ny):
		raise ParameterException("Dimensions of inmask do not match dataim.")
	if (np.shape(rc.skyvar)[0] != nx) or (np.shape(rc.skyvar)[1] != ny):
		raise ParameterException("Dimensions of skyvar do not match dataim.")

	rc.xvals1 = np.arange(rc.x1)
	rc.xvals2 = np.arange(nx - 1 - rc.x2) + rc.x2 + 1
	rc.xvals = np.array([*rc.xvals1, *rc.xvals2])

	# Subtract bias
	if rc.nobgfit:
		rc.bgim = np.ones((nx, ny)) * np.median(rc.data[:, rc.xvals])
		return rc.bgim

	# Prepare for row by row fitting
	rc.yvals = np.arange(ny)
	rc.bgim = rc.data
	rc.bgres = np.zeros((nx, ny), np.single)
	rc.bgmask = np.ones((nx, ny), np.byte)
	rc.func = "polyfunc"
	rc.berrvect = np.ones(ny)
	rc.bgim[:, rc.x1:rc.x2] = 0
	rc.allx = np.arange(nx)

	# FIT BY ROW
	# Cut up the data into rows and pass each to procvect. Tell procvect
	# to use polyfunc to estimate the background. Procvect will handle
	# bad pixel rejection, and return the polynomial over all x

	for i in range(ny):
		rc.datav = rc.data[i, :]
		rc.maskv = rc.inmask[i,:]
		rc.varv = rc.varim[i,:]
		rc.bcrv = rc.bgres[i,:]
		rc.skyvarv = rc.skyvar[i,:]
		if (sum(rc.maskv[0:rc.x1]) < 2) or (sum(rc.maskv[rc.x2 + 1:nx]) < 2):
			rc.parm = 0
		else:
			rc.parm = rc.bgdeg
		if i == rc.bgotovect:
			rc.verbose = 5


		plot_fitbg(rc.datav, rc.maskv, rc.varv, rc.skyvarv, rc.output_dir)

		if rc.verbose == 5:
			return

		rc.vectnum = i

		rc.bgim[i, :] = procvect(rc)

		if rc.errflag:
			rc.berrvect[i] = 0

		# TODO PLOTTING

		if rc.plottype == 3:

			sy1 = max((i - 10), 0)
			sy2 = min((i + 10), ny - 1)

			plot_procvect(rc.bgim, rc.yvals, rc.allx, sy1, sy2, i)
		# shade_surf, bgim[*, sy1:sy2], allx, yvals[sy1:sy2], charsize = 3,  $
		# title = "Background after # " + strtrim(i, 1), $
		# xtitle = "X pixel location", $
		# ytitle = "Y (Prefit and Postfit)"
		# wait, 0.001
		# endif

		rc.varim[i, rc.xvals] = rc.varv[rc.xvals]
		rc.bgmask[i, :] = rc.maskv
		rc.bgres[i, :] = rc.bcrv * rc.maskv

	return rc.bgim
