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


def extrspec(dataim, profim, varim, v0, q, x1, x2, **kwargs):
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
	return


def fitprof(dataim, spec, x1, x2,
            adjfunc=None,
            adjparms=None,
            adjoptions=None,
            bgim=None,
            boxcarhw=None,
            bpct=None,
            difpmask=None,
            errvect=None,
            fitboxcar=None,
            fitgauss=None,
            gotovect=None,
            inmask=None,
            noproffit=None,
            plottype=None,
            profdeg=None,
            profmask=None,
            profres=None,
            pthresh=None,
            q=None,
            skyvar=None,
            varim=None,
            verbose=None,
            v0=None):
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

	dims = np.shape(dataim)
	nx = dims[1]
	ny = dims[0]

	# FIXME: This block could be handled by setting defaults above. Test later.
	if not varim:
		varim = np.ones((ny, nx), np.double)
	if not inmask:
		inmask = np.ones((ny, nx), np.byte)
	if not pthresh:
		pthresh = 3
	if not bgim:
		bgim = np.ones((ny, nx), np.double)
	if not skyvar:
		skyvar = np.ones((ny, nx), np.double)
	if not profdeg:
		profdeg = 3
	if not boxcarhw:
		boxcarhw = 3
	if not verbose:
		verbose = 0
	if not plottype:
		plottype = 0
	if not gotovect:
		gotovect = -1

	# Check inputs
	if len(spec) != ny:
		raise VectorLengthException("spec", "ny")
	if x1 < 0 or x1 > x2:
		raise ParameterException("x1 must be greater than 0 and less than x2.")
	if x2 > nx - 1 or x2 < x1:
		raise ParameterException("x2 must be less than nx-1 and greater than x1")
	if np.shape(varim)[1] != nx or np.shape(varim)[0] != ny:
		raise ParameterException("varim must be the same dimensions as dataim.")
	if np.shape(inmask)[1] != nx or np.shape(inmask)[0] != ny:
		raise ParameterException("inmask must be the same dimensions as dataim.")
	if np.shape(bgim)[1] != nx or np.shape(bgim)[0] != ny:
		raise ParameterException("bgim must be the same dimensions as dataim.")
	if np.shape(skyvar)[1] != nx or np.shape(skyvar)[0] != ny:
		raise ParameterException("skyvar must be the same dimensions as dataim.")

	if verbose > 1:
		print("Beginning profile fitting.")

	# No Profile Smoothing
	# If profile fitting is not chosen, the data image can be divided by the spectrum at that row.
	# This is not recommended as it is sensitive to the standard extraction's bad pixels.
	# The image is then normalized and made greater than zero everywhere.

	if noproffit:
		specim = np.matmul(np.ones(nx), spec)
		profim = dataim / specim
		profim = profim * inmask > 0
		t = np.sum(profim[:, x1:x2], 1)
		t = np.matmul(np.ones(nx), t)
		profim = profim / t
		profmask = inmask
		difpmask = profmask - profmask  # FIXME: ??? what is this?
		errvect = np.ones(ny, np.byte)
		return profim

	#TODO FUNC: rest of fitprof

	return


def fitbg(dataim, x1, x2, q, v0, RunConfig):
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
	dims = np.shape(dataim)

	nx = dims[0]
	ny = dims[1]

	if RunConfig.VERBOSE > 1:
		print("Starting background fitting")

	if not RunConfig.BGDEG:
		RunConfig.BGDEG = 1
	if not RunConfig.BTHRESH:
		RunConfig.BTHRESH = 3.
	if not RunConfig.VERBOSE:
		RunConfig.VERBOSE = 0
	if not RunConfig.PLOTTYPE:
		RunConfig.PLOTTYPE = 0
	if not RunConfig.BGOTOVECT:
		RunConfig.BGOTOVECT = -1
	if not RunConfig.INMASK:
		RunConfig.INMASK = np.ones((nx, ny), np.byte)
	if not RunConfig.VARIM:
		RunConfig.VARIM = np.ones((nx, ny), np.single)
	if not RunConfig.SKYVAR:
		RunConfig.SKYVAR = np.zeros((nx, ny), np.single)

	if (x1 < 0) or (x1 > x2):
		raise ParameterException("x1 must be between 0 and x2.")
	if (x2 > nx - 1) or (x2 < x1):
		raise ParameterException("x2 cannot be greater than nx-1 or greater than x1.")
	if (np.shape(RunConfig.VARIM)[0] != nx) or (np.shape(RunConfig.VARIM)[1] != ny):
		raise ParameterException("Dimensions of varim do not match dataim.")
	if (np.shape(RunConfig.INMASK)[0] != nx) or (np.shape(RunConfig.INMASK)[1] != ny):
		raise ParameterException("Dimensions of inmask do not match dataim.")
	if (np.shape(RunConfig.SKYVAR)[0] != nx) or (np.shape(RunConfig.SKYVAR)[1] != ny):
		raise ParameterException("Dimensions of skyvar do not match dataim.")
	if RunConfig.PLOTTYPE not in range(0, 5):
		raise ParameterException("Plot type must be a value of 0-4")

	xvals1 = np.arange(x1)
	xvals2 = np.arange(nx - 1 - x2) + x2 + 1
	xvals = [xvals1, xvals2]

	# Subtract bias
	if RunConfig.NOBGFIT:
		bgim = np.ones((nx, ny)) * np.median(dataim[xvals, :])
		return bgim

	# Prepare for row by row fitting
	yvals = np.arange(ny)
	bgim = dataim
	bgres = np.zeros((nx, ny), np.single)
	bgmask = np.ones((nx, ny), np.byte)
	func = "polyfunc"
	errvect = np.ones(ny)
	bgim[x1:x2, :] = 0
	allx = np.arange(nx)

	# FIT BY ROW
	# Cut up the data into rows and pass each to procvect. Tell procvect
	# to use polyfunc to estimate the background. Procvect will handle
	# bad pixel rejection, and return the polynomial over all x

	for i in range(ny):
		datav = dataim[:, i]
		maskv = RunConfig.INMASK[:, i]
		varv = RunConfig.VARIM[:, i]
		bcrv = bgres[:, i]
		skyvarv = RunConfig.SKYVAR[:, i]
		if (sum(maskv[0:x1]) < 2) or (sum(maskv[x2 + 1:nx]) < 2):
			parm = 0
		else:
			parm = RunConfig.BGDEG
		if i == RunConfig.BGOTOVECT:
			RunConfig.VERBOSE = 5
		# TODO PLOTTING
		# if RunConfig.PLOTTYPE == 1 or RunConfig.VERBOSE == 5:
		# device, window_state = ws
		# if not ws[12] then window, 12 else wset, 12
		# !p.multi = [0,2,2,1,1]
		# plot, datav, title = 'Data Vector for BG Fitting', / ystyle
		# plot, maskv, title = 'Input Mask', yrange = [0.0, 1.1]
		# plot, varv, title = 'Variance', / ystyle
		# plot, skyvarv, title = 'Sky variance', / ystyle
		# !p.multi = 0
		# wait, 0.01

		if RunConfig.VERBOSE == 5:
			return

		bgim[:, i] = procvect(datav, RunConfig, varv=varv, maskv=maskv, crv=bcrv,
		                      xvals=xvals, q=q, v0=v0, skyvarv=skyvarv, vectnum=i, func=func,
		                      parm=parm)

		# TODO PLOTTING
		# if plottype[3] then begin
		# device, window_state = ws
		# if not ws[14] then window, 14 else wset, 14
		# sy1 = (i - 10) > 0
		# sy2 = (i + 10) < ny - 1
		# shade_surf, bgim[*, sy1:sy2], allx, yvals[sy1:sy2], charsize = 3,  $
		# title = "Background after # " + strtrim(i, 1), $
		# xtitle = "X pixel location", $
		# ytitle = "Y (Prefit and Postfit)"
		# wait, 0.001
		# endif

		RunConfig.VARIM[xvals, i] = varv[xvals]
		bgmask[:, i] = maskv
		bgres[:, i] = bcrv * maskv

	return bgim
