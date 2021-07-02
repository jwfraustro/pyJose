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


def fitprof(RunConfig):
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

	dims = np.shape(RunConfig.dataim)
	nx = dims[1]
	ny = dims[0]

	if not RunConfig.varim:
		varim = np.ones((ny, nx), np.double)
	if not RunConfig.inmask:
		inmask = np.ones((ny, nx), np.byte)
	if not RunConfig.pthresh:
		RunConfig.pthresh = 3
	if not RunConfig.bgim:
		RunConfig.bgim = np.ones((ny, nx), np.double)
	if not RunConfig.skyvar:
		RunConfig.skyvar = np.ones((ny, nx), np.double)
	if not RunConfig.profdeg:
		RunConfig.profdeg = 3
	if not RunConfig.boxcarhw:
		RunConfig.boxcarhw = 3
	if not RunConfig.verbose:
		RunConfig.verbose = 0
	if not RunConfig.plottype:
		RunConfig.plottype = 0
	if not RunConfig.Pgotovect:
		RunConfig.Pgotovect = -1

	# Check inputs
	if len(RunConfig.spec) != ny:
		raise VectorLengthException("spec", "ny")
	if RunConfig.x1 < 0 or RunConfig.x1 > RunConfig.x2:
		raise ParameterException("x1 must be greater than 0 and less than x2.")
	if RunConfig.x2 > nx - 1 or RunConfig.x2 < RunConfig.x1:
		raise ParameterException("x2 must be less than nx-1 and greater than x1")
	if np.shape(RunConfig.varim)[1] != nx or np.shape(RunConfig.varim)[0] != ny:
		raise ParameterException("varim must be the same dimensions as dataim.")
	if np.shape(RunConfig.inmask)[1] != nx or np.shape(RunConfig.inmask)[0] != ny:
		raise ParameterException("inmask must be the same dimensions as dataim.")
	if np.shape(RunConfig.bgim)[1] != nx or np.shape(RunConfig.bgim)[0] != ny:
		raise ParameterException("bgim must be the same dimensions as dataim.")
	if np.shape(RunConfig.skyvar)[1] != nx or np.shape(RunConfig.skyvar)[0] != ny:
		raise ParameterException("skyvar must be the same dimensions as dataim.")

	if RunConfig.verbose > 1:
		print("Beginning profile fitting.")

	# No Profile Smoothing
	# If profile fitting is not chosen, the data image can be divided by the spectrum at that row.
	# This is not recommended as it is sensitive to the standard extraction's bad pixels.
	# The image is then normalized and made greater than zero everywhere.

	if RunConfig.noproffit:
		RunConfig.specim = np.matmul(np.ones(nx), RunConfig.spec)
		RunConfig.profim = RunConfig.dataim / RunConfig.specim
		RunConfig.profim = RunConfig.profim * RunConfig.inmask > 0
		RunConfig.t = np.sum(RunConfig.profim[:, RunConfig.x1:RunConfig.x2], 1)
		RunConfig.t = np.matmul(np.ones(nx), RunConfig.t)
		RunConfig.profim = RunConfig.profim / RunConfig.t
		RunConfig.profmask = RunConfig.inmask
		RunConfig.difpmask = RunConfig.profmask - RunConfig.profmask
		RunConfig.Perrvect = np.ones(ny, np.byte)
		return RunConfig.profim

	#TODO FUNC: rest of fitprof

	return


def fitbg(RunConfig):
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
	dims = np.shape(RunConfig.data)

	nx = dims[0]
	ny = dims[1]

	if RunConfig.verbose > 1:
		print("Starting background fitting")

	if not RunConfig.bgdeg:
		RunConfig.bgdeg = 1
	if not RunConfig.bthresh:
		RunConfig.bthresh = 3.
	if not RunConfig.verbose:
		RunConfig.verbose = 0
	if not RunConfig.plottype:
		RunConfig.plottype = 0
	if not RunConfig.Bgotovect:
		RunConfig.Bgotovect = -1
	if not RunConfig.inmask:
		RunConfig.inmask = np.ones((nx, ny), np.byte)
	if not RunConfig.varim:
		RunConfig.varim = np.ones((nx, ny), np.single)
	if not RunConfig.skyvar:
		RunConfig.skyvar = np.zeros((nx, ny), np.single)

	if (RunConfig.x1 < 0) or (RunConfig.x1 > RunConfig.x2):
		raise ParameterException("x1 must be between 0 and x2.")
	if (RunConfig.x2 > nx - 1) or (RunConfig.x2 < RunConfig.x1):
		raise ParameterException("x2 cannot be greater than nx-1 or greater than x1.")
	if (np.shape(RunConfig.varim)[0] != nx) or (np.shape(RunConfig.varim)[1] != ny):
		raise ParameterException("Dimensions of varim do not match dataim.")
	if (np.shape(RunConfig.inmask)[0] != nx) or (np.shape(RunConfig.inmask)[1] != ny):
		raise ParameterException("Dimensions of inmask do not match dataim.")
	if (np.shape(RunConfig.skyvar)[0] != nx) or (np.shape(RunConfig.skyvar)[1] != ny):
		raise ParameterException("Dimensions of skyvar do not match dataim.")
	if RunConfig.plottype not in range(0, 5):
		raise ParameterException("Plot type must be a value of 0-4")

	RunConfig.xvals1 = np.arange(RunConfig.x1)
	RunConfig.xvals2 = np.arange(nx - 1 - RunConfig.x2) + RunConfig.x2 + 1
	RunConfig.xvals = np.array([*RunConfig.xvals1, *RunConfig.xvals2])

	# Subtract bias
	if RunConfig.nobgfit:
		RunConfig.bgim = np.ones((nx, ny)) * np.median(RunConfig.data[RunConfig.xvals, :])
		return RunConfig.bgim

	# Prepare for row by row fitting
	RunConfig.yvals = np.arange(ny)
	RunConfig.bgim = RunConfig.data
	RunConfig.bgres = np.zeros((nx, ny), np.single)
	RunConfig.bgmask = np.ones((nx, ny), np.byte)
	RunConfig.func = "polyfunc"
	RunConfig.Berrvect = np.ones(ny)
	RunConfig.bgim[RunConfig.x1:RunConfig.x2, :] = 0
	RunConfig.allx = np.arange(nx)

	# FIT BY ROW
	# Cut up the data into rows and pass each to procvect. Tell procvect
	# to use polyfunc to estimate the background. Procvect will handle
	# bad pixel rejection, and return the polynomial over all x

	for i in range(ny):
		RunConfig.datav = RunConfig.data[:, i]
		RunConfig.maskv = RunConfig.inmask[:, i]
		RunConfig.varv = RunConfig.varim[:, i]
		RunConfig.bcrv = RunConfig.bgres[:, i]
		RunConfig.skyvarv = RunConfig.skyvar[:, i]
		if (sum(RunConfig.maskv[0:RunConfig.x1]) < 2) or (sum(RunConfig.maskv[RunConfig.x2 + 1:nx]) < 2):
			RunConfig.parm = 0
		else:
			RunConfig.parm = RunConfig.bgdeg
		if i == RunConfig.Bgotovect:
			RunConfig.verbose = 5
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

		if RunConfig.verbose == 5:
			return
		#Setup pass by ref variables for procvect
		RunConfig.bcrv = RunConfig.crv
		RunConfig.bthresh = RunConfig.thresh
		RunConfig.vectnum = i

		RunConfig.bgim[:, i] = procvect(RunConfig)
		if RunConfig.errflag:
			RunConfig.Berrvect[i] = 0

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

		RunConfig.varim[RunConfig.xvals, i] = RunConfig.varv[RunConfig.xvals]
		RunConfig.bgmask[:, i] = RunConfig.maskv
		RunConfig.bgres[:, i] = RunConfig.bcrv * RunConfig.maskv

	return RunConfig.bgim
