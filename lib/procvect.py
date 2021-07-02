import numpy as np
from lib.vectfit import polyfunc
from lib.excep import *

def procvect(RunConfig):
	# TODO DOCS: procvect
	"""
	Name:
		procvect
	Purpose:
		Fit using a user defined function to the input vector,
		integrating until no more cosmic rays are found using either sigma
		or absolute comparison of the residuals rejection.
	Category:
		Optimal Spectrum Extraction Package
			- Vector Fitting driver routine
	Calling Example:

	Inputs:

	Outputs:

	History:

	Created on 4/17/2021$
	"""
	# Set defaults & check inputs
	nx = np.shape(RunConfig.datav)[0]

	if not np.any(RunConfig.xvals):
		RunConfig.xvals = np.arange(nx)
	if not np.any(RunConfig.varv):
		RunConfig.varv = np.ones(nx)
	try:
		if not np.any(RunConfig.multv):
			RunConfig.multv = np.ones(nx)
	except AttributeError:
		RunConfig.multv = np.ones(nx)
	try:
		if not np.any(RunConfig.maskv):
			RunConfig.maskv = np.ones(nx, np.byte)
	except AttributeError:
		RunConfig.maskv = np.ones(nx, np.byte)
	try:
		if not np.any(RunConfig.bgv):
			RunConfig.bgv = np.array(nx)
	except AttributeError:
		RunConfig.bgv = np.array(nx)
	if not np.any(RunConfig.skyvarv):
		RunConfig.skyvarv = np.array(nx)
	if RunConfig.bcrv.shape == 0:
		RunConfig.bcrv = np.array(nx)
	if not RunConfig.bthresh:
		RunConfig.bthresh = 5
	if not RunConfig.q:
		RunConfig.q = 1
	if not RunConfig.v0:
		RunConfig.v0 = 0
	if not RunConfig.bpct:
		RunConfig.bpct = 0.5
	if not RunConfig.func:
		RunConfig.func = 'polyfunc'
	if not RunConfig.verbose:
		RunConfig.verbose = 0
	if not RunConfig.plottype:
		RunConfig.plottype = 0

	# TODO FUNC: procvect error checking

	# str = ""
	if nx != len(RunConfig.varv):
		raise VectorLengthException("nx", "varv")
	if nx != len(RunConfig.varv):
		raise VectorLengthException("nx", "varv")
	if nx != len(RunConfig.multv) and len(RunConfig.multv) != 1:
		raise VectorLengthException("nx", "multv")
	if nx != len(RunConfig.bgv):
		raise VectorLengthException("nx", "bgv")
	if nx != len(RunConfig.skyvarv):
		raise VectorLengthException("nx", "skyvarv")
	if nx != len(RunConfig.bcrv):
		raise VectorLengthException("nx", "bcrv")
	if RunConfig.bthresh < 0:
		raise ParameterException("Threshold cannot be less than 0.")
	if RunConfig.q < 0:
		raise ParameterException("Q cannot be less than 0.")
	if RunConfig.v0 < 0:
		raise ParameterException("v0 cannot be less than 0.")
	if RunConfig.bpct > 1 or RunConfig.bpct < 0:
		raise ParameterException("bpct must be between 0 and 1.")

	# if a row or column number is given, use that for the debug plot title
	if RunConfig.vectnum:
		RunConfig.vectnum_s = " at Vector # " + str(RunConfig.vectnum)
	else:
		RunConfig.vectnum_s = ""

	# if an absolute threshold is used, use that threshold. Else, square
	# the sigma threshold so it can be used with variance calculations
	if RunConfig.absthresh:
		RunConfig.vthresh = RunConfig.bthresh
	else:
		RunConfig.vthresh = RunConfig.bthresh ** 2

	# set the error threshold to be the greatest of 6 pixels, 10% of the total pixels,
	# or the given percentage good of pixels passed in
	RunConfig.error_thresh = len(np.where(xvals=1)) * (1 - RunConfig.bpct) > len(RunConfig.xvals) * 0.10 > 6

	RunConfig.errflag = 0
	RunConfig.coeffv = 0
	RunConfig.funcdone = 1
	RunConfig.funccount = 0

	# MAIN LOOP
	# on each iteration first check to make sure there is enough good
	# pixels left to fit the data.  Next pull out the good pixels and pass
	# them to the function. Calculate the residuals of each pixel, using
	# either the difference between actual and estimated, or if sigma
	# rejection is used square the difference and divide by the variance
	# of the pixel. If the user requests a summary plot, plot the
	# estimated versus actual and the residual.  If needed, update the
	# variance to reflect the new estimation. Reject the pixel with the
	# largest residual larger then the threshold.  If no bad pixels are
	# found, exit the loop.

	while RunConfig.funcdone:
		RunConfig.funcdone = 0
		RunConfig.goodvals = np.where(RunConfig.maskv[RunConfig.xvals] == 1)
		if len(RunConfig.goodvals) < RunConfig.error_thresh:
			RunConfig.fiteval = polyfunc(np.arange(nx), RunConfig.datav, RunConfig.varv,
			                             RunConfig.multv * RunConfig.maskv, 1, RunConfig.coeffv, RunConfig.parm)
			if RunConfig.verbose > 2:
				print("Too many pixels rejected" + RunConfig.vectnum_s)
			RunConfig.errflag = 1
			return RunConfig.fiteval
		RunConfig.fitx = RunConfig.xvals[RunConfig.goodvals]  # xvals for good pixels
		RunConfig.fitdata = RunConfig.datav[RunConfig.fitx]  # data for good pixels
		RunConfig.fitvar = RunConfig.varv[RunConfig.fitx]
		RunConfig.fitmult = RunConfig.multv[RunConfig.fitx]
		RunConfig.est = polyfunc(RunConfig.fitx, RunConfig.fitdata, RunConfig.fitvar, RunConfig.fitmult, 0,
		                         RunConfig.coeffv, RunConfig.parm)

		if RunConfig.absthresh:
			RunConfig.bcrv[RunConfig.fitx] = abs(RunConfig.fitdata / RunConfig.fitmult - RunConfig.est)
		else:
			RunConfig.bcrv[RunConfig.fitx] = (RunConfig.fitdata - RunConfig.fitmult * RunConfig.est) ** 2 / (
						RunConfig.fitvar > 1E-6)

		RunConfig.badpix = np.where(RunConfig.bcrv[RunConfig.fitx] > RunConfig.vthresh)

		# TODO FUNC: procvect plotting

		# if plottype[2] or (verbose eq 5) then begin
		# device, window_state = ws
		#    if not ws[13] then window, 13 else wset, 13
		#    !p.multi = [0,1,2,1,1]
		#    plot, fitx, fitdata, $
		#         title='Actual vs. Fitted' + vectnum_s, $
		#          ytitle='Data Values (if applicable / Spec)', $
		#          xtitle='Pixel Locations', charsize=1.1
		#    oplot, fitx, est * fitmult
		# if keyword_set(absthresh) then begin
		#   plot, fitx, crv[fitx], title='Residual', $
		#     ytitle='Abs of Data - Expected', $
		#     xtitle='Pixel Locations', charsize=1.1, $
		#    yrange = [0, max([thresh, max(crv[fitx])])]
		# oplot, fitx, fitdata/fitdata*thresh
		# endif else begin
		#  plot, fitx, sqrt(crv[fitx]), title='Residuals', $
		#        ytitle='Sigma Difference of Data vs. Expected', $
		#        xtitle='Pixel Locations', charsize=1.1
		# endelse
		# !p.multi = 0
		#    wait, 0.01
		# endif

		if RunConfig.verbose == 5:
			return
		if not RunConfig.noupdate:
			RunConfig.varv[RunConfig.fitx] = (abs(RunConfig.fitmult * RunConfig.est + RunConfig.bgv[RunConfig.fitx])) / RunConfig.q + RunConfig.v0 + RunConfig.skyvarv[RunConfig.fitx]
		if RunConfig.badpix:
			RunConfig.badx = RunConfig.fitx[RunConfig.badpix]
			RunConfig.maxpos = np.where(RunConfig.crv[RunConfig.badx] == max(RunConfig.crv[RunConfig.badx]))  # only eliminate max pixel
			RunConfig.maxx = RunConfig.badx[RunConfig.maxpos]
			RunConfig.funccount = RunConfig.funccount + len(RunConfig.maxx)
			RunConfig.maskv[RunConfig.maxx] = 0
			RunConfig.funcdone = 1

	RunConfig.fiteval = polyfunc(np.arange(nx), RunConfig.datav, RunConfig.varv, RunConfig.multv * RunConfig.maskv, 1, RunConfig.coeffv, RunConfig.parm)

	if not RunConfig.noupdate:
		RunConfig.varv = (abs(RunConfig.multv * RunConfig.fiteval + RunConfig.bgv)) / RunConfig.q + RunConfig.v0 + RunConfig.skyvarv  # get var for all pixels
	return RunConfig.fiteval
