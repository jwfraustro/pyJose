import numpy as np
from lib.vectfit import polyfunc
from lib.excep import *

def procvect(rc):
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
	nx = np.shape(rc.datav)[0]

	if not np.any(rc.xvals):
		rc.xvals = np.arange(nx)
	if not np.any(rc.varv):
		rc.varv = np.ones(nx)
	try:
		if not np.any(rc.multv):
			rc.multv = np.ones(nx)
	except AttributeError:
		rc.multv = np.ones(nx)
	try:
		if not np.any(rc.maskv):
			rc.maskv = np.ones(nx, np.byte)
	except AttributeError:
		rc.maskv = np.ones(nx, np.byte)
	try:
		if not np.any(rc.bgv):
			rc.bgv = np.array(nx)
	except AttributeError:
		rc.bgv = np.array(nx)
	if not np.any(rc.skyvarv):
		rc.skyvarv = np.array(nx)
	if rc.bcrv.shape == 0:
		rc.bcrv = np.array(nx)
	if not rc.bthresh:
		rc.bthresh = 5
	if not rc.q:
		rc.q = 1
	if not rc.v0:
		rc.v0 = 0
	if not rc.bpct:
		rc.bpct = 0.5
	if not rc.func:
		rc.func = 'polyfunc'
	if not rc.verbose:
		rc.verbose = 0
	if not rc.plottype:
		rc.plottype = 0

	# TODO FUNC: procvect error checking

	# str = ""
	if nx != len(rc.varv):
		raise VectorLengthException("nx", "varv")
	if nx != len(rc.varv):
		raise VectorLengthException("nx", "varv")
	if nx != len(rc.multv) and len(rc.multv) != 1:
		raise VectorLengthException("nx", "multv")
	#if nx != len(rc.bgv):
	#	raise VectorLengthException("nx", "bgv")
	#if nx != len(rc.skyvarv):
	#	raise VectorLengthException("nx", "skyvarv")
	if nx != len(rc.bcrv):
		raise VectorLengthException("nx", "bcrv")
	if rc.bthresh < 0:
		raise ParameterException("Threshold cannot be less than 0.")
	if rc.q < 0:
		raise ParameterException("Q cannot be less than 0.")
	if rc.v0 < 0:
		raise ParameterException("v0 cannot be less than 0.")
	if rc.bpct > 1 or rc.bpct < 0:
		raise ParameterException("bpct must be between 0 and 1.")

	# if a row or column number is given, use that for the debug plot title
	if rc.vectnum:
		rc.vectnum_s = " at Vector # " + str(rc.vectnum)
	else:
		rc.vectnum_s = ""

	# if an absolute threshold is used, use that threshold. Else, square
	# the sigma threshold so it can be used with variance calculations
	if rc.absthresh:
		rc.vthresh = rc.bthresh
	else:
		rc.vthresh = rc.bthresh ** 2

	# set the error threshold to be the greatest of 6 pixels, 10% of the total pixels,
	# or the given percentage good of pixels passed in
	rc.errorthresh = len(np.where(rc.xvals==1)) * (1 - rc.bpct) > len(rc.xvals) * 0.10 > 6

	rc.errflag = 0
	rc.coeffv = 0
	rc.funcdone = 1
	rc.funccount = 0

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

	while rc.funcdone:
		rc.funcdone = 0
		rc.goodvals = np.where(rc.maskv[rc.xvals] == 1)
		if len(rc.goodvals) < rc.errorthresh:
			rc.deg = 1
			rc.fiteval = polyfunc(rc)
			if rc.verbose > 2:
				print("Too many pixels rejected" + rc.vectnum_s)
			rc.errflag = 1
			return rc.fiteval
		rc.fitx = rc.xvals[rc.goodvals]  # xvals for good pixels
		rc.fitdata = rc.datav[rc.fitx]  # data for good pixels
		rc.fitvar = rc.varv[rc.fitx]
		rc.fitmult = rc.multv[rc.fitx]
		rc.deg = 0
		rc.est = polyfunc(rc)

		if rc.absthresh:
			rc.bcrv[rc.fitx] = abs(rc.fitdata / rc.fitmult - rc.est)
		else:
			rc.bcrv[rc.fitx] = (rc.fitdata - rc.fitmult * rc.est) ** 2 / (
						rc.fitvar > 1E-6)

		rc.badpix = np.where(rc.bcrv[rc.fitx] > rc.vthresh)

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

		if rc.verbose == 5:
			return
		if not rc.noupdate:
			rc.varv[rc.fitx] = (abs(rc.fitmult * rc.est + rc.bgv[rc.fitx])) / rc.q + rc.v0 + rc.skyvarv[rc.fitx]
		if rc.badpix:
			rc.badx = rc.fitx[rc.badpix]
			rc.maxpos = np.where(rc.bcrv[rc.badx] == max(rc.bcrv[rc.badx]))  # only eliminate max pixel
			rc.maxx = rc.badx[rc.maxpos]
			rc.funccount = rc.funccount + len(rc.maxx)
			rc.maskv[rc.maxx] = 0
			rc.funcdone = 1

	rc.fiteval = polyfunc(rc)

	if not rc.noupdate:
		rc.varv = (abs(rc.multv * rc.fiteval + rc.bgv)) / rc.q + rc.v0 + rc.skyvarv  # get var for all pixels
	return rc.fiteval
