import numpy as np
from lib.vectfit import polyfunc


def procvect(datav, RunConfig, varv, maskv, crv, xvals, q, v0, skyvarv, vectnum, func, parm, multv=None, bgv=None):
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
	nx = np.shape(datav)[0]

	if not xvals:
		xvals = np.arange(nx)
	if not varv:
		varv = np.ones(nx)
	if not multv:
		multv = np.ones(nx)
	if not maskv:
		maskv = np.ones(nx, np.byte)
	if not bgv:
		bgv = np.array(nx)
	if not skyvarv:
		skyvarv = np.array(nx)
	if not crv:
		crv = np.array(nx)
	if not RunConfig.BTHRESH:
		RunConfig.BTHRESH = 5
	if not q:
		q = 1
	if not v0:
		v0 = 0
	if not RunConfig.BPCT:
		RunConfig.BPCT = 0.5
	if not func:
		func = 'polyfunc'
	if not RunConfig.VERBOSE:
		RunConfig.VERBOSE = 0
	if not RunConfig.PLOTTYPE:
		RunConfig.PLOTTYPE = 0

	# TODO FUNC: procvect error checking

	# str = ""
	# str += nx ne n_elements(varv)    ? "varv"    : ""
	# str += nx ne n_elements(maskv)   ? "maskv"   : ""
	# str += nx ne n_elements(multv) $
	#  && 1 ne n_elements(multv)   ? "multv"   : ""
	# str += nx ne n_elements(bgv)     ? "bgv"     : ""
	# str += nx ne n_elements(skyvarv) ? "skyvarv" : ""
	# str += nx ne n_elements(crv)     ? "crv"     : ""
	# str += 4 ne n_elements(plottype) ? "plottype must have 4 elements" : ""
	# str += thresh le 0               ? "thresh less than 0" : ""
	# str += Q le 0                    ? "Q is less or equal to 0" : ""
	# str += v0 lt 0                   ? "v0 is less than 0" : ""
	# str += bpct gt 1 || bpct lt 0    ? "bpct is not betwee 0 and 1" : ""
	# if str ne "" then message, "Poorly formed inputs: " + str

	# if a row or column number is given, use that for the debug plot title
	if vectnum:
		vectnum_s = " at Vector # " + str(vectnum)
	else:
		vectnum_s = ""

	# if an absolute threshold is used, use that threshold. Else, square
	# the sigma threshold so it can be used with variance calculations
	# FIXME: optspecextr_source variable name changes
	# if RunConfig.BTHRESH:
	# 	vthresh = RunConfig.BTHRESH
	# else:
	# 	thresh =

	# set the error threshold to be the greatest of 6 pixels, 10% of the total pixels,
	# or the given percentage good of pixels passed in
	error_thresh = len(np.where(xvals=1)) * (1 - RunConfig.BPCT) > len(xvals) * 0.10 > 6

	errflag = 0
	coeffv = 0
	funcdone = 1
	funccount = 0

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

	while funcdone:
		funcdone = 0
		goodvals = np.where(maskv[xvals] == 1)
		if len(goodvals) < error_thresh:
			fiteval = polyfunc(np.arange(nx), datav, varv, multv * maskv, 1, coeffv, parm)
			if RunConfig.VERBOSE > 2:
				print("Too many pixels rejected" + vectnum_s)
			errflag = 1
			return fiteval
		fitx = xvals[goodvals]  # xvals for good pixels
		fitdata = datav[fitx]  # data for good pixels
		fitvar = varv[fitx]
		fitmult = multv[fitx]
		est = polyfunc(fitx, fitdata, fitvar, fitmult, 0, coeffv, parm)
		# FIXME procvect absthresh

		# if absthresh:
		#   crv[fitx] = abs(fitdata / fitmult - est)
		# else:
		#   crv[fitx] = (fitdata - fitmult * est)**2 / (fitvar > 1E-6)

		badpix = np.where(crv[fitx] > vthresh)

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

		if RunConfig.VERBOSE == 5:
			return
		if not noupdate:
			varv[fitx] = (abs(fitmult * est + bgv[fitx])) / q + v0 + skyvarv[fitx]
		if badpix:
			badx = fitx[badpix]
			maxpos = np.where(crv[badx] == max(crv[badx]))  # only eliminate max pixel
			maxx = badx[maxpos]
			funccount = funccount + len(maxx)
			maskv[maxx] = 0
			funcdone = 1

	fiteval = polyfunc(np.arange(nx), datav, varv, multv * maskv, 1, coeffv, parm)

	if not noupdate:
		varv = (abs(multv * fiteval + bgv)) / q + v0 + skyvarv  # get var for all pixels
	return fiteval
