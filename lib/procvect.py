"""
Name:

Purpose:

Category:

Calling Example:

Inputs:

Outputs:

History:

Created on 9/19/2021$
"""

import numpy as np
from lib.excep import *
from importlib import import_module
from lib.fitting import *


def procvect(datav, **kwargs):
	# Set Defaults and Check Inputs

	nx = len(datav)

	# Check if kwarg inputs exist, if not, set their defaults
	xvals = kwargs.get("xvals", np.arange(nx))
	varv = kwargs.get("varv", np.ones(nx))
	multv = kwargs.get("multv", np.ones(nx))
	maskv = kwargs.get("maskv", np.ones(nx, np.byte))
	bgv = kwargs.get("bgv", np.zeros(nx))
	skyvarv = kwargs.get('skyvarv', np.zeros(nx))
	crv = kwargs.get("crv", np.zeros(nx))
	thresh = kwargs.get("thresh", 3)
	q = kwargs.get("q", 1)
	v0 = kwargs.get("v0", 0)
	bpct = kwargs.get("bpct", 0.5)
	func = kwargs.get('func', "polyfunc")
	verbose = kwargs.get("verbose", 0)
	plottype = kwargs.get("plottype", 0)
	parm = kwargs.get("parm", {})

	# Error checking
	if nx != len(varv): raise VectorLengthException("nx", "varv")
	if nx != len(maskv): raise VectorLengthException("nx", "maskv")
	if nx != len(multv) and len(multv) != 1:
		raise VectorLengthException("nx", "multv")
	if nx != bgv.size: raise VectorLengthException("nx", "bgv")
	if nx != len(skyvarv): raise VectorLengthException("nx", "skyvarv")
	if nx != len(crv): raise VectorLengthException("nx", "crv")
	if thresh < 0: raise ParameterException("Threshold cannot be less than 0.")
	if q < 0: raise ParameterException("Q cannot be less than 0.")
	if v0 < 0: raise ParameterException("v0 cannot be less than 0.")
	if bpct > 1 or bpct < 0: raise ParameterException(
			"bpct must be between 0 and 1.")

	if "vectnum" in kwargs:
		vectnum_s = " at Vector # " + str(kwargs.get("vectnum"))
	else:
		vectnum_s = ""

	# If an absolute threshold is used, use that threshold.  Else, square
	# the sigma threshold so it can be used with variance calculations

	absthresh = kwargs.get("absthresh", False)
	if absthresh:
		vthresh = thresh
	else:
		vthresh = thresh * thresh

	# Set the error threshold to be the greatest of 6 pixels, 10% of the total pixels,
	# or the given percentage of good pixels passed in.
	errorthresh = max(len(maskv[xvals] == 1.) * (1 - bpct), len(xvals) * 0.10, 6)

	errflag = 0
	coeffv = 0
	funcdone = 1
	funccount = 0

	noupdate = kwargs.get("noupdate", False)

	# MAIN LOOP
	# On each iteration, first check to make sure there is enough good pixels left
	# to fit the data. Next, pull out the good pixels and pass them to the function.
	# Calculate the residuals of each pixel, using either the difference between
	# actual and estimated, or if sigma rejection is used, the difference and
	# divide by the variance of the pixel. If the user requests a summary plot,
	# plot the estimated variance vs actual and the residual. If needed, update the
	# variance to reflect the new estimation. Reject the pixel with the largest
	# residual larger than the threshold. If no bad pixels are found, exit the loop.

	# Dynamically import the fitting function incase it is user-defined.
	# Function name must be the same as the module name. (i.e. myfunc.myfunc())
	fit_func = getattr(import_module("lib.fitting." + func), func)

	while funcdone:
		funcdone = 0
		goodvals = maskv[xvals] == 1

		if len(goodvals) < errorthresh:

			fiteval = fit_func(np.arange(nx), datav, varv, multv * maskv, True, coeffv, parm)

			if verbose > 2:
				print("Too many pixels rejected" + vectnum_s)
			errflag = 1
			return fiteval, errflag

		fitx = xvals[goodvals]
		fitdata = datav[fitx]
		fitvar = varv[fitx]
		fitmult = multv[fitx]

		est, coeffv = fit_func(fitx, fitdata, fitvar, fitmult, False, coeffv, parm)
		if "absthresh" == True:
			crv[fitx] = abs(fitdata / fitmult - est)
		else:
			crv[fitx] = (fitdata - fitmult * est) ** 2 / fitvar  # prevent divide by zero
		badpix = np.where(crv[fitx] > vthresh)  # get bad locations

		# TODO Procvect plotting

		# 	if plottype[2] or (verbose eq 5) then begin
		#           device, window_state = ws
		#           if not ws[13] then window, 13 else wset, 13
		#           !p.multi = [0,1,2,1,1]
		#           plot, fitx, fitdata, $
		#             title='Actual vs. Fitted' + vectnum_s, $
		#             ytitle='Data Values (if applicable / Spec)', $
		#             xtitle='Pixel Locations', charsize=1.1
		#           oplot, fitx, est * fitmult
		#   if keyword_set(absthresh) then begin
		#     plot, fitx, crv[fitx], title='Residual', $
		#           ytitle='Abs of Data - Expected', $
		#           xtitle='Pixel Locations', charsize=1.1, $
		#           yrange = [0, max([thresh, max(crv[fitx])])]
		#     oplot, fitx, fitdata/fitdata*thresh
		#   endif else begin
		#     plot, fitx, sqrt(crv[fitx]), title='Residuals', $
		#           ytitle='Sigma Difference of Data vs. Expected', $
		#           xtitle='Pixel Locations', charsize=1.1
		#   endelse
		#   !p.multi = 0
		#   wait, 0.01
		# endif

		if not noupdate:
			varv[fitx] = (abs(fitmult * est + bgv[fitx])) / q + v0 + skyvarv[fitx]
		if np.any(badpix):
			badx = fitx[badpix]
			maxpos = np.where(crv[badx] == max(crv[badx]))  # only eliminate max pixel
			maxx = badx[maxpos]
			funccount = funccount + len(maxx)  # add count for bad pixels
			maskv[maxx] = 0  # mask bad pixel
			funcdone = 1  # set so sequence loops

	fiteval, coeffv = fit_func(np.arange(nx), datav, varv, multv*maskv, True, coeffv, parm)

	if not noupdate:
		varv = (abs(multv * fiteval + bgv)) / q + v0 + skyvarv # get var for all pixels

	return fiteval, maskv, errflag, coeffv