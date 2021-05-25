# TODO fitbg docs
"""
Name:
	fitbg
Purpose:
	Creates a sky background image from a data image, interpolated across the spectrum, to create an image of sky lines only.
Category:
	Optimal Spectrum Extraction Package
		- Vector Setup functions
Calling Example:

Inputs:

Outputs:

History:

Created on 4/17/2021$
"""
import numpy as np
from lib.excep import ParameterException
from lib.procvect import procvect

def fitbg(dataim, x1, x2, q, v0, RunConfig):
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
	# TODO

	for i in range(ny):
		datav = dataim[:, i]
		maskv = RunConfig.INMASK[:, i]
		varv = RunConfig.VARIM[:, i]
		bcrv = bgres[:, i]
		skyvarv = RunConfig.SKYVAR[:, i]
		if (sum(maskv[0:x1]) < 2) or (sum(maskv[x2+1:nx]) < 2):
			parm = 0
		else:
			parm = RunConfig.BGDEG
		if i == RunConfig.BGOTOVECT:
			RunConfig.VERBOSE = 5
		# TODO PLOT INFO
		#if RunConfig.PLOTTYPE == 1 or RunConfig.VERBOSE == 5:
			#device, window_state = ws
			#if not ws[12] then window, 12 else wset, 12
			#!p.multi = [0,2,2,1,1]
			#plot, datav, title = 'Data Vector for BG Fitting', / ystyle
			#plot, maskv, title = 'Input Mask', yrange = [0.0, 1.1]
			#plot, varv, title = 'Variance', / ystyle
			#plot, skyvarv, title = 'Sky variance', / ystyle
			#!p.multi = 0
			#wait, 0.01

		if RunConfig.VERBOSE == 5:
			return

		bgim[:, i] = procvect(datav, RunConfig, varv = varv, maskv=maskv, crv=bcrv,
		                      xvals=xvals,q=q, v0=v0, skyvarv=skyvarv, vectnum=i, func=func,
		                      parm=parm)

		#TODO MORE PLOTTING
		#if plottype[3] then begin
		#device, window_state = ws
		#if not ws[14] then window, 14 else wset, 14
		#sy1 = (i - 10) > 0
		#sy2 = (i + 10) < ny - 1
		#shade_surf, bgim[*, sy1:sy2], allx, yvals[sy1:sy2], charsize = 3,  $
		#title = "Background after # " + strtrim(i, 1), $
		#xtitle = "X pixel location", $
		#ytitle = "Y (Prefit and Postfit)"
		#wait, 0.001
		#endif

		RunConfig.VARIM[xvals, i] = varv[xvals]
		bgmask[:, i] = maskv
		bgres[:, i] = bcrv * maskv

	return bgim

	return
