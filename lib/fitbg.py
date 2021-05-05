# TODO
"""
Name:

Purpose:

Category:

Calling Example:

Inputs:

Outputs:

History:

Created on 4/17/2021$
"""
import numpy as np
from lib.excep import ParameterException


def fitbg(dataim, x1, x2, q, v0, RunConfig):
	# Fill inputs and verify
	dims = np.shape(dataim)

	nx = dims[1]
	ny = dims[2]

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
	if not RunConfig.GOTOVECT:
		RunConfig.GOTOVECT = -1
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
	if (np.shape(RunConfig.VARIM)[1] != nx) or (np.shape(RunConfig.VARIM)[2] != ny):
		raise ParameterException("Dimensions of varim do not match dataim.")
	if (np.shape(RunConfig.INMASK)[1] != nx) or (np.shape(RunConfig.INMASK)[2] != ny):
		raise ParameterException("Dimensions of inmask do not match dataim.")
	if (np.shape(RunConfig.SKYVAR)[1] != nx) or (np.shape(RunConfig.SKYVAR)[2] != ny):
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

	return
