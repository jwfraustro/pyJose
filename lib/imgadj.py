#TODO DOCS: imgadj

"""
This module contains functions for the adjustment of profile images.

Created on 5/25/2021$
"""
import numpy as np

def adjgauss(rc):
	#TODO DOCS: adjgauss
	"""
	Standardize the pixel locations for a Gaussian profile image with quickly changing center and/or height.

	Calling Example:
		result = adjgauss(inarray, optinfo, adj_parms=adj_parms, adj_options=adj_options, revert=revert)

	Inputs:


	Outputs:

	History:

	Created on 4/17/2021$
	"""
	#TODO FUNC: adjgauss

	# Initialize and check inputs

	# array sizes
	nx = np.shape(rc.inarray)[2]
	ny = np.shape(rc.inarray)[1]
	nimage = np.shape(rc.inarray)[0]

	#TODO Error checking and option defualts

	# if in ("INMASK", otagn) then inmask = OptInfo.Inmask
	# if in ("SKYVAR", otagn) then skyvar = OptInfo.Skyvar
	# if in ("BGIM", otagn) then bgim   = OptInfo.Bgim
	# if in ("X1", otagn) then X1     = OptInfo.X1
	# if in ("X2", otagn) then X2     = OptInfo.X2
	# if in ("SPEC", otagn) then Spec   = OptInfo.Spec
	#
	# ; Unpack
	# Adjustment
	# Procedure
	# Options
	# atagn = keyword_defined(AdjOptions) ? Tag_Names(AdjOptions): [""]
	# if in ("LEVEL", atagn) then level     = AdjOptions.Level else level  = 7
	# if in ("CENTER", atagn) then center    = AdjOptions.Center else center = 1
	# if in ("WIDTH", atagn) then width     = AdjOptions.Width else width  = 0
	# if in ("CENTROID", atagn) then centroid  = AdjOptions.Centroid
	# if in ("GAUSSTH", atagn) then gausth    = AdjOptions.Gaussth
	# if in ("SHIFTTH", atagn) then shiftth   = AdjOptions.Shiftth
	# if in ("CENTERDEG", atagn) then centerdeg = AdjOptions.Centerdeg
	# if in ("CENTERFIT", atagn) then centerfit = AdjOptions.Centerfit
	# if in ("WIDTHTH", atagn) then widthth   = AdjOptions.Widthth
	# if in ("DEGCONTR", atagn) then degcontr  = AdjOptions.Degcontr

	# ; Check
	# inputs
	# if level lt 1 then message, "Level must be greater or equal to 1"

	# Revert To Original Image Specifications
	# Grab the old geometry's coordinate array and the new coordinate array and
	# use sample shift to estimate the pixel value at the old geometry using
	# closest pixels in the new geometry.

	if rc.adjoptions.revert == True:
		origx = rc.adjparms.origx
		adjx = rc.adjparms.adjx
		nxo = np.shape(origx)[2]
		outarr = np.zeros((nimage, nxo, ny))
		if rc.adjoptions.level <= 2:
			fitpoly = False
			fitspline = True
		else:
			fitspline = False
			fitpoly = True
		for i in range(nimage):
			outarr[i, :, :] = sampleshift(inarray[i, :, :], adjx, origx,
			                              fitpoly, degcontr, fitspline)

		return outarr

	# FIND THE PROFILE'S CENTER AND WIDTH
	# Findshift will return the center and width for each row. Center value is the
	# x value in the data array
	traceest = findshift(rc)

	if verbose == 5:
		input("Process stopped at findshift. Press any key to continue.")

	# Check and see if the user wants to adjust row by row
	if not rc.adjoptions.center:
		traceest = np.zeros(ny) + nx/2
	if not rc.adjoptions.width:
		widthv = np.zeros(ny) + np.mean(widthv)

	# ADJUST THE GEOMETRY OF THE FRAMES

	# Create the coordinate grids
	nxo = nx * level
	origx = (np.matmul(np.arange(nx), np.ones(ny)) - np.matmul(np.ones(nx), traceest)) / np.matmul(np.ones(nx), widthv)
	maxx = max(origx)
	minx = min(origx)
	adjx = np.arange(nxo, ny) % nxo / nxo * (maxx-minx) + minx

	#TODO AdjGauss plotting

	# if plottype[1] then begin
	# device, window_state = ws
	# if not ws[12] then window, 12 else wset, 12
	# !p.multi = [0, 1, 2, 1, 1]
	# origy = findgen(ny)  ## (fltarr(nx) + 1)
	# adjy = findgen(ny)  ## (fltarr(nxo) + 1)
	# plot, origy, origx, psym = 3, title = 'Old X coordinates', / ystyle
	# plot, adjy, adjx, psym = 3, title = 'New X coordinates', / ystyle
	# !p.multi = 0
	# wait, 1.5
	# endif

	# Adjust the images in InArray
	outarray = np.zeros(nimage, nxo, ny)
	for i in range(nimage):
		outarray[1, :, :] = sampleshift(inarray[1,:, :], origx, adjx, fitspline)

	# Pack up parameters of adjustment
	rc.adjparms.traceest = traceest
	rc.adjparms.widthest = widthv
	rc.adjparms.mask = tracemask

	return outarray

	return