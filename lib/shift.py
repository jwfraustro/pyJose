#TODO DOCS: shift
"""
A collection of functions for vector shifting and sampling.

Calling Example:

Inputs:

Outputs:

History:

Created on 5/25/2021$
"""

from lib.misc import rebin_nd
from lib.vectfit import centermass, gaussfunc
from lib import procvect
from lib.excep import *
import numpy as np

def sampleshift(rc):
	# TODO DOCS: sampleshift
	"""

	Shifts an array according to an input vector, then shrinks it.

	Category:
		Optimal Spectrum Extraction Package
			- Shifting Procedures

	Calling Example:
		result = sampleshift(data, inx, outx)

	Inputs:
		data:   The array to be modified.
		inx:    Array of x positions in the reference frame for each pixel in data.
		outx:   The x position values for the output array. If defined, shift is ignored. May be modified to
				reflect the true positions.

	Keyword Parameters:
		fitsample:  Set to sample the data when changing shape.
		fitinterp:  Set to use linear interpolation.
		fitspline:  Set to use a spline when upsizing.
		fitaverage: Set to use nearest neighbor sampling when downsizing.
		fitpoly:    Set to use polynomial fitting when downsizing.
		fitcubic:   Set to use a cubic convolution.
		degcontr:   Degree of fit to use when downsizing.

	Outputs:
		Returns to data array that has been expanded or contracted and/or shifted.

	History:

	Created on 4/17/2021$
	"""
	#TODO FUNC: sampleshift
	nx, ny = np.shape(rc.data)  # (size(data))[1]
	nxo = np.shape(rc.outx)[0] # (size(outx))[1]
	nxi = np.shape(rc.inx)[0]  # (size(inx))[1]
	
	if not np.any(rc.degcontr):
		rc.degcontr = 2
	if rc.degcontr < 0:
		raise ParameterException("Degree of fit must be > 0.")
	if nx != nxi or ny != np.shape(rc.inx)[1]:
		raise ParameterException("inx must be equal in size to data")
	if ny != np.shape(rc.outx)[1]:
		raise ParameterException("outx must have equal size y dimension to data") 

	# TODO INITIALIZE
	# make a large array of differences between the input coordinate and
	# output coordinate for each row.
	if rc.fitsample or rc.fitinterp or rc.fitaverage:
		# TODO Fix orientation of array on transpose
		inxa = rebin_nd(rc.inx, (nxi, ny, nxo))
		outxa = rebin_nd(rc.outx, (nxo, ny, nxi))
		inxa = np.transpose(inxa, (2, 1, 0))
		dif = (outxa - inxa)


	# SAMPLE
	# Just return the pixel with the closest coordiante to whats needed
	if rc.fitsample:
		foo = np.amin(abs(dif), rc.rmap, axis=3)
		rc.rmap /= nxo
		mdata = rc.data[rc.rmap]


	# ; INTERPOLATE
	# ; Interpolate using the closest coordiantes above and below.  Find the
	# ; pixels above and below, and add them together weighted by relative
	# ; distance from the needed coordinate.
	if rc.fitinterp:
	#   foo = np.amin(1. / dif, lmap, dimension=3, subscript_max=hmap)
		lweight = dif[rc.hmap] / (dif[rc.hmap] - dif[rc.lmap]) # the weight
		lmap = rc.lmap / nxo             # the map of old pixels to new
		hmap = hmap / nxo
		mdata = rc.data[hmap] * (1 - lweight)
		mdata += rc.data[lmap] * lweight
	#   mdata[where(1 - lweight gt 1)] = data[hmap[where(1 - lweight gt 1)]]
	#   mdata[where(    lweight gt 1)] = data[lmap[where(    lweight gt 1)]]
	# endif

	# ; SPLINE
	# ; Use the spline function on each row, then repeat the first or last
	# ; pixel so that extrapolation does not occur.
	# if keyword_set(fitspline) then begin
	#   mdata = dblarr(nxo, ny)
	#   for i = 0, ny - 1 do begin &$
	#     mdata[*,i] = spline(inx[*,i], data[*,i], outx[*,i], 0.1) &$
	#   endfor
	#   lez = where(outx le inx[0,     *] ## (fltarr(nxo) + 1), lc)
	#   gez = where(outx ge inx[nxi-1, *] ## (fltarr(nxo) + 1), gc)
	#   if lc ne 0 then mdata[lez] = data[0,     lez/nxo]
	#   if gc ne 0 then mdata[gez] = data[nxi-1, gez/nxo]
	# endif

	# ; CUBIC
	# ; Use the interpolate function with a cubic convolution of -0.5.
	# ; First calculate the output coordinates when the input coordiante is
	# ; just the pixel value.
	# if keyword_set(fitcubic) then begin
	#   mdata = dblarr(nxo, ny)
	#   for i = 0, ny-1 do begin
	#     minx = min(inx[*, i])
	#     maxx = max(inx[*, i])
	#     xvals = (outx[*, i] - minnx) / (maxx - minx) * (nxi - 1)
	#     mdata[*, i] = interpolate(data, xvals, cubic = -0.5)
	#   endfor
	# endif

	# AVERAGE
	# ; Average the closest pixels to the needed coordiante.  Add together
	# ; those pixels with a difference less than range, and divide by the
	# ; number of images
	if rc.fitaverage:
	dataa = rebin_nd(data, (nxi, ny, nxo))
	dataa = np.transpose(dataa, (2, 1, 0))
	# range = mean(outx[:, 1:nxo - 1] - outx[:, 0:nxo - 2]) / 2
	#   dmask = (dif ge -range) * (dif le range)
	#   mdata = total(dataa*dmask, 3) / total(dmask, 3)

	# ; POLYNOMIAL FITTING
	# ; Loop through every pixel in the output array.  Grab the closest
	# ; pixels and if there is enough around the needed coordinate fit a
	# ; polynomial to them. Evaluate the polynomial and the output coordinate
	if rc.fitpoly:
		# TODO switch ny and nxo 
		mdata = outx - outx
	#	for i in range(ny - 1):
	#	  range = mean(outx[1:nxo-1, i] - outx[0:nxo-2, i]) / 2
	#     minx = min(inx[*, i])
	#     maxx = max(inx[*, i])
	#     for j in range(nxo - 1):
	#       x1 = min(where(inx[*, i] ge outx[j, i]-range))
	#       x2 = max(where(inx[*, i] le outx[j, i]+range))
	#       if minx gt outx[j, i] then begin
	#         mdata[j, i] = data[x1, i]
	#       endif else if maxx lt outx[j, i] then begin
	#         mdata[j, i] = data[x2, i]
	#       endif else if x1 eq x2 then begin
	#         mdata[j, i] = data[x1, i]
	#       endif else if x1 gt x2 then begin
	#         dif = inx[x2, i] - inx[x1, i]
	#         y1 = data[x2, i] * (outx[j, i] - inx[x1, i]) / dif
	#         y2 = data[x1, i] * (inx[x2, i] - outx[j, i]) / dif 
	#         mdata[j, i] = y1 + y2
	#       endif else if (x2-x1+1) le degcontr then begin
	#         co = linfit(inx[x1:x2, i], data[x1:x2, i], /double)
	#         mdata[j, i] = co[0] + co[1]*outx[j, i]
	#       endif else begin
	#         co = poly_fit(inx[x1:x2, i], data[x1:x2, i], degcontr, $
	#                       /double, yfit = est, status = status)
	#         mdata[j, i] = poly(outx[j, i], co)
	#       endelse

	return mdata
