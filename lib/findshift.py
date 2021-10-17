"""
Name: findshift.py

Purpose: Find the shift in the trace for each row of an input vector.

Category: Shifting Procedures

Calling Example:
		findshift(dataim, centerdeg = centerdeg, centerfit = centerfit,
				centroid = centroid, gaussth=gaussth, height=height,
				inmask=inmask, plottype = plottype, shifth=shifth,
				shmask = shmask, spec = spec, tracemask=tracemask,
				varim=varim, verbose = verbose, widthth = widthth,
				widthv = widthv, x1=x1, x2=x2)
Inputs:
		dataim: sky-subtracted processed image, the horizontal direction as the
				spatial direction, and the vertical as the wavelength
Keyword Arguments:
		General
		varim: the variance image for weighting of the polynomial fit
		inmask: the mask used for all functions
		spec: array holding extracted spectra, same dimension as vertical of dataim
		x1,x2: boundaries in x which contain the spectrum
		verbose:

		#TODO FINISH DOCS

Outputs:

History:

Created on 9/22/2021$
"""

from lib.misc import rebin_nd
from lib.fitting.centermass import centermass
from lib.fitting.gaussfunc import gaussfunc
from lib.procvect import procvect
from lib.excep import *
import numpy as np

def findshift(dataim, **kwargs):

	#TODO FUNC: findshift
	nx = np.shape(dataim)[1]
	ny = np.shape(dataim)[0]

	# Defaults and checks

	gaussth = kwargs.get('gaussth', 0.03)
	shiftth = kwargs.get('shiftth', 0.5)
	widthth = kwargs.get('widthth', 0.5)
	tracedeg = kwargs.get('tracedeg', 4)
	verbose = kwargs.get('verbose', 0)
	plottype = kwargs.get('plottype', 0)
	varim = kwargs.get('varim', abs(dataim))
	inmask = kwargs.get('inmask', np.ones((nx, ny)))
	spec = kwargs.get('spec', np.sum(dataim, 1))
	x1 = kwargs.get('x1', 0)
	x2 = kwargs.get('x2', nx-1)
	gotovect = kwargs.get('gotovect', -1)
	centroid = kwargs.get('centroid', False)


	# ; Checks
	# str = ""
	# str += n_elements(spec) ne ny                             ? "spec "   : ""
	# str += x1 lt 0 || x1 gt x2                                ? "x1 "     : ""
	# str += x2 gt nx - 1 || x2 lt x1                           ? "x2 "     : ""
	# str += (size(varim ))[1] ne nx || (size(varim ))[2] ne ny ? "varim "  : ""
	# str += (size(inmask))[1] ne nx || (size(inmask))[2] ne ny ? "inmask " : ""
	# str += 4 ne n_elements(plottype)                          ? "plottype" : ""
	# if str ne "" then message, "Poorly formed inputs: " + str

	shiftv = np.zeros(ny)   # the center of each row's profile
	widthv = np.zeros(ny)   # the width of each row's profile
	height = np.zeros(ny)   # the height of each row's profile
	xrange = np.arange(x2-x1+1) + x1 # the range of pixel values to examine
	shmask = np.ones(nx)
	tracemask = inmask

	for i in range(ny):
		mainmv = inmask[i, x1:x2]
		datav = dataim[i, x1:x2] > 0
		varv = varim[i, x1:x2] * mainmv
		specv = spec[i]

		if i == gotovect:
			verbose = 5
		if centroid:
 			shiftv[i] = centermass(datav, datav, varv, specv)
			widthv[i] = 1
			if verbose == 5:
				input("Stopping at center mass. Press enter to continue.")
		else:
			# TODO findshift plotting
			# if plottype[1] or (verbose eq 5) then begin
            #     device, window_state = ws
            #     if not ws[12] then window, 12 else wset, 12
            #     !p.multi = [0, 1, 2, 1, 1]
            #     plot, datav,   title = 'Data Vector for Finding Center', /ystyle
            #     plot, mainmv,   title = 'Input Mask', yrange = [-0.1, 1.1]
            #     !p.multi = 0
            #     wait, 0.01
            #     endif
			if verbose == 5:
				input("Stopping at gassfunc")
			func = 'gaussfunc'
			parm = 1
			foo, _, coeff = procvect(datav, varv=varv, noupdate=True, thresh=gaussth, maskv=mainmv,
			               func = func, parm=parm, verbose=verbose, plottype=plottype,
			               multv=specv, vectnum=i, absthresh=True)
			tracemask[i, x1:x2] = mainmv
			foo = gaussfunc()
			if coeff[0] == 0:
				shmask[i] = 0
			else:
				shmask[i] = 1
				shiftv[i] = coeff[1]
				widthv[i] = coeff[2]
				height[i] = coeff[0]

	varv = np.zeros(ny) + (1.0/10)**2
	func = 'polyfunc'
	parm = tracedeg

	#TODO finish findshift
	if rc.centroid:
		brows = np.where(sum(rc.inmask[:, xrange], 1) < len(xrange), rcount)

	return