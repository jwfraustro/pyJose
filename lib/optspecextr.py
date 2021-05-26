#TODO DOCS: optspecextr_source
"""
Name:
	optspecextr_source
Purpose:
	Data pipeline for image processing and optimal extraction.
Category:

Calling Example:

Inputs:

Outputs:

History:

Created on 4/17/2021$
"""
from lib.vectsetup import fitbg
from lib.utils import load_config
from lib.horne import stdextr
import numpy as np
import os, sys

def optspecextr(config_file):

	data, var, rn, q, x1, x2, RunConfig, OutputConfig = load_config(config_file)

	if not RunConfig.VERBOSE:
		RunConfig.VERBOSE = 0
	if not RunConfig.PLOTTYPE:
		RunConfig.PLOTTYPE = 0

	v0 = rn * rn
	x1 = np.round(x1)
	x2 = np.round(x2)
	varim = var

	if RunConfig.VERBOSE == 5:
		input("Stopping at fitting sky background, press enter to continue.")

	#TODO FUNC: optspecextr_source
	bgim = fitbg(data, x1, x2, q, v0, RunConfig)

	verbose = save_verbose
	plot_type = save_plottype

	dataim = data - bgim

	#FIXME: stdextr kwargs: stdvar is pass by ref, adjspec is bool, inmask array
	stdspec = stdextr(dataim, varim, x1, x2, inmask, stdvar, adjspec)

	if RunConfig.INTEGRATE:
		#FIXME: where is adjspec defined?
		spec = adjspec

	else:
		spec = stdspec

	#FIXME: for pass by reference
	save_verbose = verbose
	save_plottype = plot_type

	if RunConfig.VERBOSE == 5:
		cont = input("Stopping at profile fitting. Press any key to continue or (q) to quit.")
		if cont == 'q' or "Q":
			sys.exit(0)

	profim = fitprof(dataim, spec, x1, x2, **kwargs)
	#
	# verbose = save_verbose
	# plot_type = save_plottype
	#
	# optspec = extrspec(dataim, profim, varim, v0, q, x1, x2, **kwargs)
	#
	# varout = varim
	#
	#
	# #summaryopt, data, varout, bgim, profim, inmask, bgmask, exmask, $
    # #        difpmask, optspec, stdspec, traceest, $
    # #        berrvect, perrvect, eerrvect, $
    # #        verbose, plottype, adjparms = adjparms, debughead = debughead

	return optspec