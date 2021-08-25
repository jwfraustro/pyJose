# TODO DOCS: optspecextr_source
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
from lib.vectsetup import fitbg, fitprof
from lib.utils import load_config
from lib.horne import stdextr
import numpy as np
import os, sys


def optspecextr(config_file):
	# Set up run configuration for variables
	rc = load_config(config_file)

	rc.v0 = rc.rn ** 2
	rc.x1 = np.round(rc.x1)
	rc.x2 = np.round(rc.x2)
	rc.varim = rc.var

	if rc.verbose == 5:
		input("Stopping at fitting sky background, press enter to continue.")

	rc.bgim = fitbg(rc)

	rc.dataim = rc.data - rc.bgim
	rc.stdspec = stdextr(rc)

	if rc.integrate:
		rc.spec = rc.adjspec
	else:
		rc.spec = rc.stdspec

	if rc.verbose == 5:
		cont = input("Stopping at profile fitting. Press any key to continue or (q) to quit.")
		if cont == 'q' or "Q":
			sys.exit(0)

	rc.profim = fitprof(rc)
	#
	# verbose = save_verbose
	# plot_type = save_plottype
	#
	rc.optspec = extrspec(dataim, profim, varim, v0, q, x1, x2, **kwargs)
	#
	# varout = varim
	#
	#
	# #summaryopt, data, varout, bgim, profim, inmask, bgmask, exmask, $
	# #        difpmask, optspec, stdspec, traceest, $
	# #        berrvect, perrvect, eerrvect, $
	# #        verbose, plottype, adjparms = adjparms, debughead = debughead

	return optspec
