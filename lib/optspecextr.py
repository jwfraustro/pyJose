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
from lib.vectsetup import fitprof
from lib.fitbg import fitbg
from lib.utils import load_config, save_fits
from lib.stdextr import stdextr
from lib.display import interactive_jose
import numpy as np
import os, sys
from matplotlib import pyplot as plt


def optspecextr(config_file):
	# Set up run configuration for variables
	data, opts = load_config(config_file)

	opts['v0'] = opts['rn'] ** 2
	opts['x1'] = np.round(opts['x1'])
	opts['x2'] = np.round(opts['x2'])

	if "var" in opts:
		varim = opts['var']
	else:
		varim = np.copy(data)
		varim = abs(varim) / opts['q'] + opts['rn']**2
	if opts['verbose'] == 5:
		input("Stopping at fitting sky background, press enter to continue.")


	#if rc.plottype == 5:
	#	interactive_jose(rc)

	verbose = opts['verbose']
	plottype = opts['plottype']
	#Fit Background
	bgim, varim = fitbg(data, varim=varim, **opts)

	opts['verbose'] = verbose
	opts['plottype'] = plottype

	dataim = data - bgim

	save_fits(dataim, 'bg_subtracted')
	save_fits(data, 'original_data')
	save_fits(bgim, 'bg_img')

	stdspec, stdvar, adjspec = stdextr(dataim, varim, **opts)

	if opts['integrate'] == True:
		spec = adjspec
	else:
		spec = stdspec

	if opts['verbose'] == 5:
		cont = input("Stopping at profile fitting. Press any key to continue or (q) to quit.")
		if cont == 'q' or "Q":
			sys.exit(0)

	rc.profim = fitprof(rc)
	#
	# verbose = save_verbose
	# plot_type = save_plottype
	#
	rc.optspec = extrspec(rc)
	#
	# varout = varim
	#
	#
	# #summaryopt, data, varout, bgim, profim, inmask, bgmask, exmask, $
	# #        difpmask, optspec, stdspec, traceest, $
	# #        berrvect, perrvect, eerrvect, $
	# #        verbose, plottype, adjparms = adjparms, debughead = debughead

	return optspec
