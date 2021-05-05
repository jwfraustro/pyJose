#TODO DOCS
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
from lib.fitbg import *
from lib.utils import load_config
import numpy as np

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

	bgim = fitbg(data, x1, x2, q, v0, RunConfig) #TODO

	# verbose = save_verbose
	# plot_type = save_plottype
	#
	# dataim = data - bgim
	#
	# stdspec = stdextr(dataim, varim, x1, x2, **kwargs) #TODO
	#
	# if 'integrate' in kwargs:
	# 	if kwargs['integrate'] == True:
	# 		spec = kwargs['adjspec']
	# else:
	# 	spec = stdspec
	#
	# save_verbose = verbose
	# save_plottype = plot_type
	#
	# # if (verbose eq 5 ) then stop #TODO
	# profim = fitprof(dataim, spec, x1, x2, **kwargs)
	#
	# verbose = save_verbose
	# plot_type = save_plottype
	#
	# optspec = extrspec(dataim, profim, varim, v0, q, x1, x2, **kwargs)
	#
	# varout = varim
	#
	# #TODO
	# #summaryopt, data, varout, bgim, profim, inmask, bgmask, exmask, $
    # #        difpmask, optspec, stdspec, traceest, $
    # #        berrvect, perrvect, eerrvect, $
    # #        verbose, plottype, adjparms = adjparms, debughead = debughead

	return optspec