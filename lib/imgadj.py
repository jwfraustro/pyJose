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

	# Unpack Optimal Spectrum Extraction Info
	dataim = rc.OptInfo.dataim
	varim = rc.OptInfo.varim
	verbose = rc.OptInfo.verbose
	plottype = rc.OptInfo.plottype

	inmask = rc.OptInfo.inmask



	return