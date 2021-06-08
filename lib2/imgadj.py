#TODO DOCS: imgadj

"""
This module contains functions for the adjustment of profile images.

Created on 5/25/2021$
"""

def adjgauss():
	#TODO DOCS: adjgauss
	"""
	Standardize the pixel locations for a Gaussian profile image with quickly changing center and/or height.

	Calling Example:
		result = adjgauss(inarray, optinfo, adj_parms=adj_parms, adj_options=adj_options, revert=revert)

	Inputs:
		in_array:   An [nx, ny, nimages] array of images to adjust.
		opt_info:   A structure of images and info from the optimal extraction:
			data_im:    the image to use when finding the profile
			var_im:     the variance image for weighting of the poly_fit
			sky_var:    the corresponding sky variance image (0 array)
			bg_im:      the background image to be subtracted off the data (0 array)
			in_mask:    the mask used for all functions (1 array)
			q:          the gain of the image
			v0:         the readnoise^2 of the data (0)
			x1:         the start of the object's profile
			x2:         the end of the object's profile (nx-1)
			spec:       the standard extracted spectrum (total(dataim[x1:x2, *],1)
			bpct:       the percentage of allowable bad pixels before halting (0.5)


	Outputs:

	History:

	Created on 4/17/2021$
	"""
	#TODO FUNC: adjgauss
	return