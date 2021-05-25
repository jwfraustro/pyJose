#TODO shift docs
"""
A collection of functions for vector shifting and sampling.

Calling Example:

Inputs:

Outputs:

History:

Created on 5/25/2021$
"""


def findshift():
	"""
	Finds the shift in the trace for each row of an input vector.

	Calling Example:
		#TODO findshift docs
	Inputs:

	Outputs:

	History:

	Created on 4/17/2021$
	"""
	#TODO findshift

	return


def sampleshift():
	# TODO sampleshift docs
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
		#TODO sampleshift history
	Created on 4/17/2021$
	"""
	#TODO sampleshift
	return