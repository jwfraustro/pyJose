"""
Name:
	extrspec
Purpose:
	Optimally extract spectra using weighted profiles.
Category:
	Optimal Spectrum Extraction Package
		- Vector setup functions
Calling Example:
	#TODO extrspec calling example
Inputs:
	dataim: Sky-subtracted, processed image to extract spectrum from.
			Horizontal is pixel position, vertical is wavelength, assuming no curvature.
	profim: image of spatial profiles
	varim:  variance image from processed image
	v_0:    root(v_0) is squared readout noise in DN
	q:      effective number of photons per DN
	x1, x2: boundaries in x which contain spectrum (inclusive)
Optional Keywords:
	#TODO extrspec optional keywords
Outputs:
	Returns fopt, an array of optimally extracted spectra, and their variances.
History:
	#TODO extrspec history
Created on 4/17/2021$
"""

def extrspec(dataim, profim, varim, v0, q, x1, x2, **kwargs):
	#TODO extrspec
	return