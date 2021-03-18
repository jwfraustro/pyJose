"""
Function: box_car_func

Description:
	Function for fitting the data with a boxcar average or median fit.

Category:
	Optimal Spectrum Extraction Package
		- Vector fitting functions

Parameters:
	x_vals:     The x values for the data
	data_v:     The vector with cosmic rays to be fitted
	var_v:      The variance vector
	spec_v:     The spectrum vector. Bad pixels have a spec_v of 0
	eval:       Set to True to evaluate the inputs
	coeffv:     A filler in this routine
	boxcar_hw:  The halfwidth of the median or smooth filter

Outputs:
	Returns an estimation of the data by returning data_v ran
	through a median filter. If eval is set the routine smooths
	over the good pixels, then interpolates the bad pixels.

Restrictions:
	All vectors must be the same length. Boxcar_hw must be an int >= 1

Procedure:
	If eval == False, then do a median filter over the data_v / spec_v. Else,
	smooth over the good pixels and interpolate over the bad pixels.

Example:
	TODO

Modification History:
	TODO

	

"""