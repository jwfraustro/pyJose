from .excep import VectorLengthException

#TODO centermass docs
"""
Function: centermass

Description:
	Function which returns the center of mass of datav

Category:
	Optimal Spectrum Extraction Package
		- Vector fitting functions

Parameters:
	x_vals: The x values for the data
	data_v: the sky subtracted data
	var_v:  the variance
	spec_v: the spectrum

Outputs:
	Returns the center of mass of the data vector / spec vector

Restrictions:
	All vectors must be the same length.

Modification History:
	TODO


"""

def centermass(x_vals, data_v, var_v, spec_v):
	nx = len(x_vals)

	if len(data_v) != nx:
		raise VectorLengthException(data_v, nx)

	multv = sum(data_v / spec_v)
	return sum(x_vals*data_v/spec_v/multv)