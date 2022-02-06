from lib.excep import VectorLengthException

def centermass(xvals, datav, varv, specv):
	"""
	Function: centermass

	Description:
		Function which returns the center of mass of datav

	Category:
		Optimal Spectrum Extraction Package
			- Vector fitting functions

	Parameters:
		xvals: The x values for the data
		datav: the sky subtracted data
		varv:  the variance
		specv: the spectrum

	Outputs:
		Returns the center of mass of the data vector / spec vector

	Restrictions:
		All vectors must be the same length.

	"""
	nx = len(xvals)

	if len(datav) != nx:
		raise VectorLengthException("datav", "xvals")

	multv = sum(datav / specv)
	return sum(xvals * datav / specv / multv)