



def extractfunc(rc):
	# TODO DOCS: extractfunc
	"""
	Name:
		extractfunc
	Purpose:
		A function which extracts the optimal spectrum from a vector.
	Category:
		Optimal Spectrum Extraction Package
			- Vector Fitting functions
	Calling Example:

	Inputs:

	Outputs:

	History:

	Created on 4/17/2021$
	"""
	# Check Inputs

	nx = len(rc.datav)

	if nx != len(rc.profv):
		raise VectorLengthException("datav", "profv")
	if nx != len(rc.varv):
		raise VectorLengthException("datav", "varv")

	# Always Extract

	gl = np.where(rc.profv != 0)
	if (len(gl) == 0):
		print("No good pixels in datav.")
		return 0

	denom = np.sum((rc.profv[gl] * rc.profv[gl]) / rc.varv[gl])  # avoid recalc
	rc.opt = np.sum((rc.profv[gl] * rc.datav[gl]) / rc.varv[gl]) / denom
	rc.opvar = np.sum(rc.profv[gl]) / denom

	return rc.opt
