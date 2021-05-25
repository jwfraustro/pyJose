#TODO extractfunc docs
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
from lib.done.excep import VectorLengthException
import numpy as np

def extractfunc(xvals, datav, varv, profv, eval, coeffv, opvar):

	# Check Inputs

	nx = len(datav)

	if nx != len(profv):
		raise VectorLengthException(nx, profv)
	if nx != len(varv):
		raise VectorLengthException(nx, varv)

	# Always Extract

	gl = np.where(profv != 0)
	if (len(gl) == 0):
		print("No good pixels in datav.")
		return 0

	denom = np.sum((profv[gl] * profv[gl]) / varv[gl]) # avoid recalc
	opt = np.sum((profv[gl] * datav[gl]) / varv[gl]) / denom
	opvar = np.sum(profv[gl]) / denom

	return opt