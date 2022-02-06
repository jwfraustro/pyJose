"""
Name: tests.py

Purpose: Testing Functions for JOSE

Category:

Calling Example:

Inputs:

Outputs:

History:

Created on 9/19/2021$
"""

import numpy as np

def polyfunc_test():
	from lib.fitting import polyfunc
	from matplotlib import pyplot as plt
	print("Polyfunc Test:")
	nx = 100
	xvals = np.arange(nx)
	profv = ((-xvals/nx*2 - 1)**2 + 1)*0.75/nx*2
	sbump = 10
	specv = np.sin(xvals/nx * np.pi * 6) + sbump
	gain = 10
	basey = profv*specv
	datav = np.random.normal(basey, basey,(nx,))/gain
	varv = basey/gain
	badloc = [int(x) for x in np.random.uniform(0, 50, 4)]
	datav[badloc] = np.random.uniform((4,))*max(datav)*2
	profv[badloc] = 0
	goodloc = np.where(profv != 0)
	eval = 0
	deg = 0
	coeffv = []
	est, coeffv = polyfunc.polyfunc(xvals, datav, profv, specv, eval, coeffv, deg)
	plt.plot(est)
	plt.plot(datav/specv)
	plt.title("2 Degree plot")
	plt.show()
	plt.close()

	return

def procvect_test():
	from lib.procvect import procvect

	nx = 5

	datav = np.arange(nx)
	xvals = np.arange(nx)

	procvect(datav, xvals=xvals, func="myfunc")


if __name__ == '__main__':
	#procvect_test()
	polyfunc_test()