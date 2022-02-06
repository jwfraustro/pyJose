"""
Name: 
	gaussfunc.py

Purpose: 
	Function for estimating a Gaussian curve to the data

Category: 
	Vector Fitting Functions

Calling Example:
	result = gaussfunc(xvals, datav, varv, specv, coeff, reest)
	
Created on 10/5/2021$
"""

from lib.excep import *
import numpy as np
from numpy import exp
from scipy.ndimage.filters import uniform_filter1d as smooth
from scipy.interpolate import CubicSpline
from scipy.optimize import curve_fit, OptimizeWarning

def newton_cotes(xvals, datav):

	"""
	Name: newton_cotes
	Purpose: Implementation of Newton-Cotes integration using Boole's rule.

	:param xvals: xvals of data vector
	:param datav: data vector
	:return:
	"""

	def booles_rule(y, h):

		bl = ( 7 * y[0] +
		      32 * y[1] +
		      12 * y[2] +
		      32 * y[3] +
		      7 * y[4]) * (2 * h) / 45

		return bl

	# check for sorted array
	if not all(xvals[i] <= xvals[i + 1] for i in range(len(xvals) - 1)):

		zipped_lists = zip(xvals, datav)
		sort_pairs = sorted(zipped_lists)
		tups = zip(*sort_pairs)
		xvals, datav = [list(tup) for tup in tups]

	# Divide the xvals into 5 equally spaced points

	a = min(xvals)
	b = max(xvals)
	h = (b - a) / 4

	xvals_new = [a, a + h, a + 2 * h, a + 3 * h, b]
	print(xvals_new)

	# Use cubic spline interpolation to calculate datav points at those values
	f = CubicSpline(xvals, datav)
	datav_new = f(xvals_new)
	print(datav_new)

	# Perform Boole's Rule
	bl = booles_rule(datav_new, h)

	return bl

def simple_gauss(x,a):

	z = (x - a[1]) / a[2]

	f = a[0] * exp(-z ** 2 / 2)

	return f

def gausseval(x,a, p=None):
	
	"""
	Name: gausseval
	Purpose: Small helper function for evaluating a Gaussian.
	"""

	z = (x - a[1]) / a[2]

	f = a[0] * exp(-z ** 2 / 2)

	if not p:
		return f
	else:
		p = np.zeros((len(x), 3))
		# 
		p[:, 0] = f / a[0]
		p[:, 1] = f * z / a[2]
		p[:, 2] = p[:, 1] * z

	return f, p

def gaussfunc(xvals, datav, varv, specv, coeff, reest, eval=False):
	
	"""
	Inputs:
		xvals: The x values for the data.
		datav: The vector with cosmic rays to be fitted.
		varv: The variance vector, pixels not to fit have a 0.
		specv: The spectrum at that point, bad pixels have a 0.
		coeff: Will contain the coefficients of the fit. 
		reest: Set to reestimate the spectrum after bad pixel rejection.
		eval: Do not fit, instead, evaluate based on provided coeff.
	Outputs:
		An estimation of the data by fitting a Gaussian. If a good fit is not found quickly 
		to the data, the Gaussian is fitted to the data smoothed instead.
		If eval is set it simply evaluates the Gaussian coefficients passed in at all points.
	"""
	# Check Inputs
	nx = len(xvals)

	if nx != len(datav):
		raise VectorLengthException("x_vals", "data_v")
	if nx != len(varv):
		raise VectorLengthException("x_vals", "var_v")
	if nx != len(specv) and len(specv) != 1:
		raise VectorLengthException("x_vals", "spec_v")
	if nx <= 4:
		raise ParameterException("data_v has less than 4 elements.")
	
	# Evaluate Coefficients
	if eval:
		f = gausseval(xvals, coeff)
		return f
	
	# Prepare Fitting
	if reest:
		specv = newton_cotes(xvals, datav)  # reestimate spectrum

	hbw = 4 < nx / 2 - 1  # half box width for estimation
	sdata = smooth(datav / specv, 2 * hbw + 1)  # estimate gauss
	shi = max(sdata)
	chi = np.argmax(sdata)
	slo = min(sdata)
	cli = np.argmin(sdata)

	if abs(shi) > abs(slo):
		ci = chi
	else:
		ci = cli

	hi = (datav / specv > sdata)[chi]  # estimated gaussian height
	center = xvals[ci]  # estimated gaussian center
	top = np.where(abs(datav / specv) > abs(hi / exp(1)))[0]
	fwhm = max(xvals[top]) - min(xvals[top]) + 1  # estimated gaussian FWHM
	base = np.median(sdata)  # estimated base level

	coeff = [hi, center, fwhm / 2]  # coefficients for estimate

	nz = np.where(varv != 0)[0]  # locations to fit

	est = datav / specv
	origcoeff = coeff


	try:
		parameters, covariance = curve_fit(simple_gauss, xvals[nz], datav[nz], coeff)
		for i in nz:
			est[i] = gausseval(i, parameters)
	except OptimizeWarning:
		mdata = np.median(datav/specv, 2*hbw+1)
		parameters, covariance = curve_fit(simple_gauss, xvals[nz], datav[nz], coeff)
		for i in nz:
			est[i] = gausseval(i, parameters)

	return est, coeff