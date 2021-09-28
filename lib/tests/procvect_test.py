"""
Name: procvect_test.py

Purpose: Test the procvect package.

Category: Tests

Calling Example: test_procvect()

Created on 9/27/2021$
"""

import numpy as np
from numpy.polynomial import Polynomial
from lib.procvect import procvect
from matplotlib import pyplot as plt
from astropy.io import fits
import random
import unittest


class TestProcvect(unittest.TestCase):
	def test_bad_pixs(self):
		# Generate a basic polynomial and apply it to a vector
		p = Polynomial([200, 0, -1])

		xvals = np.linspace(-1, 1, 512)

		vect = p(xvals)

		# Generate some bad pixels of varying values
		bad_pix_count = 15
		sigma_range = np.arange(4,10)

		bad_pixels = random.sample(range(len(xvals)), bad_pix_count)

		# Apply the bad pixels
		for i in bad_pixels:
			sigma = random.choice(sigma_range)
			vect[i] = vect[i]*sigma

		# Save FITS file of the initial vector
		hdu = fits.PrimaryHDU(vect)
		hdu.writeto("./test_images/procvect_1d_badpix.fits", overwrite=True)

		# Run procvect
		fiteval, maskv, _, _ = procvect(vect, func="polyfunc", deg=2)

		# Apply bad pixel mask
		vect = vect * maskv

		# Subtract fiteval
		if np.all(np.isclose(vect-fiteval, np.zeros_like(vect))):
			subtracted = np.zeros_like(vect)
		else:
			subtracted = vect - fiteval*maskv

		# Save FITS of subtracted vector
		hdu = fits.PrimaryHDU(subtracted)
		hdu.writeto("./test_images/procvect_1d_badpix_sub.fits", overwrite=True)

		self.assertTrue(np.all(np.isclose(subtracted, np.zeros_like(subtracted))))

	def test_procvect_2d(self):

		# Make a gradient
		ny, nx = 512, 512

		arr = np.zeros((ny, nx))

		xvals = np.linspace(-1, 1, 512)
		yvals = np.linspace(200, 0, 512)

		for i in range(ny):
			p = Polynomial([200 - yvals[i], 0, -1])
			vect = p(xvals)
			arr[i, :] = vect

		# Save FITS file of the gradient
		hdu = fits.PrimaryHDU(arr)
		hdu.writeto("./test_images/procvect_2d_b4.fits", overwrite=True)

		# Run procvect on each row
		# We'll check to see if the subtraction is within 1e-5 of 0 and if so, set it to zero.
		for i in range(ny):
			fiteval, _, _, _ = procvect(arr[i, :], func="polyfunc", deg=2)
			if np.all(np.isclose(arr[i,:]-fiteval, np.zeros_like(arr[i,:]))):
				arr[i,:] = np.zeros_like(arr[i,:])
			else:
				arr[i, :] = arr[i, :] - fiteval

		arr = abs(arr)

		hdu = fits.PrimaryHDU(arr)
		hdu.writeto("./test_images/procvect_2d_sub.fits", overwrite=True)

		self.assertTrue(np.all(np.isclose(arr, np.zeros_like(arr))))

	def test_procvect_1d(self):
		# Generate a basic polynomial and apply it to a vector
		p = Polynomial([200, 0, -1])

		xvals = np.linspace(-1, 1, 512)

		vect = p(xvals)

		# Save FITS file of 1D vector
		hdu = fits.PrimaryHDU(vect)
		hdu.writeto("./test_images/procvect_1d_b4.fits", overwrite=True)

		# Run procvect
		fiteval, _, _, _ = procvect(vect, func="polyfunc", deg=2)

		# Output Fit Overlay
		plt.plot(vect)
		plt.plot(fiteval)
		plt.savefig("./test_images/procvect_1d_fit.png")

		# Subtract fit from original vector and make positive
		# We'll check to see if the subtraction is within 1e-5 of 0.
		if np.all(np.isclose(vect-fiteval, np.zeros_like(vect))):
			subtracted = np.zeros_like(vect)
		else:
			subtracted = vect - fiteval

		# Save FITS of subtracted vector
		hdu = fits.PrimaryHDU(subtracted)
		hdu.writeto("./test_images/procvect_1d_sub.fits", overwrite=True)

		self.assertTrue(np.all(np.isclose(subtracted, np.zeros_like(subtracted))))


if __name__ == '__main__':
	unittest.main()
