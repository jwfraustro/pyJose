import numpy as np

"""
Function: ADJGAUSS

Description:
	Standardize the pixel locations when the profile image is a
	Gaussian with quickly changing center and/or height.

Category:
	Optimal Spectrum Extraction Package
		- Image Adjustment

Parameters:
	in_array: an [nx, ny, nimages] array of images to adjust.
	opt_info: a structure of images and info from the optimal extraction:
		data_im:    the image to use when finding the profile
		var_im:     the variance image for weighting of the poly_fit
		sky_var:    the corresponding sky variance image (0 array)
		bg_im:      the background image to be subtracted off the data (0 array)
		in_mask:    the mask used for all functions (1 array)
		q:          the gain of the image (1)
		v0:         the readnoise ^2 of the data (0)
		x1:         the start of the object's profile (0)
		x2:         the end of the object's profile (nx-1)
		spec:       the standard extracted spectrum (total(dataim[x1:x2, *], 1)
		bpct:       the percentage of allowable bad pixels before halting (.5)
		verbose:    the level of output to screen ([0,0,0,0])
		plot_type:  the type of plot to output while running (0)

Optional Inputs:
	Profile Fitting - Straightening:
	straighten:     Set to straighten the trace
	expand:         level of expansion in straightening routine (11)
	centroid:       set to straighten with centroid instead of Gaussian
	trace_deg:      degree of fitting for estimating shift in trace (2)
	no_trace_fit:   set to not fit the estimated trace
	gauss_th:       error allowed during row's center calculation (0.05)
	shift_th:       error allowed between the center and the fitted trace (0.5)
	deg_contr:      degree of fit to use when contracting the array (2)

Keyword Parameters:
	adj_params: A structure for parameters of the adjustment width
		orig_x:     the x values for each pixel in the original images
		adj_x:      the x values for each pixel in the adjusted images
		trace_est:  the estimated center for each row
		width_est:  the estimated width for each row
	adj_options: A structure with options for the gaussian adjustment
		level:      the amount to expand the array
		center:     set to True to adjust each column so the centers line up
		width:      adjust to reflect each pixel's relative distance
					from center according to that row's Gaussian width
		centroid:   set to use centroid fitting instead
		gauss_th:   the absolute value of allowable error when fitting
					a Gaussian to the data
		shift_th:   The absolute allowable error between a row's center
					and the fitted trace
		deg_contr:  degree of fit to use when contracting the array
		center_deg: degree of fit for estimating shift
		center_fit: set to not fit the estimated centers and not accept
					the fitting to be correct
	revert: set to use the data in adj_params to convert the images
			in the input array back into the original shape

Outputs:
	An array (( x2 - x2 + 1) * level by ny by n_images) where [*,*,i]
	standardized image from the corresponding image in the input array

Procedure:
	Unpack the structures first. Then find the width and centers using findshift.
	Finally adjust each image in in_array so that the profiles line up.


"""
class OptInfo(object):
	def __init__(self, opt_info):

		self.data_im = None
		self.var_im = None
		self.sky_var = None
		self.bg_im = None
		self.in_mask = None
		self.q = None
		self.v0 = None
		self.x1 = None
		self.x2 = None
		self.

		for key in opt_info:
			setattr(self, key, opt_info[key])

		if not hasattr(self, 'verbose'):
			self.verbose = False
		if not hasattr(self, 'plot_type'):
			self.plot_type = [0,0,0,0]

class AdjOptions(object):
	def __init__(self, adj_opt):
		for key in adj_opt:
			setattr(self, key, adj_opt[key])

		if not hasattr(self, 'level'):
			self.level = 7
		if not hasattr(self, 'center'):
			self.center = True
		if not hasattr(self, 'width'):
			self.width = 0

class AdjParams(object):
	def __init__(self, adj_params):
		for key in adj_params:
			setattr(self, key, adj_params[key])

def adj_gauss(in_array, opt_info: dict, adj_params: dict, adj_options: dict, revert=False, **kwargs):
	# Initialize and check inputs

	# Sizes:
	n_image = in_array.shape[0]
	nx = in_array.shape[1]
	ny = in_array.shape[2]

	# Unpack Optimal Spectrum Extraction Info
	opt_info = OptInfo(opt_info)

	# Unpack Adjustment Procedure Options
	adj_options = AdjOptions(adj_options)

	# Unpack Adjustment parameters
	adj_params = AdjParams(adj_params)

	# Check inputs
	if adj_options.level < 1:
		print("Level must be greater or equal to 1.")

	#REVERT TO ORIGINAL IMAGE SPECIFICATIONS
	#Grab the old geometry's coordinate array and the new coordinate array and
	#use sampleshift to estimate the pixel value at the old geometry using
	#closest pixels in the new geometry

	if revert:
		orig_x = adj_params.origx


	return
