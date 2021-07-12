# TODO DOCS: misc
"""
Name:
	misc
Purpose:
	Collection of functions for plotting and summarizing process output.
Category:

Calling Example:

Inputs:

Outputs:

History:

Created on 5/25/2021$
"""

import matplotlib.pyplot as plt
import os

def rebin_nd(ndarray, new_shape, operation='sum'):
    """
    Bins an ndarray in all axes based on the target shape, by summing or
    averaging.
	
    Number of output dimensions must match number of input dimensions.
    
	Notes
	-----
	This is pulled from a stack exchange question and subsequent edit, both 
	linked below:
	https://stackoverflow.com/questions/8090229/resize-with-averaging-or-rebin-a-numpy-2d-array
	https://gist.github.com/derricw/95eab740e1b08b78c03f
	
	Example
    -------
    >>> m = np.arange(0,100,1).reshape((10,10))
    >>> n = rebin_nd(m, new_shape=(5,5), operation='sum')
    >>> print(n)
    [[ 22  30  38  46  54]
     [102 110 118 126 134]
     [182 190 198 206 214]
     [262 270 278 286 294]
     [342 350 358 366 374]]
    """
    if not operation.lower() in ['sum', 'mean', 'average', 'avg']:
        raise ValueError("Operation {} not supported.".format(operation))
    if ndarray.ndim != len(new_shape):
        raise ValueError("Shape mismatch: {} -> {}".format(ndarray.shape,
                                                           new_shape))
    compression_pairs = [(d, c // d) for d, c in zip(new_shape, ndarray.shape)]
    ndarray = ndarray.reshape([l for p in compression_pairs for l in p])

    for i in range(len(new_shape)):
        if operation.lower() == "sum":
            ndarray = ndarray.sum(-1 * (i + 1))
		# This needs to be edited in case there are NaNs.
        elif operation.lower() in ["mean", "average", "avg"]:
            ndarray = ndarray.mean(-1 * (i + 1))
    return ndarray

def summaryopt():
	# TODO DOCS: summaryopt
	"""
	Name:
		summaryopt
	Purpose:
		Plot summary information from the optimal extraction to the screen so the user
		can make sure the process is running cleanly.
	Category:
		Optimal Spectrum Extraction Package
			- Miscellaneous
	Calling Example:

	Inputs:

	Outputs:

	History:

	Created on 4/17/2021$
	"""
	#TODO FUNC: summaryopt
	return

def plot_fitbg(datav, maskv, varv, skyvarv, output_dir):

	fig, ((datav_plot, maskv_plot),(varv_plot, skyvarv_plot)) = plt.subplots(2,2)
	datav_plot.plot(datav)
	datav_plot.ylabel = "Data Vector for BG Fitting"

	maskv_plot.plot(maskv)
	maskv_plot.ylabel = "Input Mask"

	varv_plot.plot(varv)
	varv_plot.ylabel = "Variance"

	skyvarv_plot.plot(skyvarv)
	skyvarv_plot.ylabel = "Sky Variance"

	plt.show()
	#TODO Save plot
	#plt.imsave(fig, os.path.join(output_dir, "fitbg_plot.png"))

def plot_procvect(bgim, yvals, allx, sy1, sy2, i):

	fig, ax = plt.subplots(subplot_kw={'projection':'3d'})

	surf = ax.plot_surface(allx, yvals[sy1:sy2], bgim[:, sy1:sy2], antialiased=False)

	fig.suptitle("Background after #"+str(i))
	ax.xlabel = "X pixel location"
	ax.ylabel = "Y (Prefit and Postfit)"

	plt.show()
	#TODO Save plot
	#plt.imsave(fig, os.path.join(output_dir, "procvect_plot.png"))

	return
def plotting():
	#TODO DOCS: plotting
	#TODO FUNC: plotting
	return
