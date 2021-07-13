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
	datav_plot.set_title("Data Vector for BG Fitting")

	maskv_plot.plot(maskv)
	maskv_plot.set_title("Input Mask")

	varv_plot.plot(varv)
	varv_plot.set_title("Variance")

	skyvarv_plot.plot(skyvarv)
	skyvarv_plot.set_title("Sky Variance")

	plt.tight_layout()
	plt.show()
	#TODO Save plot
	fig.savefig(os.path.join(output_dir, "fitbg_plot.png"))

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
