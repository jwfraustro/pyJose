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

def plot_procvect(bgim, yvals, allx, sy1, sy2, title, xtitle, ytitle):

	return
def plotting():
	#TODO DOCS: plotting
	#TODO FUNC: plotting
	return
