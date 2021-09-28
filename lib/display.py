"""
Name: Graphing

Purpose: Collection of functions for handling the real-time display and interactive plotting functions.

Category:

History:

Created on 8/25/2021$
"""

from pyqtgraph.flowchart import Flowchart
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph.metaarray as marray

def interactive_jose(data):

	app = pg.mkQApp("pyJOSE")

	win=QtGui.QMainWindow()
	pg.setConfigOptions(imageAxisOrder='row-major')
	win.resize(800,800)
	imv = pg.ImageView()
	win.setCentralWidget(imv)
	win.show()
	win.setWindowTitle('pyJOSE: example.fits')
	imv.setImage(data)
	imv.ui.roiBtn.setChecked(True)
	imv.roiClicked()
	pg.exec()

	return