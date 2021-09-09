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

def interactive_jose(rc):

	pg.setConfigOptions(imageAxisOrder='row-major')

	app = pg.mkQApp("ImageView Example")

	win = QtGui.QMainWindow()
	win.resize(800, 800)
	imv = pg.ImageView()
	win.setCentralWidget(imv)
	win.show()
	win.setWindowTitle('pyqtgraph example: ImageView')

	imv.setImage(rc.data)

	pg.exec()

	# app = pg.mkQApp("Interactive JOSE")
	#
	# # Create window with grid
	# win = QtGui.QMainWindow()
	# win.setWindowTitle('JOSE: '+rc.data_file)
	# cw = QtGui.QWidget()
	# win.setCentralWidget(cw)
	# layout = QtGui.QGridLayout()
	# cw.setLayout(layout)
	#
	# # Create flowchart
	# fc = Flowchart(
	# 		terminals={
	# 			'data_file':{'io':'in'},
	# 			'result_img':{'io':'out'}
	# 		}
	# )
	# w = fc.widget()
	#
	# # Add flowchart controls
	# layout.addWidget(fc.widget(), 0,0,2,1)
	# pw1 = pg.PlotWidget()
	# pw2 = pg.PlotWidget()
	# layout.addWidget(pw1, 0, 1)
	# layout.addWidget(pw2, 1, 1)
	#
	# win.show()
	#
	# fc.setInput(data_file=rc.data)
	# ## populate the flowchart with a basic set of processing nodes.
	# ## (usually we let the user do this)
	# plotList = {'Top Plot': pw1, 'Bottom Plot': pw2}
	#
	# pw1Node = fc.createNode('PlotWidget', pos=(0, -150))
	# pw1Node.setPlotList(plotList)
	# pw1Node.setPlot(pw1)
	#
	# pw2Node = fc.createNode('PlotWidget', pos=(150, -150))
	# pw2Node.setPlot(pw2)
	# pw2Node.setPlotList(plotList)
	#
	# pg.exec()

	return