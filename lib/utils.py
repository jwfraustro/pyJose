import yaml

"""
Name:

Purpose:

Category:

Calling Example:

Inputs:

Outputs:

History:

Created on 4/17/2021$
"""

from astropy.io import fits

class AdjGaussVars():
	def __init__(self):
		self.inarray = None

class OptInfo():
	def __init__(self):
		self.dataim = None
		self.varim = None

class AdjParms():
	def __init__(self):
		self.origx = None
		self.adjx = None
		self.traceest = None
		self.widthest = None

class AdjOptions():
	def __init__(self):
		self.level = None
		self.center = True
		self.width = False
		self.centroid = False
		self.gaussth = None
		self.shiftth = None
		self.degcontr = None
		self.centerdeg = None
		self.centerfit = None


class Vars():
	def __init__(self):
		self.data = None
		self.var = None
		self.rn = None
		self.q = None
		self.x1 = None
		self.x2 = None
		self.verbose = None
		self.plottype = None
		self.v0 = None
		self.varim = None
		self.bgdeg = None
		self.inmask = None
		self.bthresh = None
		self.bgres = None
		self.errvect = None
		self.bgmask = None
		self.skyvar = None
		self.nobgfit = None
		self.gotovect = None
		self.bgotovect = None
		self.bpcst = None
		self.dataim = None
		self.stdspec = None
		self.stdvar = None
		self.adjspec = None
		self.spec = None
		self.noproffit = None
		self.profdeg = None
		self.profmask = None
		self.pthresh = None
		self.bgim = None
		self.profres = None
		self.fitgauss = None
		self.fitboxcar = None
		self.boxcarhw = None
		self.Pgotovect = None
		self.Perrvect = None
		self.bpct = None
		self.difpmask = None
		self.adjfunc = None
		self.adjparms = None
		self.adjoptions = None
		self.opvar = None
		self.ethresh = None
		self.exres = None
		self.Egotovect = None
		self.Eerrvect = None
		self.exmask = None
		self.varout = None
		self.profim = None
		self.optspec = None
		self.stdspec = None
		self.traceest = None
		self.Berrvect = None


class Config:
	def __init__(self, **config_entries):
		self.__dict__.update(config_entries)


def load_config(config_file):
	with open(config_file) as r:
		config_map = yaml.safe_load(r)

	RunConfig = Config(**config_map)

	with fits.open(RunConfig.data_file) as r:
		data = r[0].data

	RunConfig.data = data

	return RunConfig
