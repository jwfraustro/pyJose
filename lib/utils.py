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
from lib.excep import ParameterException
import os

class Config:
	def __init__(self, **config_entries):
		self.__dict__.update(config_entries)

def check_defaults(rc):

	with open("./lib/defaults.yaml") as r:
		defaults = yaml.safe_load(r)

	for key, value in rc.__dict__.items():
		#TODO This is bad, gotta fix it
		if key == 'data':
			continue
		if value == None:
			rc.__dict__[key] = defaults.get(key)

	if (rc.x1 < 0) or (rc.x1 > rc.x2):
		raise ParameterException("x1 must be between 0 and x2.")

	if rc.plottype not in range(0, 5):
		raise ParameterException("Plot type must be a value of 0-4")

	if not os.path.exists(rc.output_dir):
		try:
			os.mkdir(rc.output_dir)
		except:
			pass


def load_config(config_file):
	with open(config_file) as r:
		config_map = yaml.safe_load(r)

	rc = Config(**config_map)

	with fits.open(rc.data_file) as r:
		fits_data = r[0].data

	rc.data = fits_data

	check_defaults(rc)

	return rc
