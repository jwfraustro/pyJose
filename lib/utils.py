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

class Config:
	def __init__(self, **config_entries):
		self.__dict__.update(config_entries)

def load_config(config_file):
	with open(config_file) as r:
		config_map = yaml.safe_load(r)

	data_file = config_map['DATA']
	var = config_map['VAR']
	rn = config_map['RN']
	q = config_map['Q']
	x1 = config_map['X1']
	x2 = config_map['X2']

	with fits.open(data_file) as r:
		data = r[0].data

	RunConfig = Config(**config_map['OPTIONS']['INPUT'])
	OutputConfig = Config(**config_map['OPTIONS']['OUTPUT'])

	return data, var, rn, q, x1, x2, RunConfig, OutputConfig
