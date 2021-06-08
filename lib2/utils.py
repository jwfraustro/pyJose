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

	RunConfig = Config(**config_map)

	with fits.open(RunConfig.data_file) as r:
		data = r[0].data

	RunConfig.data = data

	return RunConfig
