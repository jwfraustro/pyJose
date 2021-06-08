#TODO DOCS: excep

class ParameterException(Exception):
	def __init__(self, message):
		self.message = message
		super().__init__(self.message)

class VectorLengthException(Exception):
	"""Exception raised if vector lengths are not equal.

	Attributes:
		vx, nx      -- vectors being compared.
	"""

	def __init__(self, ax, bx):
		self.message = "Vector of %s does not match vector length of %s." % (ax, bx)
		super().__init__(self.message)