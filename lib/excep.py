class ParameterException(Exception):
	def __init__(self, message):
		self.message = message
		super().__init__(self.message)

class VectorLengthException(Exception):
	"""Exception raised if vector lengths are not equal.

	Attributes:
		vx      -- vector which caused the error
		nx      -- number of elements expected
	"""

	def __init__(self, ax, bx):
		self.message = "Vector of length %i does not match vector length of %i." % (len(ax), bx)
		super().__init__(self.message)