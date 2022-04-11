
"""
This module defines the classes
	curve
which is used to construct and plot diagrams.
"""

import numpy as np

class curve:

	"""
	A diagram is fundamentally a collection of curves.
	A curve is sequence of points in spacetime, which can be defined in many coordinate systems.
	Every curve must at least have diagram coordinates UV for plotting, but not all coordinate systems must be filled.
	Each point has a well-defined radius independent of coordinate system.
	Curve attributes:
		r = radius
		tr = t,r = schwarzschild block coordinates
		uv = u,v = double null block coordinates
		uvdl = \tilde{u},\tilde{v} = region coordinates (standard penrose coords on the sss region)
		UV = U,V = diagram coordinates
		sty = dictionary of style parameters
	"""

	def __init__(self):
		## radius
		self.r = np.array([])
		## block coordinates
		self.tr = np.array([[],[]])
		## block coordinates
		self.uv = np.array([[],[]])
		## region coordinates
		self.uvdl = np.array([[],[]])
		## diagram coordinates
		self.UV = np.array([[],[]])
		## style
		self.sty = dict(c='k', ls='-', lw=0.5)


