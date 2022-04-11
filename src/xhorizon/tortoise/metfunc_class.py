
"""
This module defines the class
	metfunc
which is used to encode information about the SSS metric function f(r).
"""

import numpy as np



class metfunc:

	"""
	The metfunc class contains information derived from the SSS metric function f(r).
	The default values correspond to a Minkowski spacetime with f(r)=1.
	"""

	def __init__(self):
		## metric function parameters
		self.fparams = dict()
		## tortoise function parameters
		self.Fparams = dict(eps=1e-9, inf=25.)
		## radii r_j
		self.rj = np.array([0., self.Fparams['inf']])
		## slopes k_j
		self.kj = np.array([0., 0.])
		## printable info dict
		self.info = {
		'Type' : "Minkowski",
		'Metric Function' : r"$f(r) = 1$",
		'Parameters' : "",
		}
		## interpolation reference values for F(r)
		self.r_ref = np.array([])
		self.rstar_ref = np.array([])
		## first trapped interval (0 or 1)
		self.first_trapped_j = 1

	## metric function f(r)
	def f(self, r):
		return 1. + 0. * r

	## tortoise function F(r)
	def F(self, r):
		rstar = 1. * r
		return rstar

	## inverse tortoise function Finv_j(rstar)
	def Finv(self, j, rstar):
		## j and rstar are equal length arrays
		## if interval number is invalid return nan
		r = np.nan * rstar
		## interval I_0
		mask = (j==0)
		r[mask] = rstar[mask]
		## return
		return r

	## which interval based on rj
	def winterval(self, r):
		j = 0. * r
		return j

	## sign of f in interval I_j
	def sgnf(self, j):
		sgn = np.nan
		N = len(self.rj) - 2
		## valid index?
		if j in range(N+1):
			## is sgnf positive?
			plus = ( int(j) + int(self.first_trapped_j) ) % 2
			## fill value
			if plus==0:
				sgn = -1.
			if plus==1:
				sgn = 1.
		## return
		return sgn








