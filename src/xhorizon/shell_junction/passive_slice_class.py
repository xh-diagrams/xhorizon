
"""
This module defines a class used for evaluating coordinate transformations at null shell junctions.
"""


import numpy as np
import xhorizon.shell_junction.interpolators as interp
from xhorizon.shell_junction.helpers import *



class passive_slice:

	"""
	Similar to active slice, except instead of setting new U(udl) and V(vdl) function, this one merely reads the
	existing ones. This is used to obtain functions to act as input for active slice U0,V0 params.
	"""

	def __init__(self, reg, ublocks=[], vblocks=[], r0=np.nan, u0=np.nan, v0=np.nan, mu=0., r_refs=[]):

		## process and store input values
		self.reg = reg
		self.ublocks, self.vblocks = bcon(reg,ublocks), bcon(reg,vblocks)
		self.r0, self.u0, self.v0 = set_ruv0(reg, r0=r0, u0=u0, v0=v0)
		self.mu = 1.*float(mu)

		## get r, uvdl ref arrays
		self.r = get_r_ref(self.reg, r_refs, self.r0)
		self.r, self.uvdl_u0, self.uvdl_v0 = uvdl_of_r_at_uv0(self.r, self.reg, ublocks=self.ublocks, vblocks=self.vblocks, u0=self.u0, v0=self.v0)
		
		## get UV ref arrays
		self.UV_u0 = reg.UV_of_uvdl(self.uvdl_u0)
		self.UV_v0 = reg.UV_of_uvdl(self.uvdl_v0)

	"""
	Coordinate transformation methods.
	"""

	def U_of_udl_at_v0(self, udl):
		"""
		Evaluate the function U(udl,vdl0) by interpolating from stored reference values.
		Smooth interpolation and extrapolation handled by interpolators module.
		"""
		return interp.interp_with_smooth_extrap(udl, self.uvdl_v0[0], self.UV_v0[0], mu=self.mu)

	def V_of_vdl_at_u0(self, vdl):
		"""
		Evaluate the function V(udl0,vdl) by interpolating from stored reference values.
		Smooth interpolation and extrapolation handled by interpolators module.
		"""
		return interp.interp_with_smooth_extrap(vdl, self.uvdl_u0[1], self.UV_u0[1], mu=self.mu)


	def U_of_r_at_v0(self, r):
		"""
		Evaluate the function U(r) at v0 by interpolating from stored reference values.
		Smooth interpolation and extrapolation handled by interpolators module.
		Setting mu to nan because this should not be extrapolated.
		"""
		return interp.interp_with_smooth_extrap(r, self.r, self.UV_v0[0], mu=np.nan)

	def V_of_r_at_u0(self, r):
		"""
		Evaluate the function V(r) at u0 by interpolating from stored reference values.
		Smooth interpolation and extrapolation handled by interpolators module.
		Setting mu to nan because this should not be extrapolated.
		"""
		return interp.interp_with_smooth_extrap(r, self.r, self.UV_u0[1], mu=np.nan)


	def udl_of_U_at_v0(self, U):
		"""
		Untested.
		Inverse of U(udl).
		"""
		return interp.interp_with_smooth_extrap(U, self.UV_v0[0], self.uvdl_v0[0], mu=self.mu)

	def vdl_of_V_at_u0(self, V):
		"""
		Untested.
		Inverse of V(vdl).
		"""
		return interp.interp_with_smooth_extrap(V, self.UV_u0[1], self.uvdl_u0[1], mu=self.mu)


