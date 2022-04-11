
"""
This module defines a class used for evaluating coordinate transformations at null shell junctions.
"""


import numpy as np
from xhorizon.shell_junction import interpolators as interp
from xhorizon.shell_junction.helpers import *




class active_slice:

	"""
	Class for handling shell and corner slicing of SSS regions. Given the region and the slice parameters,
	reference arrays are created for all desired transformations, which are then defined by interpolation
	and extrapolation of the reference arrays.

	Unlike passive_slice, this will actively obtain new functions U(udl),V(vdl) based on the inputs U0,V0, rather
	than just reading the existing ones.

	Designed to simultaneously incorporate both shell and corner junctions. Input determines which behavior will
	take effect. Offers protection from bad inputs in most cases. The exception is that the correct ublocks and 
	vblocks have to be given or everything will come out whacky.

	Slice Location:
		The location of the slice is determined by the values r0,u0,v0. These describe a point with coordinates
		(u0,v0) in some block, and with radius r0. The values are forced to be self-consistent by the algorithm.
		In corner junction mode, this is the corner point. Otherwise, it just describes a point on the shell.
		In shell mode, r0 and (u0 xor v0) can be nan, in which case the parameters just specify the shell.
		To fully specify the point we need to know what block it's in. So correct ublock and vblock args are needed.

	Slice Location Input Requirements:
		u0==finite xor v0==finite                   # avoids overspecifying a nonphysical point
		(optional in shell mode) r0==finite         # if finite, used to determine the other one of u0,v0

	Corner Junction Mode Input Requirements:
		r0==finite                  # corner junction point radius
		u0==finite xor v0==finite   # corner junction coordinate value
		U0 not None                 # function U0(r) at v=v0
		V0 not None                 # function V0(r) at u=u0
		ublocks==[correct]          # the correct list of block indices must be provided
		vblocks==[correct]          # the correct list of block indices must be provided

	Shell Junction Mode Input Requirements:
		For example, suppose we want a shell junction at v=v0.
			(v0==finite) xor (r0==finite and u0==finite)  # either of these uniquely determines the value v=v0 for the shell
			U0 not None                                   # function U0(r) at v=v0 shell
			vblocks==[correct]                            # the correct list of block indices must be provided
		For shell at u=u0, just switch u's and v's.

	Inputs:
		reg = The region being sliced.
		ublocks = List of indices (corresponding to reg.blocks) for the blocks containing the u=u0 slice.
		vblocks = List of indices (corresponding to reg.blocks) for the blocks containing the v=v0 slice.
		r0 = Radius of the junction corner point. Or, radius of the specified point on the shell in shell mode. See above.
		u0 = Slice location coordinate. See above.
		v0 = Slice location coordinate. See above.
		U0 = Function U0(r) at v=v0, or None. When None, U(udl) = udl.
		V0 = Function V0(r) at u=u0, or None. When None, V(vdl) = vdl.
		mu = Extrapolation parameter. When mu=0, extrapolation is just translation. When mu->inf, extrapolation is smooth linear. 
				(See interpolators module).

	Methods:
		Provides various coordinate transformation methods, which can be used by a region.

	Attributes:
		Input parameters as well as various reference arrays.

	"""

	def __init__(self, reg, ublocks=[], vblocks=[], r0=np.nan, u0=np.nan, v0=np.nan, U0=None, V0=None, mu=0., r_refs=[]):

		## process and store input values
		self.reg = reg
		self.ublocks, self.vblocks = bcon(reg,ublocks), bcon(reg,vblocks)
		self.r0, self.u0, self.v0 = set_ruv0(reg, r0=r0, u0=u0, v0=v0)
		self.U0, self.V0 = U0, V0
		self.mu = 1.*float(mu)

		## get ref arrays
		self.r = get_r_ref(self.reg, r_refs, self.r0)
		self.r, self.uvdl_u0, self.uvdl_v0 = uvdl_of_r_at_uv0(self.r, self.reg, ublocks=self.ublocks, vblocks=self.vblocks, u0=self.u0, v0=self.v0)
		self.U_v0 = U_of_udl_at_v0(self.r, self.uvdl_v0[0], self.U0)
		self.V_u0 = V_of_vdl_at_u0(self.r, self.uvdl_u0[1], self.V0)

		## use interpolated functions to get UV of ref arrays
		if np.isfinite(r0) and len(ublocks)>0 and len(vblocks)>0:
			self.UV_u0 = self.UV_of_uvdl(self.uvdl_u0)
			self.UV_v0 = self.UV_of_uvdl(self.uvdl_v0)

	"""
	Coordinate transformation methods.
	"""

	def U_of_udl_at_v0(self, udl):
		"""
		Evaluate the function U(udl) = U(r(udl,vdl0)) by interpolating from stored reference values.
		Smooth interpolation and extrapolation handled by interpolators module.
		"""
		return interp.interp_with_smooth_extrap(udl, self.uvdl_v0[0], self.U_v0, mu=self.mu)

	def V_of_vdl_at_u0(self, vdl):
		"""
		Evaluate the function V(vdl) = V(r(udl0,vdl)) by interpolating from stored reference values.
		Smooth interpolation and extrapolation handled by interpolators module.
		"""
		return interp.interp_with_smooth_extrap(vdl, self.uvdl_u0[1], self.V_u0, mu=self.mu)

	def UV_of_uvdl(self, uvdl):
		"""
		Combine U(udl) and V(vdl) into UV(uvdl).
		"""
		U_temp = self.U_of_udl_at_v0(uvdl[0])
		V_temp = self.V_of_vdl_at_u0(uvdl[1])
		UV_temp = np.array([U_temp,V_temp])
		return UV_temp

	def U_of_r_at_v0(self, r):
		"""
		Evaluate the function U(r) at v0 by interpolating from stored reference values.
		Smooth interpolation and extrapolation handled by interpolators module.
		Setting mu to nan because this should not be extrapolated.
		"""
		return interp.interp_with_smooth_extrap(r, self.r, self.U_v0, mu=np.nan)

	def V_of_r_at_u0(self, r):
		"""
		Evaluate the function V(r) at u0 by interpolating from stored reference values.
		Smooth interpolation and extrapolation handled by interpolators module.
		Setting mu to nan because this should not be extrapolated.
		"""
		return interp.interp_with_smooth_extrap(r, self.r, self.V_u0, mu=np.nan)






