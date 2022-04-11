
"""
Helpers for evaluating coordinate transformations at null shell junctions.
"""


import numpy as np




"""
Helpers.
"""

def set_ruv0(reg, r0=np.nan, u0=np.nan, v0=np.nan):
	"""
	Self-consistently fill all three values based on input.
	For a shell junction, r0 should be nan and (u0 xor v0) should be finite.
	For a corner junction, r0 should be finite and (u0 xor v0) should be finite.
	Region arg is needed to provide c and F(r).
	If invalid, warn and return all nan values.
	"""
	## initialize
	ruv0 = [np.nan, np.nan, np.nan]
	F = reg.metfunc.F
	c = float(reg.rparams['c'])
	r0, u0, v0 = 1.*float(r0), 1.*float(u0), 1.*float(v0)
	## only go if (u0 xor v0) is finite
	if np.isfinite(u0)^np.isfinite(v0):
		## shell
		if not np.isfinite(r0):
			ruv0 = [np.nan, u0, v0]
		## corner
		if np.isfinite(r0):
			## u0 given
			if np.isfinite(u0):
				v0 = u0 + 2.*(F(r0)-c)
			## v0 given
			if np.isfinite(v0):
				u0 = v0 - 2.*(F(r0)-c)	
			## set values
			ruv0 = [r0, u0, v0]
	## return
	return ruv0




def uvdl_of_r_at_uv0(r, reg, ublocks=[], vblocks=[], u0=np.nan, v0=np.nan):
	"""
	Calculate uvdl of r at fixed u0 and v0 at block points. Slices in both directions are calculated.
		r = An array of r values at which to evaluate. Non block points return nan.
		reg = A region in which to evaluate.
		ublocks =  A list of integer indices for blocks in reg.blocks. The u0 slice is evaluated at block points in these blocks only.
		vblocks =  A list of integer indices for blocks in reg.blocks. The v0 slice is evaluated at block points in these blocks only.
			Alternately, block args can be provided directly as a list of blocks (bcon converter handles this).
			Block args necessary bc radius alone can be ambiguous in extended regions. 
		u0 = Fixed value of u0 at which the slice is evaluated.
		v0 = Fixed value of v0 at which the slice is evaluated.
	Returns (r, uvdl_u0, uvdl_v0).
	"""
	## initialize
	u0, v0 = 1.*float(u0), 1.*float(v0)
	uvdl_u0 = np.nan*np.ones((2,len(r)))
	uvdl_v0 = np.nan*np.ones((2,len(r)))
	u,v = np.nan*r, np.nan*r
	rstar = reg.metfunc.F(1.*r)
	c = 1.*float(reg.rparams['c'])
	## convert block args
	ublocks = bcon(reg, ublocks)
	vblocks = bcon(reg, vblocks)
	## u0 slice
	for b in ublocks:
		## mask points outside block and calculate
		mask = np.logical_and.reduce([b.rj[0]<r, r<b.rj[1]])
		u[mask] = u0 + 0.*rstar[mask]
		v[mask] = u0 + 2.*(rstar[mask]-c)
		uvdl_u0[:,mask] = b.uvdl_of_uv(np.array([u[mask],v[mask]]))
	## v0 slice
	for b in vblocks:
		## mask points outside block and calculate
		mask = np.logical_and.reduce([b.rj[0]<r, r<b.rj[1]])
		v[mask] = v0 + 0.*rstar[mask]
		u[mask] = v0 - 2.*(rstar[mask]-c)
		uvdl_v0[:,mask] = b.uvdl_of_uv(np.array([u[mask],v[mask]]))
	## return
	return r, uvdl_u0, uvdl_v0



def U_of_udl_at_v0(r, udl, U0):
	"""
	Implements the function U(udl) = U0( r(udl,vdl0) ) at v=v0, where U0(r) is a function of radius.
	If U0==None, shell mode is activated, and the function reduces to U(udl) = udl.
	Inputs:
		r = An array of radius values at points along v=v0.
		udl = An array of udl values corresponding to the points in r.
		U0 = For corner mode, a monotonic function U0(r). If U0==None, shell mode activates.
	Returns:
		U_v0 = An array of U values corresponding the points in r at v=v0.
	"""
	## initialize
	U_v0 = np.nan*r
	## corner mode
	if not U0==None:
		U_v0 = U0(1.*r)
	## shell mode
	if U0==None:
		U_v0 = 1.*udl
	## return
	return 1.*U_v0


def V_of_vdl_at_u0(r, vdl, V0):
	"""
	Implements the function V(vdl) = V0( r(udl0,vdl) ) at u=u0, where V0(r) is a function of radius.
	If V0==None, shell mode is activated, and the function reduces to V(vdl) = vdl.
	Inputs:
		r = An array of radius values at points along u=u0.
		vdl = An array of vdl values corresponding to the points in r.
		V0 = For corner mode, a monotonic function V0(r). If V0==None, shell mode activates.
	Returns:
		V_u0 = An array of V values corresponding the points in r at u=u0.
	
	This function is identical to U_of_udl_at_v0. To avoid redundancy, this one just calls the other.
	They are kept separate only to make the code where these are called more readable.
	"""
	## return
	return U_of_udl_at_v0(r, vdl, V0)



def bcon(reg, blocks):
	"""
	Block list converter.
	If blocks is a list of integers, convert to a list of blocks using reg.blocks.
	Else, do nothing.
	Return list of blocks.
	"""
	## go
	if len(blocks)>0 and type(blocks[0])==int:
		blocks = [ reg.blocks[blocks[i]] for i in range(len(blocks)) ]
	## return
	return blocks



def dtr_calc(f0):
	"""
	Given the values f(r0) in three of the regions, return the last one.
	The one to be calculated is specified by nan. 
	Order is [A,B,C,D].
	"""
	f0 = 1.*f0
	if not np.isfinite(f0[0]):
		f0[0] = f0[2]*f0[3] / f0[1]
	if not np.isfinite(f0[1]):
		f0[1] = f0[2]*f0[3] / f0[0]
	if not np.isfinite(f0[2]):
		f0[2] = f0[0]*f0[1] / f0[3]
	if not np.isfinite(f0[3]):
		f0[3] = f0[0]*f0[1] / f0[2]
	return f0


def schwarz_dtr_masses(m, r0=1.):
	"""
	Given three of four masses, return dtr compatible remaining mass.
	"""
	f0 = 1. - 2.*m/r0
	f0 = dtr_calc(f0)
	m = r0*(1.-f0)/2.
	return m



def slicecheck(sl,reg):
	if (np.abs(sl.u0)>=2.*reg.rparams['s0']) or (np.abs(sl.v0)>=2.*reg.rparams['s0']):
		print("Error: u0 or v0 out of range.")
		print(sl.u0, sl.v0)
		sl, reg = None, None
	return sl, reg


def get_r_ref(reg, r_refs, r0):
	"""
	"""
	## init
	rrlist = []
	## add reg metfunc r_ref
	rrlist += [1.*reg.metfunc.r_ref]
	## add r_refs values
	rrlist += [1.*rrr for rrr in r_refs]
	## add values near r0
	a, npoints = 1e-21, 500
	ss = np.linspace(0., np.log(a), npoints)
	#rrlist += [r0 + np.exp(1.*ss), r0 - np.exp(1.*ss)]
	## concatenate
	rr = np.concatenate(rrlist)
	## return
	return 1.*rr









