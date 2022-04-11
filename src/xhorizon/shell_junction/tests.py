
"""
Tests for shell_junction modules.
"""


import numpy as np

from xhorizon.shell_junction.active_slice_class import active_slice
from xhorizon.shell_junction.passive_slice_class import passive_slice
from xhorizon.shell_junction.helpers import *




"""
Tests. Run when __name__=='__main__'.
"""


def test1():
	"""
	Test uvdl_of_r_at_uv0.
	"""
	##
	print("\nTEST 1\n")
	## params
	func = xh.mf.hayward()
	ublocks, vblocks = [2,5,3], [0,1,2]
	u0, v0 = 3.*np.random.randn(2)
	## go
	reg = xh.reg.MAXreg(func)
	r = reg.metfunc.r_ref
	#r = r[r<10.]
	r, uvdl_u0, uvdl_v0 = uvdl_of_r_at_uv0(r, reg, ublocks=ublocks, vblocks=vblocks, v0=v0, u0=u0)
	## print
	print("\n")
	print("u0, v0 = %s, %s"%(u0,v0))
	print("r, uvdl_u0, uvdl_v0 = ...")
	print(r)
	print(uvdl_u0)
	print(uvdl_v0)
	print("\n")
	## plot
	plt.figure(1, figsize=(12,5))
	## left subfig
	plt.subplot(121)
	plt.gca().set_aspect('equal')
	## plot region
	reg.rplot()
	## fill blocks
	cols = ['b','y']
	for b in bcon(reg,ublocks):
		sty = dict(fc=cols[0], alpha=0.3)
		b.fill(sty=sty)
	for b in bcon(reg,vblocks):
		sty = dict(fc=cols[1], alpha=0.3)
		b.fill(sty=sty)
	## plot slices
	i = 0
	cmaps = [plt.cm.Blues, plt.cm.Reds]
	for uvdl in [uvdl_u0, uvdl_v0]:
		x, y = uvdl[1] - uvdl[0], uvdl[1] + uvdl[0]
		plt.scatter(x, y, marker='x', s=5, c=-r, cmap=cmaps[i%len(cmaps)], zorder=5000)
		i+=1
	## right subfig
	plt.subplot(122)
	plt.xlabel('r')
	plt.grid()
	## plot slices
	i = 0
	cmaps = [plt.cm.Blues, plt.cm.Reds]
	for uvdl in [uvdl_u0, uvdl_v0]:
		plt.scatter(r, uvdl[(i+0)%2], marker='x', s=1, c=-r, cmap=cmaps[i%2], zorder=5000)
		plt.scatter(r, uvdl[(i+1)%2], marker='x', s=10, c=-r, cmap=cmaps[i%2], zorder=5000)
		i+=1
	##
	plt.show()
	##
	print("\nEND TEST 1\n")



def test2():
	"""
	Test set_ruv0.
	"""
	##
	print("\nTEST 2\n")
	##
	reg = xh.reg.MAXreg(xh.mf.minkowski())
	r0, u0, v0 = set_ruv0(reg, r0=5, u0=8)
	print('r0=%s, u0=%s, v0=%s'%(r0,u0,v0))
	##
	print("\nEND TEST 2\n")





def test3():
	"""
	Test active_slice plotting uvdl ref values.
	"""
	##
	print("\nTEST 3\n")
	## define func and region
	func = xh.mf.schwarzschild()
	reg = xh.reg.MAXreg(func,rlines=False)
	## set input params and evaluate
	r0, u0, v0 = 1.1, np.nan, -2.
	ublocks, vblocks = [1,2], [0,1]
	slc = active_slice(reg, ublocks=ublocks, vblocks=vblocks, r0=r0, u0=u0, v0=v0)
	##
	## add ref lines to reg
	sty = sty=dict(c='lime', lw=1, alpha=0.9, zorder=6000)
	for b in slc.ublocks:
		if np.isfinite(r0):
			b.add_curves_tr(xh.cm.rlines([r0], sty=sty))
		if np.isfinite(u0):
			b.add_curves_uv(xh.cm.uvlines([u0], uv='u', sty=sty))
	for b in slc.vblocks:
		if np.isfinite(r0):
			b.add_curves_tr(xh.cm.rlines([r0], sty=sty))
		if np.isfinite(v0):
			b.add_curves_uv(xh.cm.uvlines([v0], uv='v', sty=sty))
	## plot
	plt.figure(1, figsize=(12,5))
	## left subfig
	plt.subplot(121)
	plt.gca().set_aspect('equal')
	## plot region
	reg.rplot()
	## fill blocks
	cols = ['b','y']
	for b in bcon(reg,ublocks):
		sty = dict(fc=cols[0], alpha=0.3)
		b.fill(sty=sty)
	for b in bcon(reg,vblocks):
		sty = dict(fc=cols[1], alpha=0.3)
		b.fill(sty=sty)
	## plot ref array slices
	if True:
		cvals = -np.arange(len(slc.r))
		cmaps = [plt.cm.Blues, plt.cm.Reds]
		plt.scatter(slc.uvdl_u0[1]-slc.uvdl_u0[0],slc.uvdl_u0[1]+slc.uvdl_u0[0],marker='x',s=10,c=cvals,cmap=cmaps[0],zorder=5000)
		plt.scatter(slc.uvdl_v0[1]-slc.uvdl_v0[0],slc.uvdl_v0[1]+slc.uvdl_v0[0],marker='x',s=10,c=cvals,cmap=cmaps[1],zorder=5000)
	## right subfig
	plt.subplot(122)
	plt.xlabel('r')
	plt.grid()
	## plot slices
	if True:
		i = 0
		cvals = -np.arange(len(slc.r))
		cmaps = [plt.cm.Blues, plt.cm.Reds]
		for uvdl in [slc.uvdl_u0, slc.uvdl_v0]:
			plt.scatter(slc.r, uvdl[(i+0)%2], marker='x', s=1 , c=cvals, cmap=cmaps[i%2], zorder=5000)
			plt.scatter(slc.r, uvdl[(i+1)%2], marker='x', s=10, c=cvals, cmap=cmaps[i%2], zorder=5000)
			i+=1
	##
	plt.show()
	print("\nEND TEST 3\n")





def test4():
	"""
	Test active_slice UV ref values.
	"""
	##
	print("\nTEST 4\n")
	## define func and region
	func = xh.mf.schwarzschild()
	reg = xh.reg.MAXreg(func,rlines=False)
	## set input params and evaluate
	U0 = lambda r: - np.pi**(-1) * np.arctan(r)
	V0 = lambda r:   np.pi**(-1) * np.arctan(r)
	r0, u0, v0 = 1.5, np.nan, -2.
	ublocks, vblocks = [1,2], [0,1]
	slc = active_slice(reg, ublocks=ublocks, vblocks=vblocks, r0=r0, u0=u0, v0=v0, U0=U0, V0=V0)
	## add ref lines to reg
	sty = sty=dict(c='lime', lw=1, alpha=0.9, zorder=6000)
	for b in slc.ublocks:
		if np.isfinite(r0):
			b.add_curves_tr(xh.cm.rlines([r0], sty=sty))
		if np.isfinite(u0):
			b.add_curves_uv(xh.cm.uvlines([u0], uv='u', sty=sty))
	for b in slc.vblocks:
		if np.isfinite(r0):
			b.add_curves_tr(xh.cm.rlines([r0], sty=sty))
		if np.isfinite(v0):
			b.add_curves_uv(xh.cm.uvlines([v0], uv='v', sty=sty))
	##
	## plot
	if True:
		plt.figure(1, figsize=(12,5))
		## left subfig
		plt.subplot(121)
		plt.gca().set_aspect('equal')
		## plot region
		reg.rplot()
		## fill blocks
		cols = ['b','y']
		for b in bcon(reg,ublocks):
			sty = dict(fc=cols[0], alpha=0.3)
			b.fill(sty=sty)
		for b in bcon(reg,vblocks):
			sty = dict(fc=cols[1], alpha=0.3)
			b.fill(sty=sty)
		## plot ref array slices
		if True:
			cvals = -np.arange(len(slc.r))
			cmaps = [plt.cm.Blues, plt.cm.Reds]
			## u0,V0 slice by replacing vdl_u0=uvdl_u0[1] by V_u0
			plt.scatter(slc.V_u0-slc.uvdl_u0[0],slc.V_u0+slc.uvdl_u0[0],marker='x',s=10,c=cvals,cmap=cmaps[0],zorder=5000)
			## v0,U0 slice by replacing udl_v0=uvdl_v0[0] by U_v0
			plt.scatter(slc.uvdl_v0[1]-slc.U_v0,slc.uvdl_v0[1]+slc.U_v0,marker='x',s=10,c=cvals,cmap=cmaps[1],zorder=5000)
		## right subfig
		plt.subplot(122)
		plt.xlabel('r')
		plt.grid()
		## plot slices
		if False:
			i = 0
			cvals = -np.arange(len(slc.r))
			cmaps = [plt.cm.Blues, plt.cm.Reds]
			for uvdl in [slc.uvdl_u0, slc.uvdl_v0]:
				plt.scatter(slc.r, uvdl[(i+0)%2], marker='x', s=1 , c=cvals, cmap=cmaps[i%2], zorder=5000)
				plt.scatter(slc.r, uvdl[(i+1)%2], marker='x', s=10, c=cvals, cmap=cmaps[i%2], zorder=5000)
				i+=1
		##
		plt.show()
	print("\nEND TEST 4\n")



def test5():
	"""
	Test active_slice.
	"""
	## import
	## start
	print("\nTEST 6\n")
	## define func and region
	func = xh.mf.schwarzschild()
	reg = xh.reg.MAXreg(func,rlines=True)
	## set input params and evaluate
	U0 = lambda r: - np.pi**(-1) * np.arctan(r)
	V0 = lambda r:   np.pi**(-1) * np.arctan(r)
	r0, u0, v0 = 2.1, np.nan, 2.
	ublocks, vblocks = [1,2], [0,1]
	slc = active_slice(reg, ublocks=ublocks, vblocks=vblocks, r0=r0, u0=u0, v0=v0, U0=U0, V0=V0)
	## play with slice object
	if False:
		plt.figure(1)
		sq = 1.5
		plt.xlim(-sq,sq)
		plt.ylim(-sq,sq)
		plt.gca().set_aspect('equal')
		plt.grid()
		plt.plot(slc.uvdl_v0[0], slc.UV_v0[0], 'rx')
		plt.plot(slc.uvdl_u0[1], slc.UV_u0[1], 'bx')
		plt.show()
	## done
	print("\nEND TEST 5\n")




def test6():
	"""
	Test active_slice.
	Schwarzschild with slice outside horizon.
	"""
	##
	print("\nTEST 6\n")
	## define func and region
	func = xh.mf.schwarzschild()
	reg = xh.reg.MAXreg(func,rlines=True)
	## set input params
	U0 = lambda r: - np.pi**(-1) * np.arctan(r)
	V0 = lambda r:   np.pi**(-1) * np.arctan(r)
	## r>1
	r0, u0, v0 = 2.1, 4., np.nan
	ublocks, vblocks = [1,2], [0,1]
	## evaluate slice
	slc = active_slice(reg, ublocks=ublocks, vblocks=vblocks, r0=r0, u0=u0, v0=v0, U0=U0, V0=V0, mu=0.)
	## update region UV transformations
	if True:
		reg.U_of_udl = slc.U_of_udl_at_v0
		reg.V_of_vdl = slc.V_of_vdl_at_u0
	## add ref lines to reg
	sty = sty=dict(c='lime', lw=1, alpha=0.9, zorder=6000)
	for b in slc.ublocks:
		if np.isfinite(r0):
			b.add_curves_tr(xh.cm.rlines([r0], sty=sty))
		if np.isfinite(u0):
			b.add_curves_uv(xh.cm.uvlines([u0], uv='u', sty=sty))
	for b in slc.vblocks:
		if np.isfinite(r0):
			b.add_curves_tr(xh.cm.rlines([r0], sty=sty))
		if np.isfinite(v0):
			b.add_curves_uv(xh.cm.uvlines([v0], uv='v', sty=sty))
	##
	## plot
	if True:
		plt.figure(1, figsize=(12,5))
		## left subfig
		plt.subplot(121)
		plt.gca().set_aspect('equal')
		## plot region
		reg.rplot()
		## fill blocks
		cols = ['b','y']
		for b in bcon(reg,ublocks):
			sty = dict(fc=cols[0], alpha=0.3)
			b.fill(sty=sty)
		for b in bcon(reg,vblocks):
			sty = dict(fc=cols[1], alpha=0.3)
			b.fill(sty=sty)
		## plot ref array slices
		if True:
			cvals = -np.arange(len(slc.r))
			cmaps = [plt.cm.Blues, plt.cm.Reds]
			## plot region in UV coords
			plt.scatter(slc.UV_u0[1]-slc.UV_u0[0],slc.UV_u0[1]+slc.UV_u0[0],marker='x',s=10,c=cvals,cmap=cmaps[0],zorder=5000)
			plt.scatter(slc.UV_v0[1]-slc.UV_v0[0],slc.UV_v0[1]+slc.UV_v0[0],marker='x',s=10,c=cvals,cmap=cmaps[1],zorder=5000)
		## right subfig
		plt.subplot(122)
		plt.xlabel('uvdl')
		plt.ylabel('UV')
		plt.grid()
		## plot interpolated U(udl) and V(vdl)
		if True:
			xdl = 1.1*np.linspace(-1,1,5001)
			plt.scatter(xdl, slc.V_of_vdl_at_u0(xdl), marker='x', s=5, linestyle='-', c='b', zorder=5000)
			plt.scatter(xdl, slc.U_of_udl_at_v0(xdl), marker='x', s=5, linestyle='-', c='r', zorder=5000)
			
		## plot reference values U_v0(udl_v0) and V_v0(vdl_v0)
		if True:
			cvals = -np.arange(len(slc.r))
			cmaps = [plt.cm.Blues, plt.cm.Reds]
			plt.scatter(slc.uvdl_u0[1], slc.UV_u0[1], marker='x', s=10 , c=cvals, cmap=cmaps[0], zorder=5100)
			plt.scatter(slc.uvdl_v0[0], slc.UV_v0[0], marker='x', s=10 , c=cvals, cmap=cmaps[1], zorder=5100)
		##
		plt.show()
	print("\nEND TEST 6\n")




def test7():
	"""
	Test active_slice.
	Schwarzschild with slice inside horizon.
	"""
	##
	print("\nTEST 7\n")
	## define func and region
	func = xh.mf.schwarzschild()
	reg = xh.reg.MAXreg(func,rlines=True)
	## set input params
	U0 = lambda r: - np.pi**(-1) * np.arctan(r)
	V0 = lambda r: - np.pi**(-1) * np.arctan(r)
	## r>1
	r0, u0, v0 = 0.5, 4., np.nan
	ublocks, vblocks = [0,3], [0,1]
	## evaluate slice
	slc = active_slice(reg, ublocks=ublocks, vblocks=vblocks, r0=r0, u0=u0, v0=v0, U0=U0, V0=V0, mu=0.)
	## update region UV transformations
	if True:
		reg.U_of_udl = slc.U_of_udl_at_v0
		reg.V_of_vdl = slc.V_of_vdl_at_u0
	## add ref lines to reg
	sty = sty=dict(c='lime', lw=1, alpha=0.9, zorder=6000)
	for b in slc.ublocks:
		if np.isfinite(r0):
			b.add_curves_tr(xh.cm.rlines([r0], sty=sty))
		if np.isfinite(u0):
			b.add_curves_uv(xh.cm.uvlines([u0], uv='u', sty=sty))
	for b in slc.vblocks:
		if np.isfinite(r0):
			b.add_curves_tr(xh.cm.rlines([r0], sty=sty))
		if np.isfinite(v0):
			b.add_curves_uv(xh.cm.uvlines([v0], uv='v', sty=sty))
	##
	## plot
	if True:
		plt.figure(1, figsize=(12,5))
		## left subfig
		plt.subplot(121)
		plt.xlim(-1,1)
		plt.ylim(-1,1)
		plt.gca().set_aspect('equal')
		plt.grid()
		## plot region
		reg.rplot()
		## fill blocks
		cols = ['b','y']
		for b in bcon(reg,ublocks):
			sty = dict(fc=cols[0], alpha=0.3)
			b.fill(sty=sty)
		for b in bcon(reg,vblocks):
			sty = dict(fc=cols[1], alpha=0.3)
			b.fill(sty=sty)
		## plot ref array slices
		if True:
			cvals = -np.arange(len(slc.r))
			cmaps = [plt.cm.Blues, plt.cm.Reds]
			## plot region in UV coords
			plt.scatter(slc.UV_u0[1]-slc.UV_u0[0],slc.UV_u0[1]+slc.UV_u0[0],marker='x',s=10,c=cvals,cmap=cmaps[0],zorder=5000)
			plt.scatter(slc.UV_v0[1]-slc.UV_v0[0],slc.UV_v0[1]+slc.UV_v0[0],marker='x',s=10,c=cvals,cmap=cmaps[1],zorder=5000)
		## right subfig
		plt.subplot(122)
		plt.xlabel('uvdl')
		plt.ylabel('UV')
		plt.grid()
		## plot interpolated U(udl) and V(vdl)
		if True:
			xdl = 1.1*np.linspace(-1,1,5001)
			plt.scatter(xdl, slc.V_of_vdl_at_u0(xdl), marker='x', s=5, linestyle='-', c='b', zorder=5000)
			plt.scatter(xdl, slc.U_of_udl_at_v0(xdl), marker='x', s=5, linestyle='-', c='r', zorder=5000)
			
		## plot reference values U_v0(udl_v0) and V_v0(vdl_v0)
		if True:
			cvals = -np.arange(len(slc.r))
			cmaps = [plt.cm.Blues, plt.cm.Reds]
			plt.scatter(slc.uvdl_u0[1], slc.UV_u0[1], marker='x', s=10 , c=cvals, cmap=cmaps[0], zorder=5100)
			plt.scatter(slc.uvdl_v0[0], slc.UV_v0[0], marker='x', s=10 , c=cvals, cmap=cmaps[1], zorder=5100)
		##
		plt.show()
	print("\nEND TEST 7\n")



def test8():
	"""
	Test active_slice.
	Hayward with slice outside horizons.
	"""
	##
	print("\nTEST 8\n")
	## define func and region
	func = xh.mf.hayward()
	reg = xh.reg.MAXreg(func,rlines=True)
	## set input params
	U0 = lambda r: - np.pi**(-1) * np.arctan(r)
	V0 = lambda r:   np.pi**(-1) * np.arctan(r)
	## r>1
	r0, u0, v0 = 2.5, 1., np.nan
	ublocks, vblocks = [2,5,3], [0,1,2]
	## evaluate slice
	slc = active_slice(reg, ublocks=ublocks, vblocks=vblocks, r0=r0, u0=u0, v0=v0, U0=U0, V0=V0, mu=0.1)
	## update region UV transformations
	if True:
		reg.U_of_udl = slc.U_of_udl_at_v0
		reg.V_of_vdl = slc.V_of_vdl_at_u0
	## add ref lines to reg
	sty = sty=dict(c='lime', lw=1, alpha=0.9, zorder=6000)
	for b in slc.ublocks:
		if np.isfinite(r0):
			b.add_curves_tr(xh.cm.rlines([r0], sty=sty))
		if np.isfinite(u0):
			b.add_curves_uv(xh.cm.uvlines([u0], uv='u', sty=sty))
	for b in slc.vblocks:
		if np.isfinite(r0):
			b.add_curves_tr(xh.cm.rlines([r0], sty=sty))
		if np.isfinite(v0):
			b.add_curves_uv(xh.cm.uvlines([v0], uv='v', sty=sty))
	##
	## plot
	if True:
		plt.figure(1, figsize=(12,5))
		## left subfig
		plt.subplot(121)
		plt.xlim(-2,2)
		plt.ylim(-2,2)
		plt.gca().set_aspect('equal')
		plt.grid()
		## plot region
		reg.rplot()
		## fill blocks
		cols = ['b','y']
		for b in bcon(reg,ublocks):
			sty = dict(fc=cols[0], alpha=0.3)
			b.fill(sty=sty)
		for b in bcon(reg,vblocks):
			sty = dict(fc=cols[1], alpha=0.3)
			b.fill(sty=sty)
		## plot ref array slices
		if True:
			cvals = -np.arange(len(slc.r))
			cmaps = [plt.cm.Blues, plt.cm.Reds]
			## plot region in UV coords
			plt.scatter(slc.UV_u0[1]-slc.UV_u0[0],slc.UV_u0[1]+slc.UV_u0[0],marker='x',s=10,c=cvals,cmap=cmaps[0],zorder=5000)
			plt.scatter(slc.UV_v0[1]-slc.UV_v0[0],slc.UV_v0[1]+slc.UV_v0[0],marker='x',s=10,c=cvals,cmap=cmaps[1],zorder=5000)
		## right subfig
		plt.subplot(122)
		plt.xlabel('uvdl')
		plt.ylabel('UV')
		plt.grid()
		## plot interpolated U(udl) and V(vdl)
		if True:
			xdl = 2.5*np.linspace(-1,1,5001)
			plt.scatter(xdl, slc.V_of_vdl_at_u0(xdl), marker='x', s=5, linestyle='-', c='b', zorder=5000)
			plt.scatter(xdl, slc.U_of_udl_at_v0(xdl), marker='x', s=5, linestyle='-', c='r', zorder=5000)
			
		## plot reference values U_v0(udl_v0) and V_v0(vdl_v0)
		if True:
			cvals = -np.arange(len(slc.r))
			cmaps = [plt.cm.Blues, plt.cm.Reds]
			plt.scatter(slc.uvdl_u0[1], slc.UV_u0[1], marker='x', s=10 , c=cvals, cmap=cmaps[0], zorder=5100)
			plt.scatter(slc.uvdl_v0[0], slc.UV_v0[0], marker='x', s=10 , c=cvals, cmap=cmaps[1], zorder=5100)
		##
		plt.show()
	print("\nEND TEST 8\n")



def test9():
	"""
	Test active_slice.
	Check shell functionality when U0 or V0 is None.
	Seems to be working perfectly.
	Both None leaves UV=uvdl. One or the other does the shell transformation. Both does corner.
	"""
	##
	print("\nTEST 9\n")
	## define func and region
	func = xh.mf.schwarzschild()
	reg = xh.reg.MAXreg(func,rlines=True)
	## set input params
	U0 = lambda r: - np.pi**(-1) * np.arctan(r)
	V0 = lambda r:   np.pi**(-1) * np.arctan(r)
	## r>1
	r0, u0, v0 = 2.053, np.nan, -1.5
	ublocks, vblocks = [1,2], [0,1]
	## evaluate slice
	slc = active_slice(reg, ublocks=ublocks, vblocks=vblocks, r0=r0, u0=u0, v0=v0, U0=None, V0=V0, mu=0.)
	## update region UV transformations
	if True:
		reg.U_of_udl = slc.U_of_udl_at_v0
		reg.V_of_vdl = slc.V_of_vdl_at_u0
	## add ref lines to reg
	sty = sty=dict(c='lime', lw=1, alpha=0.9, zorder=6000)
	for b in slc.ublocks:
		if np.isfinite(r0):
			b.add_curves_tr(xh.cm.rlines([r0], sty=sty))
		if np.isfinite(u0):
			b.add_curves_uv(xh.cm.uvlines([u0], uv='u', sty=sty))
	for b in slc.vblocks:
		if np.isfinite(r0):
			b.add_curves_tr(xh.cm.rlines([r0], sty=sty))
		if np.isfinite(v0):
			b.add_curves_uv(xh.cm.uvlines([v0], uv='v', sty=sty))
	##
	## plot
	if True:
		plt.figure(1, figsize=(12,5))
		## left subfig
		plt.subplot(121)
		plt.gca().set_aspect('equal')
		## plot region
		reg.rplot()
		## fill blocks
		cols = ['b','y']
		for b in bcon(reg,ublocks):
			sty = dict(fc=cols[0], alpha=0.3)
			b.fill(sty=sty)
		for b in bcon(reg,vblocks):
			sty = dict(fc=cols[1], alpha=0.3)
			b.fill(sty=sty)
		## plot ref array slices
		if True:
			cvals = -np.arange(len(slc.r))
			cmaps = [plt.cm.Blues, plt.cm.Reds]
			## plot region in UV coords
			plt.scatter(slc.UV_u0[1]-slc.UV_u0[0],slc.UV_u0[1]+slc.UV_u0[0],marker='x',s=10,c=cvals,cmap=cmaps[0],zorder=5000)
			plt.scatter(slc.UV_v0[1]-slc.UV_v0[0],slc.UV_v0[1]+slc.UV_v0[0],marker='x',s=10,c=cvals,cmap=cmaps[1],zorder=5000)
		## right subfig
		plt.subplot(122)
		plt.xlabel('uvdl')
		plt.ylabel('UV')
		plt.grid()
		## plot interpolated U(udl) and V(vdl)
		if True:
			xdl = 1.1*np.linspace(-1,1,5001)
			plt.scatter(xdl, slc.V_of_vdl_at_u0(xdl), marker='x', s=5, linestyle='-', c='b', zorder=5000)
			plt.scatter(xdl, slc.U_of_udl_at_v0(xdl), marker='x', s=5, linestyle='-', c='r', zorder=5000)
			
		## plot reference values U_v0(udl_v0) and V_v0(vdl_v0)
		if True:
			cvals = -np.arange(len(slc.r))
			cmaps = [plt.cm.Blues, plt.cm.Reds]
			plt.scatter(slc.uvdl_u0[1], slc.UV_u0[1], marker='x', s=10 , c=cvals, cmap=cmaps[0], zorder=5100)
			plt.scatter(slc.uvdl_v0[0], slc.UV_v0[0], marker='x', s=10 , c=cvals, cmap=cmaps[1], zorder=5100)
		##
		plt.show()
	print("\nEND TEST 9\n")



def test10():
	"""
	Test passive_slice.
	Schwarz, outside horizon: working fine.
	Schwarz, inside horizon: working fine.
	Hayward, outside horizon: working fine.
	"""
	##
	print("\nTEST 10\n")
	## define func and region
	func = xh.mf.hayward()   ##schwarzschild()
	reg = xh.reg.MAXreg(func,rlines=True)
	## r>1
	r0, u0, v0 = 1.3, np.nan, 0.5
	ublocks, vblocks =  [2,5,3], [0,1,2]  ##[0,3], [0,1]
	##
	if True:
		reg.U_of_udl = lambda udl: np.arctan(udl)
		reg.V_of_vdl = lambda vdl: np.arctan(vdl)
	## evaluate slice
	slc = passive_slice(reg, ublocks=ublocks, vblocks=vblocks, r0=r0, u0=u0, v0=v0, mu=0.)
	## add ref lines to reg
	sty = sty=dict(c='lime', lw=1, alpha=0.9, zorder=6000)
	for b in slc.ublocks:
		if np.isfinite(r0):
			b.add_curves_tr(xh.cm.rlines([r0], sty=sty))
		if np.isfinite(u0):
			b.add_curves_uv(xh.cm.uvlines([u0], uv='u', sty=sty))
	for b in slc.vblocks:
		if np.isfinite(r0):
			b.add_curves_tr(xh.cm.rlines([r0], sty=sty))
		if np.isfinite(v0):
			b.add_curves_uv(xh.cm.uvlines([v0], uv='v', sty=sty))
	##
	## plot
	if True:
		plt.figure(1, figsize=(12,5))
		## left subfig
		plt.subplot(121)
		plt.gca().set_aspect('equal')
		## plot region
		reg.rplot()
		## fill blocks
		cols = ['b','y']
		for b in bcon(reg,ublocks):
			sty = dict(fc=cols[0], alpha=0.3)
			b.fill(sty=sty)
		for b in bcon(reg,vblocks):
			sty = dict(fc=cols[1], alpha=0.3)
			b.fill(sty=sty)
		## plot ref array slices
		if True:
			cvals = -np.arange(len(slc.r))
			cmaps = [plt.cm.Blues, plt.cm.Reds]
			## plot region in UV coords
			plt.scatter(slc.UV_u0[1]-slc.UV_u0[0],slc.UV_u0[1]+slc.UV_u0[0],marker='x',s=10,c=cvals,cmap=cmaps[0],zorder=5000)
			plt.scatter(slc.UV_v0[1]-slc.UV_v0[0],slc.UV_v0[1]+slc.UV_v0[0],marker='x',s=10,c=cvals,cmap=cmaps[1],zorder=5000)
		## right subfig
		plt.subplot(122)
		plt.xlabel('r')
		plt.ylabel('UV')
		#plt.xlim(0,2)
		plt.grid()
		## plot interpolated U(r) and V(r)
		if True:
			rr = 5.*np.linspace(0,1,5001)
			plt.scatter(rr, slc.V_of_r_at_u0(rr), marker='x', s=5, linestyle='-', c='b', zorder=5000)
			plt.scatter(rr, slc.U_of_r_at_v0(rr), marker='x', s=5, linestyle='-', c='r', zorder=5000)
			
		## plot reference values U_v0(r) and V_v0(r)
		if True:
			cvals = -np.arange(len(slc.r))
			cmaps = [plt.cm.Blues, plt.cm.Reds]
			plt.scatter(slc.r, slc.UV_u0[1], marker='x', s=10 , c=cvals, cmap=cmaps[0], zorder=5100)
			plt.scatter(slc.r, slc.UV_v0[0], marker='x', s=10 , c=cvals, cmap=cmaps[1], zorder=5100)
		##
		plt.show()
	print("\nEND TEST 10\n")




if __name__=="__main__":

	import matplotlib.pyplot as plt
	import xhorizon as xh

	print("\n\nTESTS\n")
	test9()
	print("\nEND TESTS\n")


