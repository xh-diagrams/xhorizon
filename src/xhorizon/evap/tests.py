


"""
Tests for evap module.
"""

import numpy as np
import matplotlib.pyplot as plt
import xhorizon as xh
import pprint
import datetime

from xhorizon.evap.evap import *
from xhorizon.evap.helpers import *


def main():
	#test1()
	#test2()
	#test3()
	#test4()
	#test5()
	#test6()
	#test7()
	#test8()
	#test9()
	#test10()
	#test11()
	#test12()
	#test13()
	#test14()
	#test15()
	#test16()
	pass



def test1():
	"""
	Test functionality of get_rinf_uv0.
	"""
	##
	print("\nTEST 1\n")
	## regions
	reg1 = xh.reg.EFreg(xh.mf.schwarzschild())
	reg2 = xh.reg.EFreg(xh.mf.hayward())
	reg3 = xh.reg.EFreg(xh.mf.minkowski())
	reglist = [reg1,reg2,reg3]
	## v values
	v0 = [-1.,1.,0.]
	## get
	rinf = get_rinf_uv0(reglist,v0=v0)
	## print
	print(rinf)
	##
	print("\nEND TEST 1\n")



def test2():
	"""
	Test functionality of accrete.
	"""
	##
	print("\nTEST 2\n")
	## create initial region
	reglist = [xh.reg.EFreg(xh.mf.minkowski(),rlines=False,boundary=False)]
	## create accreted regions
	reglist += xh.evap.accrete(reglist.pop(), v1=0., v2=0., R2=0.8)
	reglist += xh.evap.accrete(reglist.pop(), v1=.5, v2=.5, R2=.9)
	reglist += xh.evap.accrete(reglist.pop(), v1=1., v2=1., R2=1.)
	## add lines
	xh.evap.colorlines(reglist)
	xh.evap.boundarylines(reglist)
	## draw
	plt.figure()
	plt.gca().set_aspect('equal')
	for reg in reglist:
		reg.rplot()
	plt.show()
	##
	print("\nEND TEST 2\n")



def test3():
	"""
	Test functionality of accretion.
	"""
	##
	print("\nTEST 3\n")
	## params
	v = np.linspace(0,1,5)
	m = .2 + .8*v
	## create initial region
	reglist = accretion(m,v)
	## add lines
	xh.evap.colorlines(reglist)
	xh.evap.boundarylines(reglist)
	## draw
	plt.figure()
	plt.gca().set_aspect('equal')
	for reg in reglist:
		reg.rplot()
	plt.show()
	##
	print("\nEND TEST 3\n")





def test4():
	"""
	Test functionality of evap.
	"""
	##
	print("\nTEST 4\n")
	## create initial region
	reg0 = xh.reg.EFreg(xh.mf.hayward(R=1.,l=0.1),rlines=False,boundary=False)
	reglist = [reg0]
	## create evaporated regions
	reglist += xh.evap.evap(reglist.pop(), u1=10., v1=-10., u2=10., R2=0.9)
	#reglist += xh.evap.evap(reglist.pop(), u1=2., v1=2., u2=0., R2=0.98)
	## check
	check_uvr(reglist)
	## draw diagram?
	if True:
		## add lines
		xh.evap.colorlines(reglist, sty=dict(marker='.', markersize=.3, ls='none'), rmin=0.8, rmax=1.1, dr=.02, npoints=2001, inf=25.)
		xh.evap.boundarylines(reglist)
		xh.evap.s0_lines(reglist)
		## draw
		xh.newfig(tex=False,sqaxis=3)
		plt.title("Test 4")
		for reg in reglist:
			reg.rplot()
		## fill
		fill_by_R(reglist)
		## show plot
		plt.savefig("temp-figs/test4.png", dpi=400)
		##
	print("\nEND TEST 4\n")


def test5():
	"""
	Test functionality of evaporation().
	"""
	##
	print("\nTEST 5\n")
	## create initial region
	reg0 = xh.reg.EFreg(xh.mf.hayward(R=1.,l=0.1),rlines=False,boundary=False)
	reglist = [reg0]
	## init params
	u0 = 10.
	v0 = 0.
	## subsequent params
	Nreg = 1
	tt = np.linspace(0.,1.,Nreg+2)[1:-1]
	R = 1. * (1.-tt)**(1./3.)
	du = 0.4 + 0.*tt
	dv = 0.2 + 0.*tt
	## create evaporated regions
	reglist += xh.evap.evaporation(reglist.pop(), R=R, du=du, dv=dv, u0=u0, v0=v0, l=0.1, rparams={})
	## check uvr
	check_uvr(reglist)
	## draw regions?
	if True:
		## add lines
		xh.evap.colorlines(reglist)
		xh.evap.boundarylines(reglist)
		## draw
		xh.newfig(tex=False,sqaxis=3)
		plt.title('Test 5')
		for reg in reglist:
			reg.rplot()
		## fill
		fill_by_R(reglist, cm=plt.cm.prism)
		## show plot
		plt.savefig("temp-figs/test5.png", dpi=200)
		##
	print("\nEND TEST 5\n")
	## return
	return reglist


def test6():
	"""
	Test functionality of formevap().
	"""
	##
	print("\nTEST 6\n")
	## accrete params
	N = 4.
	R0, R = .8, 1.
	accrete_R = R0 + (R-R0)*np.linspace(0,1,N)
	v0, v = -1., -.5
	accrete_v = np.linspace(v0,v,N)
	## evap params
	evap_R = [.9 ,.0]
	evap_u = [6.,9.]
	## create evaporated regions
	reglist = xh.evap.formevap(accrete_R=accrete_R, accrete_v=accrete_v, evap_R=evap_R, evap_u=evap_u, l=0.01, rparams=dict(c=0., s0=10.))
	## resquish
	fU = None
	fV = None
	reglist = UVcompose(reglist, fU=fU, fV=fV)
	## add lines
	xh.evap.colorlines(reglist)
	xh.evap.boundarylines(reglist)
	## draw
	xh.newfig(tex=False,sqaxis=3)
	#plt.figure()
	plt.title('Test 6')
	plt.xlim(0.5,2.5)
	plt.ylim(-1,1)
	for reg in reglist:
		reg.rplot()
	## fill
	fill_by_R(reglist)
	## mink compare
	#xh.reg.EFreg(xh.mf.minkowski()).rplot()
	## show plot
	plt.savefig("temp-figs/test6.png", dpi=600)
	##
	print("\nEND TEST 6\n")


def test7():
	"""
	Test functionality of get_uvdl_of_UV().
	"""
	##
	print("\nTEST 7\n")
	## create evaporated regions
	reglist = xh.evap.formevap(rparams=dict(s0=10.))
	## resquish
	fU, fV = get_uvdl_of_UV(reglist[-1])
	## plot
	xx = np.linspace(-2,2,50001)
	plt.figure()
	plt.plot(xx, fU(xx), 'rx-')
	plt.plot(xx, fV(xx), 'bx-')
	plt.grid()
	plt.show()
	## edit regions
	UVcompose(reglist, fU=fU, fV=fV)
	## add lines
	xh.evap.colorlines(reglist)
	xh.evap.boundarylines(reglist)
	## draw
	xh.newfig(tex=False,sqaxis=3)
	plt.title('Test 7')
	for reg in reglist:
		reg.rplot()
	fill_by_R(reglist)
	## show plot
	plt.show()
	plt.savefig("temp-figs/test7.png", dpi=400)
	##
	print("\nEND TEST 7\n")



def test8():
	"""
	Test functionality of EF cornermask.
	"""
	##
	print("\nTEST 8\n")
	## make region
	func = xh.mf.schwarzschild()
	reg = xh.reg.EFreg(func)
	##
	if True:
		abcd = 'a'
		reg = xh.cornermask.EFreg(reg, abcd=abcd, u0=0., v0=0.)
	##
	plt.figure()
	plt.xlim(-3,3)
	plt.gca().set_aspect('equal')
	reg.rplot()
	plt.show()
	##
	print("\nEND TEST 8\n")




def test9():
	"""
	Test functionality of split abcd helper.
	"""
	##
	print("\nTEST 9\n")
	## make region
	func = xh.mf.schwarzschild()
	reg = xh.reg.EFreg(func)
	##
	if True:
		abcd = 'bcd'
		reglist = xh.evap.split_reg_abcd(reg, abcd=abcd, u0=0., v0=0.)
	##
	print(len(reglist))
	print(reglist)
	##
	plt.figure()
	plt.xlim(-3,3)
	plt.gca().set_aspect('equal')
	for reg in reglist:
		reg.rplot()
	plt.show()
	##
	print("\nEND TEST 9\n")



def test10():
	"""
	Test direct formation of joint. Formerly was evap.debug.go().
	"""
	##
	print("\nTEST 10\n")
	for ux in [15.]:
		## params
		R1 = 1.
		R2 = .5
		l = .1
		rparams = dict(s0=10.)
		u1 = ux
		v1 = 0.
		u2 = 1.*u1
		## funcs
		if False:
			func1 = xh.mf.hayward(R=R1,l=l)
			func2 = xh.mf.hayward(R=R2,l=l)
		else:
			func1 = xh.mf.schwarzschild(R=R1)
			func2 = xh.mf.schwarzschild(R=R2)
		## regs
		reg1 = xh.reg.EFreg(func1,rparams=rparams,rlines=False,boundary=False)
		reg2 = xh.reg.EFreg(func2,rparams=rparams,rlines=False,boundary=False)
		## horizon radii
		Rh1, Rh2 = reg1.blocks[-1].rj[0], reg2.blocks[-1].rj[0]
		## passive slice of reg1
		r1 = reg1.blocks[-1].r_of_uv(np.array([[u1],[v1]]))[0]
		pslice = xh.junc.pslice(reg1, ublocks=[-1], vblocks=range(len(reg1.blocks)), r0=r1, u0=u1)
		## active slice of reg2
		aslice = xh.junc.aslice(reg2, ublocks=[-1], vblocks=range(len(reg2.blocks)), r0=pslice.r0, u0=u2, U0=pslice.U_of_r_at_v0, V0=pslice.V_of_r_at_u0, r_refs=[pslice.reg.metfunc.r_ref])
		## plot slice
		if False:
			print('REG1 PSLICE U(r,v0)=red, V(r,u0)=blue')
			print('REG2 ASLICE U(r,v0)=red, V(r,u0)=blue')
			plt.figure()
			plt.title('REG1 PSLICE U(r,v0)=red, V(r,u0)=blue\nREG2 ASLICE U(r,v0)=mag, V(r,u0)=cyan')
			plt.xlabel('r')
			plt.ylabel('U(r),V(r) on slice')
			plt.grid()
			plt.xlim(.5,1.5)
			## passive slice
			style1 = dict(marker='x', markersize=8, lw=1, ls='none')
			plt.plot(pslice.r, pslice.UV_v0[0], c='r', **style1)
			plt.plot(pslice.r, pslice.UV_u0[1], c='b', **style1)
			## active slice
			style2 = dict(marker='o', markersize=6, lw=1, ls='none')
			plt.plot(aslice.r, aslice.UV_v0[0], c='m', **style2)
			plt.plot(aslice.r, aslice.UV_u0[1], c='c', **style2)
			## junction radius
			plt.plot([pslice.r0, pslice.r0], [-1.,1.], 'k-' )
		## print slice params
		print("\n")
		print("                          %22r, %22r, %22r, %22r, %22r"%('Rh', 'r', 'r/Rh', 'u', 'v'))
		print("Pslice Inputs:")
		print("Rh1, r1, r1/Rh1, u1, v1 = %22r, %22r, %22r, %22r, %22r"%(Rh1, r1, r1/Rh1, u1, v1))
		print("From PSlice:")
		print("Rh1, r1, r1/Rh1, u1, v1 = %22r, %22r, %22r, %22r, %22r"%(Rh1, pslice.r0, pslice.r0/Rh1, pslice.u0, pslice.v0))
		print("Aslice Inputs:")
		print("Rh2, r2, r2/Rh2, u2, v2 = %22r, %22r, %22r, %22r, %22r"%(Rh2, pslice.r0, pslice.r0/Rh2, u2, '?'))
		print("From ASlice:")
		print("Rh2, r2, r2/Rh2, u2, v2 = %22r, %22r, %22r, %22r, %22r"%(Rh2, aslice.r0, aslice.r0/Rh2, aslice.u0, aslice.v0))
		print("\n")
		## reglist
		reglist = []
		## split regions
		reglist1 = xh.evap.split_reg_abcd(reg1,  abcd='bcd', u0=pslice.u0, v0=pslice.v0)
		reglist2 = xh.evap.split_reg_abcd(reg2,  abcd='a',   u0=aslice.u0, v0=aslice.v0)
		## update transformations for reg2
		for i in range(len(reglist2)):
			reglist2[i].U_of_udl = aslice.U_of_udl_at_v0
			reglist2[i].V_of_vdl = aslice.V_of_vdl_at_u0
		## add to reglist
		reglist += reglist1
		reglist += reglist2
		## print block uvbounds
		print("\n")
		print('%22r, %22r, %22r, %22r, %22r, %22r'%('R', 'b.rj', 'vmin', 'vmax', 'umin', 'umax'))
		for i in range(len(reglist)):
			reg = reglist[i]
			for b in reg.blocks:
					print('%22s, %22s, %22r, %22r, %22r, %22r'%(reg.metfunc.fparams['R'], b.rj, b.uvbounds['vmin'], b.uvbounds['vmax'], b.uvbounds['umin'], b.uvbounds['umax']))
		print("\n")
		## draw diagram
		if True:
			## add lines
			xh.evap.colorlines(reglist, sty=dict())
			xh.evap.boundarylines(reglist, sty=dict(), npoints=5001)
			xh.evap.s0_lines(reglist, sty=dict(lw=3))
			## draw
			plt.figure()
			plt.xlim(-3,3)
			plt.gca().set_aspect('equal')
			for reg in reglist:
				reg.rplot()
			## fill
			xh.evap.fill_by_R(reglist, cm=plt.cm.prism)
	## show plots
	plt.show()
	##
	print("\nEND TEST 10\n")






def test11(u0=-12., v0=0., du=.5, dv=.3, R=np.array([1.,.98,.96])):
	"""
	One intermediate region.
	Build from bottom.
	"""
	##
	print("\nTEST 11\n")
	## reglist
	reglist = []
	## funcs
	funclist = [xh.mf.schwarzschild(R=Rx) for Rx in R]
	## first region
	reglist += [xh.reg.EFreg(funclist.pop(0), boundary=False, rlines=False)]
	## second region
	reg1 = reglist.pop()
	reg2 = xh.reg.EFreg(funclist.pop(0), boundary=False, rlines=False)
	u0 = 1.*u0
	v0 = 1.*v0
	r0 = reg1.blocks[-1].r_of_uv(np.array([[u0],[v0]]))[0]
	pslice = xh.junc.pslice(reg1, ublocks=[-1], vblocks=range(len(reg1.blocks)), r0=1.*r0, u0=1.*u0)
	aslice = xh.junc.aslice(reg2, ublocks=[-1], vblocks=range(len(reg2.blocks)), r0=1.*pslice.r0, v0=1.*pslice.v0, U0=pslice.U_of_r_at_v0, V0=pslice.V_of_r_at_u0, r_refs=[pslice.reg.metfunc.r_ref])
	reg2.U_of_udl = aslice.U_of_udl_at_v0
	reg2.V_of_vdl = aslice.V_of_vdl_at_u0
	reglist += xh.evap.split_reg_abcd(reg1,  abcd='dbc', u0=1.*pslice.u0, v0=1.*pslice.v0)
	reglist += xh.evap.split_reg_abcd(reg2,  abcd='a'  , u0=1.*aslice.u0, v0=1.*aslice.v0)
	sliceprint(pslice, aslice, reg1.blocks[-1].rj[0], reg2.blocks[-1].rj[0])
	## third region
	reg1 = reglist.pop()
	reg2 = xh.reg.EFreg(funclist.pop(0), boundary=False, rlines=False)
	u0 = 1.*aslice.u0 + du
	v0 = 1.*aslice.v0 + dv
	r0 = reg1.blocks[-1].r_of_uv(np.array([[u0],[v0]]))[0]
	pslice = xh.junc.pslice(reg1, ublocks=[-1], vblocks=range(len(reg1.blocks)), r0=1.*r0, u0=1.*u0)
	aslice = xh.junc.aslice(reg2, ublocks=[-1], vblocks=range(len(reg2.blocks)), r0=1.*pslice.r0, v0=1.*pslice.v0, U0=pslice.U_of_r_at_v0, V0=pslice.V_of_r_at_u0, r_refs=[pslice.reg.metfunc.r_ref])
	reg2.U_of_udl = aslice.U_of_udl_at_v0
	reg2.V_of_vdl = aslice.V_of_vdl_at_u0
	reglist += xh.evap.split_reg_abcd(reg1,  abcd='dbc', u0=1.*pslice.u0, v0=1.*pslice.v0)
	reglist += xh.evap.split_reg_abcd(reg2,  abcd='a'  , u0=1.*aslice.u0, v0=1.*aslice.v0)
	sliceprint(pslice, aslice, reg1.blocks[-1].rj[0], reg2.blocks[-1].rj[0])
	## plot
	rgp(reglist)
	## show
	plt.show()
	##
	print("\nEND TEST 11\n")




def test12(u0=12., dr0=.1, du=.5, R=np.array([1.,.98,.96])):
	"""
	One intermediate region.
	Build from top.
	"""
	##
	print("\nTEST 12\n")
	## reglist
	reglist = []
	## funcs
	funclist = [xh.mf.schwarzschild(R=Rx) for Rx in R]
	## first region
	reglist += [xh.reg.EFreg(funclist.pop(-1), boundary=False, rlines=False)]
	## second region
	reg1 = reglist.pop()
	reg2 = xh.reg.EFreg(funclist.pop(-1), boundary=False, rlines=False)
	u0 = 1.*u0
	r0 = 1.*reg2.blocks[-1].rj[0] + 1.*dr0
	pslice = xh.junc.pslice(reg1, ublocks=[-1], vblocks=range(len(reg1.blocks)), r0=1.*r0, u0=1.*u0)
	aslice = xh.junc.aslice(reg2, ublocks=[-1], vblocks=range(len(reg2.blocks)), r0=1.*pslice.r0, u0=1.*pslice.u0, U0=pslice.U_of_r_at_v0, V0=pslice.V_of_r_at_u0, r_refs=[pslice.reg.metfunc.r_ref])
	reg2.U_of_udl = aslice.U_of_udl_at_v0
	reg2.V_of_vdl = aslice.V_of_vdl_at_u0
	reglist += xh.evap.split_reg_abcd(reg1,  abcd='a', u0=1.*pslice.u0, v0=1.*pslice.v0)
	reglist += xh.evap.split_reg_abcd(reg2,  abcd='dbc'  , u0=1.*aslice.u0, v0=1.*aslice.v0)
	sliceprint(pslice, aslice, reg1.blocks[-1].rj[0], reg2.blocks[-1].rj[0])
	## third region
	reg3 = xh.reg.EFreg(funclist.pop(-1), boundary=False, rlines=False)
	u0 = 1.*aslice.u0 - 1.*du
	r0 = 1.*reg3.blocks[-1].rj[0] + 1.*dr0
	pslice = xh.junc.pslice(reg2, ublocks=[-1], vblocks=range(len(reg1.blocks)), r0=1.*r0, u0=1.*u0)
	aslice = xh.junc.aslice(reg3, ublocks=[-1], vblocks=range(len(reg2.blocks)), r0=1.*pslice.r0, u0=1.*pslice.u0, U0=pslice.U_of_r_at_v0, V0=pslice.V_of_r_at_u0, r_refs=[pslice.reg.metfunc.r_ref])
	reg3.U_of_udl = aslice.U_of_udl_at_v0
	reg3.V_of_vdl = aslice.V_of_vdl_at_u0
	for reg in reglist[-3:-1]:
		for b in reg.blocks:
			b.uvbounds.update(dict(vmin=1.*pslice.v0))
	for reg in reglist[-2:]:
		for b in reg.blocks:
			b.uvbounds.update(dict(umin=1.*pslice.u0))
	reglist += xh.evap.split_reg_abcd(reg3,  abcd='dbc' , u0=1.*aslice.u0, v0=1.*aslice.v0)
	sliceprint(pslice, aslice, reg2.blocks[-1].rj[0], reg3.blocks[-1].rj[0])
	## plot
	rgp(reglist)
	## show
	plt.show()
	##
	print("\nEND TEST 12\n")





def test13():
	"""
	Test evap.evap.funclist_chain().
	"""
	##
	print("\nTEST 13\n")
	## make funclist
	R=np.array([0.,.5,1.,.7,0.])
	funclist  = xh.evap.formevap_funclist(R=R)
	## init
	t0 = np.zeros(len(funclist))
	x0 = np.zeros(len(funclist))
	## params
	t0 = -2. + t0
	x0 = 2. + x0
	## params
	## params
	ss = np.ones(len(funclist))
	du  = 1.*ss
	dv  = .9*ss
	##
	reglist, chainparams = xh.evap.funclist_chain(funclist, seed=-2, u0=t0-x0, v0=t0+x0, du=1.*du, dv=1.*dv, eta=1., matchmode='rv')
	pprint.pprint(chainparams)
	##
	reglist, chainparams = xh.evap.chain_masker(reglist, chainparams)
	pprint.pprint(chainparams)
	##
	rgp(reglist)
	plt.show()
	##
	print("\nEND TEST 13\n")




def test14():
	"""
	Test evap.evap.shellparams_list().
	"""
	##
	print("\nTEST 14\n")
	##
	## shellparams
	sp = xh.evap.shellparams_list(Rmin=.2, Rmax=1., dv=1., l=.1, A=10., functype=xh.mf.schwarzschild, fparams=dict())
	spt = xh.evap.sp_transpose(sp)
	## outs
	RR = spt['R']
	du = spt['du']
	dv = spt['dv']
	uu = spt['uu']
	vv = spt['vv']
	A = spt['A'][0]
	## smooth uu,vv
	ss = np.linspace(0,1.1*np.max(uu),5001)
	## plot
	if True:
		plt.figure()
		plt.title('R vs u')
		plt.xlabel('u')
		plt.ylabel('A*R^3')
		plt.xlim(0,10)
		plt.ylim(0,10)
		plt.plot(uu,A*RR**3,'k.')
		plt.plot(ss,ss,'b-')
		plt.grid(1)
	## plot
	if False:
		plt.figure()
		plt.title('R vs v')
		plt.xlabel('v')
		plt.ylabel('A*R^3')
		plt.ylim(0,A)
		plt.plot(vv,A*RR**3,'k.')
		plt.grid(1)
	## plot
	if False:
		plt.figure()
		plt.title('du/dv vs v')
		plt.xlabel('v')
		plt.ylabel('du/dv')
		plt.plot(vv,du/dv,'k.')
		plt.grid(1)
	## show
	plt.show()
	##
	##
	print("\nEND TEST 14\n")




def test15():
	"""
	Test evap.evap.shellparams_list().
	"""
	##
	print("\nTEST 15\n")
	##
	## shellparams
	sp = xh.evap.shellparams_list(Rmin=.2, Rmax=1., dv=1., l=.1, A=1., functype=xh.mf.schwarzschild, fparams=dict())
	spt = xh.evap.sp_transpose(sp)
	## outs
	funclist = spt['funclist']
	du = spt['du']
	dv = spt['dv']
	r0f = spt['R'] + spt['l']
	r0p = np.roll(r0f,1)
	v0 = 0.*np.linspace(-4, 4, len(funclist))
	##
	print("chain")
	reglist, chainparams = xh.evap.funclist_chain(funclist, seed=0, du=1.*du, dv=1.*dv, r0p=1.*r0p, r0f=1.*r0f, v0=1.*v0)
	pprint.pprint(chainparams)
	##
	print("mask")
	reglist, chainparams = xh.evap.chain_masker(reglist, chainparams)
	pprint.pprint(chainparams)
	##
	print("plot")
	rgp(reglist)
	plt.savefig('temp-figs/test15.png', dpi=200)
	##
	print("\nEND TEST 15\n")



def test16():
	"""
	Test evap.evap.formevap_input() with funclist_chain().
	"""
	##
	print("\nTEST 16\n")
	## params
	params = dict()
	## funcs
	params.update(dict(functype0=xh.mf.minkowski, fparams0=dict(), functype1=xh.mf.hayward, fparams1=dict(l=.01)))
	## evap
	params.update(dict(Rmin=1., Rmax=1., dv_evap=.5, l=.1, A=.2))
	## accrete
	params.update(dict(B=.5, Naccrete=1))
	## offset
	params.update(dict(voff=0., veta=1., uoff=0., ueta=0.))
	## seed
	seed = 0
	##
	print("inputs")
	funclist, cp = formevap_input(**params)
	##
	print("chain")
	reglist, chainparams = funclist_chain(funclist, seed=seed, **cp)
	##
	print("mask")
	reglist, chainparams = xh.evap.chain_masker(reglist, chainparams)
	pprint.pprint(chainparams)
	##
	print("plot")
	rgp(reglist[:])
	plt.xlim(0,3)
	plt.ylim(-1.5,1.5)
	##
	print("save fig")
	fname =  datetime.datetime.now().strftime("Test16_%Y-%m-%d_%H-%M-%S-%f")
	plt.title(fname)
	plt.savefig("temp-figs/%s.png"%(fname), dpi=400)
	##
	print("save txt")
	ff = open("temp-figs/%s.txt"%(fname), 'w')
	ff.write("%s\n"%(fname))
	ff.write('\n')
	ff.write('Input:\nparams=\n%s\nseed=\n%s\n'%(pprint.pformat(params),seed))
	ff.write('\n')
	ff.write('Output:\nchainparams=\n%s\n'%(pprint.pformat(chainparams)))
	##
	print("\nEND TEST 16\n")


main()
