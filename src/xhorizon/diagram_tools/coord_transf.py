
"""
This module specifies standard transformations between
the double null block coordinates $(u,v)$ and region
coordinates $(\tilde{u},\tilde{v})$ (the standard
Penrose coordinates on an SSS region), along with the
associated helper functions.
"""

"""
The coordinate hierarchy is as follows:

"Schwarzschild Coordinates"
tr = $ (t,r) $
Schwarzschild coordinates covering a single block with metric
$ ds^2 = -f(r) dt^2 + f(r)^{-1} dr^2 + r^2 d\Omega^2 $.

"Double Null Block Coordinates"
uv = $ (u,v) $
Double-null block coordinates covering a single block with metric
$ ds^2 = - du dv + r^2 d\Omega^2 $.

"Region Coordinates"
uvdl = $ (\tilde{u},\tilde{v}) $
Standard Penrose coordinates for an SSS spacetime.

"Diagram Coordinates"
UV = $ (U,V) $
Coordinate system to be directly plotted in Penrose diagram.
Obtained from region coords by transformations U(udl) and V(vdl).
"""


import numpy as np


## block to block transformations

def uv_of_tr(tr, F, c):
	t, r = tr
	u = t - F(r) + c
	v = t + F(r) - c
	uv = np.array([u,v])
	return uv

def tr_of_uv(uv, Finv, c):
	u, v = uv
	t = 0.5*(v+u)
	r = Finv( 0.5*(v-u) + c)
	tr = np.array([t,r])
	return tr

def r_of_uv(uv, Finv, c):
	r = (tr_of_uv(uv, Finv, c))[1]
	return r



## block to region transformations

def uvdl_of_uv(uv, cdlu, cdlv, epsu, epsv, kpm, s0):
	h = lambda s: hks(s, kpm, s0)
	u, v = uv
	udl = np.pi**(-1) * np.arctan(  epsu * h( u/2.) ) + cdlu
	vdl = np.pi**(-1) * np.arctan( -epsv * h(-v/2.) ) + cdlv
	uvdl = np.array([udl,vdl])
	return uvdl

def uv_of_uvdl(uvdl, cdlu, cdlv, epsu, epsv, kpm, s0):
	hinv = lambda s: hksinv(s, kpm, s0)
	udl, vdl = uvdl
	u =  2. * hinv(  epsu**(-1) * np.tan( np.pi * (udl - cdlu) ) )
	v = -2. * hinv( -epsv**(-1) * np.tan( np.pi * (vdl - cdlv) ) )
	uv = np.array([u,v])
	return uv











## supplementary functions

def hks(s, kpm, s0):
	## get k_+ and k_- values
	km, kp = np.min(kpm), np.max(kpm)
	## initialize array w nan
	h = np.nan * s
	## implement piecewise function
	mask = s < -s0
	h[mask] = -s0 + H(s[mask]+s0, km)
	mask = np.abs(s) <= s0
	h[mask] = s[mask]
	mask = s > s0
	h[mask] =  s0 + H(s[mask]-s0, kp)
	## return
	return h

def H(s,k):
	## initialize array w nan
	H = np.nan * s
	## implement piecewise function
	if k!=0.:
		H = k**(-1) * ( np.exp(k*s) - 1. )
	if k==0.:
		H = s
	## return
	return H


def Hinv(s,k):
	## initialize array w nan
	Hinv = np.nan * s
	## implement piecewise function
	if k!=0.:
		Hinv = k**(-1) * np.log( k*s + 1 )
	if k==0.:
		Hinv = s
	## return
	return Hinv


def hksinv(s, kpm, s0):
	## get k_+ and k_- values
	km, kp = np.min(kpm), np.max(kpm)
	## initialize array w nan
	hinv = np.nan * s
	## implement piecewise function
	mask = s < -s0
	hinv[mask] = -s0 + Hinv(s[mask]+s0, km)
	mask = np.abs(s) <= s0
	hinv[mask] = s[mask]
	mask = s > s0
	hinv[mask] =  s0 + Hinv(s[mask]-s0, kp)
	## return
	return hinv













"""
Tests. Tests run if __name__=="__main__".
"""

import matplotlib.pyplot as plt


def test1():
	"""
	Test that uv_of_tr and tr_of_uv are stable numerical inverses.
	"""
	print( "\nTEST 1")
	## input ##
	npts = 100
	t = 10. * ( 2.*np.random.rand(npts) - 1. )
	r = 10. * np.random.rand(npts)
	c = 10. * ( 2.*np.random.random() - 1. )
	F = lambda rr: np.arctan(rr)
	Finv = lambda rrstar: np.tan(rrstar)
	##
	## test block to block transformations
	tr = np.array([t,r])
	uv = uv_of_tr(tr, F, c)
	uvnew = uv_of_tr(tr, F, c)
	trnew = tr_of_uv(uv, Finv, c)
	for i in range(4):
		uvnew = uv_of_tr(trnew, F, c)
		trnew = tr_of_uv(uvnew, Finv, c)
		rnew = r_of_uv(uvnew, Finv, c)
	## results of test
	if False:
		print( "uv = ", uv)
		print( "tr = ", tr)
		print( "r = ", r)
	if True:
		print( "max diff uv = ", np.max(np.abs(uvnew-uv)))
		print( "max diff tr = ", np.max(np.abs(trnew-tr)))
		print( "max diff r = ", np.max(np.abs(rnew-r)))
	## end
	print( "END TEST 1\n")


def test2():
	"""
	Test that H and Hinv are stable numerical inverses.
	"""
	print( "\nTEST 2")
	## input ##
	s = np.linspace(-20.,20.,3001)
	k = 1.
	##
	## test H
	x = s
	y = H(s, k)
	xnew = Hinv(y, k)
	ynew = H(xnew, k)
	for i in range(4):
		xnew = Hinv(ynew, k)
		ynew = H(xnew, k)
	if False:
		print( "x = ", x)
		print( "y = ", y)
	if True:
		print( "max diff x = ", np.max(np.abs(xnew-x)))
		print( "max diff y = ", np.max(np.abs(ynew-y)))
	## end
	print( "END TEST 2\n")


def test3():
	"""
	Test that h and hinv are stable numerical inverses.
	"""
	print( "\nTEST 3")
	## input ##
	s = np.linspace(-20.,20.,2001)
	kpm = np.array([-1.,2.])
	s0 = 5.
	## setup h
	h = lambda s: hks(s, kpm, s0)
	hinv = lambda s: hksinv(s, kpm, s0)
	## test h
	x = s
	y = h(s)
	xnew = hinv(y)
	ynew = h(xnew)
	for i in range(4):
		xnew = hinv(ynew)
		ynew = h(xnew)
	if False:
		print( "x = ", x)
		print( "y = ", y)
	if True:
		print( "max diff x = ", np.max(np.abs(xnew-x)))
		print( "max diff y = ", np.max(np.abs(ynew-y)))
	## end
	print( "END TEST 3\n")


def test4():
	"""
	Plot h and hinv.
	"""
	print( "\nTEST 4")
	## input ##
	s = np.linspace(-20.,20.,2001)
	kpm = np.array([-0.6,0.2])
	kpm = np.sort(kpm)
	s0 = 5.
	## setup h
	h = lambda s: hks(s, kpm, s0)
	hinv = lambda s: hksinv(s, kpm, s0)
	## test plot h
	plt.figure()
	plt.plot(s,s,'0.8')
	plt.plot(s,h(s),'k-', label='hks(s,kpm,s0)')
	plt.plot(s,hinv(s),'k--', label='hksinv(s,kpm,s0)')
	sq = 2.*s0
	plt.xlim(-sq,sq)
	plt.ylim(-sq,sq)
	plt.xticks([-s0,0,s0])
	plt.yticks([-s0,0,s0])
	plt.grid()
	plt.gca().set_aspect('equal')
	plt.xlabel('s')
	plt.title('TEST 4 (k-=%s, k+=%s, s0=%s)'%(kpm[0],kpm[1],s0))
	plt.legend(loc='lower right', fontsize=10)
	plt.show()
	## end
	print( "END TEST 4\n")


def test5():
	"""
	Plot H and Hinv.
	"""
	print( "\nTEST 5")
	## input ##
	s = np.linspace(-20.,20.,10001)
	k = -10.
	## test plot H
	plt.figure()
	plt.plot(s,s,'0.8')
	plt.plot(s, H(s,k),'k-', label='H(s,k)')
	## test plot Hinv
	s = np.linspace(-.999999/k, 10.*k, 10001)
	plt.plot(s, Hinv(s,k), 'k--', label='Hinv(s,k)')
	## format
	sq = 2.
	plt.xlim(-sq,sq)
	plt.ylim(-sq,sq)
	plt.grid()
	plt.gca().set_aspect('equal')
	plt.xlabel('s')
	plt.title('TEST 5 (k=%r)'%(k))
	plt.legend(loc='lower right', fontsize=10)
	plt.show()
	## end
	print( "END TEST 5\n")

## run test routines if this module is main

if __name__=='__main__':
	test1()
	test2()
	test3()
	test4()
	test5()

