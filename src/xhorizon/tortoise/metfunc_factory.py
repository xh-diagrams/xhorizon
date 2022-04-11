
"""
This module generates
	metfunc
objects from templates.
"""

import numpy as np


from xhorizon.tortoise.metfunc_class import metfunc
from xhorizon.tortoise import tortoise
from xhorizon.tortoise import math_util


"""
Templates.
Each template calls build_metfunc(params) in order to build the object from specified parameters.
"""

def minkowski():
	############################# INPUT #############################
	## parameters
	fparams = dict()
	Fparams = dict(eps=1e-9, inf=25., npoints_interp=500)
	## metric function
	f = lambda r: 1. + 0.*r
	## zeroes and slopes of f(r)
	ri = np.array([])
	ki = np.array([])
	## info
	info = {
		'Type' : "Minkowski",
		'Metric Function' : r"$f(r) = 1$",
		'Parameters' : '',
		}
	#################################################################
	## store params in dict
	params = dict(fparams=fparams, Fparams=Fparams, f=f, ri=ri, ki=ki, info=info)
	## build
	func = build_metfunc(params)
	## return
	return func





def schwarzschild(R=1.):
	############################# INPUT #############################
	## parameters
	R = float(R)
	fparams = dict(R=R)
	Fparams = dict(eps=1e-6*R, inf=25., npoints_interp=500)
	## metric function
	f = lambda r: 1. - R / r
	## zeroes and slopes of f(r)
	ri = np.array([R])
	ki = np.array([1./R])
	## info
	info = {
		'Type' : "Schwarzschild",
		'Metric Function' : r"$f(r) = 1 - R/r$",
		'Parameters' : ', '.join([r"$R=%r$"%R])
		}
	#################################################################
	## store params in dict
	params = dict(fparams=fparams, Fparams=Fparams, f=f, ri=ri, ki=ki, info=info)
	## build
	func = build_metfunc(params)
	## return
	return func


def reissner_nordstrom(M=0.5, Q=0.3):
	############################# INPUT #############################
	## parameters
	M, Q = float(M), float(Q)
	fparams = dict(M=M,Q=Q)
	Fparams = dict(eps=1e-9, inf=25., npoints_interp=500)
	## metric function
	f = lambda r: 1. - 2.*M/r + Q**2/r**2
	## zeroes and slopes of f(r)
	ri = math_util.pos_real_roots(np.array([1., -2.*M, Q**2]))
	ki = 2. * ri**(-3) * ( M*ri - Q**2)
	## info
	info = {
		'Type' : "Reissner Nordstrom",
		'Metric Function' : r"$f(r) = 1 - 2M/r + Q^2/r^2$",
		'Parameters' : ', '.join([r"$M=%r$"%M,r"$Q=%r$"%Q])
		}
	#################################################################
	## store params in dict
	params = dict(fparams=fparams, Fparams=Fparams, f=f, ri=ri, ki=ki, info=info)
	## build
	func = build_metfunc(params)
	## return
	return func


def hayward(R=1.,l=0.1):
	############################# INPUT #############################
	## parameters
	R, l = float(R), float(l)
	fparams = dict(R=R,l=l)
	Fparams = dict(eps=1e-9, inf=25., npoints_interp=500)
	## metric function
	f = lambda r: 1. - R*r**2 / ( R*l**2 + r**3 )
	## zeroes and slopes of f(r)
	ri = math_util.pos_real_roots(np.array([1.,-R,0.,R*l**2]))
	ki = ( R * ri**4 - 2. * R**2 * l**2 * ri ) / ( R * l**2 + ri**3 )**2
	## info
	info = {
		'Type' : "Hayward",
		'Metric Function' : r"$f(r) = 1 - R r^2 / (R l^2 + r^3)$",
		'Parameters' : ', '.join([r"$R=%r$"%R,r"$l=%r$"%l])
		}
	#################################################################
	## store params in dict
	params = dict(fparams=fparams, Fparams=Fparams, f=f, ri=ri, ki=ki, info=info)
	## build
	func = build_metfunc(params)
	## return
	return func


def dS(L=10.):
	############################# INPUT #############################
	## parameters
	L = float(L)
	fparams = dict(L=L)
	Fparams = dict(eps=1e-9, inf=20.*L, npoints_interp=500)
	## metric function
	f = lambda r: 1. - (r/L)**2
	## zeroes and slopes of f(r)
	ri = np.array([L])
	ki = np.array([-2./L])
	## info
	info = {
		'Type' : "de Sitter",
		'Metric Function' : r"$f(r) = 1 - (r/L)^2$",
		'Parameters' : r"$L=%r$"%L
		}
	#################################################################
	## store params in dict
	params = dict(fparams=fparams, Fparams=Fparams, f=f, ri=ri, ki=ki, info=info)
	## build
	func = build_metfunc(params)
	## return
	return func


def AdS(L=10.):
	############################# INPUT #############################
	## parameters
	L = float(L)
	fparams = dict(L=L)
	Fparams = dict(eps=1e-9, inf=20.*L, npoints_interp=500)
	## metric function
	f = lambda r: 1. + (r/L)**2
	## zeroes and slopes of f(r)
	ri = np.array([])
	ki = np.array([])
	## info
	info = {
		'Type' : "Anti de Sitter",
		'Metric Function' : r"$f(r) = 1 + (r/L)^2$",
		'Parameters' : r"$L=%r$"%L
		}
	#################################################################
	## store params in dict
	params = dict(fparams=fparams, Fparams=Fparams, f=f, ri=ri, ki=ki, info=info)
	## build
	func = build_metfunc(params)
	## return
	return func



def Hay_dS(l=.1, R=1., L=10.):
	############################# INPUT #############################
	## parameters
	l, R, L = float(l), float(R), float(L)
	fparams = dict(l=l,R=R,L=L)
	Fparams = dict(eps=1e-9, inf=20.*L, npoints_interp=500)
	## metric function
	f = lambda r: 1. - R*r**2/(R*l**2+r**3) - r**2/L**2
	## zeroes and slopes of f(r)
	ri = math_util.pos_real_roots(np.array([-1.,0.,L**2,-R*l**2-R*L**2,0.,R*(L**2)*(l**2)]))
	ki = ( R * ri**4 - 2. * R**2 * l**2 * ri ) / ( R * l**2 + ri**3 )**2 - 2.*ri / L**2
	## info
	info = {
		'Type' : "Hayward - de Sitter",
		'Metric Function' : r"$f(r) = 1 - Rr^2/(Rl^2+r^3) - (r/L)^2$",
		'Parameters' : ', '.join([r"$l=%r$"%l,r"$R=%r$"%R,r"$L=%r$"%L])
		}
	#################################################################
	## store params in dict
	params = dict(fparams=fparams, Fparams=Fparams, f=f, ri=ri, ki=ki, info=info)
	## build
	func = build_metfunc(params)
	## return
	return func


def Hay_AdS(l=.1, R=1., L=10.):
	############################# INPUT #############################
	## parameters
	l, R, L = float(l), float(R), float(L)
	fparams = dict(l=l,R=R,L=L)
	Fparams = dict(eps=1e-9, inf=20.*L, npoints_interp=500)
	## metric function
	f = lambda r: 1. - R*r**2/(R*l**2+r**3) + r**2/L**2
	## zeroes and slopes of f(r)
	ri = math_util.pos_real_roots(np.array([ 1.,0.,L**2, R*l**2-R*L**2,0.,R*(L**2)*(l**2)]))
	ki = ( R * ri**4 - 2. * R**2 * l**2 * ri ) / ( R * l**2 + ri**3 )**2 + 2.*ri / L**2
	## info
	info = {
		'Type' : "Hayward - Anti de Sitter",
		'Metric Function' : r"$f(r) = 1 - Rr^2/(Rl^2+r^3) + (r/L)^2$",
		'Parameters' : ', '.join([r"$l=%r$"%l,r"$R=%r$"%R,r"$L=%r$"%L])
		}
	#################################################################
	## store params in dict
	params = dict(fparams=fparams, Fparams=Fparams, f=f, ri=ri, ki=ki, info=info)
	## build
	func = build_metfunc(params)
	## return
	return func


def Hay_dS_cutoff(l=.1, R=1., L=10.):
	############################# INPUT #############################
	## parameters
	l, R, L = float(l), float(R), float(L)
	fparams = dict(l=l,R=R,L=L)
	Fparams = dict(eps=1e-9, inf=20.*L, npoints_interp=500)
	## L0
	ri = math_util.pos_real_roots(np.array([-1.,0.,L**2,-R*l**2-R*L**2,0.,R*(L**2)*(l**2)]))
	L0 = ri[-1] - 0.1*(ri[-1] - ri[-2])
	## metric function
	def f(r):
		if type(r) in [float, np.float64]:
			if r<=L0:
				return 1. - r**2/L**2
			if r>L0:
				return 1. - R*r**2/(R*l**2+r**3) - r**2/L**2
		if type(r)==np.ndarray:
			ff = np.nan*r
			mask = r> L0
			ff[mask] = 1. - R*r[mask]**2/(R*l**2+r[mask]**3) - r[mask]**2/L**2
			mask = r<=L0
			ff[mask] = 1. - r[mask]**2/L**2
			return 1.*ff
	## zeroes and slopes of f(r)
	ri = math_util.pos_real_roots(np.array([-1.,0.,L**2,-R*l**2-R*L**2,0.,R*(L**2)*(l**2)]))
	ki = ( R * ri**4 - 2. * R**2 * l**2 * ri ) / ( R * l**2 + ri**3 )**2 - 2.*ri / L**2
	## should be just one though because of cutoff
	ri = ri[-1:]
	ki = ki[-1:]
	## info
	info = {
		'Type' : "Hayward - de Sitter with Cutoff",
		'Metric Function' : r"De Sitter inside cutoff length, Hay-dS outside cutoff length.",
		'Parameters' : ', '.join([r"$l=%r$"%l,r"$R=%r$"%R,r"$L=%r$"%L])
		}
	#################################################################
	## store params in dict
	params = dict(fparams=fparams, Fparams=Fparams, f=f, ri=ri, ki=ki, info=info)
	## build
	func = build_metfunc(params)
	## return
	return func



def negativeAdS(L=1.):
	############################# INPUT #############################
	## parameters
	L = float(L)
	fparams = dict(L=L)
	Fparams = dict(eps=1e-9, inf=50., npoints_interp=500)
	## metric function
	f = lambda r: -( 1. + (r/L)**2 )
	## zeroes and slopes of f(r)
	ri = np.array([])
	ki = np.array([])
	## info
	info = {
		'Type' : "Negative Anti de Sitter",
		'Metric Function' : r"$f(r) = -(1 + (r/L)^2)$",
		'Parameters' : r"$L=%r$"%L
		}
	#################################################################
	## store params in dict
	params = dict(fparams=fparams, Fparams=Fparams, f=f, ri=ri, ki=ki, info=info)
	## build
	func = build_metfunc(params)
	## return
	return func



def schwarzschild_dS(R=1.,L=10.):
	############################# INPUT #############################
	## parameters
	R, L = float(R), float(L)
	fparams = dict(R=R,L=L)
	Fparams = dict(eps=1e-9, inf=25., npoints_interp=500)
	## metric function
	f = lambda r: 1. - R/r - r**2/L**2
	## zeroes and slopes of f(r)
	ri = math_util.pos_real_roots(np.array([1.,0.,-(L**2),R*L**2]))
	ki = R/ri**2 - 2.*ri/L**2
	## info
	info = {
		'Type' : "Schwarzschild - de Sitter",
		'Metric Function' : r"$f(r) = 1 - R/r - (r/L)^2$",
		'Parameters' : ', '.join([r"$R=%r$"%R,r"$L=%r$"%L])
		}
	#################################################################
	## store params in dict
	params = dict(fparams=fparams, Fparams=Fparams, f=f, ri=ri, ki=ki, info=info)
	## build
	func = build_metfunc(params)
	## return
	return func


"""
Generic build from params.
"""

def build_metfunc(params):
	## get build parameters from input dictionary
	fparams = params['fparams']
	Fparams = params['Fparams']
	f       = params['f']
	ri      = params['ri']
	ki      = params['ki']
	info    = params['info']
	## create object
	func = metfunc()
	## fill fields using specified values
	func.f = f
	func.fparams = fparams
	func.Fparams = Fparams
	func.rj = np.concatenate( [ np.array([0.]), ri, np.array([Fparams['inf']]) ] )
	func.kj = np.concatenate( [ np.array([0.]), ki, np.array([0.]) ] )
	func.info = info
	## generate tortoise function and inverse
	print("\nThis part is slow, reference F(r) arrays being calculated by integration...")
	tort_by_interp = tortoise.F_by_interp(func.f, func.rj, func.Fparams['eps'], npoints=func.Fparams['npoints_interp'])
	print("...done.\n")
	func.F = tort_by_interp.make_F()
	func.Finv = tort_by_interp.make_Finv()
	## define interval getter, first trapped interval
	func.winterval = lambda r,asint=True: tortoise.which_interval(r, func.rj, asint)
	func.first_trapped_j = tortoise.first_trapped_j(func.f,func.rj)
	## define interpolation reference values
	func.r_ref = tort_by_interp.r_ref
	func.rstar_ref = tort_by_interp.rstar_ref
	## return
	return func










"""
Tests.
Tests run if __name__='__main__'.
"""


def test1():
	"""
	"""
	######### input ########
	func  = hayward(R=.05, l=.01)
	########################
	print("\nTEST 9")
	print(func)
	print(func.rj)
	print(func.info['Type'])
	##
	import matplotlib.pyplot as plt
	plt.plot(func.r_ref, func.rstar_ref, 'kx')
	plt.xlim(0,25)
	plt.grid()
	plt.show()
	## end
	print("END TEST 9\n")



## run tests if __name__="__main__"

if __name__=="__main__":
	print("RUN TESTS")
	test1()

