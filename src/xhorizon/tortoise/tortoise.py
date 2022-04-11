
"""
This module contains functions for numerically evaluating the global tortoise function F(r),
and its inverse F^{-1}_{j}(rstar).

After the main block of code, there are supporting utilities and tests.
"""


import numpy as np
import scipy.integrate as spint



class F_by_interp:

	"""
	Our numerical implementation of the global tortoise function is always defined by linear interpolation
	of densely spaced reference arrays. This method is used even when an analytically defined alternative is
	possible (note that analytical alternatives are rarely available for Finv). The interpolation method is 
	generally valid, and allows F and Finv to be implemented symmetrically.

	The F_by_interp class is necessary in order to avoid accidentally overwriting the objects containing the
	reference arrays. Since a new F_by_interp instance is created when the class is called, each function
	F (and Finv) will have its own copy of the reference arrays.

	The class methods make_F() and make_Finv() generate functions F(r) and Finv(j, rstar) respectively:
		The inputs r and rstar should be arrays.
		The input j is the interval index label. It may be either a scalar (float or int), or an array of equal length to rstar.

	In order to create the reference arrays, F_by_interp calls create_reference, which in turn calls other subroutines.
	The call hierarchy is:
		F_by_interp < create_reference < F_by_integration < F_refpoints.

	The heavy lifting is done in F_by_integration, which implements the integral definition of the tortoise function.
	The function F_refpoints supports this task by preliminarily determining the reference values F(r_j+eps) in each interval.
	"""

	def __init__(self, f, rj, eps, npoints=500):

		## store input data
		self.f = f
		self.rj = rj
		self.eps = eps
		self.npoints = npoints

		## create reference r and rstar arrays
		ref = create_reference(f, rj, eps, npoints)

		## store reference arrays
		self.r_ref = ref[0]
		self.rstar_ref = ref[1]

		## initialize empty piecewise subarrays (filled by make_Finv)
		self.r_ref_j = []
		self.rstar_ref_j = []

	def make_F(self):

		"""
		Using the reference arrays, returns a function object F(r).
		Values outside the valid range return nan.

		Input r to F(r) should be an array.
		"""

		## define function by linear interpolation of reference array
		def F(r):
			rstar = np.interp(r, self.r_ref, self.rstar_ref, left=np.nan, right=np.nan)
			return rstar

		## return a function object
		return F

	def make_Finv(self):

		"""
		Using the reference arrays, returns a function object Finv(j, rstar).
		Values outside the valid range return nan.

		Inputs rstar and j to Finv(j,r) should be equal length arrays.

		If input j is a scalar (float or int), it is converted to an array of correct length with uniform value.
		"""

		## determine intervals and initialize reference subarrays
		rj = self.rj
		intervals = range(len(rj)-1)
		r_ref = self.r_ref
		rstar_ref = self.rstar_ref
		r_ref_j = []
		rstar_ref_j = []

		## make reference subarrays for each interval
		for k in intervals:
			mask = np.logical_and( r_ref>rj[k] , r_ref<rj[k+1] )
			r_ref_j += [r_ref[mask]]
			rstar_ref_j += [rstar_ref[mask]]

		## sort reference subarrays for each interval such that rstar is increasing
		for k in intervals:
			idx = np.argsort( rstar_ref_j[k] )
			rstar_ref_j[k] = rstar_ref_j[k][idx]
			r_ref_j[k] = r_ref_j[k][idx]

		## save reference subarrays
		self.r_ref_j = r_ref_j
		self.rstar_ref_j = rstar_ref_j

		## define function by linearly interpolating reference subarrays
		def Finv(j, rstar):

			## convert j to an array of integers
			if np.isscalar(j):
				j = j + np.zeros_like(rstar)
			j = j.astype(int)

			## initialize output as nan
			r = np.nan * rstar

			## cycle through intervals, interpolating relevant points in each one
			for k in intervals:
				mask = (j==k)
				r[mask] = np.interp(rstar[mask], self.rstar_ref_j[k], self.r_ref_j[k], left=np.nan, right=np.nan)

			## return
			return r

		## return a function object
		return Finv



def create_reference(f, rj, eps, npoints=500):
	"""
	Create a well-spaced array of r values and calculate the corresponding values rstar = F(r).
	Return both arrays to use as interpolation reference for forward and inverse tortoise functions.
	"""
	r = rspace1(rj, eps, npoints)
	rstar = F_by_integration(f, rj, eps, r)
	## return
	return r, rstar



def F_by_integration(f, rj, eps, r):
	"""
	Given the metric function f(r) and the values of r_j and epsilon,
	returns an array rstar=F(r) from an array r.
	The values are determined by performing the defining definite integral.
	"""
	## initialize output to nan
	rstar = np.nan * r
	## get reference values
	Fj = F_refpoints(f, rj, eps)
	## determine appropriate interval for each point in r
	jj = which_interval(r, rj)
	## define integrand for calculation and loop over each point in r
	integrand = lambda r: f(r)**(-1)
	for m in range(len(r)):
		j = jj[m]
		a = rj[j] + eps
		b = r[m]
		rstar[m] = Fj[j] + integrate(integrand, a, b)
	## return
	return rstar



def F_refpoints(f, rj, eps):
	"""
	Given the metric function f(r) and the values of r_j and epsilon,
	returns an array of of the values F(r_j+eps) == F_j.
	Special Values:
		F(eps) = F(r_0 + eps) = 0
		F(inf + eps) = F(r_{N+1} + eps) = nan
		F(r_i + eps) = F(r_i - eps)
	Recursive Formula:
		F_0 = 0
		F_k = F_{k-1} + int_{r_{k-1}+eps}^{r_{k}-eps}
	"""
	## initialize output array to nan
	Fj = np.nan * rj
	## first value is necessarily zero
	Fj[0] = 0.
	## calculate remaining values by integrating
	integrand = lambda r: f(r)**(-1)
	for k in range(1, len(rj)-1 ):
		a = rj[k-1] + eps
		b = rj[k  ] - eps
		Fj[k] = Fj[k-1] + integrate(integrand, a, b)
	## return
	return Fj











"""
Utilities.
"""

def integrate(func,a,b):
	"""
	Wrapper for definite integration using scipy.integrate.quad().
	"""
	## calculate integral and error
	val, abserr = spint.quad(func,a,b, limit=100)
	## get relative error
	relerr = np.nan
	if val != 0:
		relerr = np.abs(abserr/val)
	## warning parameters
	maxrel, maxabs = 1e-4, 1e-6
	## execute warnings
	if (val!=0. and relerr > maxrel):
		print("\nIntegration Warning: Relative error exceeds %e"%(maxrel))
		print("%s from %s to %s"%(func,a,b))
		print("val %s, abserr %e, relerr %e"%(val,abserr,relerr))
	if (np.abs(val) < 1. and abserr > maxabs):
		print("\nIntegration Warning: Absolute error exceeds %e"%(maxabs))
		print("%s from %s to %s"%(func,a,b))
		print("val %s, abserr %e, relerr %e"%(val,abserr,relerr))
	## return
	return val


def which_interval(r,rj,asint=True):
	"""
	Given the radii r_j, determine which interval I_j a radius belongs to. Return nan if error.
	Inputs are the array rj and an arbitrary length array r.
	Output is an array j of equal length as r giving the corresponding interval
	for each radius in array r.
	"""
	j = np.nan * r
	for jj in range(len(rj)-1):
		mask = np.logical_and(r>rj[jj], r<rj[jj+1])
		j[mask] = jj
	if np.any(np.isnan(j)):
		print("\nWARNING: Function which_interval() encountered invalid values. Returning nan.\n")
	if asint==True:
		j = j.astype(int)
	return j


def rspace0(r0, r1, eps, npoints=500):
	"""
	Provide a dense array of r values for interpolating F and Finv.
	Must have extra values near exponential tails to smoothly interpolate there.
	"""
	## interval
	a = r0 + eps
	b = r1 - eps
	## linear spacing
	rr0 = np.linspace(a, b, 2*npoints)
	## exponential near tails
	B = np.log(eps**2)
	rr1 = a + 0.5 * (b-a) * np.exp( np.linspace(0, B, npoints) )
	rr2 = b - 0.5 * (b-a) * np.exp( np.linspace(0, B, npoints) )
	## concatenate and sort
	rr = np.sort(np.concatenate( [rr0,rr1,rr2] ))
	## return
	return rr


def rspace1(rj, eps, npoints=500):
	"""
	Like rspace0(), but gives values over all intervals I_j simultaneously.
	"""
	rx = []
	for j in range(len(rj)-1):
		rx += [rspace0(rj[j], rj[j+1], eps, npoints)]
	rvals = np.sort(np.concatenate(rx))
	## return
	return rvals


def first_trapped_j(f,rj):
	r0 = 0.5 * ( rj[0] + rj[1] )
	f0 = f(r0)
	out = np.nan
	if f0>0.:
		out = 1
	if f0<0.:
		out = 0
	return out









"""
Tests. 
Tests run if __name__ == "__main__".
"""

import matplotlib.pyplot as plt

def test1():
	"""
	Tests whether the function scipy.quad(h,a,b) obeys int_a^b = - int_b^a, even when endpoints
	approach singular points of the integral.
	This would be unlikely for numerical integration of an ODE, but for this area-based 
	quadrature it appears to be true.
	In fact, it appears to be true even when the integral method provides a warning.
	Since this is true, integrals of the form int_{rj+eps}^b are admissible, and using a
	less extreme integration point in the middle of the interval will not increase accuracy.
	"""
	######### input ########
	f = lambda r: 1. - 1. / r
	h = lambda r: f(r)**(-1)
	a = 1.5
	b = 1. + 1e-11
	########################
	print("\nTEST 1")
	forward      =  integrate(h,a,b)
	neg_backward = -integrate(h,b,a)
	abs_diff = forward - neg_backward
	sumtot = forward + neg_backward
	print(" forward  = % .18e"%forward)
	print("-backward = % .18e"%neg_backward)
	print(" abs diff = % .18e"%abs_diff)
	if sumtot != 0.:
		rel_diff = 2. * abs_diff / sumtot
		print(" rel diff = % .18e"%rel_diff)
	## end
	print("END TEST 1\n")




def test2():
	"""
	Demonstrate the functionality of F_refpoints.
	"""
	######### input ########
	f = lambda r: ( 1. - 1. / r )
	rj = np.array([0.,1.,25.])
	eps = 1e-11
	########################
	print("\nTEST 2")
	Fj = F_refpoints(f,rj,eps)
	print("Fj = ", Fj)
	## end
	print("END TEST 2\n")


def test3():
	"""
	Test of which_interval().
	"""
	######### input ########
	r = np.linspace(1e-12, 15, 20)
	rj = np.array([0.,1.,2.,8.,25.])
	########################
	print("\nTEST 3")
	j = which_interval(r, rj)
	print("j = ", j)
	## end
	print("END TEST 3\n")


def test4():
	"""
	Demonstrate the functionality of F_by_integration.
	"""
	######### input ########
	f = lambda r:  1. - 1. / r 
	rj = np.array([0.,1.,25.])
	eps = 1e-9
	r = np.linspace(1e-9,7,1000)
	########################
	print("\nTEST 4")
	## numerically calculate values
	rstar = F_by_integration(f,rj,eps,r)
	## know analytical solution
	rstar2 = r + np.log(np.abs(r-1.))
	## plot
	plt.plot(r, rstar,  'bx',  zorder=10)
	plt.plot(r, rstar2, 'r', lw=1.3, zorder=20)
	plt.grid()
	plt.show()
	## end
	print("END TEST 4\n")



def test5():
	"""
	Demonstrate the functionality of rspace0().
	"""
	######### input ########
	r0 = 0.
	r1 = 1.
	eps = 1e-9
	npoints = 500
	########################
	print("\nTEST 5")
	rvals = rspace0(r0,r1,eps,npoints)
	print("rvals = ", rvals)
	plt.plot(rvals,'kx')
	plt.show()
	## end
	print("END TEST 5\n")


def test6():
	"""
	Demonstrate the functionality of rspace1().
	"""
	######### input ########
	rj = np.array([0.,1.,10.,25.])
	eps = 1e-9
	npoints = 500
	########################
	print( "\nTEST 6")
	rvals = rspace1(rj,eps,npoints)
	print( "rvals = ", rvals)
	plt.plot(rvals,'kx')
	plt.show()
	## end
	print( "END TEST 6\n")


def test7():
	"""
	Demonstrate the functionality of create_reference().
	"""
	######### input ########
	f = lambda r: 1. -  1./r
	rj = np.array([0.,1.,25.])
	eps = 1e-9
	npoints=500
	########################
	print("\nTEST 7")
	## go
	r, rstar = create_reference(f, rj, eps, npoints)
	## plot
	plt.plot(r,rstar,'kx')
	plt.grid()
	plt.show()
	## end
	print("END TEST 7\n")


def test8():
	"""
	Demonstrate the functionality of F_by_interp.make_F().
	"""
	######### input ########
	f = lambda r: 1. -  1./r
	rj = np.array([0.,1.,25.])
	eps = 1e-9
	npoints=500
	r = np.concatenate([rspace1(rj,eps,npoints), np.linspace(50,-50,1000)])
	########################
	print("\nTEST 8")
	## build F_by_interp
	print("only this part is slow....")
	tort_by_interp = F_by_interp(f, rj, eps, npoints)
	print("see?")
	## make F
	F = tort_by_interp.make_F()
	## plot
	plt.plot(r,F(r),'kx')
	plt.grid()
	plt.show()
	## end
	print("END TEST 8\n")


def test9():
	"""
	Demonstrate the functionality of F_by_interp.make_Finv().
	"""
	######### input ########
	f = lambda r: 1. -  1./r
	rj = np.array([0.,1.,2.,25.])
	eps = 1e-9
	npoints=100
	rstar = np.linspace(-10,50,5000)
	jvals = [0,1,2]
	########################
	print("\nTEST 9")
	## build F_by_interp
	print("only this part is slow....")
	tort_by_interp = F_by_interp(f, rj, eps, npoints)
	print("see?")
	## make F
	Finv = tort_by_interp.make_Finv()
	## plot
	for j in jvals:
		plt.plot(Finv(j,rstar), rstar,'x')
	plt.grid()
	plt.show()
	## end
	print("END TEST 9\n")




## run tests if __name__="__main__"

if __name__=="__main__":
	print("RUN TESTS")
	test1()
	test2()
	test3()
	test4()
	test5()
	test6()
	test7()
	test8()
	test9()

