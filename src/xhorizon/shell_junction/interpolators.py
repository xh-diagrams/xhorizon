
"""
Functions used to interpolate and extrapolate coord transformations from reference arrays.
"""

import numpy as np




def interp_with_smooth_extrap(x, x_ref, y_ref, mu=0.):
	"""
	Interpolate with smooth extrapolation from the reference endpoints,
	approaching a line of unit slope as mu*x goes to infinity.

	A wrapper for numpy.interp allowing:
		a. custom handling of x values outside the x_ref domain,
		b. scrubbing nan values and sorting ref arrays.
	
	The smooth continuation and decay to unit slope line uses the function 
		f1(x) = x + (mu/b)*(s0-1.)*(1.-np.exp(-x*b/mu)).

	We are specifically working with monotonic functions, which allows us to know whether to
	approach a slope of +1 or -1. The algorithm is basically tailored to monotonic increasing
	functions. Monotonic decreasing functions are handled by vertically flipping the reference
	y_ref array, then proceeding with the monotonic increasing routine. The final result y is
	then flipped back at the end.

	Output y(x) is continuous and once differentiable.

	Parameter mu multiplies the decay length for exponentally approaching unit slope line.
	When mu=0 the extrapolation is not differentiable, the line of extrapolation is just y=x.
	As mu goes to infinity, the line of extrapolation becomes linear, and this function 
	approaches interp_with_linear_extrap.
	"""
	## init
	x, x_ref, y_ref = 1.*x, 1.*x_ref, 1.*y_ref
	mu = np.abs(1.*float(mu))
	y = np.nan*x
	## scrub refs
	mask = np.logical_and.reduce( [np.isfinite(x_ref), np.isfinite(y_ref)] )
	x_ref, y_ref = x_ref[mask], y_ref[mask]
	## sort refs
	idx = np.argsort(x_ref)
	x_ref, y_ref = x_ref[idx], y_ref[idx]
	## flip if decreasing
	flip = False
	if y_ref[-1]<y_ref[0]:
		flip = True
		y_ref = -y_ref
	## get first and last ref points
	x0, y0 = x_ref[0],  y_ref[0]
	x1, y1 = x_ref[-1], y_ref[-1]
	## get slopes at first and last ref points
	dy0 = get_dy0(x_ref, y_ref)
	dy1 = get_dy1(x_ref, y_ref)
	## interpolate
	mask = np.logical_and.reduce( [x0<x, x<x1] )
	y[mask] = np.interp(x[mask], x_ref, y_ref)
	## extrapolate
	## left
	mask = x<=x0
	y[mask] = y0 + f1((x[mask]-x0), s0=dy0, mu=-mu, b=0.1+np.abs(np.log(dy0)))
	## right
	mask = x>=x1
	y[mask] = y1 + f1((x[mask]-x1), s0=dy1, mu= mu, b=0.1+np.abs(np.log(dy1)))
	## unflip if yref was flipped
	if flip==True:
		y = -y
	## return
	return y



def f1(x, s0=1., mu=0., b=1.):
	"""
	Monotonic function f(x) passing through origin with slope s0, 
	such that slope decays to 1 as mu*x goes to infinity.
	Inputs:
		x =  An array of x values to evaluate.
		s0 = Slope at the origin, used for differentiable matching.
		mu = Exponential decay length is mu/b. Typically use mu as overall user controlled multiplier.
		b  = Exponential decay length is mu/b. Typically use b to set an algorithmic length scale control.
			e.g. setting b = |ln(s0)| causes a stronger exponential decay when s0>>1 or s0<<1.
		Parameters mu and b are separate only for convenience and implementation.
	Returns:
		y = Array of values f(x).
	"""
	if np.any(np.logical_not(np.isfinite(np.array([b,s0,1./b])))):
		print('mu, s0, b = %s, %s, %s'%(mu, s0, b))

	## init
	x  = 1.*x
	mu = 1.*float(mu)
	b = 1.*float(b)
	s0 = 1.*float(s0)
	## case mu=0
	if mu==0.:
		y = 1.*x
	## go
	else:
		y = x + (mu/b)*(s0-1.)*(1.-np.exp(-x*b/mu))
	## return
	return y


def get_dy0(x,y):
	"""
	Find initial slope by iterating until finite.
	Avoids issues with repeated initial point.
	"""
	## loop until finite
	i, dy = 0, np.nan
	while (not np.isfinite(dy)) or dy==0.:
		i += 1
		dy  = (y[i] - y[0]) / (x[i] - x[0])
	## return
	return dy

def get_dy1(x,y):
	"""
	Find final slope by iterating until finite.
	Avoids issues with repeated final point.
	"""
	## loop until finite
	i, dy = 1, np.nan
	while (not np.isfinite(dy)) or dy==0.:
		i += 1
		dy  = (y[-1] - y[-i]) / (x[-1] - x[-i])
	## return
	return dy


















"""
Tests. Run if __name__=="__main__".
"""


def test2():
	"""
	Test interp_with_smooth_extrap.
	"""
	##
	print("\nTEST 2\n")
	## params
	A = 1e1
	mux = [.1,1.,10.]
	## ref values
	x_ref = np.linspace(-.8,1.2,1001)
	y_ref = A * x_ref**(3)
	## go
	## plot ref
	plt.plot(x_ref,y_ref,'rx', zorder=5)
	for i in range(len(mux)):
	## interp values
		mu = mux[i]
		x = np.linspace(-5,5,5001)
		y = interp_with_smooth_extrap(x,x_ref,y_ref,mu=mu)
		## plot interp
		c = 1.- float(i) / float(len(mux)) 
		plt.plot(x,y, 'b-', lw=2, alpha=c, zorder=100-i)
	##
	plt.gca().set_aspect('equal')
	sq = 5.
	ymult = 5.
	plt.xlim(-sq,sq)
	plt.ylim(-ymult*sq,ymult*sq)
	plt.grid()
	plt.show()
	##
	print("\nEND TEST 2\n")



def test3():
	"""
	Test interp_with_smooth_extrap.
	Decreasing functions now work.
	"""
	##
	print("\nTEST 3\n")
	## params
	A = -10.
	mu = 1.
	## ref values
	x_ref = np.linspace(-.95,1.05,1001)
	y_ref = A * x_ref**(3)
	## interp values
	x = np.linspace(-5,5,5001)
	y = interp_with_smooth_extrap(x, x_ref, y_ref, mu=mu)
	## plot ref
	plt.plot(x_ref,y_ref,'rx', zorder=5)
	## plot interp
	plt.plot(x,y, 'b-', lw=2, zorder=100)
	##
	plt.gca().set_aspect('equal')
	sq = 5.
	ymult = 2.
	plt.xlim(-sq,sq)
	#plt.ylim(-ymult*sq,ymult*sq)
	plt.grid()
	plt.show()
	##
	print("\nEND TEST 3\n")





if __name__=="__main__":

	import matplotlib.pyplot as plt

	print("\n\nTESTS\n")
	#test2()
	test3()
	print("\nEND TESTS\n")


