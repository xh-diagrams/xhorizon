
"""

"""

import numpy as np

from xhorizon.diagram_tools.curve_class import curve


"""
Specific masks.
Each of these subs returns
"""

def uv_range(blk, crvlist):
	"""
	"""
	## get uvbounds from block
	uv=dict(umin=-np.inf, umax=np.inf, vmin=-np.inf, vmax=np.inf)
	uv.update(blk.uvbounds)
	umin, umax = uv['umin'], uv['umax']
	vmin, vmax = uv['vmin'], uv['vmax']
	## go
	out = []
	for crv in crvlist:
		u,v = crv.uv
		mask = np.logical_and.reduce( [umin<=u, u<=umax, vmin<=v, v<=vmax] )
		out += mask_split_curve(crv, mask)
	return out

def rstar_minmax(blk, crvlist):
	"""
	Allow only valid values of rstar in the given block.
	"""
	## get min and max rstar values
	j = blk.j
	lims = rstar_limits(blk.master.metfunc)
	left, right = lims[j], lims[j+1]
	rstar_min, rstar_max = np.sort(np.array([left, right]))
	## proceed to mask curves
	out = []
	for crv in crvlist:
		u, v = crv.uv
		rstar = 0.5*(v-u) + blk.master.rparams['c']
		mask = np.logical_and.reduce( [rstar_min<=rstar, rstar<=rstar_max] )
		out += mask_split_curve(crv, mask)
	## return
	return out


def rstar_limits(func):
	"""
	Find limiting values of rstar at each rj given the metfunc.
	"""
	## get N
	N = len(func.rj) - 2
	## first N-1 blocks always follow same pattern
	rstar_limits = - np.inf * np.sign(func.kj)
	rstar_limits[0] = 0.
	## final block is trickier
	inf = 1e3 * func.rj[-1]
	finf = func.f(inf)
	## default to +- inf
	rstar_limits[N+1] = finf * np.inf
	## if |f(r)|>r as r to inf then |1./f(r)| < 1./r as r to inf, so F(r) converges
	## if F(r) converges use last computed value
	if np.abs(finf) > inf:
		rstar_limits[N+1] = (func.rstar_ref[np.argsort(func.r_ref)])[-1]
	## return
	return rstar_limits





"""
General tools.
"""

def mask_split_array(array, mask):
	"""
	Given an array and a mask of proper length, split the array into valid subarrays.
	Return list of valid subarrays. Discard invalid subarrays.
	Apply mask to itself to get length info.
	If array is empty, split into the correct number of empty arrays.
	"""
	## mask to int
	mask = mask.astype(int)
	## prepare output list of arrays
	out = []
	## get split indices for subarrays
	idx = np.nonzero(mask[1:] - mask[:-1])[0] + 1
	idx0 = np.concatenate([np.array([0]),idx])
	Nsubs = len(idx0)
	## if array is filled, split and return valid subs
	if len(array)>0:
		## split input array into subarrays
		subs = np.split(array, idx)
		## choose valid arrays and add to output
		for i in range(Nsubs):
			if mask[idx0[i]]==True:
				out += [subs[i]]
	## if array is empty, return correct number of empty subarrays
	if len(array)==0:
		for i in range(Nsubs):
			if mask[idx0[i]]==True:
				out += [np.array([])]
	## return
	return out

def mask_split_2darray(array, mask):
	"""
	Applies mask_split_array to each of two columns of a 2d array.
	Input is an array of the form xy = np.array([x,y]) and a boolean mask.
	Output is a list of valid 2d subarrays of xy.
	"""
	## prep output list
	out = []
	## get lists for each axis
	out0 = mask_split_array(array[0], mask)
	out1 = mask_split_array(array[1], mask)
	## recombine
	for i in range(len(out0)):
		out += [ np.array([out0[i],out1[i]]) ]
	## return
	return out


def mask_split_curve(crv, mask):
	"""
	"""
	## prep output list
	out = []
	## get unmasked array values
	r, tr, uv, uvdl, UV  = crv.r, crv.tr, crv.uv, crv.uvdl, crv.UV
	## create list of masked array values
	r  = mask_split_array(r, mask)
	tr = mask_split_2darray(tr, mask)
	uv = mask_split_2darray(uv, mask)
	uvdl = mask_split_2darray(uvdl, mask)
	UV = mask_split_2darray(UV, mask)
	## for each valid subcurve create new curve and fill values
	Nsubs = len(r)
	for i in range(Nsubs):
		newcrv = curve()
		newcrv.sty = crv.sty
		newcrv.r, newcrv.tr, newcrv.uv, newcrv.uvdl, newcrv.UV = r[i], tr[i], uv[i], uvdl[i], UV[i]
		out += [newcrv]
	## return
	return out






















"""
Tests.
Run if __name__='__main__'.
"""



def test1():
	"""
	Demonstrate the functionality of mask_split_array(array, mask).
	"""
	## define input arrays and mask
	x = np.linspace(-10,10,1001)
	y = np.sin(x)
	z = np.array([])
	mask = np.logical_or( np.abs(x)>5., y>.5)
	## get 
	xout = mask_split_array(x,mask)
	yout = mask_split_array(y,mask)
	zout = mask_split_array(z,mask)
	## print outputs
	print(xout)
	print(yout)
	print(zout)
	## plot
	plt.plot(x,y,'b-',lw=15, label='unmasked')
	plt.plot(x[mask],y[mask],'y-',lw=10, label='naive numpy mask')
	label='split mask with mask_split_array'
	for i in range(len(xout)):
		plt.plot(xout[i],yout[i],'r-', lw=5, label=label)
		label=None
	## format
	plt.title('test of mask_split_array(array, mask)')
	plt.legend(loc='lower right', fontsize=11)
	plt.grid()
	plt.show()



def test2():
	"""
	Demonstrate the functionality of mask_split_2darray(array, mask).
	Input is an array of the form xy = np.array([x,y]) and a boolean mask.
	Output is a list of valid 2d subarrays of xy.
	"""
	## define input arrays and mask
	x = np.linspace(-10,10,1001)
	y = np.sin(x)
	z = np.array([])
	mask = np.logical_or( np.abs(x)>5., y>.5)
	## define input 2d arrays
	xy = np.array([x,y])
	xz = np.array([x,z])
	## get 
	xyout = mask_split_2darray(xy,mask)
	xzout = mask_split_2darray(xz,mask)
	## print outputs
	print(xyout)
	print(xzout)
	## plot
	plt.plot(xy[0],xy[1],'b-',lw=15, label='unmasked')
	plt.plot(xy[0][mask],xy[1][mask],'y-',lw=10, label='naive numpy mask')
	label='split mask with mask_split_2darray'
	for i in range(len(xyout)):
		plt.plot(xyout[i][0],xyout[i][1],'r-', lw=5, label=label)
		label=None
	## format
	plt.title('test of mask_split_2darray(array, mask)')
	plt.legend(loc='lower right', fontsize=11)
	plt.grid()
	plt.show()



def test3():
	"""
	Demonstrate the functionality of mask_split_curve(crv, mask).
	"""
	## create input curve
	crv = curve()
	## input curve coordinates
	s = np.linspace(-5,3,1001)
	crv.r = s**2
	crv.tr = np.array([s,s**2])
	crv.uv = np.array([s,s])
	crv.uvdl = np.array([s**3,s**2])
	## define mask
	mask = np.logical_or(np.sin(5.*crv.uv[0]**2) > .5, np.abs(crv.uv[1]) > 3.5)
	## run
	out = mask_split_curve(crv, mask)
	# plot
	for s in ['tr','uv','uvdl','UV']:
		x = crv.__dict__[s]
		## plot unmasked
		plt.plot(x[0], x[1], 'b-', lw=15, label='unmasked')
		## plot naive numpy mask application
		if len(x[0])>0:
			plt.plot(x[0][mask], x[1][mask], 'y-', lw=10, label='naive numpy mask')
		## plot valid subcurves
		label = 'split mask with mask_split_curve'
		for i in range(len(out)):
			x = out[i].__dict__[s]
			plt.plot(x[0], x[1], 'r-', lw=5, label=label)
			label = None
		plt.title(s)
		plt.legend(loc='lower left', fontsize=10)
		plt.show()






## run tests if __name__='__main__'.
if __name__=='__main__':
	## import
	import matplotlib.pyplot as plt
	## test
	# test1()
	# test2()
	# test3()

