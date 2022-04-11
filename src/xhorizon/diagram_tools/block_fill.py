
"""
This module contains subroutines for filling parts of the diagram with a background color.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from xhorizon.diagram_tools import curve_class as cc
from xhorizon.diagram_tools.curvemakers import rstarlines_special_2, block_boundary, block_boundary_2, rlines, rstarlines




def fill_between_curves_uv(blk, crvlist, tr=False, reverse=True, sty={}, eps=1e-24):
	"""
	Fill in a patch using a list of curves to provide the boundary of polygon.
	The curves are read from their uv coords.
	Every other input curve is reversed by default.
		block = block in which to fill curves
		crvlist = list of two curves with uv coords
		tr = set to true if tr coords are primary
		reverse = reverse point order of every other curve
	Curves are clipped to given rectangular size using snap_crvlist_to_bounds.
	"""
	## default style
	style = dict(fc='c', ec='none', lw=0.25, zorder=100)
	## get uvbounds from block
	uv=dict(umin=-np.inf, umax=np.inf, vmin=-np.inf, vmax=np.inf)
	uv.update(blk.uvbounds)
	umin, umax = uv['umin'], uv['umax']
	vmin, vmax = uv['vmin'], uv['vmax']
	## make sure limits are float
	umin, umax, vmin, vmax = float(umin), float(umax), float(vmin), float(vmax)
	## reverse odd index input curves
	if reverse==True:
		for i in range(len(crvlist)):
			if (i % 2) == 1:
				crvlist[i] = reverse_curve(crvlist[i])
	## start from tr values if requested
	if tr==True:
		crvlist = blk.update_curves_from_tr(crvlist)
	## snap uv values to within desired rectangle
	crvlist = snap_crvlist_to_bounds(crvlist, umin=umin, umax=umax, vmin=vmin, vmax=vmax, eps=eps)
	## update crvlist from new uv coords and apply masks
	crvlist = blk.apply_masks( blk.update_curves_from_uv( crvlist ) )
	## concatenate UV coords to create boundary
	U = np.concatenate([ crv.UV[0] for crv in crvlist ])
	V = np.concatenate([ crv.UV[1] for crv in crvlist ])
	## create xy array
	x, y = V-U, V+U
	xy = np.dstack([x,y])[0]
	## plot polygon
	style.update(sty)
	p = mpl.patches.Polygon(xy, **style)
	plt.gca().add_patch(p)



def fill_block(blk, sty={}):
	"""
	Fill in an entire block.
	"""
	## get uvbounds from block
	uv=dict(umin=-np.inf, umax=np.inf, vmin=-np.inf, vmax=np.inf)
	uv.update(blk.uvbounds)
	umin, umax = uv['umin'], uv['umax']
	vmin, vmax = uv['vmin'], uv['vmax']
	## default style
	style = dict(fc='k', ec='none', zorder=10)
	style.update(sty)
	## generate block boundaries
	crvlist = block_boundary_2(blk)
	## fill between boundaries
	fill_between_curves_uv(blk, crvlist, sty=sty)



def fill_between_r(blk, rvals=np.array([0.,0.]), sty={}, inf=100., npoints=1000):
	"""
	"""
	## default to invalid params
	valid = False
	## check that there is something to fill
	eps = blk.master.metfunc.Fparams['eps']
	rsnap = snap_to_bounds(rvals, ymin=blk.rj[0]+1.1*eps, ymax=blk.rj[1]-1.1*eps)
	if np.abs(rsnap[1]-rsnap[0]) > 5.*eps:
		valid = True
	## get block boundary curves
	boundary_curves = block_boundary_2(blk, inf=inf, npoints=npoints)
	## get params and initialize boundary rstar values
	boundary_rstar = np.nan * rvals
	c = float(blk.master.rparams['c'])
	## determine block boundary rstar values in order of increasing radius
	for i in range(2):
		u,v = boundary_curves[i].uv
		boundary_rstar[i] = 0.5*(v[0]-u[0])+c
	## for each rval either use F(r) or boundary value as rstar
	rvals = np.sort(rvals)
	rstar_vals = np.nan * rvals
	for i in range(len(rvals)):
		rval = rvals[i]
		## check if radius value is finite
		if not np.isfinite(rval):
			valid = False
		## if radius too small use left boundary
		elif rval <= blk.rj[0]+1.1*eps:
			rstar_vals[i] = boundary_rstar[0]
		## if radius too big use right boundary
		elif rval >= blk.rj[1]-1.1*eps:
			rstar_vals[i] = boundary_rstar[1]
		## if radius is valid use F(r)
		else:
			rstar_vals[i] = blk.master.metfunc.F(rval)
	## check that new values are finite
	if not np.all(np.isfinite(rstar_vals)):
		valid = False
	## go if valid
	if valid==True:
		## can be buggy if inf or npoints are too low
		#crvlist = rstarlines(rstar_vals, c=c, inf=2.*inf, npoints=npoints, sty=dict(marker='x'))
		crvlist = rstarlines_special_2(rstar_vals, blk.uvbounds, c=1.*c, inf=inf, npoints=npoints, eps=1e-15)
		fill_between_curves_uv(blk, crvlist, sty=sty, eps=eps)






"""
Useful utilities for filling.
"""



def reverse_curve(crv):
	"""
	Reverse the order of points in a curve.
	"""
	## reverse each 1d subarray
	crv.r    =    crv.r[::-1]
	crv.tr   =   crv.tr[:,::-1]
	crv.uv   =   crv.uv[:,::-1]
	crv.uvdl = crv.uvdl[:,::-1]
	crv.UV   =   crv.UV[:,::-1]
	## return
	return crv




def snap_to_bounds(y, ymin=-np.inf, ymax=np.inf):
	"""
	Limit an array y to the range ymin <= y <= ymax by replacing values 
	outside the range with the cap values.
	"""
	## copy y into z
	z = 1.*y
	## snap z to bounds
	z[y<ymin] = ymin + 0. * z[y<ymin] 
	z[y>ymax] = ymax + 0. * z[y>ymax] 
	## return
	return z


def snap_crvlist_to_bounds(crvlist, umin=-np.inf, umax=np.inf, vmin=-np.inf, vmax=np.inf, eps=1e-15):
	"""
	Apply snap_to_bounds to uv coordinates of a curve.
	Displaces curve into a rectangle of specified size.
	Move inward by eps to avoid edge cases.
	"""
	## define new curve list
	crvlist2 = []
	## alter min and max values
	umin, umax = umin+eps, umax-eps
	vmin, vmax = vmin+eps, vmax-eps
	## snap to bounds
	for crv in crvlist:
		crv.uv[0] = snap_to_bounds(crv.uv[0], ymin=umin, ymax=umax)
		crv.uv[1] = snap_to_bounds(crv.uv[1], ymin=vmin, ymax=vmax)
		crvlist2 += [crv]
	## return
	return crvlist2








"""
Test routines. Will run if __name__=='__main__'.
"""

def test1():
	"""
	Test matplotlib polygon patch functionality.
	"""
	##
	x1 = np.linspace(-3,3,10000)
	y1 = np.sin(2.*np.pi*x1)
	##
	x = np.concatenate([x1,y1])
	y = np.concatenate([y1,x1])
	xy = np.dstack([x,y])[0]
	##
	plt.figure()
	plt.xlim(-4,4)
	plt.ylim(-4,4)
	plt.gca().set_aspect('equal')
	##
	sty = dict(fc='c', ec='none', alpha=.5)
	p = mpl.patches.Polygon(xy, **sty)
	plt.gca().add_patch(p)
	##
	plt.show()


def test2():
	"""
	Test of snap_to_bounds.
	"""
	## define arrays to clip
	x = np.linspace(0,10,1000)
	y = np.sin(x)
	## clip
	y1 = snap_to_bounds(y)
	y2 = snap_to_bounds(y, ymin=-.8, ymax=.7)
	## plot
	plt.plot(x,y1,'b-',lw=15)
	plt.plot(x,y1,'r-',lw=10)
	plt.plot(x,y2,'y-',lw=5)
	## show
	plt.show()



if __name__=='__main__':
	print("TESTS")
	#test1()
	#test2()







