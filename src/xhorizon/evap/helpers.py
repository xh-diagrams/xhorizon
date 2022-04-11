"""
This module contains helpers to help with evap code.
"""

import numpy as np
import matplotlib.pyplot as plt
import copy

import xhorizon as xh
from xhorizon.shell_junction import interpolators as interp

### constants ###
irr = .5 * ((np.pi+np.e)%1.)

################ tools for choosing junction corner radius ###############



def get_rinf_uv0(reglist, v0=[]):
	"""
	Given a list of regions, assume that `outermost' region has index [-1].
	Assume r->inf as u->-inf.
	Find the r0 value corresponding to (v=v0 and u=2.*s0) in each region.
	Take minimum r0 over reglist.
	Return a value slightly less than r0.
	"""
	r0 = np.nan * np.ones(len(reglist))
	for i in range(len(reglist)):
		reg = reglist[i]
		for b in [reg.blocks[-1]]:
			u = np.array([-2.*reg.rparams['s0']])
			v = np.array([v0[i]])
			uv = np.array([u,v])
			r0[i] = b.tr_of_uv(uv)[1]
	r0 = .85*np.min(r0)
	## bug fix for AdS
	if not np.isfinite(r0):
		if 'L' in reglist[0].metfunc.fparams.keys():
			r0=1.9*reglist[0].metfunc.fparams['L']
	## return
	return r0


###############################################################################
				

############# tools to rescale all coordinates at the end ###################

def UVcompose1(reg, fU=None, fV=None):
	"""
	Take the final UV coords for region and compose with a function in each null direction.
	"""
	reg2 = copy.deepcopy(reg)
	if not fU==None:
		reg.U_of_udl = lambda x: fU(reg2.U_of_udl(x))
	if not fV==None:
		reg.V_of_vdl = lambda x: fU(reg2.V_of_vdl(x))
	return reg

def UVcompose(reglist, fU=None, fV=None):
	"""
	Take the final UV coords for region and compose with a function in each null direction.
	Takes list input.
	"""
	for reg in reglist:
		reg = UVcompose1(reg, fU=fU, fV=fV)
	return reglist

def get_uvdl_of_UV(reg):
	"""
	Get current udl_of_U and vdl_of_V for a region by evaluation and interpolation.
	"""
	## copy region
	reg2 = copy.deepcopy(reg)
	## get x_ref arrays
	ss = np.linspace(-4,4,50001)
	uudl = 1.*ss
	vvdl = 1.*ss
	## get y_ref arrays
	UU = reg.U_of_udl(uudl)
	VV = reg.V_of_vdl(vvdl)
	## make functions
	udl_of_U = lambda U: interp.interp_with_smooth_extrap(U, UU, uudl, mu=0.)
	vdl_of_V = lambda V: interp.interp_with_smooth_extrap(V, VV, vvdl, mu=0.)
	## return
	return udl_of_U, vdl_of_V











############################################################################


############## plot reglist ####################

def rgp(reglist):
	"""
	Plot all regions in reglist.
	"""
	## lines
	xh.evap.colorlines(reglist, sty=dict(alpha=0.5, c='k'))
	#xh.evap.uv_lines(reglist, uv='u', sty=dict())
	xh.evap.boundarylines(reglist, sty=dict(), npoints=5001)
	xh.evap.s0_lines(reglist, sty=dict(lw=3))
	## figure
	plt.figure(figsize=(6,6))
	plt.axes([.1, .1, .8, .8])
	plt.xlim(-3,3)
	plt.ylim(-3,3)
	plt.gca().set_aspect('equal')
	## plot
	for reg in reglist:
		reg.rplot()
	## fill
	xh.evap.fill_by_R_2(reglist, col='r', amax=.2)

#############################################################


################# tools for drawing ########################################

def boundarylines(reglist, npoints=5001, sty={}):
	"""
	Add boundary lines to regions in reglist.
	"""
	for reg in reglist:
		for b in reg.blocks:
			style = dict(lw=0.9, c='0.5', zorder=2000)
			style.update(sty)
			b.add_curves_uv(xh.cm.block_boundary_2(b, sty=style, npoints=npoints))


def colorlines(reglist, rmin=0.05, rmax=5., dr=.2, sty={}, cm=plt.cm.hsv_r, npoints=2001, inf=25.):
	"""
	Add colorscaled lines of constant radius to region.
	Useful to check for matching.
	"""
	rvals = np.arange(rmin,rmax,dr)
	colvals = (rvals - np.min(rvals)) / (np.max(rvals) - np.min(rvals))
	cols = cm(colvals)
	for reg in reglist:
		for b in reg.blocks:
			for i in range(len(rvals)):
				if (b.rj[0]<rvals[i]) and (rvals[i]<b.rj[1]):
					style = dict(c=cols[i], ls='-', markersize=1, lw=.5, zorder=1500)
					style.update(sty)
					rstarvals = 1.*b.master.metfunc.F(rvals[i:i+1])
					b.add_curves_uv(xh.cm.rstarlines_special_2(rstarvals, b.uvbounds, c=b.master.rparams['c'], sty=style, inf=2.*inf, npoints=npoints))
	return reglist

def s0_lines(reglist, sty={}, npoints=1001, inf=50., eps=1e-24):
	"""
	"""
	for reg in reglist:
		for b in reg.blocks:
			x = np.array([(2.*reg.rparams['s0'] - eps)])
			style1 = dict(c='m', alpha=.5, lw=3, ls=':', zorder=5000)
			style2 = dict(c='c', alpha=.5, lw=3, ls=':', zorder=5000)
			style1.update(sty)
			style2.update(sty)
			b.add_curves_uv(xh.cm.uvlines( x, uv='uv', uvbounds=b.uvbounds, sty=style1, c=0., inf=inf, npoints=npoints))
			b.add_curves_uv(xh.cm.uvlines(-x, uv='uv', uvbounds=b.uvbounds, sty=style2, c=0., inf=inf, npoints=npoints))
	return reglist


def uv_lines(reglist, uv='uv', sty={}, npoints=1001, inf=50., eps=1e-24):
	"""
	"""
	for reg in reglist:
		for b in reg.blocks:
			smin, smax, ds = -5., 25., 1.
			vals = np.arange(smin, smax, ds)
			cm = plt.cm.gist_rainbow
			cv = np.linspace(0,1,len(vals))
			for i in range(len(vals)):
				style1 = dict(c=cm(cv[i]), lw=.8, ls='-', zorder=6000)
				style1.update(sty)
				b.add_curves_uv(xh.cm.uvlines([vals[i]], uv=uv, uvbounds=b.uvbounds, sty=style1, c=0., inf=inf, npoints=npoints))
	return reglist

##########################################################################################




######################### tools for filling #############################################

def fillcols_by_R(reglist):
	"""
	Get fill color values based on radius.
	"""
	## get rvals
	Rvals = np.zeros(len(reglist))
	for i, reg in enumerate(reglist):
		if 'R' in reg.metfunc.fparams.keys():
			x = reg.metfunc.fparams['R']
			Rvals[i] = x
	## get colvals
	norm = np.max(Rvals)
	colvals = 0.2 + 0.7 * Rvals / norm
	#print("Rvals = %s"%(Rvals))
	#print("colvals = %s"%(colvals))
	## return
	return colvals

def alphas_by_R_2(reglist, amax=.9):
	"""
	Get fill color values based on radius.
	"""
	## rvals
	Rvals = np.zeros(len(reglist))
	for i, reg in enumerate(reglist):
		if 'R' in reg.metfunc.fparams.keys():
			x = reg.metfunc.fparams['R']
			Rvals[i] = 1.*x
	## log
	Rvals = np.log(1.+Rvals)
	## alphas
	amin, amax = 0., amax
	avals = amin + (amax-amin)*(Rvals/np.max(Rvals))
	## return
	return 1.*avals


def fill_by_R(reglist, cm=plt.cm.prism):
	colvals = fillcols_by_R(reglist)
	for i in range(len(reglist)):
		reg = reglist[i]
		col = cm(colvals[i])
		for b in reg.blocks:
			sty = dict(fc=col, alpha=.2)
			b.fill(sty=sty)

def fill_by_R_2(reglist, col='r', amax=.9):
	##
	alphas = alphas_by_R_2(reglist, amax=amax)
	for i in range(len(reglist)):
		reg = reglist[i]
		for b in reg.blocks:
			sty = dict(fc=col, alpha=alphas[i])
			b.fill(sty=sty)

############################################################################################


################################ checkers ################################

def evaporation_input_good(R, du, dv):
	"""
	"""
	## initialize
	good = True
	## check
	if not (len(R)==len(du) and len(R)==len(dv)):
		print("MISMATCHED INPUT ARRAY LENGTHS")
		good = False
	if not np.all(R>0.):
		print("RADIUS MUST BE POSITIVE")
		good = False
	## return
	return good

def check_uvr(reglist):
	"""
	"""
	## first values
	iu = [0,0]
	iv = [1,1]
	uu = [reglist[0].blocks[-1].uvbounds['umin'], reglist[0].blocks[-1].uvbounds['umax']]
	vv = [reglist[1].blocks[-1].uvbounds['vmin'], reglist[1].blocks[-1].uvbounds['vmax']]
	## middle values
	for i in range(len(reglist)):
		reg = reglist[i]
		## uu
		a, b = reg.blocks[-1].uvbounds['umin'], reg.blocks[-1].uvbounds['umax']
		if np.isfinite(a) and np.isfinite(b):
			uu += [a,b]
			iu += [i,i]
		## vv
		a, b = reg.blocks[-1].uvbounds['vmin'], reg.blocks[-1].uvbounds['vmax']
		if np.isfinite(a) and np.isfinite(b):
			vv += [a,b]
			iv += [i,i]
	## last values
	iu += [-1,-1]
	iv += [-1,-1]
	uu += [reglist[-1].blocks[-1].uvbounds['umin'], reglist[-1].blocks[-1].uvbounds['umax']]
	vv += [reglist[-1].blocks[-1].uvbounds['vmin'], reglist[-1].blocks[-1].uvbounds['vmax']]
	## radius
	RRu = []
	for i in iu:
		if 'R' in reglist[i].metfunc.fparams.keys():
			RRu += [reglist[i].metfunc.fparams['R']]
		else:
			RRu += [0.]
	## radius
	RRv = []
	for i in iv:
		if 'R' in reglist[i].metfunc.fparams.keys():
			RRv += [reglist[i].metfunc.fparams['R']]
		else:
			RRv += [0.]
	## radius
	rru = []
	for k in range(len(iu)):
		rru += [reglist[iu[k]].blocks[-1].r_of_uv(np.array([[uu[k]],[vv[k]]]))[0]]
	## radius
	rrv = []
	for k in range(len(iv)):
		rrv += [reglist[iv[k]].blocks[-1].r_of_uv(np.array([[uu[k]],[vv[k]]]))[0]]
	## array
	iu,iv,uu,vv,rru,rrv,RRu,RRv = np.array(iu),np.array(iv),np.array(uu),np.array(vv),np.array(rru),np.array(rrv),np.array(RRu),np.array(RRv)
	## print
	print("\n")
	print("CHECK_UVR:")
	print("iu  = %r"%(iu))
	print("iv  = %r"%(iv))
	print("uu  = %r"%(uu))
	print("vv  = %r"%(vv))
	print("rru = %r"%(rru))
	print("rrv = %r"%(rrv))
	print("RRu = %r"%(RRu))
	print("RRv = %r"%(RRv))
	print("rru/RRu = " + repr(rru/RRu))
	print("\n")


############################################################################


######### get mass parameters #######
def getmm(reglist):
	## initialize masses
	m = np.zeros(len(reglist))
	## loop
	for i in range(len(reglist)):
		fp = reglist[i].metfunc.fparams
		if 'R' in fp.keys():
			m[i] = 0.5*fp['R']
	## return
	return 1.*m
#############################


############## split region into abcd ####################

def split_reg_abcd(reg, abcd='abcd', u0=0., v0=0.):
	"""
	"""
	##
	reglist = []
	##
	for x in abcd:
		reglist += [xh.cornermask.EFreg(copy.deepcopy(reg), abcd=x, u0=u0, v0=v0)]
	## return
	return reglist

#############################################################


#################### other #########################################


def sliceprint(pslice, aslice, Rh1, Rh2):
	print("\n")
	print("                          %22r, %22r, %22r, %22r, %22r"%('Rh', 'r', 'r-Rh', 'u', 'v'))
	print("PSlice:")
	print("Rh1, r1, r1/Rh1, u1, v1 = %22r, %22r, %22r, %22r, %22r"%(Rh1, pslice.r0, pslice.r0-Rh1, pslice.u0, pslice.v0))
	print("ASlice:")
	print("Rh2, r2, r2/Rh2, u2, v2 = %22r, %22r, %22r, %22r, %22r"%(Rh2, aslice.r0, aslice.r0-Rh2, aslice.u0, aslice.v0))
	print("\n")


def mp(matchmode):
	matchpop = 'uvr0'
	for s in matchmode:
		matchpop = matchpop.replace(s,'')
	return matchpop

#############################################################

