"""
This module contains helpers to draw with evap code.
"""

import numpy as np
import matplotlib.pyplot as plt
import copy
import matplotlib.hatch
import matplotlib as mpl
import matplotlib.colors

import xhorizon as xh
from xhorizon.shell_junction import interpolators as interp
from xhorizon.evap.helpers import irr


############## plot reglist ####################


## styles
rline_zero_sty      = dict(lw=1.5, c='0.48', zorder=11000)
singularity_sty     = dict(ls='dashed', dashes=(2,2))
rline_inf_sty       = dict(lw=1.5, c='0.48', zorder=10000)
rline_hor_sty       = dict(lw=0.4, c='k', zorder=9999)
rline_sty           = dict(ls='-', zorder=10, lw=.6, alpha=.5)
acc_shells_sty      = dict(lw=0.6, ls='dashed', dashes=(3,3), zorder=2000)
evap_shells_out_sty = acc_shells_sty
evap_shells_in_sty  = dict(lw=1., ls='dashed', dashes=(1,1.2), zorder=2000)
fill_horizons_sty   = dict(fc='none', ec='k', lw=0, hatch='c', zorder=9990)
fill_density_sty    = dict(zorder=100)
tick_sty            = dict(markersize=10, markeredgecolor='k', alpha=.3, zorder=100)


def drawreg(reglist, chainparams, fparams=dict()):
	"""
	Plot all regions in reglist.
	"""       
	l=fparams['l']
	R=fparams['R']

	## lines
	rline_zero(reglist, sty=rline_zero_sty, sty2=singularity_sty, npoints=1001)
	rline_inf_col(reglist, sty=rline_inf_sty, npoints=5001, inf=100., M=0.5*R)
	rline_hor(reglist, sty=rline_hor_sty)
	acc_shells(reglist, chainparams, sty=acc_shells_sty, inf=100., npoints=1001)
	evap_shells_out(reglist, chainparams, sty=evap_shells_out_sty, inf=100., npoints=1001)
	evap_shells_in(reglist, chainparams, sty=evap_shells_in_sty, inf=100., npoints=1001)
	make_rlines(reglist, chainparams, l=l, R=R, Tevap=5., sty=rline_sty)
	vticks(reglist, dv=0.5, sty=tick_sty)
	uticks(reglist, du=0.5, sty=tick_sty)
	tticks(reglist, dt=0.5, sty=tick_sty)
	## plot
	for reg in reglist:
		reg.rplot()
	## fill
	fill_horizons(reglist, sty=fill_horizons_sty)
	fill_density(reglist, sty=fill_density_sty)




#############################################################


## being actually used ############

def rline_zero(reglist, npoints=5001, inf=100., sty={}, sty2={}):
	"""
	Add boundary lines to regions in reglist.
	"""
	for reg in reglist:
		for b in reg.blocks:
			if b.j==0:
				## style
				style = dict(lw=0.9, c='0.5', zorder=11000)
				style.update(sty)
				## singularity style
				if b.sgnf<0.:
					style.update(sty2)
				## cutoff style
				if reg.metfunc.info['Type'] in ['Hayward - de Sitter with Cutoff']:
					lw = 1.5*(1.+ 2.)
					style.update(lw=lw)
				## curve
				cv = xh.cm.rstarlines_special_2([0.], b.uvbounds, c=1.*reg.rparams['c'], sty=style, inf=1.*inf, npoints=1.*npoints, eps=1e-12)
				b.add_curves_uv(cv)
				




def rline_inf_col(reglist, npoints=5001, inf=100., sty={}, M=.5):
	"""
	Add boundary lines to regions in reglist.
	"""
	for reg in reglist:
		for b in reg.blocks:
			## outermost blocks only
			if b.j==len(b.master.metfunc.rj)-2:
				## style
				style = dict(c='0.5', zorder=10000)
				style.update(sty)
				## rstar value
				rstar = b.master.metfunc.rstar_ref[-1]
				## mass
				m = 0.
				if 'R' in b.master.metfunc.fparams.keys():
					m = 0.5 * b.master.metfunc.fparams['R']
				## linewidth
				lw0 = 1.*style['lw']
				lw = lw0 * (1.+2.*m/M)
				style.update(dict(lw=lw))
				## curve
				style.update(solid_capstyle='projecting')
				b.add_curves_uv(xh.cm.block_boundary_2(b, sty=style, inf=2., npoints=npoints)[-1:])
				style.update(solid_capstyle='butt')
				b.add_curves_uv(xh.cm.block_boundary_2(b, sty=style, inf=inf, npoints=npoints)[-1:])


def rline_hor(reglist, npoints=5001, sty={}):
	"""
	Add horizon lines to regions in reglist.
	"""
	for reg in reglist:
		for b in reg.blocks:
			## trapped blocks only
			if b.sgnf<0.:
				## style
				style = dict(lw=0.2, c='k', zorder=9999)
				style.update(sty)
				## rstar value
				rstar = b.master.metfunc.rstar_ref[-1]
				## curve
				cvs = xh.cm.block_boundary_2(b, sty=style, inf=100., npoints=npoints)
				if b.j==0:
					cvs = cvs[-1:]
				## add
				b.add_curves_uv(cvs)

def shell_col(x):
	Cmax = .8
	return np.array([1.,1.,1.]) - 0.33*(1.+2.*x)*Cmax

def acc_shells(reglist, chainparams, sty={}, inf=50., npoints=5001):
	"""
	Add boundary lines to regions in reglist.
	"""
	Rh = chainparams['Rh']
	v0 = chainparams['fs_v0']-1e-15
	m = chainparams['m']
	M = np.max(m)
	for i in range(len(reglist[:-1])):
		reg = reglist[i]
		## accretion only
		if Rh[i+1]>Rh[i]:
			##
			for b in reg.blocks:
				## mass
				dm = m[i+1] - m[i]
				print("ACC SHELL COLOR dm/M=%s"%(dm/M))
				col = shell_col(dm/M)
				## style
				style = dict(lw=0.2, ls='dashed', dashes=(4,4), c=1.*col, zorder=2000)
				style.update(sty)
				## v value
				vv = v0[i:i+1]
				## curve
				b.add_curves_uv(xh.cm.uvlines(vv, uv='v', uvbounds=b.uvbounds, sty=style, c=0., inf=inf, npoints=npoints))


def evap_shells_out(reglist, chainparams, sty={}, inf=100., npoints=5001):
	"""
	Add boundary lines to regions in reglist.
	"""
	Rh = chainparams['Rh']
	u0 = chainparams['fs_u0']-1e-15
	m = chainparams['m']
	M = np.max(m)
	for i in range(len(reglist[:-1])):
		reg = reglist[i]
		## evap only
		if Rh[i+1]<Rh[i]:
			## outgoing
			for b in reg.blocks:
				if b.bparams['epsu']>0.:
					if not np.isfinite(b.uvbounds['vmax']):
						## mass
						dm = - (m[i+1] - m[i])
						print("EVAP OUT SHELL COLOR dm/M=%s"%(dm/M))
						col = shell_col(dm/M)
						## style
						style = dict(lw=0.2, ls='dashed', dashes=(4,4), c=1.*col, zorder=2000)
						style.update(sty)
						## u value
						uu = u0[i:i+1]
						## curve
						b.add_curves_uv(xh.cm.uvlines(uu, uv='u', uvbounds=b.uvbounds, sty=style, c=0., inf=inf, npoints=npoints))


def evap_shells_in(reglist, chainparams, sty={}, inf=5., npoints=5001):
	"""
	Add boundary lines to regions in reglist.
	"""
	Rh = chainparams['Rh']
	u0 = chainparams['fs_u0']
	v0 = chainparams['fs_v0']-1e-15
	r0 = chainparams['fs_r0']
	m = chainparams['m']
	M = np.max(m)
	for i in range(len(reglist[:-1])):
		reg = reglist[i]
		## evap only
		if Rh[i+1]<Rh[i]:
			## uv value
			uu = u0[i:i+1]
			vv = v0[i:i+1]
			rr = r0[i:i+1]
			## ingoing
			for b in reg.blocks:
				## mass
				dm = - (m[i+1] - m[i])
				print("EVAP IN SHELL COLOR dm/M=%s"%(dm/M))
				col = shell_col(dm/M)
				## style
				style = dict(lw=0.9, ls=':', c=1.*col, zorder=2000)
				style.update(sty)
				## inner blocks
				if b.rj[1]<rr:
					## values
					u = np.sort(np.concatenate([np.linspace(-1e-3,-1.,npoints//2),np.linspace(1e-3,1.,npoints//2),np.linspace(-inf,inf,npoints)]))
					v = vv+0.*u
					## make curve
					cv = xh.curve()
					cv.uv = np.array([u,v])
					cv.sty.update(style)
					## add
					b.add_curves_uv([cv])
				## outer block
				elif b.rj[0]<rr and not np.isfinite(b.uvbounds['umax']):
					## values
					u = np.linspace(uu,inf,npoints)
					v = vv+0.*u
					## make curve
					cv = xh.curve()
					cv.uv = np.array([u,v])
					cv.sty.update(style)
					## add
					b.add_curves_uv([cv])




def fill_horizons(reglist, sty={}):
	##
	for i in range(len(reglist)):
		reg = reglist[i]
		for b in reg.blocks:
			if b.sgnf < 0.:
				style = dict(fc='none', ec='k', lw=0, hatch='.', zorder=9990)
				style.update(sty)
				b.fill(sty=style)


def fill_density(reglist, sty={}):
	## regions
	for reg in reglist:
		## linespacing
		x = np.arange(0.,10.1,.5)
		## default to zero
		density = 0.*x
		## vac spacetimes
		vac = ['Schwarzschild', 'Minkowski', 'Anti de Sitter', 'de Sitter', 'Hayward - de Sitter with Cutoff']
		if reg.metfunc.info['Type'] in vac:
			for b in reg.blocks:
				b.fill(dstyle(0.))
		## hayward
		if reg.metfunc.info['Type'] in ['Hayward', 'Hayward - Anti de Sitter', 'Hayward - de Sitter']:
			## get params
			l = reg.metfunc.fparams['l']
			R = reg.metfunc.fparams['R']
			## rcore and r
			rcore = (R*l**2.)**(1./3.)
			r = rcore * x
			rmax = 100.
			r = np.concatenate([r, np.array([rmax])])
			## density 
			density = 1./(1.+x**3.)**2.   # rho/rho0 = 1/(1+x^3)^2
			## blocks
			for b in reg.blocks:
				## loop over radii
				for n in range(len(x)):
					## rr and rho
					rr = r[n:n+2]
					rho = density[n]
					## style
					style = dstyle(rho)
					## fill
					rvals=b.fill_between_r(rr, sty=style, npoints=1001, inf=50.)



def dstyle(rho, cm=plt.cm.Oranges, sty={}, ovr={}):
	"""
	Density rho in (0,1), output fill style.
	Since all Hayward have same l, all regions same scale.
	"""
	## style
	style = dict(fc='none', ec='none', lw=0, zorder=900)
	style.update(sty)
	## cnorm
	cmin, cmax = .1, .6
	cnorm = cmin + (cmax-cmin)*rho
	## color and alpha
	style.update(fc=cm(cnorm), alpha=.75)
	## override
	style.update(ovr)
	## return
	return style.copy()
		
			

def make_rlines(reglist, chainparams, l=.05, R=1., Tevap=None, sty={}):
	"""
	"""
	## params
	R = 1.*R
	l = 1.*l
	rcore = (R*l**2)**(1./3.)
	## general
	style = dict(ls='-', zorder=10, lw=.4, alpha=.2)
	style.update(sty)
	## l scale
	x = np.arange(0.,.5*R/l,.5)
	scale = l
	style.update(dict(c='c'))
	rlines(reglist, scale*x, sty=style, inf=10., npoints=2001.)
	## R scale
	x = np.arange(0.,20.01,.5)
	scale = R
	style.update(dict(c='m'))
	rlines(reglist, scale*x, sty=style, inf=25., npoints=1001.)
	## Tevap scale
	if Tevap is not None:
		x = np.arange(.8,25.01,.5)
		scale = Tevap
		style.update(dict(c='#0000f0'))
		rlines(reglist, scale*x, sty=style, inf=25., npoints=1001.)


def vticks(reglist, dv=1., inf1=100., inf2=50.+irr, sty={}):
	"""
	Remainder lets dv carry across regions.
	"""
	## init
	remainder = 0.
	## regions
	for reg in reglist:
		## blocks
		for b in reg.blocks:
			## last blocks only
			if b.j == len(reg.metfunc.rj)-2:
				## outer last blocks only
				if not np.isfinite(b.uvbounds['umin']):
					## get min and max
					vmin = np.max([b.uvbounds['vmin'], -inf2])
					vmax = np.min([b.uvbounds['vmax'],  inf2])
					## is it too small
					toosmall = False
					if vmax-vmin <= remainder:
						toosmall = True
					## if too small adjust remainder
					if toosmall==True:
						remainder = remainder - (vmax-vmin)
					## if big enough go
					if toosmall==False:
						## make array
						vv = remainder + np.arange(vmin, vmax, dv)
						## get new remainder
						remainder = dv - (vmax - vv[-1])
						## make and add curve
						crv = xh.curve()
						crv.uv = np.array([-inf1+0.*vv, 1.*vv])
						style = dict(markeredgecolor='k', alpha=0.3, ls='none', marker=(2,0,45), markersize=7, zorder=-100)
						style.update(sty)
						crv.sty.update(style)
						b.add_curves_uv([crv])




def uticks(reglist, du=1., inf1=100., inf2=50.+irr, sty={}):
	"""
	Remainder lets dv carry across regions.
	"""
	## init
	remainder = 0.
	## regions
	for reg in reglist:
		## blocks
		for b in reg.blocks:
			## last blocks only
			if b.j == len(reg.metfunc.rj)-2:
				## outer last blocks only
				if not np.isfinite(b.uvbounds['vmax']):
					## get min and max
					umin = np.max([b.uvbounds['umin'], -inf2])
					umax = np.min([b.uvbounds['umax'],  inf2])
					## is it too small
					toosmall = False
					if umax-umin <= remainder:
						toosmall = True
					## if too small adjust remainder
					if toosmall==True:
						remainder = remainder - (umax-umin)
					## if big enough go
					if toosmall==False:
						## make array
						uu = remainder + np.arange(umin, umax, du)
						## get new remainder
						remainder = du - (umax - uu[-1])
						## make and add curve
						crv = xh.curve()
						crv.uv = np.array([1.*uu, inf1+0.*uu])
						style = dict(markeredgecolor='k', alpha=.3, marker=(2,0,-45), ls='none', markersize=7, zorder=100)
						style.update(sty)
						crv.sty.update(style)
						b.add_curves_uv([crv])



def tticks(reglist, dt=1., inf1=1., inf2=50.+irr, sty={}):
	"""
	Remainder lets dt carry across regions.
	"""
	if reglist[0].metfunc.info['Type']=='Anti de Sitter':
		inf1 = 9.*reglist[0].metfunc.fparams['L']
		## init
		remainder = 0.
		## regions
		for reg in reglist:
			## blocks
			for b in reg.blocks:
				## last blocks only
				if b.j == len(reg.metfunc.rj)-2:
					## outer last blocks only
					if (not np.isfinite(b.uvbounds['vmax'])) or (not np.isfinite(b.uvbounds['umin'])):
						## get min and max
						umin = np.max([b.uvbounds['umin'], -inf2])
						umax = np.min([b.uvbounds['umax'],  inf2])
						vmin = np.max([b.uvbounds['vmin'], -inf2])
						vmax = np.min([b.uvbounds['vmax'],  inf2])
						tmin = 0.5*(vmin+umin)
						tmax = 0.5*(vmax+umax)
						## is it too small
						toosmall = False
						if tmax-tmin <= remainder:
							toosmall = True
						## if too small adjust remainder
						if toosmall==True:
							remainder = remainder - (tmax-tmin)
						## if big enough go
						if toosmall==False:
							print("HELLLOOOOOOOOO tticks")
							print(tmin, tmax, dt)
							## make array
							tt = remainder + np.arange(tmin, tmax, dt)
							## get new remainder
							remainder = dt - (tmax - tt[-1])
							## make and add curve
							crv = xh.curve()
							crv.tr = np.array([1.*tt,1.*inf1])
							style = dict(markeredgecolor='k', alpha=.3, marker=(2,0,90), ls='none', markersize=7, zorder=100)
							style.update(sty)
							crv.sty.update(style)
							b.add_curves_tr([crv])


## custom hatch
class CustomHatch(matplotlib.hatch.SmallFilledCircles):
    size = 0.1
    filled = True

    def __init__(self, hatch, density):
        density = 2
        self.num_rows = (hatch.count('c')) * density
        matplotlib.hatch.Circles.__init__(self, hatch, density)

## assign custom hatch
matplotlib.hatch._hatch_types.append(CustomHatch)



##############################

################# tools for drawing ########################################



def rlines(reglist, rvals, sty={}, npoints=2001, inf=25.):
	"""
	Add colorscaled lines of constant radius to region.
	Useful to check for matching.
	"""
	rvals = 1.*rvals
	for reg in reglist:
		for b in reg.blocks:
			for i in range(len(rvals)):
				if (b.rj[0]<rvals[i]) and (rvals[i]<b.rj[1]):
					style = dict(c='c', ls='-', lw=.5, zorder=1500)
					style.update(sty)
					rstarvals = 1.*b.master.metfunc.F(rvals[i:i+1])
					b.add_curves_uv(xh.cm.rstarlines_special_2(rstarvals, b.uvbounds, c=b.master.rparams['c'], sty=style, inf=1.*inf, npoints=1.*npoints))
	return reglist


def colorlines(reglist, rmin=1.1, rmax=25., dr=.5, sty={}, cm=plt.cm.hsv_r, npoints=2001, inf=25.):
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

##########################################################################################



