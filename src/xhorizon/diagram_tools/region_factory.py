
"""
"""

import numpy as np

from xhorizon.diagram_tools.curve_class import curve
from xhorizon.diagram_tools.diagram_classes import region, block
from xhorizon.diagram_tools import curvemakers as cm






def EFreg(func, rparams={}, io=0, lr=0, basepoint=np.array([0.,0.]), boundary=True, uvgrid=False, rlines=True):
	"""
	Create an EF region with given metfunc.
		basepoint_uvdl gives coordinates of the outermost horizon vertex.
		io makes the region ingoing (0) or outgoing (1)
		lr makes the region left-directed or right-directed
	"""
	## get func params
	N = len(func.rj)-2
	eps0 = func.sgnf(0)
	## get region params
	rparam = dict(c=0.,s0=10.)
	rparam.update(rparams)
	## create region
	reg = region(func, rparams=rparam)
	## get basepoint
	x0, y0 = basepoint
	udl0, vdl0 = 0.5*(y0-x0), 0.5*(y0+x0)
	## add blocks
	## ingoing io=0 lr=0
	if io==0 and lr==0:
		for j in range(N+1):
			cdlu = udl0 - 0.5 + (N-j)
			cdlv = vdl0 + 0.5
			reg.add_block(j, bparams=dict(cdlu=cdlu, cdlv=cdlv, epsu=func.sgnf(j), epsv=1.))
	## ingoing io=0 lr=1
	if io==0 and lr==1:
		for j in range(N+1):
			cdlu = udl0 - 0.5
			cdlv = vdl0 - 1.5 + (j)
			reg.add_block(j, bparams=dict(cdlu=cdlu, cdlv=cdlv, epsu=1., epsv=func.sgnf(j) ))
	## outgoing io=1 lr=0
	if io==1 and lr==0:
		for j in range(N+1):
			cdlu = udl0 - 1.5 + (j)
			cdlv = vdl0 - 0.5
			reg.add_block(j, bparams=dict(cdlu=cdlu, cdlv=cdlv, epsu=-func.sgnf(j), epsv=-1.))
	## outgoing io=1 lr=1
	if io==1 and lr==1:
		for j in range(N+1):
			cdlu = udl0 + 0.5
			cdlv = vdl0 - 0.5 + (N-j)
			reg.add_block(j, bparams=dict(cdlu=cdlu, cdlv=cdlv, epsu=-1., epsv=-func.sgnf(j) ))
	## add boundary
	if boundary==True:
		for b in reg.blocks:
			b.add_curves_uv(cm.block_boundary(b))
	## add grid
	if uvgrid==True:
		for b in reg.blocks:
			uvvals = np.arange(-30.,30., 0.5)
			b.add_curves_uv(cm.uvlines(uvvals, uv='uv', sty=dict(c='0.5')))
	## add default rlines
	if rlines==True:
		for b in reg.blocks:
			b.add_curves_tr(cm.default_rlines(reg.metfunc))
	## return
	return reg


def MAXreg(func, rparams={}, boundary=True, uvgrid=False, rlines=True):
	"""
	Create a maximally extended region with given metfunc.
	"""
	## get N and eps0
	N = len(func.rj)-2
	eps0 = func.sgnf(0)
	## send to relevant subroutine
	if N==0:
		reg = MAXreg0(func, rparams=rparams, boundary=boundary, uvgrid=uvgrid, rlines=rlines)
	elif N==1 and eps0==-1.:
		reg = MAXreg1a(func, rparams=rparams, boundary=boundary, uvgrid=uvgrid, rlines=rlines)
	elif N==1 and eps0== 1.:
		reg = MAXreg1b(func, rparams=rparams, boundary=boundary, uvgrid=uvgrid, rlines=rlines)
	elif N == 2 and eps0== 1.:
		reg = MAXreg2a(func, rparams=rparams, boundary=boundary, uvgrid=uvgrid, rlines=rlines)
	elif N == 2 and eps0==-1.:
		reg = MAXreg2b(func, rparams=rparams, boundary=boundary, uvgrid=uvgrid, rlines=rlines)
	else:
		reg = None
		print("\nNo existing subroutine for MAXreg at this value of N and eps0.\n")
	## format
	## add boundary
	if boundary==True:
		for b in reg.blocks:
			b.add_curves_uv(cm.block_boundary(b))
	## add grid
	if uvgrid==True:
		for b in reg.blocks:
			uvvals = np.arange(-30.,30., 0.5)
			b.add_curves_uv(cm.uvlines(uvvals, uv='uv', sty=dict(c='0.5')))
	## add default rlines
	if rlines==True:
		for b in reg.blocks:
			b.add_curves_tr(cm.default_rlines(reg.metfunc))
	## return
	return reg


def MAXreg2a(func, rparams={}, boundary=True, uvgrid=False, rlines=True):
	"""
	Create a maximally extended region with given metfunc with two horizons.
	"""
	## get func params
	N = len(func.rj)-2
	eps0 = func.sgnf(0)
	## make sure there are exactly zero horizons
	if N != 2:
		return None
	## get region params
	rparam = dict(c=0.,s0=10.)
	rparam.update(rparams)
	## create region
	reg = region(func, rparams=rparam)
	## add blocks
	## set basepoint
	udl0, vdl0 = 0., 0.
	## ingoing io=0 lr=0
	for j in [0,1,2]:
		cdlu = udl0 - 0.5 + (N-j)
		cdlv = vdl0 + 0.5
		reg.add_block(j, bparams=dict(cdlu=cdlu, cdlv=cdlv, epsu=func.sgnf(j), epsv=1.))
	## ingoing io=0 lr=1
	for j in [0]:
		cdlu = udl0 - 0.5
		cdlv = vdl0 - 1.5 + (j)
		reg.add_block(j, bparams=dict(cdlu=cdlu, cdlv=cdlv, epsu=1., epsv=func.sgnf(j) ))
	## outgoing io=1 lr=0
	for j in [0,1,2]:
		cdlu = udl0 - 1.5 + (j)
		cdlv = vdl0 - 0.5
		reg.add_block(j, bparams=dict(cdlu=cdlu, cdlv=cdlv, epsu=-func.sgnf(j), epsv=-1.))
	## outgoing io=1 lr=1
	for j in [0]:
		cdlu = udl0 + 0.5
		cdlv = vdl0 - 0.5 + (N-j)
		reg.add_block(j, bparams=dict(cdlu=cdlu, cdlv=cdlv, epsu=-1., epsv=-func.sgnf(j) ))
	## return
	return reg

def MAXreg2b(func, rparams={}, boundary=True, uvgrid=False, rlines=True):
	"""
	Create a maximally extended region with given metfunc with two horizons.
	"""
	## get func params
	N = len(func.rj)-2
	eps0 = func.sgnf(0)
	## make sure there are exactly zero horizons
	if N != 2:
		return None
	## get region params
	rparam = dict(c=0.,s0=10.)
	rparam.update(rparams)
	## create region
	reg = region(func, rparams=rparam)
	## add blocks
	## set basepoint
	udl0, vdl0 = 0., 0.
	## ingoing io=0 lr=0
	for j in [0,1,2]:
		cdlu = udl0 - 0.5 + (N-j)
		cdlv = vdl0 - 0.5
		reg.add_block(j, bparams=dict(cdlu=cdlu, cdlv=cdlv, epsu=func.sgnf(j), epsv=1.))
	## ingoing io=0 lr=1
	for j in [0]:
		cdlu = udl0 + 0.5
		cdlv = vdl0 - 1.5 + (j)
		reg.add_block(j, bparams=dict(cdlu=cdlu, cdlv=cdlv, epsu=1., epsv=func.sgnf(j) ))
	## outgoing io=1 lr=0
	for j in [0,1,2]:
		cdlu = udl0 - 1.5 + (j)
		cdlv = vdl0 + 0.5
		reg.add_block(j, bparams=dict(cdlu=cdlu, cdlv=cdlv, epsu=-func.sgnf(j), epsv=-1.))
	## outgoing io=1 lr=1
	for j in [0]:
		cdlu = udl0 - 0.5
		cdlv = vdl0 - 0.5 + (N-j)
		reg.add_block(j, bparams=dict(cdlu=cdlu, cdlv=cdlv, epsu=-1., epsv=-func.sgnf(j) ))
	## return
	return reg

def MAXreg1a(func, rparams={}, boundary=True, uvgrid=False, rlines=True):
	"""
	Create a maximally extended region with given metfunc with two horizons.
	"""
	## get func params
	N = len(func.rj)-2
	eps0 = func.sgnf(0)
	## make sure there are exactly zero horizons
	if N != 1:
		return None
	## get region params
	rparam = dict(c=0.,s0=10.)
	rparam.update(rparams)
	## create region
	reg = region(func, rparams=rparam)
	## add blocks
	## set basepoint
	udl0, vdl0 = 0., 0.
	## ingoing io=0 lr=0
	for j in [0,1]:
		cdlu = udl0 - 0.5 + (N-j)
		cdlv = vdl0 + 0.5
		reg.add_block(j, bparams=dict(cdlu=cdlu, cdlv=cdlv, epsu=func.sgnf(j), epsv=1.))
	## outgoing io=1 lr=0
	for j in [0,1]:
		cdlu = udl0 - 0.5 + (j)
		cdlv = vdl0 - 0.5
		reg.add_block(j, bparams=dict(cdlu=cdlu, cdlv=cdlv, epsu=-func.sgnf(j), epsv=-1.))
	## return
	return reg




def MAXreg1b(func, rparams={}, boundary=True, uvgrid=False, rlines=True):
	"""
	Create a maximally extended region with given metfunc with two horizons.
	"""
	## get func params
	N = len(func.rj)-2
	eps0 = func.sgnf(0)
	## make sure there are exactly zero horizons
	if N != 1:
		return None
	## get region params
	rparam = dict(c=0.,s0=10.)
	rparam.update(rparams)
	## create region
	reg = region(func, rparams=rparam)
	## add blocks
	## set basepoint
	udl0, vdl0 = 0., 0.
	## ingoing io=0 lr=0
	for j in [0,1]:
		cdlu = udl0 - 0.5 + (N-j)
		cdlv = vdl0 - 0.5
		reg.add_block(j, bparams=dict(cdlu=cdlu, cdlv=cdlv, epsu=func.sgnf(j), epsv=1.))
	## outgoing io=1 lr=0
	for j in [0,1]:
		cdlu = udl0 - 0.5 + (j)
		cdlv = vdl0 + 0.5
		reg.add_block(j, bparams=dict(cdlu=cdlu, cdlv=cdlv, epsu=-func.sgnf(j), epsv=-1.))
	## return
	return reg





def MAXreg0(func, rparams={}, boundary=True, uvgrid=False, rlines=True):
	"""
	Create a maximally extended region with given metfunc with two horizons.
	"""
	## get func params
	N = len(func.rj)-2
	eps0 = func.sgnf(0)
	## make sure there are exactly zero horizons
	if N != 0:
		return None
	## get region params
	rparam = dict(c=0.,s0=10.)
	rparam.update(rparams)
	## create region
	reg = region(func, rparams=rparam)
	## add blocks
	## set basepoint
	udl0, vdl0 = 0., 0.
	## ingoing io=0 lr=0
	for j in [0]:
		cdlu = udl0 - 0.5 + (N-j)
		cdlv = vdl0 + 0.5
		reg.add_block(j, bparams=dict(cdlu=cdlu, cdlv=cdlv, epsu=func.sgnf(j), epsv=1.))
	## return
	return reg

