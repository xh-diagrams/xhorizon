
"""
This module provides method for making forming and evaporation BH diagrams.
This module imports the entire xhorizon package. 
It is meant for a higher level usage than the other subpackages, none of the
guts of xhorizon rely on this.
"""

import numpy as np
import matplotlib.pyplot as plt
import copy, pprint
import scipy.optimize as opt

import xhorizon as xh
from xhorizon.evap.helpers import *



###############################################################################################################3



def funclist_chain(funclist, seed=0, du=None, dv=None, r0p=None, r0f=None, u0=None, v0=None, ps_matchmode=None, fs_matchmode=None):
	"""
	Create a chain of matched regions, starting at seed region which is unmodified.
	Each region except ends has two slices through it, a future slice fslice and past slice pslice.
	Each fslice and pslice can be either active or passive, but there can only be one active slice per region.
	The index i refers to each region in the sequence for all variables.

	Inputs:
		funclist = list of func objects, in order, to chain together
		seed = index value for seed region (seed region has trivial transforms to target coords)
		du = list of du values so that du[i] will always be size of region[i]
		dv = list of du values so that du[i] will always be size of region[i]
		r0p = list of r0 values for past slice so that r0p will always be ps_r0 when pslice is active
		r0f = list of r0 values for future slice so that r0f will always be fs_r0 when fslice is active
		u0 = list of offset values for range of u values in slice, defaults to zero
		v0 = list of offset values for range of v values in slice, defaults to zero
		ps_matchmode = list of strings, each either 'ru' or 'rv', to determine how past slice is sliced when pslice is active
		ps_matchmode = list of strings, each either 'ru' or 'rv', to determine how future slice is sliced when fslice is active
	"""
	print("du funclist_chain")
	print(repr(du))
	print("dv funclist_chain")
	print(repr(dv))
	## init default values
	if u0 is None:
		u0 = np.zeros(len(funclist))
	if v0 is None:
		v0 = np.zeros(len(funclist))
	if ps_matchmode==None:
		ps_matchmode = ['rv' for func in funclist]
	if fs_matchmode==None:
		fs_matchmode = ['rv' for func in funclist]
	## set irrelevant first and last du and dv values to zero
	du[0], du[-1] = 0., 0.
	dv[0], dv[-1] = 0., 0.	
	## init internal variables
	reglist = [xh.reg.EFreg(funcx, boundary=False, rlines=False) for funcx in funclist]
	pslice  = [None for funcx in funclist]
	fslice  = [None for funcx in funclist]
	Rh      = [funcx.rj[-2] for funcx in funclist]
	ps_r0   = [np.nan for funcx in funclist]
	ps_u0   = [np.nan for funcx in funclist]
	ps_v0   = [np.nan for funcx in funclist]
	fs_r0   = [np.nan for funcx in funclist]
	fs_u0   = [np.nan for funcx in funclist]
	fs_v0   = [np.nan for funcx in funclist]
	i0 = range(len(funclist))[1*seed]
	ps_matchpop = [mp(mmm) for mmm in ps_matchmode]
	fs_matchpop = [mp(mmm) for mmm in fs_matchmode]
	## seed region
	i = 1*i0
	for i in [1*i0]:
		###### past passive slice
		## past passive slice input params (not mutually consistent)
		ps_u0[i]  = u0[i] - 0.5*du[i]
		ps_v0[i]  = v0[i] - 0.5*dv[i]
		ps_r0[i]  = 1.*r0p[i]
		## get past passive slice location from inputs and matchpop
		sliceloc = dict(u0=ps_u0[i], v0=ps_v0[i], r0=ps_r0[i])
		sliceloc.pop(ps_matchpop[i])
		print("i=%s pslice loc: %s"%(i,sliceloc))
		## execute past passive slice at sliceloc
		pslice[i] = xh.junc.pslice(reglist[i], ublocks=[-1], vblocks=range(len(reglist[i].blocks)), **sliceloc)
		## update past passive slice location to true values
		ps_u0[i], ps_v0[i], ps_r0[i]  = 1.*pslice[i].u0, 1.*pslice[i].v0, 1.*pslice[i].r0
		#### future passive slice
		## future passive slice input params (not mutually consistent)
		fs_u0[i]  = 1.*ps_u0[i] + 1.*du[i]
		fs_v0[i]  = 1.*ps_v0[i] + 1.*dv[i]
		fs_r0[i]  = 1.*r0f[i]
		## get future passive slice location from inputs and matchpop
		sliceloc = dict(u0=fs_u0[i], v0=fs_v0[i], r0=fs_r0[i])
		sliceloc.pop(fs_matchpop[i])
		print("i=%s fslice loc: %s"%(i,sliceloc))
		## execute future passive slice at sliceloc
		fslice[i] = xh.junc.pslice(reglist[i], ublocks=[-1], vblocks=range(len(reglist[i].blocks)), **sliceloc)
		## update future passive slice location to true values
		fs_u0[i], fs_v0[i], fs_r0[i]  = 1.*fslice[i].u0, 1.*fslice[i].v0, 1.*fslice[i].r0
	## forward regions
	i = 1*i0 + 1
	while i < len(reglist):
		###### past active slice
		## past active slice input params (not mutually consistent)
		ps_u0[i]  = u0[i] - 0.5*du[i]
		ps_v0[i]  = v0[i] - 0.5*dv[i]
		ps_r0[i]  = 1.*fs_r0[i-1]
		## get past active slice location from inputs and matchpop
		sliceloc = dict(u0=ps_u0[i], v0=ps_v0[i], r0=ps_r0[i])
		sliceloc.pop(ps_matchpop[i])
		print("i=%s pslice loc: %s"%(i,sliceloc))
		## execute past active slice at sliceloc
		pslice[i] = xh.junc.aslice(reglist[i], ublocks=[-1], vblocks=range(len(reglist[i].blocks)), U0=fslice[i-1].U_of_r_at_v0, V0=fslice[i-1].V_of_r_at_u0, r_refs=[fslice[i-1].reg.metfunc.r_ref], **sliceloc)
		## update past active slice location to true values
		ps_u0[i], ps_v0[i], ps_r0[i]  = 1.*pslice[i].u0, 1.*pslice[i].v0, 1.*pslice[i].r0
		#### modify transformations
		## adjust transformations
		reglist[i].U_of_udl = pslice[i].U_of_udl_at_v0
		reglist[i].V_of_vdl = pslice[i].V_of_vdl_at_u0
		#### future passive slice
		## future passive slice input params (not mutually consistent)
		fs_u0[i]  = 1.*ps_u0[i] + 1.*du[i]
		fs_v0[i]  = 1.*ps_v0[i] + 1.*dv[i]
		fs_r0[i]  = 1.*r0f[i]
		## get past active slice location from inputs and matchpop
		sliceloc = dict(u0=fs_u0[i], v0=fs_v0[i], r0=fs_r0[i])
		sliceloc.pop(fs_matchpop[i])
		print("i=%s fslice loc: %s"%(i,sliceloc))
		## execute future passive slice at sliceloc
		fslice[i] = xh.junc.pslice(reglist[i], ublocks=[-1], vblocks=range(len(reglist[i].blocks)), **sliceloc)
		## update future passive slice location to true values
		fs_u0[i], fs_v0[i], fs_r0[i]  = 1.*fslice[i].u0, 1.*fslice[i].v0, 1.*fslice[i].r0
		##### iterate
		## iterate
		i += 1
	## backward regions
	i = 1*i0 - 1
	while i>=0:
		###### future active slice
		## past active slice input params (not mutually consistent)
		fs_u0[i]  = u0[i] - 0.5*du[i]
		fs_v0[i]  = v0[i] - 0.5*dv[i]
		fs_r0[i]  = 1.*ps_r0[i+1]
		## get future active slice location from inputs and matchpop
		sliceloc = dict(u0=fs_u0[i], v0=fs_v0[i], r0=fs_r0[i])
		sliceloc.pop(fs_matchpop[i])
		print("i=%s fslice loc: %s"%(i,sliceloc))
		## execute future active slice at sliceloc
		fslice[i] = xh.junc.aslice(reglist[i], ublocks=[-1], vblocks=range(len(reglist[i].blocks)), U0=pslice[i+1].U_of_r_at_v0, V0=pslice[i+1].V_of_r_at_u0, r_refs=[pslice[i+1].reg.metfunc.r_ref], **sliceloc)
		## update future active slice location to true values
		fs_u0[i], fs_v0[i], fs_r0[i]  = 1.*fslice[i].u0, 1.*fslice[i].v0, 1.*fslice[i].r0
		#### modify transformations
		## adjust transformations
		reglist[i].U_of_udl = fslice[i].U_of_udl_at_v0
		reglist[i].V_of_vdl = fslice[i].V_of_vdl_at_u0
		#### past passive slice
		## past passive slice input params (not mutually consistent)
		ps_u0[i]  = 1.*fs_u0[i] - 1.*du[i]
		ps_v0[i]  = 1.*fs_v0[i] - 1.*dv[i]
		ps_r0[i]  = 1.*r0p[i]
		## get past passive slice location from inputs and matchpop
		sliceloc = dict(u0=ps_u0[i], v0=ps_v0[i], r0=ps_r0[i])
		sliceloc.pop(ps_matchpop[i])
		print("i=%s pslice loc: %s"%(i,sliceloc))
		## execute past passive slice at sliceloc
		pslice[i] = xh.junc.pslice(reglist[i], ublocks=[-1], vblocks=range(len(reglist[i].blocks)), **sliceloc)
		## update future passive slice location to true values
		ps_u0[i], ps_v0[i], ps_r0[i]  = 1.*pslice[i].u0, 1.*pslice[i].v0, 1.*pslice[i].r0
		##### iterate
		## iterate
		i -= 1
	## make sliceparams dict
	chainparams = dict(Rh=1.*np.array(Rh), ps_u0=1.*np.array(ps_u0), ps_v0=1.*np.array(ps_v0), ps_r0=1.*np.array(ps_r0), fs_u0=1.*np.array(fs_u0), fs_v0=1.*np.array(fs_v0), fs_r0=1.*np.array(fs_r0), i0=1*i0, ps_matchmode=ps_matchmode, fs_matchmode=fs_matchmode, funclist=funclist)
	##
	print("\n")
	pprint.pprint(chainparams)
	print("\n")
	## return
	return reglist, chainparams






def chain_masker(reglist, chainparams):
	"""
	"""
	##
	for i in range(len(reglist)):
		## mask interior blocks
		for b in reglist[i].blocks[:-1]:
			## past
			if i>0:
				b.uvbounds.update(dict(vmin=chainparams['ps_v0'][i]))
			## future
			if i<len(reglist)-1:
				b.uvbounds.update(dict(vmax=chainparams['fs_v0'][i]))
		## mask final blocks for part that is always there
		for b in reglist[i].blocks[-1:]:
			## past
			if i>0:
				b.uvbounds.update(dict(vmin=chainparams['ps_v0'][i], umin=chainparams['ps_u0'][i]))
			## future
			if i<len(reglist)-1:
				b.uvbounds.update(dict(vmax=chainparams['fs_v0'][i]))
		## copy final block for parts which depend on radius change values
		for b in reglist[i].blocks[-1:]:
			## copies
			ba = xh.block(b.master, b.j, b.bparams)
			bb = xh.block(b.master, b.j, b.bparams)
			bc = xh.block(b.master, b.j, b.bparams)
			## mask a=top b=bottom c=right
			ba.uvbounds.update(dict(vmin=chainparams['fs_v0'][i], vmax= np.inf, umin=chainparams['ps_u0'][i], umax=chainparams['fs_u0'][i]))
			bb.uvbounds.update(dict(vmin=chainparams['ps_v0'][i], vmax=chainparams['fs_v0'][i], umin=-np.inf, umax=chainparams['ps_u0'][i]))
			bc.uvbounds.update(dict(vmin=chainparams['fs_v0'][i], vmax= np.inf, umin=-np.inf, umax=chainparams['ps_u0'][i]))
		## add bottom if increasing from past
		if i>0 and chainparams['Rh'][i-1]<chainparams['Rh'][i]:
			reglist[i].blocks += [bb]
		## add top if decreasing to future
		if i<len(reglist)-1 and chainparams['Rh'][i+1]<chainparams['Rh'][i]:
			reglist[i].blocks += [ba]
		## add right if both
		if i>0 and i<len(reglist)-1 and chainparams['Rh'][i-1]<chainparams['Rh'][i] and chainparams['Rh'][i+1]<chainparams['Rh'][i]:
			reglist[i].blocks += [bc]
	## add masses to chainparams
	chainparams.update(dict(m=getmm(reglist)))
	## return
	return reglist, chainparams


def shellparams_list(Rmax=1., le=.1, Nevap=5, Tevap=10., Tacc=1., Naccrete=1, functype=xh.mf.schwarzschild, fparams=dict()):
	"""
	"""
	## init
	m, du, dv = xh.evap.SSp.SSduvm(Nevap=1*Nevap, Tevap=1.*Tevap, M=0.5*Rmax, le=1.*le)
	m, du, dv = m[::-1], du[::-1], dv[::-1]
	mdudv = [m, du, dv]
	## get shellparams
	sp = []
	for i in range(len(m)):
		func = functype(R=2.*m[i], **fparams)
		sp += [dict(func=copy.deepcopy(func), Rself=1.*func.fparams['R'], du=1.*du[i], dv=1.*dv[i], le=1.*le, Tevap=1.*Tevap, Nevap=1*Nevap, mdudv=mdudv)]
	## edit final one
	sp[-1]['dv'] = 1.*Tacc/float(max(Naccrete-1,1))
	## print
	pprint.pprint(sp)
	## return
	return sp



def cp_from_fdudv(funclist, du=None, dv=None, le=None, uoff=0., voff=0., ueta=1., veta=1.):
	"""
	"""
	## init
	funclist = funclist
	reglist = [xh.reg.EFreg(funcx, boundary=None, rlines=None) for funcx in funclist]
	Rh = np.array([funclist[i].rj[-2] for i in range(len(funclist))])
	du  = 1.*du
	dv  = 1.*dv
	r0f = 1.*Rh + 1.*le
	r0p = 1.*np.roll(r0f,1)
	u0  = 1.*ueta*np.cumsum(du-du[0]) + 1.*uoff
	v0  = 1.*veta*np.cumsum(dv-dv[0]) + 1.*voff
	ps_matchmode = None #['ru' for i in range(len(funclist))]         ### Edit matchmode here
	fs_matchmode = None #['ru' for i in range(len(funclist))]         ### Edit matchmode here
	## iterator
	ii = range(len(funclist))
	## get rinf
	rinf = np.nan * Rh
	for i in ii:
		ia, ib = max(0, i-1), min(i+2, len(ii))
		rinf[i] = get_rinf_uv0(reglist[ia:ib], v0=1.*v0)
	print(rinf)
	## correct first and last r0 values
	r0p[0]  = 1.*rinf[0]
	r0f[-1] = 1.*rinf[-1]
	## correct r0 values for formation and evaporation
	for i in ii:
		## past
		if i>0:
			## accretion
			if Rh[i]>=Rh[i-1]:
				r0p[i] = 1.*rinf[i]
			## evaporation
			if Rh[i]< Rh[i-1]:
				r0p[i] = 1.*Rh[i-1] + 1.*le
		## future
		if i<len(ii)-1:
			## accretion
			if Rh[i]<=Rh[i+1]:
				r0f[i] = 1.*rinf[i]
			## evaporation
			if Rh[i]> Rh[i+1]:
				r0f[i] = 1.*Rh[i] + 1.*le
	## make cp
	cp = dict(du=1.*du, dv=1.*dv, r0p=1.*r0p, r0f=1.*r0f, u0=1.*u0, v0=1.*v0, ps_matchmode=ps_matchmode, fs_matchmode=fs_matchmode)
	# ## return
	return cp.copy()




def formevap_input(Rmax=1., le=.01, Tevap=1., Tacc=1., Nevap=5, Naccrete=5, uoff=0., voff=0., ueta=1., veta=1., functype0=xh.mf.minkowski, fparams0=dict(), functype1=xh.mf.schwarzschild, fparams1=dict()):
	"""
	Build inputs in reverse order starting from far future.

	funclist, seed=0, du=None, dv=None, r0p=None, r0f=None, u0=None, v0=None, ps_matchmode=None, fs_matchmode=None
	"""
	## init
	funclist = []
	du  = []
	dv  = []
	## final region
	funclist += [functype0(**fparams0)]
	du  += [0.]
	dv  += [0.]
	## evap
	sp = shellparams_list(Rmax=1.*Rmax, Nevap=Nevap, le=1.*le, Tevap=1.*Tevap, Naccrete=1*Naccrete, Tacc=1.*Tacc, functype=functype1, fparams=fparams1)
	for i in range(len(sp)):
		funclist += [sp[i]['func']]
		du  += [sp[i]['du']]
		dv  += [sp[i]['dv']]
	## max radius
	Rmax = sp[-1]['Rself']
	## accrete params
	RR = np.linspace(Rmax,0.5*Rmax, Naccrete)[1:]
	for R in RR:
		funclist += [functype1(R=1.*R, **fparams1)]
		du += [0.]
		dv += [Tacc/float(Naccrete-1)]
	## first region
	funclist += [functype0(**fparams0)]
	du  += [0.]
	dv  += [0.]
	## prep for output
	funclist = funclist[::-1]
	du = np.array(du[::-1])
	dv = np.array(dv[::-1])
	le = 1.*le
	## get chain params
	cp = cp_from_fdudv(funclist, du=1.*du, dv=1.*dv, le=1.*le, uoff=1.*uoff, voff=1.*voff, ueta=1.*ueta, veta=1.*veta)
	##
	pprint.pprint(cp)
	## return
	return funclist, cp



def create_evap(params, seed=0):
	"""
	Takes input parameters of the form:

	"""
	## 
	import pprint
	## print
	pprint.pprint("params = %s"%(params))
	pprint.pprint("seed = %s"%(seed))
	## formevap_input
	print("inputs")
	funclist, cp = xh.evap.formevap_input(**params)
	## funclist_chain
	print("chain")
	reglist, chainparams = xh.evap.funclist_chain(funclist, seed=seed, **cp)
	## chain_masker
	print("mask")
	reglist, chainparams = xh.evap.chain_masker(reglist, chainparams)
	## print
	pprint.pprint(chainparams)
	## return
	return reglist, chainparams


def evapsave(path="temp/temp", params=None, chainparams=None, seed=None, sfp=dict(), temp_only=False, massplot=False):
	"""
	Save figure with timestamp and txt notes.
	"""
	##
	import shutil
	import time
	import pprint
	import matplotlib.pyplot as plt
	## get path with timestamp
	ts = str(time.time()).replace(".","")
	## save figure
	print( "save...")
	plt.figure(1)
	sfpp = dict(dpi=400)
	sfpp.update(sfp)
	plt.savefig("%s_%s.png"%(path,ts), **sfpp)
	print( "save done")
	##save text
	print( "save txt")
	ff = open("%s_%s.txt"%(path,ts), 'w')
	ff.write("%s_%s\n"%(path,ts))
	ff.write('\n')
	ff.write('Input:\nparams=\n%s\nseed=\n%s\n'%(pprint.pformat(params),seed))
	ff.write('\n')
	ff.write('Output:\nchainparams=\n%s\n'%(pprint.pformat(chainparams)))
	ff.close()
	##save massplot
	if massplot==True:
		print( "save massplot...")
		xh.evap.massplot.massplotrc()
		plt.figure(99)
		plt.savefig("%s_%s_mass.png"%(path,ts), **sfpp)
		print( "save done")
	## copy to temp
	print( "copy...")
	## copy normally
	if temp_only==False:
		tempsave = shutil.copy
	if temp_only==True:
		tempsave = shutil.move
	## copy or move
	tempsave("%s_%s.png"%(path,ts), path+"_temp.png")
	tempsave("%s_%s.txt"%(path,ts), path+"_temp.txt")
	tempsave("%s_%s_mass.png"%(path,ts), path+"_temp_mass.png")
	## print
	print( "copy done")




if __name__=='__main__':
	pass








##################################################################################################################



