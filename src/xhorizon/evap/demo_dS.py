
def demo():
	import numpy as np
	import matplotlib.pyplot as plt
	import xhorizon as xh
	from xhorizon.evap.massplot import massplot
	from xhorizon.evap import helpers

	"""
	"""

	## setup figure #################
	## rc
	plt.rcParams['font.size'] = 7
	plt.rcParams['axes.linewidth'] = .4
	plt.rcParams['xtick.major.size'] = 2
	plt.rcParams['ytick.major.size'] = 2
	plt.rcParams['lines.solid_capstyle'] = 'butt'
	plt.rcParams['lines.dash_capstyle'] = 'butt'
	## figure size
	fig_width  = 2.5  ## inches
	fig_aspect = 1.68   ## height/width
	fig_height = fig_width * fig_aspect
	## axis size
	width  = 0.86   ## x fraction
	aspect = 1.84   ## absolute height/width
	height = aspect * width / fig_aspect    ## y fraction
	## axis loc
	left, bottom = .12, .05     ## x fraction, y fraction
	## define figure and axes
	plt.figure(1,figsize=(fig_width,fig_height))
	ax1 = plt.axes([left,bottom,width,height], aspect=1.)
	ax = [ax1]
	## format axes
	for axx in ax:
		## set axis
		plt.sca(axx)
		## labels
		plt.xlabel('$V-U$', labelpad=-7)
		plt.ylabel('$V+U$', labelpad=-4)
		## ticks
		plt.xticks([-1,1])
		plt.yticks([-1,1])
		## lims
		sq = 2.4
		x0 = -1.2
		y0 = 0.
		plt.xlim(x0,x0+sq)
		plt.ylim(y0-0.5*aspect*sq,y0+0.5*aspect*sq)
	################################


	## save
	path = "temp-figs/demo"
	sfp = dict(dpi=800)
	temp_only = False

	## draw
	draw = True

	## label
	label = '$(b)$'

	## input	
	l  = .01
	le = .01

	## input
	R  = .2

	## input
	L = 5.

	## accrete
	Tacc = 1.25


	## vv
	voff = -.4*Tacc - 0.

	## params
	seed = 0
	mpgo=True
	params = dict()
	chainparams = dict()

##################################

	## func
	func0 = xh.mf.dS(L=1.*L)
	func1 = xh.mf.Hay_dS(l=1.*l, R=1.*R, L=1.*L)
	func2 = xh.mf.Hay_dS_cutoff(l=1.*l, R=1.*R, L=1.*L)
	func3 = xh.mf.dS(L=1.*L)
	funclist = [func0, func1, func2, func3]

	## need to match 1 and 2 at infinty
	c12 = func2.F(10.*L) - func1.F(10.*L)

	## regions
	deflines = False
	reg0 = xh.reg.MAXreg(func0, boundary=deflines, rlines=deflines)
	reg1 = xh.reg.EFreg( func1, boundary=deflines, rlines=deflines)
	reg2 = xh.reg.MAXreg(func2, rparams=dict(c=1.*c12), boundary=deflines, rlines=deflines)
	reg3 = xh.reg.MAXreg(func3, boundary=deflines, rlines=deflines)
	reglist = [reg0, reg1, reg2, reg3]

	## modify Hay-AdS region to add extra blocks
	## top
	b  = reg1.blocks[3]
	cdlu = 1.*b.bparams['cdlu'] + 1.
	cdlv = 1.*b.bparams['cdlv'] + 1.
	epsu = -1.*b.bparams['epsu']
	epsv = -1.*b.bparams['epsv']
	bparams = dict(cdlu=1.*cdlu, cdlv=1.*cdlv, epsu=1.*epsu, epsv=1.*epsv)
	reg1.blocks += [xh.block(b.master, b.j, dict(cdlu=1.*cdlu, cdlv=1.*cdlv, epsu=1.*epsu, epsv=1.*epsv))]
	## right
	b  = reg1.blocks[2]
	cdlu = 1.*b.bparams['cdlu'] - 1.
	cdlv = 1.*b.bparams['cdlv'] + 1.
	epsu = -1.*b.bparams['epsu']
	epsv = -1.*b.bparams['epsv']
	bparams = dict(cdlu=1.*cdlu, cdlv=1.*cdlv, epsu=1.*epsu, epsv=1.*epsv)
	reg1.blocks += [xh.block(b.master, b.j, dict(cdlu=1.*cdlu, cdlv=1.*cdlv, epsu=1.*epsu, epsv=1.*epsv))]


	## chain params
	m  = np.array([0., 0.5*R, 0., 0.])
	Rh = np.array([0., func1.rj[2], 0., 0.])
	va = 1.*voff
	vb = 1.*va + 1.*Tacc

	#### initial region
	## passive slice
	sliceloc0 = dict(v0=1.*va, r0=3.*L)
	pslice0 = xh.junc.pslice(reg0, ublocks=[1,2], vblocks=[0,1], **sliceloc0)

	#### Hay-AdS region
	## active slice
	aslice1 = xh.junc.aslice(reg1, ublocks=[3,5], vblocks=[0,1,2,3], U0=pslice0.U_of_r_at_v0, V0=pslice0.V_of_r_at_u0, r_refs=[pslice0.reg.metfunc.r_ref], **sliceloc0)
	## modify transformations
	reg1.U_of_udl = aslice1.U_of_udl_at_v0
	reg1.V_of_vdl = aslice1.V_of_vdl_at_u0
	## passive slice
	sliceloc1 = dict(v0=1.*vb, r0=Rh[1]+le)
	pslice1 = xh.junc.pslice(reg1, ublocks=[2,4], vblocks=[0,1,2,3], **sliceloc1)

	#### Cutoff region
	## modify transformations
	reg2.U_of_udl = reg1.U_of_udl
	reg2.V_of_vdl = reg1.V_of_vdl

	#### Final region
	## active slice
	aslice3 = xh.junc.aslice(reg3, ublocks=[0,3], vblocks=[0,1], U0=pslice1.U_of_r_at_v0, V0=pslice1.V_of_r_at_u0, r_refs=[pslice1.reg.metfunc.r_ref], **sliceloc1)
	ub = 1.*aslice3.u0
	## modify transformations
	reg3.U_of_udl = aslice3.U_of_udl_at_v0
	reg3.V_of_vdl = aslice3.V_of_vdl_at_u0


	#### mask
	## initial region
	reg0.blocks = reg0.blocks[0:2]
	for b in reg0.blocks:
		b.uvbounds.update(dict(vmax=1.*va))
	# Cutoff region
	reg2.blocks = reg2.blocks[2:3]
	reg2.blocks[0].bparams.update(reg1.blocks[-1].bparams)
	reg2.blocks[0].update()
	# Hay-AdS region past slice
	reg1.blocks = reg1.blocks[0:5]
	for b in reg1.blocks[0:4]:
		b.uvbounds.update(dict(vmin=1.*va))
	# Hay-AdS future slice
	newblocks = []
	for i in [0,1]:
		newblocks += [xh.block(reg1.blocks[i].master, reg1.blocks[i].j, reg1.blocks[i].bparams)]
		newblocks[-1].uvbounds.update(reg1.blocks[i].uvbounds)
		newblocks[-1].uvbounds.update(dict(vmax=1.*vb))
	for i in [2]:
		newblocks += [xh.block(reg1.blocks[i].master, reg1.blocks[i].j, reg1.blocks[i].bparams)]
		newblocks[-1].uvbounds.update(reg1.blocks[i].uvbounds)
		newblocks[-1].uvbounds.update(dict(umax=1.*ub, vmax=1.*vb))
	for i in [2]:
		newblocks += [xh.block(reg1.blocks[i].master, reg1.blocks[i].j, reg1.blocks[i].bparams)]
		newblocks[-1].uvbounds.update(reg1.blocks[i].uvbounds)
		newblocks[-1].uvbounds.update(dict(umax=1.*ub, vmin=1.*vb))
	for i in [4]:
		newblocks += [xh.block(reg1.blocks[i].master, reg1.blocks[i].j, reg1.blocks[i].bparams)]
		newblocks[-1].uvbounds.update(reg1.blocks[i].uvbounds)
		newblocks[-1].uvbounds.update(dict(umax=1.*ub))
	for i in [3]:
		newblocks += [xh.block(reg1.blocks[i].master, reg1.blocks[i].j, reg1.blocks[i].bparams)]
		newblocks[-1].uvbounds.update(reg1.blocks[i].uvbounds)
		newblocks[-1].uvbounds.update(dict())	
	## reset Hay-AdS blocks to masked ones
	reg1.blocks = newblocks
	## final region
	reg3.blocks = [reg3.blocks[0],reg3.blocks[3]]
	reg3.blocks[0].uvbounds.update(dict(vmin=1.*vb, umin=1.*ub))
	reg3.blocks[1].uvbounds.update(dict(umin=1.*ub))

	## chainparams
	m = 1.*m
	Rh = 1.*Rh
	fs_v0 = np.array([1.*va, 1.*vb, 0., 0.])
	fs_u0 = np.array([0., 1.*ub, 0., 0.])
	fs_r0 = 1.*Rh+le
	chainparams.update(dict(m=1.*m, Rh=1.*Rh, fs_r0=1.*fs_r0, fs_u0=1.*fs_u0, fs_v0=1.*fs_v0))

	print(chainparams)

######################################

	## draw
	if draw:
		print("plot")
		pp = dict(l=1.*l, R=1.*R)
		xh.evap.draw_dS.drawreg(reglist, chainparams, fparams=pp)

	## label
	if True:
		plt.annotate(s=label, xy=(.95,.97), xycoords='axes fraction', ha='right', va='top', size=8)


	## auto param label
	if True:
		#plabel = [r"", r"$\tau_{acc}=%3s$"%(Tacc), r"$\tau_{ev}=%3s$"%(Tevap)]
		plabel = [r"$l_{ev}=%3s$"%(le), r"", r"$l=%3s$"%(l), r"$2M=%3s$"%(R), r"$L=%3s$"%(L)]
		plabel = "\n".join(plabel)
		plt.annotate(s=plabel, xy=(.95,.03), xycoords='axes fraction', ha='right', va='bottom', size=8)

	## manual param label
	if False:
		plabel = [r"$\tau_{ev}=%s$"%(Tevap), r"$\tau_{acc}=%s$"%(Tacc)]
		plabel += [r"", r"$M=10^{-3}$", r"$l_{ev}=l=10^{-4}$"]
		plabel = "\n".join(plabel)
		plt.annotate(s=plabel, xy=(.95,.03), xycoords='axes fraction', ha='right', va='bottom', size=8)

	## r=infty label
	if True:
		plt.annotate(s=r"$r=\infty$", xy=(0., 1.), ha='center', va='center', size=7)
		plt.annotate(s=r"$r=\infty$", xy=(0.,-1.05), ha='center', va='center', size=7)
		plt.annotate(s=r"$r=0$", xy=(-1.1, 0.), ha='center', va='center', size=7, rotation=90)
		plt.annotate(s=r"$r=0$", xy=(1.02, 0.), ha='center', va='center', size=7, rotation=-90)

	## save
	if True:
		xh.evap.evapsave(path=path, params=params, chainparams=chainparams, seed=seed, sfp=sfp, temp_only=temp_only, massplot=mpgo)


	## show
	if False:
		plt.show()


	## return
	return reglist, chainparams




if __name__=="__main__":
	demo()


