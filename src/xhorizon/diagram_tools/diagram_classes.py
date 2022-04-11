
"""
This module defines the classes
	diagram
	region
	block
which are used to construct and plot diagrams.
"""

import numpy as np
import matplotlib.pyplot as plt

from xhorizon.diagram_tools import coord_transf as coord
from xhorizon.diagram_tools.block_masks import rstar_minmax, uv_range
from xhorizon.diagram_tools import block_fill



class diagram:

	"""
	A diagram is fundamentally a collection of curves.
	A diagram object consists of:
		(1) a list of curves in UV coordinates, and 
		(2) a list of sss regions, each of which bring their own curves.
	"""

	def __init__(self):
		self.curves = []
		self.regions = []

	def add_curves_UV(self, crvlist):
		pass

	def dplot(self):
		### "diagram plot"
		for reg in self.regions:
			reg.rplot()




class region:

	"""
	A region object corresponds to a single piece of an sss spacetime with metric function f(r).
	This will often, but not always, be a single Eddington-Finklestein region of the spacetime.
	A region object consists of:
		(1) a metfunc object which determines the metric function f(r) and its associated transformations
		(2) a list of curve objects explicitly belonging to the region
		(3) a list of block objects belonging to the region
		(4) methods providing transformations from region to diagram coordinates
	"""

	def __init__(self, metfunc, rparams={}):
		self.metfunc = metfunc
		self.rparams = dict(c=0., s0=10.)
		self.rparams.update(rparams)
		self.curves = []
		self.blocks = []

	"""
	UV to uvdl transformations.
	"""

	def UV_of_uvdl(self, uvdl):
		udl, vdl = uvdl
		U = self.U_of_udl(udl)
		V = self.V_of_vdl(vdl)
		UV = np.array([U,V])
		return UV

	def uvdl_of_UV(self, UV):
		U, V = UV
		udl = self.udl_of_U(U)
		vdl = self.vdl_of_V(V)
		uvdl = np.array([udl,vdl])
		return uvdl

	"""
	Null direction transformations.
	"""

	## forward

	def U_of_udl(self,udl):
		U = 1. * udl
		return U

	def V_of_vdl(self,vdl):
		V = 1. * vdl
		return V

	## backward

	def udl_of_U(self,U):
		udl = 1. * U
		return udl

	def vdl_of_V(self,V):
		vdl = 1. * V
		return vdl

	"""
	Utilities.
	"""

	def add_block(self, j, bparams):
		"""
		Add a new block to the region, given index j and parameters bparams.
		"""
		b = block(self, j, bparams)
		if b.consistency_check()==True:
			self.blocks += [b]
		return b

	def add_curves_uvdl(self, crvlist):
		pass
	
	"""
	Plotting.
	"""

	def rplot(self):
		### "region plot"
		for b in self.blocks:
			b.bplot()




class block:

	"""
	A block object corresponds to a single sss block.
	It is owned by a region, its "master", which determines its metric function and region params c,s0.

	The block object does most of the heavy lifting for creating a diagram.
	"""

	def __init__(self, master, j, bparams={}):
		## primary attributes
		self.master = master
		self.j = int(j)
		self.bparams = dict(cdlu=0., cdlv=0., epsu=1., epsv=1.)
		self.bparams.update(bparams)
		self.curves = []
		self.masks = []
		self.uvbounds = dict(umin=-np.inf, umax=np.inf, vmin=-np.inf, vmax=np.inf)
		## derived attributes
		self.rj = None
		self.kj = None
		self.kpm = None
		self.sgnf = None
		self.F    = None
		self.Finv = None
		self.s0 = None
		self.c = None
		## coordinate transformations
		self.uv_of_tr = None
		self.tr_of_uv = None
		self.r_of_uv  = None
		self.uvdl_of_uv = None
		self.uv_of_uvdl = None
		## supplementary functions
		self.h = None
		self.hinv = None
		## set up block
		self.masks += [ rstar_minmax, uv_range ]
		self.update()

	"""
	Methods for use by the user.
	"""

	def update(self):
		"""
		Use current block primary attributes to set derived attributes and coordinate transformations.
		Check for consistency.
		"""
		self.update_derived_attributes()
		self.consistency_check()
		self.update_coord_transf()

	def add_curves_uv(self, crvlist):
		"""
		Add a curve to the block, using uv coords to fill all other coordinate fields.
		"""
		self.curves += self.apply_masks(self.update_curves_from_uv(crvlist))

	def refresh_curves_from_uv(self):
		"""
		Recalculate coordinates of all curves based on the current uv coords.
		"""
		self.curves = self.apply_masks(self.update_curves_from_uv(self.curves))

	def add_curves_tr(self, crvlist):
		"""
		Add a curve to the block, using tr coords to fill all other coordinate fields.
		Only curves lying entirely within the block can be added.
		If any points have an invalid radius, whole curve is discarded.
		"""
		self.curves += self.apply_masks(self.update_curves_from_tr(crvlist))

	def update_curves_from_tr(self, crvlist):
		"""
		Update curves to the block, using tr coords to fill all other coordinate fields.
		Only curves lying entirely within the block can be added.
		If any points have an invalid radius, whole curve is discarded.
		"""
		crvlist2 = []
		for crv in crvlist:
			t, r = crv.tr
			mask = np.logical_and( self.rj[0]<r , r<self.rj[1] )
			if not np.any(mask==False):
				crv.uv = self.uv_of_tr(crv.tr)
				crvlist2 += [crv]
		crvlist2 = self.update_curves_from_uv(crvlist2)
		return crvlist2

	def fill(self, sty={}):
		"""
		Fill in the block with a background patch.
		"""
		block_fill.fill_block(self, sty=sty)

	def fill_between_r(self, rvals=np.array([0.,0.]), sty={}, npoints=1000, inf=100):
		"""
		Fill block between r values if valid.
		"""
		block_fill.fill_between_r(self, rvals=rvals, sty=sty, npoints=npoints, inf=inf)


	def bplot(self):
		"""
		Plot all curves belonging to the block.
		"""
		self.refresh_curves_from_uv()
		for crv in self.curves:
			U, V = crv.UV[0], crv.UV[1]
			plt.plot(V-U, V+U, **crv.sty)

	"""
	Methods called by other methods.
	"""

	def update_derived_attributes(self):
		"""
		Use the current block primary attributes to set the derived attributes.
		"""
		## from metfunc
		self.rj = self.master.metfunc.rj[self.j:self.j+2]
		self.kj = self.master.metfunc.kj[self.j:self.j+2]
		self.kpm = np.sort( self.kj )
		self.sgnf = self.master.metfunc.sgnf(self.j)
		self.F    = lambda r: self.master.metfunc.F(r)
		self.Finv = lambda rstar: self.master.metfunc.Finv(self.j, rstar)
		## from region
		self.s0 = float(self.master.rparams['s0'])
		self.c  = float(self.master.rparams['c'])

	def update_coord_transf(self):
		"""
		Use the current block attributes to set the coordinate transformations.
		"""
		## extract block parameters from dict
		cdlu, cdlv = float(self.bparams['cdlu']), float(self.bparams['cdlv'])
		epsu, epsv = float(self.bparams['epsu']), float(self.bparams['epsv'])
		## define coord transformations with current parameters
		self.uv_of_tr   = lambda tr: coord.uv_of_tr(tr, self.F, self.c)
		self.tr_of_uv   = lambda uv: coord.tr_of_uv(uv, self.Finv, self.c)
		self.r_of_uv    = lambda uv: coord.r_of_uv(uv, self.Finv, self.c)
		self.uvdl_of_uv = lambda uv:   coord.uvdl_of_uv(uv  , cdlu, cdlv, epsu, epsv, self.kpm, self.s0)
		self.uv_of_uvdl = lambda uvdl: coord.uv_of_uvdl(uvdl, cdlu, cdlv, epsu, epsv, self.kpm, self.s0)
		## define supplementary functions for reference
		self.h    = lambda s: coord.hks(s, self.kpm, self.s0)
		self.hinv = lambda s: coord.hksinv(s, self.kpm, self.s0)

	def consistency_check(self):
		"""
		Check if current block attributes are self-consistent.
		True if consistent. False if not.
		"""
		## initialize
		consistent = True
		messages = [
			"WARNING: Block attributes are not self-consistent.",
			"Block Info: %s"%(self.master.metfunc.info),
			"Block Info: j = %r"%(self.j),
			"Issues:",
			]
		## # valid interval label?
		if not self.j in range(len(self.master.metfunc.rj)-1):
			consistent = False
			messages += ["Check block index.", "j=%r, jmax=%r"%( self.j, (len(self.master.metfunc.rj)-2) )]
		## # valid block orientation?
		epsu, epsv = float(self.bparams['epsu']), float(self.bparams['epsv'])
		if not epsu*epsv == self.sgnf:
			consistent = False
			messages += ["Check block orientation.", "epsu*epsv.sgnf = %r"%( epsu*epsv/self.sgnf )]
		## # warn if not consistent
		if consistent == False:
			messages += ["END WARNING"]
			print('\n' + '\n'.join(messages) + '\n')
		## return
		return consistent

	def update_curves_from_uv(self, crvlist):
		"""
		Update curves in the block, using uv coords to fill all other coordinate fields.
		"""
		for crv in crvlist:
			crv.r    = self.r_of_uv(crv.uv)
			crv.tr   = self.tr_of_uv(crv.uv)
			crv.uvdl = self.uvdl_of_uv(crv.uv)
			crv.UV   = self.master.UV_of_uvdl(crv.uvdl)
		return crvlist

	def apply_masks(self, crvlist):
		"""
		Apply all masks to self.curves.
		Masks are functions of the form 
			mask(block, crvlist): return new_crvlist
		"""
		for mask in self.masks:
			crvlist = mask(self, crvlist)
		return crvlist








