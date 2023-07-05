
import numpy as np
import matplotlib.pyplot as plt
from pprint import pprint
import matplotlib.patches

def massplotrc():
	## tex
	#plt.rcParams['text.usetex'] = True
	plt.rcParams['font.family'] = 'serif'
	plt.rcParams['font.serif'] = []
	#plt.rcParams['text.latex.preamble'] = []
	## rc
	plt.rcParams['font.size'] = 7
	plt.rcParams['axes.linewidth'] = .4
	plt.rcParams['xtick.major.size'] = 2
	plt.rcParams['ytick.major.size'] = 2


def massplot(chainparams, params):
	"""
	"""
	## setup figure #################
	## figure size
	fig_width  = 2.5  ## inches
	fig_aspect = 1.68   ## height/width
	fig_height = fig_width * fig_aspect
	## define figure and axes
	plt.figure(99, figsize=(fig_width,fig_height))
	## axes
	h, hh = .08, .33
	ax1 = plt.axes([.08,1.-h-hh,.89,hh]) # [left,bottom,width,height]
	ax2 = plt.axes([.08,h      ,.89,hh]) # [left,bottom,width,height]
	ax = [ax1,ax2]
	## format axes
	for axx in [ax1,ax2]:
		plt.sca(axx)
		## axes
		plt.xticks([0,1])
		plt.yticks([0,1])
		plt.ylabel(r"$m/M$", labelpad=-4)
		plt.xlim(-.3,1.3)
		plt.ylim(-.1,1.1)
		plt.grid(1, alpha=.15, zorder=0)
	## labels
	ax1.set_xlabel(r"$v/\tau_{acc}$", labelpad=-6)
	ax2.set_xlabel(r"$u/\tau_{ev}$", labelpad=-6)
	## labels
	ax1.set_title("$(d)$ Formation Dynamics")
	ax2.set_title("$(e)$ Evaporation Dynamics")
	################################

	"""
	Parameters.
	"""

	## ideal total mass
	M = 0.5 * params['Rmax']
	
	## ideal lifetimes
	Tacc  = params['Tacc']
	Tevap = params['Tevap']

	## horizon radii
	Rh = chainparams['Rh']

	## initialize masses
	m = chainparams['m']

	## icrit=nacc and mcrit
	Nacc = params['Naccrete']
	icrit = np.argmax(m)
	mcrit = m[icrit]

	## junction locations
	r0, u0, v0 = chainparams['ps_r0'], chainparams['ps_u0'], chainparams['ps_v0']
	r1, u1, v1 = chainparams['fs_r0'], chainparams['fs_u0'], chainparams['fs_v0']

	## inf
	npoints = 5001
	inf = 100.

	"""
	Accretion.
	v is split into Nacc+1 intervals, each of length dv
	    vv = past starting point of each interval
	    mm = mass in each interval
	"""

	## init
	vacc  = 1.*np.zeros(Nacc+1)
	macc  = 1.*np.zeros(Nacc+1)

	## fill vacc values
	vacc[0] = -1.*inf
	dvacc = (v1-v0)[1:Nacc]
	vacc[1:] = np.cumsum(np.concatenate([np.array([0.]),dvacc]))

	## fill macc values
	macc = 1.*m[:Nacc+1]

	## reference vacc np.arrays
	vv = Tacc * np.linspace(-1, 2, npoints)

	## mass step function
	mmacc = 0.*vv
	for i in range(Nacc+1):
		mmacc[vv>vacc[i]] = macc[i]

	## ideal mass function
	x = 1.*vv/Tacc
	macc_ideal = 0.*vv
	macc_ideal[x>0.] = M*(0.5 + 0.5*vv[x>0.]/Tacc)
	macc_ideal[x>1.] = M + 0.*vv[x>1.]

	"""
	Evap.
	u is split into Nevap+1 intervals, each of length du
	    uu = past ending point of each interval
	    mm = mass in each interval
	"""

	## Nevap
	Nevap = len(m[icrit:])-1

	## init
	uevap  = 1.*np.zeros(Nevap+1)
	mevap  = 1.*np.zeros(Nevap+1)

	## fill uevap values
	uevap[0] = -1.*inf
	duevap = (u1-u0)[-Nevap:-1]
	uevap[1:] = np.cumsum(np.concatenate([ np.array([0.]), duevap]))

	## dvevap
	dvevap = (v1-v0)[-Nevap:-1]

	## fill mevap values
	mevap = 1.*m[-Nevap-1:]


	## reference uevap np.arrays
	uu = Tevap * np.linspace(-1, 2, npoints)

	## mass step function
	mmevap = 0.*uu
	for i in range(Nevap+1):
		mmevap[uu>uevap[i]] = mevap[i]

	## ideal mass function
	y = 1.*uu/Tevap
	mevap_ideal = 0.*uu
	mevap_ideal[y<1.] = M * (1. - (uu[y<1.]/Tevap))**(1./3.)
	mevap_ideal[y<0.] = 1.*M + 0.*uu[y<0.]


	"""
	Checks.
	"""


	## more checks
	print("r1-Rh")
	print(repr(r1-Rh))
	print("Rh/2m")
	print(repr(Rh/2*m))
	print("mcrit/M")
	print(repr(mcrit/M))

	## print checks
	print("m evap")
	print(repr(mevap))
	print("du evap")
	print(repr(duevap))
	print("dv evap")
	print(repr(dvevap))




	"""
	Plot.
	"""



	## style
	fill_sty = dict(alpha=.15, edgecolor='none', zorder=100)
	border_sty = dict(alpha=.3, lw=1, zorder=150)
	ideal_sty = dict(alpha=1., lw=1.5, ls='--')

	## ax1
	plt.sca(ax1)

	## plot accretion
	plt.fill_between(vv/Tacc, 0., mmacc/M, color='b', **fill_sty)
	plt.plot(vv/Tacc, mmacc/M, color='b', **border_sty)
	plt.plot(vv/Tacc, macc_ideal/M, color='b', **ideal_sty)

	## legend
	blue1, blue2 = (0.,0.,1.,.15), (0.,0.,1.,.4)
	a1, = plt.plot([], [], label=r"$m_v\!(v)$ (ideal)", color='b', linestyle="dashed", dashes=(3,3), lw=1.5)
	a2 = matplotlib.patches.Patch(label=r"$m_v\!(v)$ (numerical)", color=blue1, ec=blue2, lw=.7)

	## make legend
	leg1 = plt.legend(handles=[a1,a2], loc='lower right', numpoints=10, fontsize=6)
	leg1.zorder = 900
	leg1.get_frame().set_linewidth(.4)

	## ax2
	plt.sca(ax2)

	## plot evaportion
	plt.fill_between(uu/Tevap, 0., mmevap/M, color='r', **fill_sty)
	plt.plot(uu/Tevap, mmevap/M, color='r', **border_sty)
	plt.plot(uu/Tevap, mevap_ideal/M, color='r', **ideal_sty)

	## legend
	red1, red2 = (1.,0.,0.,.15), (1.,0.,0.,.4)
	a1, = plt.plot([], [], label=r"$m_u\!(u)$ (ideal)", color='r', linestyle="dashed", dashes=(3,3), lw=1.5)
	a2 = matplotlib.patches.Patch(label=r"$m_u\!(u)$ (numerical)", color=red1, ec=red2, lw=.7)

	## make legend
	leg1 = plt.legend(handles=[a1,a2], loc='lower left', numpoints=10, fontsize=6)
	leg1.zorder = 900
	leg1.get_frame().set_linewidth(.4)


	## save
	if False:
		plt.savefig("temp.png", dpi=800)
		plt.show()


	## show
	if False:
		plt.show()



####################### fake input #######################
def fakeinputs():
	"""
	"""
	## params
	params={'Tevap': 2.0,
		 'Tacc': 0.4,
		 'Naccrete': 5,
		 'Rmax': 0.5,
		 'Rmin': 0.2,
		 'dv_evap': 0.5,
		 'fparams0': {},
		 'fparams1': {'l': 0.01},
		 'l': 0.1,
		 'ueta': 0.0,
		 'uoff': 0.0,
		 'veta': 1.0,
		 'voff': -0.4}

	## chainparams
	chainparams= {'Rh': np.array([ 0.        ,  0.23932156,  0.29934074,  0.3593317 ,  0.41930663,
			        0.47927155,  0.45431214,  0.42702186,  0.39668451,  0.3622525 ,
			        0.32190157,  0.27169069,  0.19949748,  0.        ]),
			 'fs_r0': np.array([ 7.6411929 ,  7.56156025,  7.45288371,  7.35467021,  7.26557663,
			        0.57927155,  0.55431214,  0.52702186,  0.49668451,  0.4622525 ,
			        0.42190157,  0.37169069,  0.29949748,  7.68776392]),
			 'fs_u0': np.array([-15.6823858 , -17.00840717, -16.95371909, -16.89077914,
			       -16.82124943,   1.09869074,   1.52189764,   1.94109687,
			         2.35540545,   2.76383322,   3.16496844,   3.5563797 ,
			         3.93435595, -11.37552785]),
			 'fs_v0': np.array([-0.4 , -0.25, -0.15, -0.05,  0.05,  0.75,  1.25,  1.75,  2.25,
			        2.75,  3.25,  3.75,  4.25,  4.  ]),
			 'm': np.array([ 0.        ,  0.11987007,  0.14983759,  0.17980511,  0.20977263,
			        0.23974014,  0.22726618,  0.21362808,  0.19846838,  0.18126438,
			        0.16110626,  0.13602963,  0.1       ,  0.        ]),
			 'ps_r0': np.array([ 7.6411929 ,  7.6411929 ,  7.56156025,  7.45288371,  7.35467021,
			        7.26557663,  0.57927155,  0.55431214,  0.52702186,  0.49668451,
			        0.4622525 ,  0.42190157,  0.37169069,  0.29949748]),
			 'ps_u0': np.array([-15.6823858 , -17.27285923, -17.28010953, -17.19723407,
			       -17.11028635, -16.81971891,   0.76917269,   1.17997904,
			         1.58456074,   1.98070854,   2.36565531,   2.73424739,
			         3.07213381,   3.40100505]),
			 'ps_v0': np.array([-0.4 , -0.35, -0.25, -0.15, -0.05,  0.25,  0.75,  1.25,  1.75,
			        2.25,  2.75,  3.25,  3.75,  4.  ])}

	## return
	return chainparams, params
#######################################################################################


###########################################

if __name__=='__main__':
	## fake
	fake = True
	## go
	if fake:
		cp, p = fakeinputs()
		massplot(cp, p)
	## demo
	else:
		import xhorizon as xh
		xh.evap.demo()


###########################################
