
"""
This module provides test routines for verifying that metfunc object properties have been
properly calculated.
"""

import numpy as np
import matplotlib.pyplot as plt




def func_test(func, saveplots=False):
	"""
	Runs a series of tests to verify that metfunc object properties have been properly calculated.
	"""
	print('\nRUN FUNC TEST')
	print('\n' + '\n'.join(['%s: %s'%(key,func.info[key]) for key in func.info.keys()]) +'\n')
	f_test(func, saveplots)
	F_test_0(func, saveplots)
	F_test_1(func, saveplots)
	F_test_2(func, saveplots)
	#print('\n' + '\n'.join(['%s: %s'%(key,func.info[key]) for key in func.info.keys()]) +'\n')
	print('END FUNC TEST\n')





def nice_plot_background(func):
	"""
	Fill the intervals I_j, with blue as trapped and green as untrapped.
	Separate intervals with vertical grey lines.
	Plot a horizontal black dotted line at y=0.
	"""
	## styles
	fillparams = [dict(facecolor=r'b', alpha=.15, edgecolor='none', zorder=1),
				  dict(facecolor=r'g', alpha=.2 , edgecolor='none', zorder=1)]
	vert_sty = dict(ls='-', c='0.8', lw=1, zorder=5)
	zero_sty = dict(c='k',lw=.5,ls=':',zorder=3)
	##
	rj = func.rj
	ri = rj[1:-1]
	## fill intervals
	for i in range(len(rj)-1):
		plt.fill_betweenx([-100,100],rj[i],rj[i+1], **fillparams[(i+func.first_trapped_j) % 2])
	## plot interval separators
	for i in range(len(ri)):
		plt.plot([ri[i],ri[i]],[-100.,100.], **vert_sty)
	## plot zero line
	plt.plot([0,100],[0,0], **zero_sty)






def f_test(func, saveout=False):
	"""
	Plot metfunc.f(r) along with metfunc.rj and metfunc.kj to check that zeroes and slopes are accurate.
	"""
	## styles
	f_sty = dict(c='k', ls='-', lw=2, zorder=50, label=r'f(r)')
	tang_sty = dict(ls='-', c='0.8', lw=1, zorder=5, label='tangent at rj')
	## make figure and background
	plt.figure()
	nice_plot_background(func)
	## plot metric function
	rr = func.r_ref
	ff = func.f(rr)
	plt.plot(rr,ff,**f_sty)
	## plot tangents
	ri = func.rj[1:-1]
	ki = func.kj[1:-1]
	for i in range(len(ri)):
		plt.plot(rr, ki[i]*(rr-ri[i]), **tang_sty)
		if i==0:
			tang_sty.pop('label')
	## format
	plt.title('f_test: check zeroes and slopes of f(r) at r_i')
	plt.xlabel('r')
	plt.ylabel('f(r)')
	plt.xlim(0, np.max(rr))
	plt.ylim(-10,10 )
	plt.legend(loc='lower right', fontsize=10)
	if saveout==True:
		plt.savefig('f_test.pdf')
	plt.show()




def F_test_0(func, saveout=False):
	"""
	Plot metfunc.F(r) and metfunc.Finv(j,rstar) on top of metfunc.f(r) and metfunc.rj.
	Check that Finv is truly the inverse of F.
	"""
	## styles
	F_sty = dict(c='b', ls='-', lw=2, zorder=50, label='F(r)')
	Finv_sty = dict(ls='none', marker='x', c='r', lw=1, zorder=40, label='Finv(j,rstar)')
	f_sty = dict(c='0.8', ls='-', lw=1, zorder=10, label='f(r)')
	## make figure and background
	plt.figure()
	nice_plot_background(func)
	## plot
	rr = func.r_ref
	rrstar = func.F(rr)
	jj = func.winterval(rr)
	plt.plot(rr,rrstar,**F_sty)
	plt.plot(func.Finv(jj,rrstar),rrstar,**Finv_sty)
	plt.plot(rr,func.f(rr),**f_sty)
	## format
	plt.title('F_test_0: test inverse relation of F and Finv')
	plt.xlabel('r')
	plt.ylabel('')
	plt.xlim(0, np.max(rr))
	plt.ylim(-10,10 )
	plt.legend(loc='lower right', fontsize=10)
	if saveout==True:
		plt.savefig('F_test_0.pdf')
	plt.show()




def F_test_1(func, saveout=False):
	"""
	Plot the reciprocal numerical derivative dF/dr on top of f(r).
	Check that 1./(dF/dr) = f(r) so that tortoise function is an integral of f(r)^{-1}.
	"""
	## styles
	F_sty = dict(c='0.5', ls='-', lw=1.5, zorder=10, label='F(r)')
	dFdr_sty = dict(c='r', ls='none', marker='x', zorder=20, label='numerical 1./(dF/dr)')
	f_sty = dict(c='b', ls='-', lw=1.5, zorder=40, label='f(r)')
	## make figure and background
	plt.figure()
	nice_plot_background(func)
	## plot f(r) and F(r)
	rr = func.r_ref
	rrstar = func.F(rr)
	plt.plot(rr,rrstar,**F_sty)
	plt.plot(rr,func.f(rr),**f_sty)
	## get derivative dF/dr by right difference
	dF = np.roll(rrstar,-1) - rrstar
	dF[-1] = np.nan
	dr = np.roll(rr,-1) - rr
	dFdr = dF / dr
	## plot reciprocal derivative 1./dFdr vs r omitting final point
	plt.plot(rr, dFdr**(-1), **dFdr_sty)
	## format
	plt.title('F_test_1: test dF/dr = 1/f(r)')
	plt.xlabel('r')
	plt.ylabel('')
	plt.xlim(0, np.max(rr))
	plt.ylim(-10,10 )
	plt.legend(loc='lower right', fontsize=10)
	if saveout==True:
		plt.savefig('F_test_1.pdf')
	plt.show()


def F_test_2(func, saveout=False):
	"""
	Plot the continuous version of F(r), with logarithmic infinities subtracted out.
	Plot the exponential function appearing in the Penrose coordinate metric.
	Check continuity at ri.
	"""
	## styles
	F_sty = dict(c='0.7', ls='-', lw=1.5, zorder=10, label=r'$F(r)$')
	Fcont_sty = dict(c='k', ls='-', lw=4, zorder=50, label=r'$F_{cont}(r) = F(r) - \sum_i \, k_i^{-1} \ln(r-r_i)$')
	Fexp_sty = dict(ls='-', lw=2.5, zorder=40)
	Fexp_cols = ['m','c','g','r','b']
	## make figure and background
	plt.figure()
	nice_plot_background(func)
	## plot F(r)
	rr = func.r_ref
	rrstar = func.F(rr)
	jj = func.winterval(rr)
	plt.plot(rr,rrstar,**F_sty)
	## get ri, ki
	ri = func.rj[1:-1]
	ki = func.kj[1:-1]
	## get and plot Fcont = F with logarithmic infinities subtracted out
	Fcont = rrstar
	for i in range(len(ri)):
		Fcont = Fcont - ki[i]**(-1) * np.log( np.abs( rr-ri[i] ) )
	plt.plot(rr, Fcont, **Fcont_sty)
	## get and plot exponential analytic functions
	for i in range(len(ri)):
		Fexp = np.abs( func.f(rr) ) * np.exp( - ki[i] * rrstar )
		Fexp_sty.update(c=Fexp_cols[i%len(ri)], label=r'$F_{exp}(r) = |f(r)| e^{-k_i F(r)}, \quad (i=%d)$'%(i+1) )
		plt.plot(rr, Fexp, **Fexp_sty)
		Fexp_sty['zorder'] += -1
	## format
	plt.title('F_test_2: check continuity properties of F(r)')
	plt.xlabel('r')
	plt.ylabel('')
	plt.xlim(0, np.max(rr))
	plt.ylim(-10,10 )
	plt.legend(loc='lower right', fontsize=11)
	if saveout==True:
		plt.savefig('F_test_2.pdf')
	plt.show()

