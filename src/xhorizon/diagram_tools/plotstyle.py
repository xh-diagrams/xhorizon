
"""
This module provides routines for setting the matplotlib style and creating new figures in the proper style.
"""

import matplotlib as mpl
import matplotlib.pyplot as plt



def pubrc(tex=False):
	## if matplotlib version 2, use classic style
	if int(mpl.__version__.split('.')[0]) > 1:
		plt.style.use('classic')
	## set default rc params for diagrams
	plt.rcParams['font.family'] = 'serif'
	## latex settings
	if tex==True:
		plt.rcParams['text.usetex'] = True
		plt.rcParams['font.family'] = 'serif'
		plt.rcParams['font.serif'] = []
		plt.rcParams['text.latex.preamble'] = [r'\renewcommand{\seriesdefault}{\bfdefault}']
	## line params
	plt.rcParams['lines.linewidth'] = 0.5
	## font params
	plt.rcParams['font.size'] = 10
	## legend params
	plt.rcParams['legend.fontsize'] = 'medium'



def newfig(sqfig=4, sqaxis=4, tilde=False, tex=False):
	## set rc params
	pubrc(tex=tex)
	## new figure
	plt.figure(figsize=(sqfig,sqfig))
	## new axes
	left, bottom, right, top = .105, .105, .97, .97
	width, height = right-left, top-bottom
	plt.axes([left, bottom, width, height])
	## format axes
	xpad = 4
	ypad = 0
	plt.xlabel(r'$V-U$', labelpad=xpad)
	plt.ylabel(r'$V+U$', labelpad=ypad)
	## tilde instead of cap?
	if tilde==True:
		plt.xlabel(r'$\tilde{v}-\tilde{u}$', labelpad=xpad)
		plt.ylabel(r'$\tilde{v}+\tilde{u}$', labelpad=ypad)
	## range
	sqaxis = sqaxis * 1.1
	plt.xlim(-sqaxis,sqaxis)
	plt.ylim(-sqaxis,sqaxis)







