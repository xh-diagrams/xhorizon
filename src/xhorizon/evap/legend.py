

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from xhorizon.evap.draw import *
import matplotlib.colors


## legend
def legend():

	## setup legend #################
	## tex
	plt.rcParams['text.usetex'] = True
	plt.rcParams['font.family'] = 'serif'
	plt.rcParams['font.serif'] = []
	plt.rcParams['text.latex.preamble'] = []
	## rc
	fs = 5
	plt.rcParams['font.size'] = fs
	plt.rcParams['axes.linewidth'] = .4
	## figure size
	fig_height = 1.68 * 2.5  ## inches
	fig_width = 1.54  ## inches
	## axis loc
	## define figure and axes
	plt.figure(1, figsize=(fig_width, fig_height))
	## axes loc
	x = .31
	ax1 = plt.axes([.025,x+.047,.95,.43 ]) ## left bottom width height
	ax2 = plt.axes([.075,x,     .85,.037]) ## left bottom width height
	## format axes
	for axx in [ax1,ax2]:
		## set axis
		plt.sca(axx)
		## lims
		plt.xlim(0,1)
		plt.ylim(0,1)
		## ticks
		plt.xticks([])
		plt.yticks([])
		## off
		#plt.axis('off')
	################################

	if True:
		## legend
		plt.sca(ax1)

		## init
		labels = []
		hand   = []


		## r=0
		labels += [r"$r=0$ (no singularity)"]
		hand += [mpl.lines.Line2D([], [], **rline_zero_sty)]
		## r=0 singularity
		labels += [r"$r=0$ (singularity)"]
		rline_zero_sty.update(singularity_sty)
		hand += [mpl.lines.Line2D([], [], **rline_zero_sty)]

		## blank
		labels += [r""]
		hand += [mpl.lines.Line2D([], [], c='w')]

		## horizons
		labels += [r"horizons ($f(r)=0$)"]
		hand += [mpl.lines.Line2D([], [], **rline_hor_sty)]
		## trapped region
		labels += [r"trapped spheres region"]
		hand += [mpl.patches.Patch(fc='w', ec='k', hatch='.', lw=.3)]

		## blank
		labels += [r""]
		hand += [mpl.lines.Line2D([], [], c='w')]

		## shells
		labels += [r"positive-mass shells ($\sigma>0$)"]
		acc_shells_sty.update(dict(dashes=(2.5,2.5), lw=.8))
		hand += [mpl.lines.Line2D([], [], c='0.5', **acc_shells_sty)]
		## shells
		labels += [r"negative-mass shells ($\sigma<0$)"]
		evap_shells_in_sty.update(dict())
		hand += [mpl.lines.Line2D([], [], c='0.5', **evap_shells_in_sty)]
		## shells
		labels += [r"(shell darkness $ \sim \Delta m /M$)"]
		hand += [mpl.lines.Line2D([], [], c='w')]

		## blank
		labels += [r""]
		hand += [mpl.lines.Line2D([], [], c='w')]

		## r = const
		labels += [r"$r=\textrm{constant} \;\; (dr = l/2)$"]
		alpha = rline_sty['alpha']
		col = np.array(mpl.colors.colorConverter.to_rgb('c'))
		col2 = 0.75*np.array(mpl.colors.colorConverter.to_rgb(plt.cm.Oranges(.1))) + .25
		col = alpha*col + (1.-alpha)*col2
		rline_sty.update(c=col, alpha=1.)
		hand += [mpl.lines.Line2D([], [], **rline_sty)]
		## r = const
		labels += [r"$r=\textrm{constant} \;\; (dr = 2M/2)$"]
		col = np.array(mpl.colors.colorConverter.to_rgb('m'))
		col2 = 0.75*np.array(mpl.colors.colorConverter.to_rgb(plt.cm.Oranges(.1))) + .25
		col = alpha*col + (1.-alpha)*col2
		rline_sty.update(c=col, alpha=1.)
		hand += [mpl.lines.Line2D([], [], **rline_sty)]

		## blank
		labels += [r""]
		hand += [mpl.lines.Line2D([], [], c='w')]

		## r = inf
		labels += [r"$r=\infty \;\; (\textrm{width} \! \propto \! (1 \! + \! 2m/M) )$"]
		rline_inf_sty2 = rline_inf_sty.copy()
		rline_inf_sty2.update(dict(lw=3, ls='dashed', dashes=(6,4)))
		hand += [(mpl.lines.Line2D([], [], **rline_inf_sty ),
		          mpl.lines.Line2D([], [], **rline_inf_sty2) )]

		## ticks
		labels += [r"$d\tau_{\infty}=\textrm{constant}$"]
		tsty = dict(markeredgecolor='0.5', ls='none', marker=(2,0,0), markersize=5)
		hand += [mpl.lines.Line2D([], [], **tsty)]



		## make legend
		legend = plt.legend(hand, labels, loc='center', borderpad=1, labelspacing=.75, fontsize=fs, numpoints=4, frameon=False)
		legend.get_frame().set_linewidth(.5)





	if True:
		## colorbar
		plt.sca(ax2)

		## params
		npoints = 200
		cmin, cmax, alpha = .1, .6, .75

		## custom colors
		cc = np.linspace(cmin, cmax, npoints)
		cols = alpha * plt.cm.Oranges(cc) + (1.-alpha)
		cols[:,-1] = 1. + 0.*cols[:,-1]

		## colormap norm etc
		cmap = matplotlib.colors.ListedColormap(cols)
		norm = mpl.colors.Normalize(vmin=0., vmax=1.)
		ticks = []
		orientation = 'horizontal'

		##
		plt.annotate(s=r"proper density ($\rho/\rho_0$)", xy=(.5,-.35), xycoords='axes fraction', ha='center', va='top', size=fs)
		plt.annotate(s=r"$0$", xy=(0.01,-.22), xycoords='axes fraction', ha='center', va='top', size=fs)
		plt.annotate(s=r"$1$", xy=(1.,-.22), xycoords='axes fraction', ha='center', va='top', size=fs)


		## colorbar
		cb1 = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, ticks=ticks, orientation=orientation)



	## save
	if True:
		path = "temp-figs/demo"
		plt.savefig("%s_legend.png"%(path), dpi=800)



if __name__=="__main__":
	legend()

