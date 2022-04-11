
"""
Minimal working example to generate a Penrose diagram with xhorizon.
"""

# Imports.
import numpy as np
import matplotlib.pyplot as plt
import xhorizon as xh

# Create a "metfunc" (metric function) object, which stores info about f(r) in the metric.
f = xh.mf.schwarzschild(R=1.)

# Or try commenting out the above line and replacing it with one of these:
# f = xh.mf.minkowski()
# f = xh.mf.reissner_nordstrom(M=0.5, Q=0.3)

# Create a maximally extended "region" object from the metric function, with auto radius lines turned on.
reg = xh.reg.MAXreg(f, rlines=True)

# Tell the region object to plot itself.
reg.rplot()

# Set figure aspect ratio to one, so that null rays are at 45 degrees.
plt.gca().set_aspect(1.)

# Save the plot.
plt.savefig("out.png", dpi=300)

# Uncomment to additionally show the plot in interactive window.
#plt.show()

