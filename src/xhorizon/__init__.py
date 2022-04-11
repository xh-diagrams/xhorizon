
"""
Hello! I am the docstring heading the top level __init__.py in xhorizon.
"""

######### FROM SKELETON ###########

from importlib.metadata import version
__version__ = version("xhorizon")

###################################

## external
import numpy as np
import matplotlib.pyplot as plt
plt.style.use("classic")

## main subpackages
from xhorizon import tortoise
from xhorizon import diagram_tools
from xhorizon import shell_junction

## convenience names
from xhorizon.tortoise import metfunc_factory as mf
from xhorizon.diagram_tools import region_factory as reg
from xhorizon.diagram_tools import curvemakers as cm

## alternate sub names
from xhorizon import shell_junction as junc

## more convenience names
from xhorizon.diagram_tools.plotstyle import newfig, pubrc
from xhorizon.diagram_tools.curve_class import curve
from xhorizon.diagram_tools.diagram_classes import diagram, region, block
from xhorizon.tortoise.metfunc_tests import func_test

## more convenience names
from xhorizon.diagram_tools import block_masks as mask
from xhorizon.diagram_tools import block_fill as fill
from xhorizon.shell_junction import reg_corner_masks as cornermask


## evap subpackage needs to be loaded after the above due to circular import issue
from xhorizon import evap

