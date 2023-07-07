Getting Started: Minimal Examples
=================================

Basic syntax
############

The minimum syntax to generate a diagram with ``xhorizon`` is::
    
    import xhorizon as xh
    import matplotlib.pyplot as plt

    f = xh.mf.schwarzschild(R=1)
    reg = xh.reg.MAXreg(f)

    reg.rplot()
    plt.show()

.. image:: ../_static/tutorial/001a.png
  :width: 400
  :alt: Diagram.

Above you can think of:

* ``f`` as a "metric function" object, here (Schwarzschild example) representing the function f(r)=1-R/r along with its integral, roots, and other information.

* ``reg`` as a "region" object, as in "a maximally extended region with metric function f(r)", consisting of a collection of blocks, curves, and coordinate transformation functions.

By default the boundaries, horizons, and some constant radius lines are included, but the diagrams are completely customizeable. Later in the tutorial we will learn to turn off default features and add custom curves that can be defined in any of several standard coordinate systems. We will also see how to generate various types of regions and use built-in or custom metric functions to define the metric. But for now, back to the basics.


The output above was not very pretty. To make null rays at 45 degrees you must fix the aspect ratio. You can do this (and make a prettier plot) by adjusting matplotlib options in the usual ways::


    
    import xhorizon as xh
    import matplotlib.pyplot as plt

    f = xh.mf.minkowski()
    reg = xh.reg.MAXreg(f)

    plt.figure(1, figsize=(4,4))
    plt.xlim(.5,2.5)
    plt.ylim(-1.2,1.2)
    plt.xticks([1,2])
    plt.yticks([-1,0,1])
    plt.gca().set_aspect(1)

    reg.rplot()
    plt.show()

.. image:: ../_static/tutorial/001b.png
  :width: 400
  :alt: Diagram.


To make this easier you can use ``xh.newfig()`` to automatically "prettify" matplotlib settings and set aspect ratio::

    
    import xhorizon as xh
    import matplotlib.pyplot as plt

    f = xh.mf.reissner_nordstrom(M=0.5, Q=0.3)
    reg = xh.reg.MAXreg(f)

    xh.newfig()
    reg.rplot()
    plt.show()

.. image:: ../_static/tutorial/001c.png
  :width: 400
  :alt: Diagram.

In this Reissner-Nordstrom example, the blue and green default lines are lines of constant radius at two different length scales, corresponding to the scales of the inner and outer horizon.
