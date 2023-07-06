Getting Started: Minimal Examples
=========

Basic syntax
######

The absolute minimum syntax to generate a diagram with ``xhorizon`` is::
    
    import xhorizon as xh
    import matplotlib.pyplot as plt

    f = xh.mf.schwarzschild(R=1)
    reg = xh.reg.MAXreg(f)

    reg.rplot()
    plt.show()

.. image:: ../_static/tutorial/001a.png
  :width: 400
  :alt: Diagram.

You can think of:

* ``f`` as a "metric function" object, eg here in the Schwarzschild example f(r)=1-2M/r.
* ``reg`` as a "region" object, as in "maximally extended region with metric function f(r)".

By default the boundaries, horizons, and some constant radius lines were plotted.


The output was not very pretty. To make null rays at 45 degrees you must fix the aspect ratio. You can do this (and make a prettier plot) by adjusting matplotlib options as normal::


    
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


Or you can use ``xh.newfig()`` to automatically "prettify" matplotlib settings and set aspect ratio::

    
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
