0. Package Capabilities
=======================

This package implements the methods of the paper:

JC Schindler, A Aguirre. *Algorithms for the explicit computation of Penrose diagrams.* Class Quantum Grav 35 105019 (2018). `doi:10.1088/1361-6382/aabce2 <https://doi.org/10.1088/1361-6382/aabce2>`_. [`arxiv:1802.02263 <https://arxiv.org/abs/1802.02263>`_.]



Allowed metrics
###############

It allows one to plot causal diagrams for 3+1 dimensional spacetimes in two classes:

* Any metric of the form

.. math::
   ds^2 = -f(r) \, dt^2 + f(r)^{-1} \, dr^2 + r^2 \, d\Omega^2.


* Metrics that are piecewise-of the above form, glued along radial null shell junctions. Arbitrarily many shells and piecewise regions can be included to approximate smooth dynamical evolution.


Nature of the diagrams
######################

The diagrams are constructed by exactly constructing coordinate transformations from a defining coordinate system (e.g. Schwarzchild :math:`(t,r,\theta,\phi)` or Eddington-Finklestein :math:`(u,v,\theta,\phi)` coordinates) into a system of "diagram" coordinates :math:`(U,V,\theta,\phi)` with the following properties:

* *The diagram is simultaneously a Penrose diagram and an exact coordinate diagram.*

      - The coordinates :math:`(U,V)` cover the entire spacetime in a global compact double-null coordinate patch, so causal cones are bounded by lines of unit slope, without losing conformal information. 
      - Thus internal features such as matter content, observers, and arbitrary spherically symmetric functions, can be unambiguously plotted.


* *Each point in the diagram corresponds to a spherical symmetry surface with a well defined radius.*

      - Each point :math:`(U_0,V_0)` on the diagram corresponds to a spherical symmetry surface :math:`(U_0,V_0,\theta,\phi)`, and has a physical radius. 

      - The radius :math:`r=r(U_0,V_0)` is defined by the area of symmetry spheres.


* *The diagram extends the metric smoothly across arbitrarily many horizons, without coordinate singularity.*

      - Horizons occur wherever :math:`f(r)=0`.


* *All coordinates and parameters are unitless. Unitful values are restored by a choice of an overall length scale.*

Custumizability
###############

The package makes available functions for performing coordinate transformations, creating regions based off metrics, adding curves and contours to regions, and plotting diagrams. One can:

   * Quickly plot **basic diagrams** (by default depicting causal boundaries and lines of constant radius).

   * Create **custom diagrams depicting arbitrary curves and coordinate systems** on the diagram.

   * Calculate the spherical **radius** at each point in the diagram and plot contours of spherically symmetric functions.

   * Plot **internal features** such as densities, curvature scalars, vector fields, trapped surfaces, etc, using the above.

   * Specify **custom metrics** via the metric function :math:`f(r)`.

   * Glue together regions across properly matched **null shell junctions**.


Gallery
#######

You can see some examples of what is possible in the :ref:`gallery`.


