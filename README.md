[ewoms-numerics logo]{ewoms-numerics.svg image-align:right clear: both}

The ewoms-numerics module
=========================

This eWoms module implements the numerical code which is required to
solve a system of coupled non-linear partial differential equations
that emerges in the context of fluid flow simulations in porous
media. In particular, this means spatial and temporal discretizations,
non-linear and linear solvers and the actual flow models.

Installation
============

Like all eWoms modules, ewoms-numerics uses the DUNE build system to
handle the building and installation process. It is thus supposed to
be build using the `dunecontrol` utility. Please confer to the [INSTALL](https://github.com/ewoms/ewoms-common/blob/master/INSTALL.md)
file located in the ewoms-common module for detailed build
instructions.
