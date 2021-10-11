Burgers' Playground
===================

This repository is intended to gather different software implementations
of a numerical scheme discretizing the two-dimensional viscid Burgers' equations:

.. math::
    \dfrac{\partial u}{\partial t} + u \, \dfrac{\partial u}{\partial x} + v \, \dfrac{\partial u}{\partial y} & = \varepsilon \left( \dfrac{\partial^2 u}{\partial x^2} + \dfrac{\partial^2 u}{\partial y^2} \right) \\
    \dfrac{\partial v}{\partial t} + u \, \dfrac{\partial v}{\partial x} + v \, \dfrac{\partial v}{\partial y} & = \varepsilon \left( \dfrac{\partial^2 v}{\partial x^2} + \dfrac{\partial^2 v}{\partial y^2} \right) \, .

Here, :math:`x` and :math:`y` are the spatial coordinates, :math:`t` is the time,
:math:`u = u(x, y, t)` and :math:`v = v(x, y, t)` are the velocity components
respectively in the :math:`x`- and :math:`y`-direction, and :math:`\varepsilon` is
the diffusion coefficient. The equations are set over the spatial domain
:math:`\Omega = \left[ 0, 1 \right] \times \left[ 0, 1 \right]` and are completed
with the initial conditions

.. math::
    TODO

and the Dirichlet boundary conditions

.. math::
    TODO.

The equations are numerically solved over the Cartesian grid :math:`\Omega_h`.
The grid features :math:`n_x` equidistant points along the :math:`x`-direction
and :math:`n_y` equidistant points along the :math:`y`-direction.
The equations are discretized in space using finite differences: the terms involving
first-order spatial derivates are treated with a fifth-order upwind scheme,
while the second-order spatial derivates are computed using a fourth-order
centered scheme.
The equations are integrated in time using the three-stage Runge-Kutta method
proposed by Wicker and Skamarock (2002).
