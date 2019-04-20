.. pyatmos documentation master file, created by
   sphinx-quickstart2 on Tue Aug 14 13:26:34 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pyatmos's documentation for Master!
==============================================
The pyatmos software is a atmospheric properties library that produces reasonable
atmospheric properties up to 500,000 feet.  This goes far beyond the 1976 standard
atmosphere's limit of ~80,000 feet and is suitable for reentery problems.

Additionally, mixed units are supported.  

Altitudes may be input in:
 - m, meters
 - ft, feet
 - kft, kilofeet

Velocities are in units of:
 - m/s
 - feet/s
 - knots

Densities are in units of:
 - kg/m^3
 - slug/ft^3
 - slinch/in^3

The following properties may be calculated:
 - pressure
 - temperature
 - density
 - dynamic pressure
 - speed of sound
 - equivalent airspeed
 - Reynolds number/unit length
 - dynamic viscosity
 - kinematic viscosity

Additionally, there is some support for vectorization with numpy.

.. toctree::

   reference/pyatmos.utils
   reference/pyatmos.utils.atmosphere
   reference/pyatmos.utils.atmosphere_vectorized
   reference/pyatmos.utils.sweep
..   manual/index


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

