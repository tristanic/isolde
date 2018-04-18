*n*-Dimensional Interpolation
=============================

.. contents::
    :local:

General Principles
------------------

Many structure validation tasks boil down to checking an *n*-dimensional set of
values (e.g. the *chi* dihedral angles for a sidechain) against a map of
probability values based on what we know from very high-resolution structures.
For a given set of (x1, x2, ... xn) coordinates, the resulting P-value must
be calculated via linear interpolation from the 2^n corners of the n-orthotope
surrounding the point in the map.

Providing rotamer and Ramachandran validation in real time for a 1,000-residue
simulation requires on the order of 4-5,000 dihedral angle measurements, 1,000
2D interpolations and 1,000 1-4D interpolations (spread over 25 different maps)
for every coordinate update. Maintaining a high graphics framereate in the face
of 20+ coordinate updates per second therefore requires this pipeline to be very
highly optimised. The :py:class:`RegularGridInterpolator` described below may be
used as-is as a Python class, taking approximately :math:`(20 + 0.13*m) \mu s`
to perform *m* 3D interpolations on a single map. While this is good enough for
many situations, the Python class is not used by ISOLDE's core :class:`Rama` and
:class:`Rotamer` classes. Instead the underlying C++ class is used, allowing all
the necessary :cpp:class:`RegularGridInterpolator` instances to be stored in an
efficient C++ mapping so that the entire validation task is handled in a single
Python call.

.. automodule:: chimerax.isolde.interpolation

*n*-dimensional Linear Interpolation
------------------------------------

RegularGridInterpolator
~~~~~~~~~~~~~~~~~~~~~~~
    .. autoclass:: RegularGridInterpolator
        :members:

Testing
~~~~~~~
    test_interpolator
