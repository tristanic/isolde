Simulation Setup and Management
===============================

.. contents::
    :local:

General Principles
------------------
At the heart of ISOLDE's capabilities is the ability to quickly create an
interactive simulation from an arbitrary selection of residues. In order to
ensure the interface and graphical rendering remain smooth regardless of the
size of the construct, once created most simulation tasks are pushed away to
a separate C++ thread, communicating with ISOLDE in a semi-asynchronous manner.
To be more specific, each call to

.. code:: python

    OpenmmThreadHandler.step(n)

spawns a C++ thread tasked with running n simulation steps, during which time an
arbitrary number of graphics updates may occur. The thread will terminate after
the desired number of steps, and the simulation will not move forward further
until the next call to :func:`step`. This prevents the simulation from running
out of control no matter what happens in the GUI.

Simulation Management Classes
-----------------------------
.. automodule:: chimerax.isolde.openmm

SimParams
~~~~~~~~~
    .. autoclass:: SimParams
        :members:

SimConstruct
~~~~~~~~~~~~~
    .. autoclass:: SimConstruct
        :members:

SimManager
~~~~~~~~~~~
    .. autoclass:: SimManager
        :members:

SimHandler
~~~~~~~~~~~
    .. autoclass:: SimHandler
        :members:

OpenmmThreadHandler
~~~~~~~~~~~~~~~~~~~~~
    .. autoclass:: OpenmmThreadHandler
        :members:

Custom Forces
-------------
The actual custom OpenMM force classes used to implement MDFF and the various
interactive restraints are listed below.

.. automodule:: chimerax.isolde.openmm.custom_forces

MDFF Potentials
~~~~~~~~~~~~~~~
    .. autoclass:: CubicInterpMapForce
        :members:

    .. autoclass:: CubicInterpMapForce_Old
        :members:

    .. autoclass:: LinearInterpMapForce
        :members:



Distance Restraints
~~~~~~~~~~~~~~~~~~~
    .. autoclass:: AdaptiveDistanceRestraintForce
        :members:

    .. autoclass:: TopOutBondForce
        :members:

Position Restraints
~~~~~~~~~~~~~~~~~~~
    .. autoclass:: TopOutRestraintForce
        :members:

Dihedral Restraints
~~~~~~~~~~~~~~~~~~~
    .. autoclass:: FlatBottomTorsionRestraintForce
        :members:

Adaptive Dihedral Restraints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    .. autoclass:: TopOutTorsionForce
        :members:

CMAP Correction terms
~~~~~~~~~~~~~~~~~~~~~
    .. autoclass:: AmberCMAPForce
        :members:

Implicit Solvent
~~~~~~~~~~~~~~~~
    .. autoclass:: GBSAForce
        :members:
