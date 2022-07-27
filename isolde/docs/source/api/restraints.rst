Restraints and Steering Forces
==============================

.. contents::
    :local:

General Principles
------------------

Virtually all manipulations of the model in ISOLDE are achieved via the use of
custom forces implemented in the OpenMM API. Most of these need to be visualised
in some way - it is generally a bad idea to have a custom restraint imposed if
it's easy to forget it's there!

While not all of the manager classes described below, each is implemented as a
:py:class:`chimerax.Model` subclass for simplicity and easy housekeeping. With
the exception of :py:class:`MDFFMgr` (which is placed below the
:py:class:`chimerax.Volume` holding the map it manages), each manager will
appear as a child model to the main :py:class:`chimerax.AtomicStructure`.

.. automodule:: chimerax.isolde.restraints

General Utilities
-----------------
    .. automethod:: chimerax.isolde.restraints.restrain_torsions_to_template

    .. automethod:: chimerax.isolde.restraints.restrain_atom_distances_to_template

    .. automethod:: chimerax.isolde.restraints.restrain_small_ligands

    .. automethod:: chimerax.isolde.restraints.restrain_secondary_structure

    .. automethod:: chimerax.isolde.restraints.restrain_ca_distances_to_template

Change Tracking
---------------

RestraintChangeTracker
~~~~~~~~~~~~~~~~~~~~~~~~
    .. autoclass:: RestraintChangeTracker
        :members:

Map Steering Forces
-------------------

MDFFMgr
~~~~~~~~
    .. autoclass:: MDFFMgr
        :members:
        :inherited-members:

MDFFAtom
~~~~~~~~~
    .. autoclass:: MDFFAtom
        :members:

MDFFAtoms
~~~~~~~~~~
    .. autoclass:: MDFFAtoms
        :members:


Position Restraints
-------------------

PositionRestraintMgr
~~~~~~~~~~~~~~~~~~~~~~
    .. autoclass:: PositionRestraintMgr
        :members:
        :inherited-members:

PositionRestraint
~~~~~~~~~~~~~~~~~~
    .. autoclass:: PositionRestraint
        :members:

PositionRestraints
~~~~~~~~~~~~~~~~~~~
    .. autoclass:: PositionRestraints
        :members:

Interactive Tugging Forces
--------------------------

NOTE: The classes :py:class:`TuggableAtom` and :py:class:`TuggableAtoms`
are simply wrappers over :py:class:`PositionRestraint` and
:py:class:`PositionRestraints` respectively, with no API changes.

TuggableAtomsMgr
~~~~~~~~~~~~~~~~~~
    .. autoclass:: TuggableAtomsMgr
        :members:
        :inherited-members:

TuggableAtom
~~~~~~~~~~~~~
    See :py:class:`Position_Restraint`

TuggableAtoms
~~~~~~~~~~~~~~
    See :py:class:`PositionRestraints`


Distance Restraints
-------------------

DistanceRestraintMgr
~~~~~~~~~~~~~~~~~~~~~~
    .. autoclass:: DistanceRestraintMgr
        :members:
        :inherited-members:

DistanceRestraint
~~~~~~~~~~~~~~~~~~
    .. autoclass:: DistanceRestraint
        :members:

DistanceRestraints
~~~~~~~~~~~~~~~~~~~
    .. autoclass:: DistanceRestraints

Adaptive Distance Restraints
----------------------------

AdaptiveDistanceRestraintMgr
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    .. autoclass:: AdaptiveDistanceRestraintMgr
        :members:
        :inherited-members:

AdaptiveDistanceRestraint
~~~~~~~~~~~~~~~~~~~~~~~~~~~
    .. autoclass:: AdaptiveDistanceRestraint
        :members:

AdaptiveDistanceRestraints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    .. autoclass:: AdaptiveDistanceRestraints
        :members:

Proper Dihedral Restraints
--------------------------

ProperDihedralRestraintMgr
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    .. autoclass:: ProperDihedralRestraintMgr
        :members:
        :inherited-members:

ProperDihedralRestraint
~~~~~~~~~~~~~~~~~~~~~~~~~
    .. autoclass:: ProperDihedralRestraint
        :members:

ProperDihedralRestraints
~~~~~~~~~~~~~~~~~~~~~~~~~~
    .. autoclass:: ProperDihedralRestraints
        :members:

Adaptive Dihedral Restraints
----------------------------

AdaptiveDihedralRestraintMgr
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    .. autoclass:: AdaptiveDihedralRestraintMgr
        :members:
        :inherited-members:
    
AdaptiveDihedralRestraint
~~~~~~~~~~~~~~~~~~~~~~~~~
    .. autoclass:: AdaptiveDihedralRestraint
        :members:

AdaptiveDihedralRestraints
~~~~~~~~~~~~~~~~~~~~~~~~~~
    .. autoclass:: AdaptiveDihedralRestraints
        :members:

Chirality Restraints
--------------------

ChiralRestraintMgr
~~~~~~~~~~~~~~~~~~~~
    .. autoclass:: ChiralRestraintMgr
        :members:
        :inherited-members:

ChiralRestraint
~~~~~~~~~~~~~~~~
    .. autoclass:: ChiralRestraint
        :members:

ChiralRestraints
~~~~~~~~~~~~~~~~~
    .. autoclass:: ChiralRestraints
        :members:


Rotamer Restraints
------------------

RotamerRestraintMgr
~~~~~~~~~~~~~~~~~~~~~
    .. autoclass:: RotamerRestraintMgr
        :members:
        :inherited-members:

RotamerRestraint
~~~~~~~~~~~~~~~~~
    .. autoclass:: RotamerRestraint
        :members:

RotamerRestraints
~~~~~~~~~~~~~~~~~~
    .. autoclass:: RotamerRestraints
        :members:
