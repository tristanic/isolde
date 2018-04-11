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
the exception of :py:class:`MDFF_Mgr` (which is placed below the
:py:class:`chimerax.Volume` holding the map it manages), each manager will
appear as a child model to the main :py:class:`chimerax.AtomicStructure`.

.. automodule:: chimerax.isolde.restraints

Change Tracking
---------------

Restraint_Change_Tracker
~~~~~~~~~~~~~~~~~~~~~~~~
.. autoclass:: Restraint_Change_Tracker
    :members:

Map Steering Forces
-------------------

    MDFF_Mgr
    ~~~~~~~~
    .. autoclass:: MDFF_Mgr
        :members:

Position Restraints
-------------------

    Position_Restraint_Mgr
    ~~~~~~~~~~~~~~~~~~~~~~
    .. autoclass:: Position_Restraint_Mgr
        :members:

Interactive Tugging Forces
--------------------------

NOTE: The classes :py:class:`Tuggable_Atom` and :py:class:`Tuggable_Atoms`
are simply wrappers over :py:class:`Position_Restraint` and
:py:class:`Position_Restraints` respectively, with no API changes.

    Tuggable_Atoms_Mgr
    ~~~~~~~~~~~~~~~~~~
    .. autoclass:: Tuggable_Atoms_Mgr
        :members:

Distance Restraints
-------------------

    Distance_Restraint_Mgr
    ~~~~~~~~~~~~~~~~~~~~~~
    .. autoclass:: Distance_Restraint_Mgr
        :members:

Proper Dihedral Restraints
--------------------------

    Proper_Dihedral_Restraint_Mgr
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    .. autoclass:: Proper_Dihedral_Restraint_Mgr
        :members:

Rotamer Restraints
------------------

    Rotamer_Restraint_Mgr
    ~~~~~~~~~~~~~~~~~~~~~
    .. autoclass:: Rotamer_Restraint_Mgr
        :members:
