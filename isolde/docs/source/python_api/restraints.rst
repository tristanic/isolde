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

MDFF_Atom
~~~~~~~~~
    .. autoclass:: MDFF_Atom
        :members:

MDFF_Atoms
~~~~~~~~~~
    .. autoclass:: MDFF_Atoms
        :members:


Position Restraints
-------------------

Position_Restraint_Mgr
~~~~~~~~~~~~~~~~~~~~~~
    .. autoclass:: Position_Restraint_Mgr
        :members:

Position_Restraint
~~~~~~~~~~~~~~~~~~
    .. autoclass:: Position_Restraint
        :members:

Position_Restraints
~~~~~~~~~~~~~~~~~~~
    .. autoclass:: Position_Restraints
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

Tuggable_Atom
~~~~~~~~~~~~~
    See :py:class:`Positon_Restraint`

Tuggable_Atoms
~~~~~~~~~~~~~~
    See :py:class:`Position_Restraints`


Distance Restraints
-------------------

Distance_Restraint_Mgr
~~~~~~~~~~~~~~~~~~~~~~
    .. autoclass:: Distance_Restraint_Mgr
        :members:

    .. autoclass:: Distance_Restraint
        :members:

    .. autoclass:: Distance_Restraints

Proper Dihedral Restraints
--------------------------

Proper_Dihedral_Restraint_Mgr
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    .. autoclass:: Proper_Dihedral_Restraint_Mgr
        :members:

    .. autoclass:: Proper_Dihedral_Restraint
        :members:

    .. autoclass:: Proper_Dihedral_Restraints
        :members:

Rotamer Restraints
------------------

Rotamer_Restraint_Mgr
~~~~~~~~~~~~~~~~~~~~~
    .. autoclass:: Rotamer_Restraint_Mgr
        :members:

    .. autoclass:: Rotamer_Restraint
        :members:

    .. autoclass:: Rotamer_Restraints
        :members:
