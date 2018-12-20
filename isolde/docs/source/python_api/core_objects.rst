Core Dihedral Classes
=====================

.. toctree::
   :maxdepth: 2

.. contents::
    :local:

General Principles
------------------

Dihedrals are the key fundamental units defining many aspects of a molecule's
structure: backbone "twist" or secondary structure, peptide bond conformation,
rotamer (sidechain) conformation, planarity, chirality, etc. Fast methods for
finding, measuring and tracking them are thus vital to ISOLDE's success.

Each of the groups below are arranged in essentially the same way, so I will
explain the layout using proper dihedrals as an example. There are three key
Python classes:

    * :py:class:`Proper_Dihedral_Mgr` exists as a single instance per session
      and is responsible for creation, deletion and retrieval of C++
      :cpp:class:`Proper_Dihedral` objects. The preferred way to find the
      manager instance is to use
      :func:`session_extensions.get_proper_dihedral_manager`. This will create
      the manager if it doesn't yet exist, or simply retrieve it if it does.
    * :py:class:`Proper_Dihedral` is built on the :py:class:`chimerax.State`
      framework and behaves similarly to e.g. :py:class:`chimerax.Atom`.
    * :py:class:`Proper_Dihedrals` is built on the
      :py:class:`chimerax.Collection` framework and behaves similarly to
      e.g. :py:class:`chimerax.Atoms`.

Some general key points:
    * In general, the C++ :cpp:class:`Proper_Dihedral` objects are only created
      when needed. :py:func:`Proper_Dihedral_Mgr.get_dihedral` by default will
      work through the input list of residues, returning dihedrals that already
      exist and attempting to create those that don't.
    * If any constituent atom in a C++ :cpp:class:`Proper_Dihedral` is deleted,
      the dihedral will be deleted and automatically removed from any
      :py:class:`Proper_Dihedrals` instances. Any corresponding Python
      :py:class:`Proper_Dihedral` instances will become invalid, raising an
      exception if an attempt is made to use them. Similarly,
      :class:`Rotamer` objects will be cleaned up if any constituent
      dihedral is deleted. The exception to this rule is :class:`Rama`, where
      each :class:`Rama` will persist as long as its associated residue has at
      least a CA atom. :class:`Rama` instances with missing dihedrals will
      return :attr:`valid` = False, and will attempt to complete themselves
      on each call to the underlying C++ class.
    * The combination of the above rules means that these classes for the most
      part should "just work", with the user or developer rarely if ever needing
      to worry about object creation/deletion. If, for example, a residue is
      mutated from lysine to arginine, the old rotamer will disappear and the
      next call to :py:func:`Rota_Mgr.get_rotamers` will have the new one in
      the expected position.
    * For best performance, you should try to do most tasks using the plural
      :py:class:`chimerax.Collection` calls, which loop over their constituent
      objects in C++. Looping over individual :py:class:`State` objects is at
      least an order of magnitude slower. So, for example:

        .. code-block:: python

            angles = proper_dihedrals.angles

      ... is very fast, whereas:

        .. code-block:: python

            angles = []
            for d in proper_dihedrals:
                angles.append(d.angle)

      is very slow, since it involves creating and deleting a
      :py:class:`Proper_Dihedral` for every iteration of the loop.

.. automodule:: chimerax.isolde.atomic

Proper Dihedrals
----------------

Proper_Dihedral_Mgr
~~~~~~~~~~~~~~~~~~~
    .. autoclass:: Proper_Dihedral_Mgr
        :members:


Proper_Dihedral
~~~~~~~~~~~~~~~
    .. autoclass:: Proper_Dihedral
        :members:
        :inherited-members:

Proper_Dihedrals
~~~~~~~~~~~~~~~~
    .. autoclass:: Proper_Dihedrals
        :members:
        :inherited-members:

Chiral Centres
--------------

Chiral_Mgr
~~~~~~~~~~
    .. autoclass:: Chiral_Mgr
        :members:

Chiral_Center
~~~~~~~~~~~~~
    .. autoclass:: Chiral_Center
        :members:
        :inherited-members:

Chiral_Centers
~~~~~~~~~~~~~~
    .. autoclass:: Chiral_Centers
        :members:
        :inherited-members:

Ramachandran validation
-----------------------

Rama_Mgr
~~~~~~~~
    .. autoclass:: Rama_Mgr
        :members:

Rama
~~~~
    .. autoclass:: Rama
        :members:

Ramas
~~~~~
    .. autoclass:: Ramas
        :members:
        :inherited-members:

Amino acid rotamers
-------------------

Rota_Mgr
~~~~~~~~
    .. autoclass:: Rota_Mgr
        :members:

Rotamer
~~~~~~~
    .. autoclass:: Rotamer
        :members:

Rotamers
~~~~~~~~
    .. autoclass:: Rotamers
        :members:
        :inherited-members:
