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

    * :py:class:`ProperDihedralMgr` exists as a single instance per session
      and is responsible for creation, deletion and retrieval of C++
      :cpp:class:`ProperDihedral` objects. The preferred way to find the
      manager instance is to use
      :func:`session_extensions.get_proper_dihedral_manager`. This will create
      the manager if it doesn't yet exist, or simply retrieve it if it does.
    * :py:class:`ProperDihedral` is built on the :py:class:`chimerax.State`
      framework and behaves similarly to e.g. :py:class:`chimerax.Atom`.
    * :py:class:`ProperDihedrals` is built on the
      :py:class:`chimerax.Collection` framework and behaves similarly to
      e.g. :py:class:`chimerax.Atoms`.

Some general key points:
    * In general, the C++ :cpp:class:`ProperDihedral` objects are only created
      when needed. :py:func:`ProperDihedralMgr.get_dihedral` by default will
      work through the input list of residues, returning dihedrals that already
      exist and attempting to create those that don't.
    * If any constituent atom in a C++ :cpp:class:`ProperDihedral` is deleted,
      the dihedral will be deleted and automatically removed from any
      :py:class:`ProperDihedrals` instances. Any corresponding Python
      :py:class:`ProperDihedral` instances will become invalid, raising an
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
      next call to :py:func:`RotaMgr.get_rotamers` will have the new one in
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
      :py:class:`ProperDihedral` for every iteration of the loop.

.. automodule:: chimerax.isolde.atomic

Proper Dihedrals
----------------

ProperDihedralMgr
~~~~~~~~~~~~~~~~~~~
    .. autoclass:: ProperDihedralMgr
        :members:


ProperDihedral
~~~~~~~~~~~~~~~
    .. autoclass:: ProperDihedral
        :members:
        :inherited-members:

ProperDihedrals
~~~~~~~~~~~~~~~~
    .. autoclass:: ProperDihedrals
        :members:
        :inherited-members:

Chiral Centres
--------------

ChiralMgr
~~~~~~~~~~
    .. autoclass:: ChiralMgr
        :members:

ChiralCenter
~~~~~~~~~~~~~
    .. autoclass:: ChiralCenter
        :members:
        :inherited-members:

ChiralCenters
~~~~~~~~~~~~~~
    .. autoclass:: ChiralCenters
        :members:
        :inherited-members:

Ramachandran validation
-----------------------

RamaMgr
~~~~~~~~
    .. autoclass:: RamaMgr
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

RotaMgr
~~~~~~~~
    .. autoclass:: RotaMgr
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
