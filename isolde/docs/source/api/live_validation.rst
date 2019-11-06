Live Structure Validation
=========================

.. contents::
    :local:

General Principles
------------------

One of the core design principles underlying ISOLDE is the need for real-time,
continuous and visually-clear feedback. When you are working on a structure with
thousands of residues, having to stop and generate a table or plot to tell you
how you did, then navigate back to each outlier and attempt to fix it before
repeating the process is simply too slow, and leads people to give up while
their structure is still far from perfect. The aim in ISOLDE is to instead
continually show you how you're going *right now* by directly marking up
problems where and when they occur, telling you directly whether your
manipulations are making things better or worse.

.. automodule:: chimerax.isolde.validation

    Rama_Annotator
    ~~~~~~~~~~~~~~
    .. autoclass:: Rama_Annotator
        :members:

    Rotamer_Annotator
    ~~~~~~~~~~~~~~~~~
    .. autoclass:: Rotamer_Annotator
        :members:
