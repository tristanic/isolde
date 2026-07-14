# @Author: Tristan Croll
# @Date:   13-Jul-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 13-Jul-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
GAFF<->ff14SB boundary ("seam") parameters for covalent parameterisation.

Under Strategy A (see the design plan and ``covalent.py``) the atoms of a covalent
unit that are chemically unchanged keep their standard **ff14SB** atom types, while
the modified atoms + any ligand get **GAFF2** types from ANTECHAMBER. Wherever a
bond/angle/dihedral spans the two, the base force field has no parameter for that
atom-type combination -- that crossing is the *seam*.

Two mechanisms supply the seam parameters, in priority order:

1. **Curated overrides** (:data:`CURATED_SEAM`) -- hand-vetted values for common
   anchor chemistries, generalising the ``CYScyc`` blocks already shipped in
   ``termods.xml``. Empty to begin with; populated as chemistries are validated.
2. **Analogy re-keying** (:func:`gaff_equivalent`) -- ANTECHAMBER GAFF-types the
   *whole* capped super-residue, so ``parmchk2`` already emits a frcmod covering
   every internal term in GAFF types. At emission we look a seam term up in that
   frcmod using the **GAFF-equivalent** of each ff14SB boundary type, then emit it
   keyed by the real ff14SB type. No second ``parmchk2`` run is needed.

This is a graft, not a fit: analogy-derived seam stiffnesses are approximate and
the curated library is the place to correct the ones that matter.
'''

#: Nearest GAFF2 atom type for each ff14SB (amber14) atom-type *class*, used to
#: look seam parameters up by analogy. Approximate by construction -- correct any
#: consequential case via :data:`CURATED_SEAM`. Classes are those found in
#: ``amberff14SB.xml`` ``<AtomTypes>``.
GAFF_EQUIVALENT = {
    # sp3 carbon
    'CT': 'c3', 'CX': 'c3', '2C': 'c3', '3C': 'c3', 'CI': 'c3', 'CO': 'c3',
    'C8': 'c3', 'C4': 'c3', 'C5': 'c3',
    # carbonyl / sp2 carbon
    'C': 'c', 'CM': 'c2',
    # aromatic / conjugated carbon
    'CA': 'ca', 'CB': 'ca', 'CC': 'cc', 'CD': 'cd', 'CN': 'ca', 'CR': 'cc',
    'CV': 'cc', 'CW': 'cc', 'CK': 'cc', 'CQ': 'cc', 'CP': 'ca', 'CS': 'ca',
    'C*': 'c*',
    # nitrogen
    'N': 'n', 'N2': 'nh', 'N3': 'n3', 'N*': 'na', 'NA': 'na', 'NB': 'nb',
    'NC': 'nb',
    # oxygen
    'O': 'o', 'O2': 'o', 'OH': 'oh', 'OS': 'os',
    # sulfur / phosphorus
    'S': 'ss', 'SH': 'sh', 'P': 'p5',
    # hydrogen
    'H': 'hn', 'H1': 'h1', 'H2': 'h1', 'H4': 'h4', 'H5': 'h5', 'HA': 'ha',
    'HC': 'hc', 'HO': 'ho', 'HP': 'hp', 'HS': 'hs',
}


def gaff_equivalent(amber_type):
    '''GAFF2 type to stand in for an ff14SB ``amber_type`` when looking seam
    parameters up by analogy. Returns the input unchanged if it is already a GAFF
    type (lower-case) or unknown -- so a GAFF<->GAFF seam re-keys to itself.'''
    if amber_type in GAFF_EQUIVALENT:
        return GAFF_EQUIVALENT[amber_type]
    return amber_type


#: Curated seam-parameter overrides, keyed by ``frozenset`` of the two atom types
#: at the seam *bond*. Each value is a dict with any of ``'bonds'``/``'angles'``/
#: ``'dihedrals'`` (and ``'impropers'``), whose entries are OpenMM-ready parameter
#: dicts. **Empty to begin with**; add entries as anchor chemistries are validated
#: in the GUI. The schema mirrors the trailing force blocks of ``CYScyc`` in
#: ``termods.xml`` -- e.g. a curated Cys-thioether seam (Cys ``SG`` type ``S`` to a
#: warhead sp3 carbon ``c3``) would look like::
#:
#:     CURATED_SEAM[frozenset(('S', 'c3'))] = {
#:         'bonds':  [{'type1': 'S',  'type2': 'c3', 'length': 0.181, 'k': 189953.6}],
#:         'angles': [{'type1': 'CX', 'type2': 'S', 'type3': 'c3',
#:                     'angle': 1.726, 'k': 518.816}, ...],
#:         'dihedrals': [{'type1': 'H1', 'type2': 'CX', 'type3': 'S', 'type4': 'c3',
#:                        'periodicity1': 3, 'phase1': 0.0, 'k1': 1.3947}, ...],
#:     }
#:
#: Lengths are nm, angles radians, energies kJ/mol -- OpenMM ffXML units.
CURATED_SEAM = {}


def lookup_seam(type_a, type_b):
    '''Return the curated seam-parameter dict for the bond between ``type_a`` and
    ``type_b`` (order-independent), or ``None`` if none is curated -- in which case
    the caller falls back to analogy re-keying via :func:`gaff_equivalent`.'''
    return CURATED_SEAM.get(frozenset((type_a, type_b)))
