# @Author: Tristan Croll
# @Date:   24-Aug-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 24-Aug-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Extract the primitives garnet's featuriser needs, in-process, from a ChimeraX
model.

This mirrors ``garnet-isolde/validation/md_harness/cx_extract.py`` (which runs
as a separate ChimeraX subprocess in the reference harness) but returns the
primitives directly, since ISOLDE already *is* a ChimeraX session. ChimeraX is
the topology + idatm source (the inference-native path).

Formal charges are deliberately **not** taken from ChimeraX/ISOLDE here:
per-residue perception mis-charges polymer backbones. They are inferred from
idatm types + connectivity by ``garnet_core.idatm_charges`` (in
:mod:`params_cache`). ``estimate_net_charge`` is emitted only as an independent
per-structure cross-check.
'''


def extract_primitives(structure):
    '''
    Pull garnet's input primitives from a ChimeraX :class:`AtomicStructure`.

    Returns ``(primitives, atoms)`` where ``primitives`` is the dict garnet's
    ingest path consumes, and ``atoms`` is the ordered Python list of ChimeraX
    :class:`Atom` objects defining the 0-based index space of every field
    (``bonds`` and ``aromatic_rings`` reference it). Callers key parameters back
    to the model by that list.

    The index space is ``list(structure.atoms)`` — the same convention as
    ``cx_extract`` — so ``primitives[...][i]`` is the value for ``atoms[i]``.
    '''
    m = structure
    atoms = list(m.atoms)
    idx = {a: i for i, a in enumerate(atoms)}

    try:
        from chimerax.add_charge import estimate_net_charge
        est = estimate_net_charge(m.atoms)
    except Exception as e:                        # cross-check only; never fatal
        est = 'FAIL:{}'.format(repr(e)[:60])

    # Aromatic rings, in the same 0-based index space as ``bonds``. These unlock
    # the delocalised charge modes of ``infer_formal_charges`` (azolate -1,
    # aromatic-cation collapse, protonated amidinium, porphyrin -2), none of
    # which has a per-atom signal. Absent/failed -> [] -> those modes are simply
    # skipped, which is byte-identical to the historical ring-free behaviour.
    try:
        arom = [[idx[a] for a in r.atoms if a in idx]
                for r in m.rings() if getattr(r, 'aromatic', False)]
        arom = [r for r in arom if len(r) >= 3]
    except Exception:                             # charge-mode enabler, not required
        arom = []

    primitives = {
        'atomic_numbers': [a.element.number for a in atoms],
        'idatm_types': [a.idatm_type for a in atoms],
        'coords': [[float(c) for c in a.coord] for a in atoms],          # Angstrom
        'names': [a.name for a in atoms],
        'residues': [[a.residue.name, int(a.residue.number), str(a.residue.chain_id)]
                     for a in atoms],
        'bonds': [[idx[b.atoms[0]], idx[b.atoms[1]]] for b in m.bonds],
        'aromatic_rings': arom,
        'estimate_net_charge': est,
    }
    return primitives, atoms
