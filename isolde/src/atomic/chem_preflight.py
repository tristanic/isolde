# @Author: Tristan Croll
# @Date:   02-Sep-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 02-Sep-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Read-only chemistry preflight for ``isolde parameterise`` and the
Unparameterised-residues widget.

The guiding principle: parameterisation must **never modify** the residue it is
handed, nor silently "repair" its chemistry. Instead it inspects the residue and
either proceeds, warns, or refuses. This module supplies the inspections; the
callers decide whether a finding is a warning (widget) or a hard refusal
(command).

Three independent checks:

* :func:`residue_valence_problems` -- per-atom substituent count vs the atom's
  IDATM type. Two classes are distinguished, because the user's phrasing is
  "substituents out of step with idatm_type in a way that cannot be explained by
  ionisation":
    - **over-coordinated**: more bonds than the (already charge-state-aware)
      IDATM type permits. Ionisation can only add/remove *hydrogens*, and a
      protonated atom is assigned the higher-substituent type, so an atom with
      more bonds than its assigned type allows is a genuine bonding error.
    - **under-coordinated**: fewer bonds than the type expects -- an incomplete
      valence, overwhelmingly "the modelled hydrogens are missing" (which the
      command needs present and correct), occasionally a missing heavy-atom bond.

* :func:`ccd_comparison` -- compares the residue against the CCD component of the
  same name, for **warnings only** (never a rebuild): flag when it looks like the
  CCD molecule but is missing atoms, or suggest a new name when it looks nothing
  like it.
'''


def _heavy(residue):
    return [a for a in residue.atoms if a.element.number != 1]


def _idatm_substituents(atom):
    '''Expected number of substituents for ``atom`` from its IDATM type, or
    ``None`` if the type carries no geometry expectation (unknown type, or a
    bare-element type assigned when perception failed).'''
    from chimerax.atomic import Atom
    info = Atom.idatm_info_map.get(atom.idatm_type)
    if info is None:
        return None
    return info.substituents


def residue_valence_problems(residue):
    '''Inspect one residue's heavy atoms. Returns ``(over, under)`` where each is
    a list of ``(atom, num_bonds, expected)``:

    * ``over``  -- num_bonds > expected  (a bonding pathology ionisation cannot
      explain);
    * ``under`` -- num_bonds < expected  (incomplete valence: usually missing
      hydrogens, sometimes a missing bond).

    Metal atoms and atoms bonded to a metal are skipped: metal coordination is
    not an ordinary covalent bond and is handled by the metal-site pipeline.
    '''
    over, under = [], []
    for a in _heavy(residue):
        if a.element.is_metal:
            continue
        if any(nb.element.is_metal for nb in a.neighbors):
            continue
        expected = _idatm_substituents(a)
        if not expected:            # None or 0
            continue
        nb = a.num_bonds
        if nb > expected:
            over.append((a, nb, expected))
        elif nb < expected:
            under.append((a, nb, expected))
    return over, under


def _fmt_atom(a):
    return '/{}{}{}@{}'.format(a.residue.chain_id, a.residue.name,
                               a.residue.number, a.name)


def format_valence_problems(residue, over, under):
    '''A user-facing message describing the valence problems of one residue, or
    ``None`` if there are none.'''
    if not over and not under:
        return None
    lines = ['Residue {} /{}:{} has chemistry that ISOLDE will not parameterise '
             'as modelled:'.format(residue.name, residue.chain_id, residue.number)]
    for (a, nb, exp) in over:
        lines.append('  - {} is bonded to {} atoms, but its chemical type ({}) '
                     'allows at most {}. This is a bonding error (a spurious or '
                     'duplicated bond, or a mis-identified atom).'
                     .format(_fmt_atom(a), nb, a.idatm_type, exp))
    for (a, nb, exp) in under:
        lines.append('  - {} is bonded to {} atoms, but its chemical type ({}) '
                     'expects {} -- most likely missing hydrogen(s).'
                     .format(_fmt_atom(a), nb, a.idatm_type, exp))
    lines.append('')
    lines.append('Parameterisation only builds parameters; it does not alter your '
                 'model. Complete or correct the chemistry first (for missing '
                 'hydrogens, add and check them with e.g. "addh"; verify '
                 'protonation states), then re-run.')
    return '\n'.join(lines)


def _ccd_heavy_atoms(session, resname, allow_fetch=False):
    '''``{atom_id: ELEMENT}`` (upper-case element symbols) for the heavy atoms of
    the CCD component ``resname``, taken from the **local, offline** ChemComp
    store. Never touches the network unless ``allow_fetch`` is set. Returns
    ``None`` if the component is not available locally.'''
    rec = None
    try:
        from chimerax.chemcomp import lookup as _lookup
        rec = _lookup(session, resname)               # local store, no network
    except Exception:
        rec = None
    if rec is None and allow_fetch:
        try:
            from .rdkit_bridge import ccd_records       # local-first, may fetch
            rec = ccd_records(session, resname)
        except Exception:
            rec = None
    if rec is None:
        return None
    atoms = rec[0]                                      # (atom_id, type_symbol, charge, arom)
    heavy = {}
    for row in atoms:
        aid, sym = row[0], (row[1] or '')
        e = sym.upper()
        if e in ('', 'H', 'D'):
            continue
        heavy[aid] = e
    return heavy or None


def ccd_comparison(session, residue, allow_fetch=False):
    '''Compare ``residue`` against the CCD component of the same name, by
    **heavy-atom element composition** (naming- and protonation-independent, so a
    ligand that merely uses non-CCD atom names is not mistaken for a different
    molecule). Local/offline by default -- the command must never block on the
    network. Returns a dict with a ``verdict``:

    * ``'no_ccd'``        -- no local CCD component of that name (a novel name is
      fine, and a name absent from the local store simply can't be checked);
    * ``'match'``         -- same heavy-atom formula;
    * ``'missing_atoms'`` -- a clean subset of the CCD component (same atom
      naming), i.e. looks like it but with atoms left out;
    * ``'unlike'``        -- little compositional resemblance -> probably a
      different species reusing the code;
    * ``'ambiguous'``     -- neither clearly the same nor clearly different (no
      warning worth making).

    Never modifies anything.
    '''
    ref = _ccd_heavy_atoms(session, residue.name, allow_fetch=allow_fetch)
    if ref is None:
        return {'verdict': 'no_ccd'}
    model = {a.name: a.element.name.upper()
             for a in residue.atoms if a.element.number != 1}
    if not model:
        return {'verdict': 'no_ccd'}
    from collections import Counter
    mc, rc = Counter(model.values()), Counter(ref.values())
    inter = sum((mc & rc).values())
    union = sum((mc | rc).values())
    sim = inter / union if union else 0.0
    missing = sorted(set(ref) - set(model))
    extra = sorted(set(model) - set(ref))
    name_coverage = len(set(model) & set(ref)) / len(ref)
    submultiset = all(mc[e] <= rc.get(e, 0) for e in mc)   # model formula ⊆ CCD formula
    info = {'ccd_name': residue.name, 'missing': missing, 'extra': extra,
            'similarity': sim, 'n_ref_heavy': len(ref), 'n_model_heavy': len(model)}
    if mc == rc:
        info['verdict'] = 'match'
    elif sim < 0.5:
        info['verdict'] = 'unlike'
    elif submultiset and not extra and missing and name_coverage >= 0.5:
        info['verdict'] = 'missing_atoms'
    else:
        info['verdict'] = 'ambiguous'
    return info


def format_ccd_warning(info):
    '''Warning string for a :func:`ccd_comparison` result, or ``None`` when there
    is nothing worth saying (``no_ccd`` / ``match`` / ``ambiguous``).'''
    v = info.get('verdict')
    if v == 'missing_atoms':
        return ('Residue {name} resembles CCD component {name} but is missing '
                'heavy atom(s): {miss} ({nm} of {nr} heavy atoms present). '
                'Parameterising it exactly as modelled; if those atoms belong, '
                'add them first (a separate, upstream task) rather than expecting '
                'parameterisation to.'
                .format(name=info['ccd_name'], miss=', '.join(info['missing']),
                        nm=info['n_model_heavy'], nr=info['n_ref_heavy']))
    if v == 'unlike':
        return ('Residue {name} does not resemble CCD component {name} in heavy-'
                'atom composition (similarity {sim:.0f}%). If this is a genuinely '
                'different chemical species, give it a new residue name so it does '
                'not collide with the standard component of that name.'
                .format(name=info['ccd_name'], sim=100 * info['similarity']))
    return None
