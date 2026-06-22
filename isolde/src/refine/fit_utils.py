# @Author: Tristan Croll
# @Date:   21-Jun-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 21-Jun-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Shared scoring utilities for pose fitting (sidechain backrub, rebuild torsion
optimisation, ...). Two cheap, coordinate-only terms:

* a clash term -- the summed overlap distance of clashing atom pairs, via
  ``chimerax.clashes.find_clashes`` (lower is better);
* a density term -- the element-number-weighted mean map value at the atom
  centres, normalised by the map sigma (higher is better).

These were originally embedded in ``backrub_rotamer.py``; they live here now so
both the Backrub rotamer fitter and the template-rebuild torsion optimiser score
poses identically.
'''


def clash_score(session, atoms, surrounds):
    '''
    Sum of overlap distances for atom pairs in ``atoms`` that clash with
    ``surrounds``. Lower is better; 0.0 means clash-free.
    '''
    from chimerax.clashes import clashes
    clash_dict = clashes.find_clashes(session, atoms, restrict=surrounds)
    # Make list of unique clashing atom pairs and their distances
    seen = set()
    clash_sum = 0
    for a1, partners in clash_dict.items():
        seen.add(a1)
        for a2, dist in partners.items():
            if a2 in seen:
                continue
            clash_sum += dist
    return clash_sum


def near_atoms(atoms, model, distance):
    '''
    Atoms of ``model`` (excluding ``atoms`` themselves) within ``distance`` of
    any atom in ``atoms`` -- the candidate clash partners for a moving group.
    '''
    from chimerax.geometry import find_close_points
    all_atoms = model.atoms
    other_atoms = all_atoms.subtract(atoms)
    other_coords = other_atoms.coords
    test_coords = atoms.coords
    a, _ = find_close_points(other_coords, test_coords, distance)
    return other_atoms.filter(a)


def density_score(volume, atoms, weights=None):
    '''
    Element-number-weighted mean map value at the atom centres, normalised by the
    map sigma. Higher is better (atoms sitting in strong density).

    Returns ``None`` if any atom falls outside the map box -- callers should then
    fall back to a clash-only score rather than treat the pose as scoreless. (A
    rebuild often runs before the map box has been recentred on the residue.)
    '''
    import numpy
    if weights is None:
        weights = atoms.elements.numbers
    dvals, outside = volume.interpolated_values(atoms.coords, out_of_bounds_list=True)
    if len(outside):
        return None
    wsum = numpy.sum(weights)
    if wsum == 0:
        return None
    return float(numpy.sum(dvals * weights) / wsum / volume.sigma)


def active_mdff_map(model):
    '''
    Find the model's preferred density map for fitting and return
    ``(volume, weight)``, or ``(None, None)`` if no usable map exists.

    Preference order:
      1. an enabled MDFF manager's volume (the map ISOLDE is actively fitting to;
         a live ``MDFF potential`` map wins as it already excludes the free set);
      2. failing that -- e.g. before any simulation has been set up, which is the
         usual state during a template rebuild -- the model's best crystallographic
         map straight from the Clipper map manager, excluding difference maps and
         live MDFF-potential maps (which only exist during simulation).
    '''
    try:
        from chimerax.clipper.symmetry import get_map_mgr
    except ImportError:
        return None, None
    mgr = get_map_mgr(model)
    if mgr is None:
        return None, None
    from chimerax.isolde.session_extensions import get_mdff_mgr
    preferred = None
    fallback = None
    for v in mgr.all_maps:
        mdff_mgr = get_mdff_mgr(model, v)
        if mdff_mgr is None or not getattr(mdff_mgr, 'enabled', False):
            continue
        if 'MDFF potential' in v.name:
            preferred = mdff_mgr
            break
        if fallback is None:
            fallback = mdff_mgr
    chosen = preferred or fallback
    if chosen is not None:
        return chosen.volume, chosen.global_k
    # No enabled MDFF manager: fall back to a plain crystallographic map for
    # scoring. Skip difference maps (mFo-DFc) and the live MDFF-potential maps.
    candidates = [
        v for v in mgr.all_maps
        if not getattr(v, 'is_difference_map', False) and 'MDFF potential' not in v.name
    ]
    if candidates:
        return candidates[0], 1.0
    return None, None


def pose_score(session, moving_atoms, surrounds, volume, map_weight=1.0, clash_weight=1.0):
    '''
    Combined fit score for a moving group: ``clash_weight * clash`` minus
    ``map_weight * density``, so **lower is better**. When ``volume`` is None (or
    the group is currently outside the map box) only the clash term contributes.
    '''
    score = clash_weight * clash_score(session, moving_atoms, surrounds)
    if volume is not None:
        d = density_score(volume, moving_atoms)
        if d is not None:
            score -= map_weight * d
    return score
