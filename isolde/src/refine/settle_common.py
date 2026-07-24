# @Author: Tristan Croll
# @Date:   20-Jul-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 20-Jul-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Shared helpers for ISOLDE's in-simulation fitting engines (``rotafit`` and
``settle_poses``): the parts that are independent of *how* candidate poses are
generated or scored, so both engines call one implementation. The sim-driving
primitives (coord read, 0 K settle, MDFF decoupling, coupling constant) live on
:class:`SimHandler` instead; these are the residue/selection-level utilities.
'''

import numpy as np

SEVERE_OVERLAP = 2.0     # Angstrom: a moved heavy atom closer than this to a fixed
                         #   environment heavy atom is a severe clash (cull the pose).


def severe_clash(pose_coords, atoms, moved_mask, env_coords, cutoff=SEVERE_OVERLAP):
    '''True if any MOVED heavy atom of this pose comes within ``cutoff`` Angstrom of a
    fixed environment heavy atom -- an obvious non-starter to cull before settling.

    Args:
        * pose_coords: ``(N,3)`` candidate coordinates in ``atoms`` order (Angstrom).
        * atoms: the residue/fragment ``Atoms`` the pose describes (for element + mask).
        * moved_mask: bool array over ``atoms`` marking the atoms this pose actually
          moved -- only these are clash-tested (a static backbone can't newly clash).
        * env_coords: ``(M,3)`` fixed environment heavy-atom coordinates (Angstrom).
        * cutoff: clash distance (Angstrom).
    '''
    moved_mask = np.asarray(moved_mask, dtype=bool)
    is_heavy = atoms.element_names != 'H'
    sel = moved_mask & is_heavy
    if not np.any(sel) or not len(env_coords):
        return False
    from chimerax.geometry import find_close_points
    close, _ = find_close_points(pose_coords[sel], env_coords, cutoff)
    return len(close) > 0


def start_sim_on(session, isolde, residues, *, decouple_atoms=None, lambda_decouple=None):
    '''Start an ISOLDE simulation around ``residues`` using ISOLDE's STANDARD sim-start
    selection (its normal padding + soft-shell). The generous buffer gives the local
    environment room to relax around a re-fitted pose -- important for fixes that need the
    neighbours to accommodate (an over-contracted region can leave no such room). Raises
    (via the sim-start path) if a target residue cannot be parameterised.

    If ``decouple_atoms`` is given, that (mobile) selection is softened to
    ``lambda_decouple`` (default ghost-like 0.01) the instant the simulation starts -- BEFORE
    its first minimisation -- so a freshly-rebuilt or clashing group is dropped in as a
    near-ghost and never drives a start-up "jump" that would molest the environment. This is
    exactly the ``decouple`` option of ``isolde sim start`` (same primitive), applied here so
    the caller need not build an atom spec.'''
    from chimerax.atomic import Residues
    from chimerax.core.commands import run
    session.selection.clear()
    Residues(residues).atoms.selected = True
    run(session, 'isolde sim start sel')
    if decouple_atoms is not None and len(decouple_atoms):
        # Reuse the command's primitive (validates subset-of-mobile + needs group soft-core).
        from chimerax.isolde.cmd.cmd import _apply_start_decouple
        _apply_start_decouple(session, isolde, decouple_atoms, lambda_decouple)
