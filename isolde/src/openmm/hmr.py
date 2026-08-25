# @Author: Tristan Croll
# @Date:   25-Aug-2026
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tcroll
# @Last modified time: 25-Aug-2026
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2026 Tristan Croll
'''
Hydrogen mass repartitioning (HMR).

Raising each hydrogen's mass toward a target (taking the added mass from the
heavy atom it is bonded to) lowers the X-H vibrational frequency. That lets the
variable integrator take larger steps and, in ISOLDE, lets flexible (un-constrained)
hydrogens run without ringing at a gentler friction.

**The floor (the gotcha this exists for).** A naive "add mass to every H, subtract
from the parent" can bleed a heavy atom dry -- e.g. a light heavy atom carrying
several hydrogens can end up *lighter* than its own hydrogens, which makes the
dynamics worse, not better. So the redistribution is floored **per centre**: a
heavy atom's mass is never taken below the (new) hydrogen mass. When the requested
target would breach that, the centre's total mass is split equally over the heavy
atom and its hydrogens instead (heavy == each H). Total mass is conserved at every
centre.
'''


def repartition_hydrogen_masses(system, bonds, hydrogen_indices, target_h_mass):
    '''
    Apply HMR to ``system`` in place.

    Args:
        system:           the ``openmm.System`` (masses read/written in daltons).
        bonds:            iterable of ``(i, j)`` particle-index pairs (the
                          intramolecular bond graph, same indexing as the system).
        hydrogen_indices: iterable of particle indices that are hydrogens.
        target_h_mass:    desired hydrogen mass in daltons (amu). Values ``<= 0``
                          are a no-op.

    Returns ``(n_hydrogens, n_centres, n_floored)``: how many hydrogens were
    repartitioned, over how many heavy-atom centres, and how many of those centres
    hit the equal-split floor (the requested target was too heavy and was clamped).
    '''
    from openmm import unit
    from collections import defaultdict
    da = unit.dalton
    if target_h_mass is None or target_h_mass <= 0:
        return (0, 0, 0)

    h = set(hydrogen_indices)
    # Map each heavy atom to the hydrogens bonded to it. (A hydrogen bonded to more
    # than one heavy atom is pathological; assign it to the first seen.)
    hs_of = defaultdict(list)
    claimed = set()
    for i, j in bonds:
        i = int(i); j = int(j)
        if i in h and j not in h and i not in claimed:
            hs_of[j].append(i); claimed.add(i)
        elif j in h and i not in h and j not in claimed:
            hs_of[i].append(j); claimed.add(j)

    n_h = n_centres = n_floored = 0
    for heavy, hlist in hs_of.items():
        mX0 = system.getParticleMass(heavy).value_in_unit(da)
        mH0 = [system.getParticleMass(k).value_in_unit(da) for k in hlist]
        if mX0 <= 0:
            continue                                   # already massless (fixed) -- skip
        if target_h_mass <= min(mH0):
            continue                                   # target below current H mass -- nothing to do
        k = len(hlist)
        total = mX0 + sum(mH0)
        equal_split = total / (k + 1)                  # heavy == each H at this value
        mH = min(target_h_mass, equal_split)           # per-centre floor
        floored = (mH < target_h_mass)
        donated = sum(mH - m0 for m0 in mH0)
        mX_new = mX0 - donated
        for kk in hlist:
            system.setParticleMass(kk, mH * da)
        system.setParticleMass(heavy, mX_new * da)
        n_h += k
        n_centres += 1
        if floored:
            n_floored += 1
    return (n_h, n_centres, n_floored)
