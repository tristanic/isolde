# @Author: Tristan Croll <tic20>
# @Date:   19-Sep-2020
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tic20
# @Last modified time: 19-Sep-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll

from math import radians

def test_compare_torsions(session):
    from chimerax.atomic import AtomicStructure
    m1, m2 = session.models.list(type=AtomicStructure)
    return compare_torsions(session, m1.residues, m2.residues)

def chord_formula(angle_deltas):
    from numpy import cos
    return (1-cos(angle_deltas))/2

def backbone_score(phi_deltas, psi_deltas, chi_deltas):
    import numpy
    combined = numpy.column_stack([chord_formula(d) for d in (phi_deltas, psi_deltas, chi_deltas)])
    return numpy.nanmean(combined, axis=1)

TWISTED_THRESHOLD = radians(30)

def sidechain_buried_score(residue):
    '''
    Defines how "buried" a sidechain is by counting the number of heavy atoms
    from other residues coming within 4A of any heavy atom from the sidechain.
    The returned score is a value ranging from 0 to 1, where 0 indicates no
    contact with other atoms, and 1 indicates 3 or more other atoms per
    sidechain atoms. The score scales linearly with the number of contacting
    atoms in between these values.
    '''
    from chimerax.geometry import find_close_points
    from chimerax.atomic import Residues
    import numpy
    r = residue
    m = r.structure
    other_residues = m.residues.subtract(Residues([r]))
    sidechain_atoms = r.atoms[numpy.logical_not(numpy.in1d(r.atoms.names, ['N', 'C', 'CA', 'O']))]
    if not len(sidechain_atoms):
        return 0
    other_atoms = other_residues.atoms
    cp = find_close_points(sidechain_atoms.coords, other_atoms.coords, 4.0)[1]
    score = (len(cp)/len(sidechain_atoms))/3
    if score > 1:
        score = 1
    return score

def sidechain_score(burial_scores, n_chis, chi1_deltas, chi2_deltas, chi1_cutoff = radians(30)):
    import numpy
    no_chi1 = numpy.isnan(chi1_deltas)
    no_chi2 = numpy.isnan(chi2_deltas)
    scores = numpy.zeros(chi1_deltas.shape)
    # Sidechains with no rotamers do not contribute to the score
    scores[n_chis==0] = numpy.nan
    # If a rotameric sidechain is present in the target but not in the model, it
    # gets the maximum possible score
    scores[numpy.logical_and(n_chis>0, no_chi1)] = 1
    chi1_scores = (1+numpy.cos(chi1_deltas))/2
    chi2_scores = (1+numpy.cos(chi2_deltas))/2
    single_chi_filter = numpy.logical_and(~no_chi1, no_chi2)
    two_chi_filter = numpy.logical_and(~no_chi1, ~no_chi2)

    # print('Scores: {}\nBurial_scores: {}\nChi1_scores: {}'.format(scores, burial_scores, chi1_scores))
    scores[single_chi_filter] = (burial_scores[single_chi_filter] *
        (1-chi1_scores[single_chi_filter]))

    scores[two_chi_filter] = (burial_scores[two_chi_filter] *
        (1-2/3*chi1_scores[two_chi_filter]
          - 1/3 * numpy.exp(-(chi1_deltas[two_chi_filter]/chi1_cutoff)**2)
                * chi2_scores[two_chi_filter]))
    return scores

_torsion_adjustments = {
    'THR': ('chi1', radians(-120)),
    'TRP': ('chi2', radians(180)),
}

def minimum_difference(angle1, angle2):
    from math import pi
    a = angle1-angle2
    return abs(a+pi)%(2*pi)-pi

def compare_torsions(session, model_residues, ref_residues, alignment_cutoff=2):
    '''
    Residue-by-residue comparison of key torsion angles for two models.
    model_residues and ref_residues should be from two models of the same
    protein, each containing a single chain.
    '''
    import numpy
    from math import pi
    from chimerax.atomic import Residues
    from chimerax.isolde import session_extensions as sx

    proper_dihedral_mgr = sx.get_proper_dihedral_mgr(session)
    protein_dihedral_dict = proper_dihedral_mgr.dihedral_dict['residues']['protein']
    backbone_dihedral_names = ('phi', 'psi', 'omega')
    sidechain_dihedral_names = ('chi1', 'chi2', 'chi3', 'chi4')

    from chimerax.isolde.restraints.restraint_utils import (
        sequence_align_all_residues,
        paired_principal_atoms
        )
    rrs, mrs = sequence_align_all_residues(session, [ref_residues], [model_residues])
    rpa, mpa = paired_principal_atoms([rrs, mrs])
    rpa = rpa
    mpa = mpa
    from chimerax.std_commands import align
    ra, ma, _, _, tf = align.align(session, rpa, mpa, move=False, cutoff_distance=alignment_cutoff)
    moved_rpa_coords = tf*rpa.coords

    ca_ca_distances = numpy.linalg.norm(moved_rpa_coords-mpa.coords, axis=1)

    r_res = rpa.residues
    m_res = mpa.residues

    residue_numbers = m_res.numbers
    n_res = len(r_res)

    backbone_mean_error = numpy.empty(n_res)
    sidechain_mean_error = numpy.empty(n_res)
    twists = numpy.zeros(n_res, bool)
    cis_trans_flips = numpy.zeros(n_res, bool)
    target_num_chi = numpy.zeros(n_res, int)

    deltas = {
    'phi': numpy.ones(n_res)*numpy.nan,
    'psi': numpy.ones(n_res)*numpy.nan,
    'omega': numpy.ones(n_res)*numpy.nan,
    'chi1': numpy.ones(n_res)*numpy.nan,
    'chi2': numpy.ones(n_res)*numpy.nan,
    'chi3': numpy.ones(n_res)*numpy.nan,
    'chi4': numpy.ones(n_res)*numpy.nan
    }

    sidechain_burial_scores = numpy.ones(n_res)*numpy.nan

    for i, (tr, dr) in enumerate(zip(r_res, m_res)):
        nchi = protein_dihedral_dict[tr.name]['nchi']
        target_num_chi[i] = nchi
        symm = protein_dihedral_dict[tr.name]['symm']
        t_name = tr.name
        d_name = dr.name
        chi_adjustment = None
        if t_name != d_name:
            if t_name in _torsion_adjustments.keys():
                chi_adjustment = _torsion_adjustments[t_name]
            elif d_name in _torsion_adjustments.keys():
                chi_adjustment = list(_torsion_adjustments[d_name])
                chi_adjustment[1] = -chi_adjustment[1]

        abs_b_mean = 0
        for bd_name in backbone_dihedral_names:
            td = proper_dihedral_mgr.get_dihedral(tr, bd_name)
            dd = proper_dihedral_mgr.get_dihedral(dr, bd_name)
            if bd_name == 'omega' and dd is not None:
                dd_angle = dd.angle
                if abs(dd_angle) > TWISTED_THRESHOLD and abs(dd_angle) < pi-TWISTED_THRESHOLD:
                    twists[i] = True
                if td is not None:
                    if abs(minimum_difference(dd_angle, td.angle)) > pi/2:
                        cis_trans_flips[i] = True
            if td is None or dd is None:
                # No penalty for missing backbone dihedrals (only occur at
                # chain breaks)
                abs_b_mean += 0
            else:
                diff = minimum_difference(td.angle, dd.angle)
                deltas[bd_name][i] = diff
                abs_b_mean += abs(diff)
        abs_b_mean = abs_b_mean/3
        backbone_mean_error[i] = abs_b_mean

        abs_s_mean = 0
        sum_weights = 0

        # If a sidechain has no rotamers, it doesn't contribute to the score
        buried_score = sidechain_buried_score(tr)
        sidechain_burial_scores[i] = buried_score

        if nchi==0:
            sidechain_mean_error[i] = numpy.nan
            continue

        if buried_score == 0:
            sidechain_mean_error[i] = 0
            continue

        for j in range(nchi):
            sd_name =sidechain_dihedral_names[j]
            td = proper_dihedral_mgr.get_dihedral(tr, sd_name)
            dd = proper_dihedral_mgr.get_dihedral(dr, sd_name)

            if td is not None:
                # downweight angles further from the backbone
                weight = 1/2**j
                sum_weights += weight
                if dd is None:
                    # Big penalty for missing sidechain in decoy
                    abs_s_mean += pi * weight
                else:
                    td_angle = td.angle
                    if chi_adjustment is not None:
                        if sd_name == chi_adjustment[0]:
                            td_angle += chi_adjustment[1]
                    dd_angle = dd.angle
                    if symm and j == nchi-1:
                        if td_angle < 0:
                            td_angle += pi
                        if dd_angle < 0:
                            dd_angle += pi
                    diff = minimum_difference(td_angle, dd_angle)
                    deltas[sd_name][i] = diff
                    abs_s_mean += abs(diff) * weight
        if sum_weights > 0:
            abs_s_mean = abs_s_mean/sum_weights
            sidechain_mean_error[i]=abs_s_mean*buried_score
        else:
            sidechain_mean_error[i] = numpy.nan
    backbone_mean_error = numpy.degrees(backbone_mean_error)
    sidechain_mean_error = numpy.degrees(sidechain_mean_error)
    existing_sidechain_mask = numpy.logical_not(numpy.isnan(sidechain_mean_error))
    results = {
        'Target residue numbers': residue_numbers,
        'CA-CA distances': ca_ca_distances,
        'Backbone angle error': backbone_mean_error,
        'Overall backbone angle error': numpy.mean(backbone_mean_error),
        'Weighted chi angle error': sidechain_mean_error,
        'Overall weighted chi angle error': numpy.mean(sidechain_mean_error[existing_sidechain_mask]),
        'Twisted peptide bonds': twists,
        'Peptide cis/trans flips': cis_trans_flips,
        'Aligned residues': (rrs, mrs),
        'Torsion deltas': deltas,
        'Backbone score': backbone_score(deltas['phi'], deltas['psi'], deltas['omega']),
        'Sidechain score': sidechain_score(sidechain_burial_scores, target_num_chi, deltas['chi1'], deltas['chi2'])
        }
    return results
