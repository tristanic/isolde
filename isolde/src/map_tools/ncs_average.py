# @Author: Tristan Croll <tic20>
# @Date:   24-May-2020
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tic20
# @Last modified time: 01-Sep-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll

def ncs_average_map(session, asu_map, full_map, reference_chain, ncs_chains):
    from chimerax.match_maker.match import (
        defaults, align
        )
    grid_points = asu_map.grid_points(asu_map.model_transform().inverse())
    from chimerax.atomic import Residues, Atoms
    interpolated_maps = []
    interpolated_maps.append(asu_map.data.matrix())
    original_map_position = full_map.position
    dssp_cache={}
    for ncs in ncs_chains:
        score, s1, s2 = align(session, ncs, reference_chain, defaults['matrix'],
            'nw', defaults['gap_open'], defaults['gap_extend'], dssp_cache)
        ref_residues = []
        ncs_residues = []
        for i, (mr, rr) in enumerate(zip(s1, s2)):
            if mr=='.' or rr=='.':
                continue
            ref_res = s1.residues[s1.gapped_to_ungapped(i)]
            match_res = s2.residues[s2.gapped_to_ungapped(i)]
            if not ref_res or not match_res:
                continue
            ref_residues.append(ref_res)
            ncs_residues.append(match_res)
        ref_residues = Residues(ref_residues)
        ncs_residues = Residues(ncs_residues)

        (ref_pa, ncs_pa) = paired_principal_atoms([ref_residues, ncs_residues])

        from chimerax.geometry import align_points
        tf, rms = align_points(ncs_pa.coords, ref_pa.coords)
        #full_map.position = tf * full_map.position
        from chimerax.core.commands import run
        #run(session, "fitmap #{} in #{}".format(full_map.id_string, asu_map.id_string))
        interpolated_maps.append(
            full_map.interpolated_values(grid_points, point_xform=tf.inverse()).reshape(asu_map.data.size))
        #full_map.position = original_map_position

    asu_map.data.matrix()[:] = sum(interpolated_maps)/len(interpolated_maps)
    asu_map.data.values_changed()

def paired_principal_atoms(aligned_residues):
    from chimerax.atomic import Atoms
    return [Atoms(aa) for aa in zip(*((a1, a2) for a1, a2 in zip(
        aligned_residues[0].principal_atoms, aligned_residues[1].principal_atoms
    ) if a1 is not None and a2 is not None
    ))]
