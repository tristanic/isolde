from chimerax.isolde.atomic.building.build_utils import next_chain_id

tooltip=('Merge all models with any atoms selected into the first selected model.')


def run_script(session):
    from chimerax.atomic import selected_residues
    from chimerax.isolde.atomic.building.merge import merge_fragment
    from chimerax.geometry import Place
    selres = selected_residues(session)
    us = selres.unique_structures
    if not len(us) > 1:
        from chimerax.core.errors import UserError
        raise UserError('Must have at least two atomic models selected!')
    us = list(sorted(us, key=lambda m: m.id_string))
    target = us[0]
    target.atoms.coords = target.atoms.scene_coords
    target.position=Place()
    import numpy
    with target.triggers.block_trigger('changes'):
        for m in us[1:]:
            for cid in m.residues.unique_chain_ids:
                new_cid = cid
                cres = m.residues[m.residues.chain_ids==cid]
                if cid in target.residues.unique_chain_ids:
                    tres = target.residues[target.residues.chain_ids==cid]
                    if any(numpy.isin(cres.numbers, tres.numbers)):
                        new_cid = next_chain_id(target)
                merge_fragment(target, m.residues[m.residues.chain_ids==cid], chain_id=new_cid, transform=m.position)
