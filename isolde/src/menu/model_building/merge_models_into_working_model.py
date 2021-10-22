from chimerax.isolde.atomic.building.build_utils import next_chain_id

tooltip=("Merge all models with any atoms selected into ISOLDE's current working model.")


def run_script(session):
    from chimerax.atomic import selected_residues
    from chimerax.core.errors import UserError
    from chimerax.atomic import next_chain_id
    #from chimerax.isolde.atomic.building.merge import merge_fragment
    #from chimerax.geometry import Place
    session.logger.info('ISOLDE: merge models')
    if not hasattr(session, 'isolde'):
        raise UserError('ISOLDE must be running')
    target = session.isolde.selected_model
    if target is None:
        raise UserError('ISOLDE has no model initialised!')
    selres = selected_residues(session)
    us = selres.unique_structures
    others = [s for s in us if s != target]
    if not len(others):
        from chimerax.core.errors import UserError
        raise UserError('Must have at least two atomic models selected!')
    session.logger.info(f'Merging models {",".join([f"#{m.id_string}" for m in others])} into #{target.id_string}.')
    seen_ids = set(target.residues.unique_chain_ids)

    for s in others:
        chain_id_mapping = {}
        chain_ids = sorted(s.residues.unique_chain_ids)
        for cid in chain_ids:
            if cid in seen_ids:
                new_id = next_chain_id(cid)
                while new_id in seen_ids or new_id in chain_ids:
                    new_id = next_chain_id(new_id)
                session.logger.info(f"Remapping chain ID {cid} in #{s.id_string} to {new_id}")
                chain_id_mapping[cid] = new_id
                seen_ids.add(new_id)
            else:
                seen_ids.add(cid)
        target.combine(s, chain_id_mapping, target.scene_position)