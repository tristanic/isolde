
tooltip=('Remove all altlocs in the currently selected model(s)')

def run_script(session):
    from chimerax.atomic import selected_residues
    models = selected_residues(session).unique_structures
    for m in models:
        atoms_with_alt_locs = m.atoms[m.atoms.num_alt_locs>0]
        m.delete_alt_locs()
        atoms_with_alt_locs.occupancies = 1
    session.logger.info((
        'Removed all altlocs in #{}.'
        'Occupancies for all affected atoms have been reset to 1.'
    ).format(', '.join(m.id_string for m in models)))