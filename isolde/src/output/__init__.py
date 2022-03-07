def save_with_glycam_names(session, model, filename):
    '''
    Temporarily renames sugar residues to their GLYCAM names, saves the coordinates to file,
    then renames them back.
    '''
    from chimerax.isolde.openmm.openmm_interface import find_residue_templates

    atom_name_dict = {}
    res_name_dict = {}

    ff = session.isolde.forcefield_mgr['amber14']

    template_dict = find_residue_templates(m.residues, ff)

    from chimerax.isolde.graph import make_graph_from_residue

    for i, template_name in template_dict.items():
        if 'GLYCAM' in template_name:
            tmpl = ff._templates[template_name]
            r = m.residues[i]
            rgraph = make_graph_from_residue(r)
            tgraph = tmpl.graph
            ris, tis, _ = rgraph.maximum_common_subgraph(tgraph)
            for ri, ti in zip(ris, tis):
                a = r.atoms[ri]
                atom_name_dict[a] = a.name
                a.name = tmpl.atoms[ti].name
            res_name_dict[r] = r.name
            r.name = template_name.split('_')[1]
        elif '_' not in template_name:
            r = m.residues[i]
            res_name_dict[r] = r.name
            r.name = template_name
    from chimerax.save_command.cmd import provider_save
    provider_save(session, filename, models=[model])
    for a, name in atom_name_dict.items():
        a.name = name
    for r, name in res_name_dict.items():
        r.name = name