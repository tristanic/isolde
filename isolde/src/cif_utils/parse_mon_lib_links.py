

def get_all_defined_links(restraint_lib):
    from chimerax.mmcif import get_cif_tables
    table_list = get_cif_tables(restraint_lib, ('chem_link_bond',), all_data_blocks=True)
    defined_links = {}
    for entry in table_list:
        t = entry[1][0]
        try:
            fields = t.fields(['link_id','atom_id_1','atom_id_2'])[0]
            link_name, atom1, atom2 = fields
            defined_links[link_name]=(atom1, atom2)
        except ValueError:
            print(fields)
    return defined_links

