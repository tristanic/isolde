
def _get_template(session, name):
    """Get Chemical Component Dictionary (CCD) entry"""
    from chimerax.core.fetch import fetch_file
    filename = '%s.cif' % name
    url = "http://ligand-expo.rcsb.org/reports/%s/%s/%s.cif" % (name[0], name,
                                                                name)
    try:
        return fetch_file(session, url, 'CCD %s' % name, filename, 'CCD')
    except UserError:
        session.logger.warning(
            "Unable to fetch template for '%s': might be missing bonds"
            % name)
        return None


def residue_info_from_ccd_cif(session, residue):
    template = _get_template(session, residue.name)
    from chimerax.atomic.mmcif import get_mmcif_tables
    res_info_raw = get_mmcif_tables(template, (
        'chem_comp',
        'chem_comp_atom',
        'chem_comp_bond',
        ))

class Residue_Info:
    from chimerax.atomic import Residue
    def __init__(self, raw_info):
        comp_labels, comp_data = raw_info['chem_comp']
        self.name = comp_data[comp_labels.index('id')]
        self.long_name = comp_data[comp_labels.index('name')]
        res_type = comp_data[comp_labels.index('type')].upper()
        if "PEPTIDE LINKING" in res_type:
            self.polymer_type = Residue.PT_AMINO
        elif "NA LINKING" in res_type:
            self.polymer_type = Residue.PT_NUCLEIC
        else:
            self.polymer_type = Residue.PT_NONE
        self.type = res_type

        atom_labels, atom_data = raw_info['chem_comp_atom']
        nlabels = len(atom_labels)
        natoms = len(atom_data)/nlabels
