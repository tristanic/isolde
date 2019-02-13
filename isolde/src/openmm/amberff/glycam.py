_residue_name_to_glycam_code = {
    '64K':  'AA',    # alpha-D-arabinopyranose
    'ARA':  'aA',    # alpha-L-arabinopyranose
    'ARB':  'aB',    # beta-L-arabinopyranose
                     # NOTE: beta-D-arabinopyranose doesn't seem to appear in the CCD

    'LDY':  'DA',    # alpha-D-lyxopyranose

    'RIP':  'RB',    # beta-D-ribopyranose
    '0MK':  'rB',    # beta-L-ribopyranose

    'XYP':  'XB',    # beta-D-xylopyranose
    'XYS':  'XA',    # alpha-D-xylopyranose
    'LXC':  'xB',    # beta-L-xylopyranose
    'HSY':  'xA',    # alpha-L-xylopyranose


    'ALL':  'NB',    # beta-D-allopyranose
    'AFD':  'NA',    # alpha-D-allopyranose
    'WOO':  'nB',    # beta-L-allopyranose

    'SHD':  'EA',    # alpha-D-altropyranose

    'GAL':  'LB',    # beta-D-galactopyranose
    'GLA':  'LA',    # alpha-D-galactopyranose
    'GXL':  'lA',    # alpha-L-galactopyranose
    'GIV':  'lB',    # beta-L-galactopyranose

    'BGC':  'GB',    # beta-D-glucopyranose
    'GLC':  'GA',    # alpha-D-glucopyranose

    'GL0':  'KB',    # beta-D-gulopyranose
    'GUP':  'kA',    # alpha-L-gulopyranose

    '4N2':  'IB',    # beta-L-idopyranose

    'BMA':  'MB',    # beta-D-mannopyranose
    'MAN':  'MA',    # alpha-D-mannose

    'BDF':  'CB',    # beta-D-fructopyranose

    'SOE':  'bA',    # alpha-L-sorbopyranose

    'T6T':  'JA',    # alpha-D-tagatopyranose

    'FCA':  'FA',    # alpha-D-fucose
    'FCB':  'FB',    # beta-D-fucose
    'FUC':  'fA',    # alpha-L-fucose
    'FUL':  'fB',    # beta-L-fucose

    'G6D':  'QA',    # alpha-D-quinovose

    'RAM':  'hA',    # alpha-L-rhamnose

    'ADA':  'OA',    # alpha-D-galacturonic acid
    'GTR':  'OB',    # beta-D-galacturonic acid

    'BDP':  'ZB',    # beta-D-glucuronic acid
    'GCU':  'ZA',    # alpha-D-glucuronic acid

    'IDR':  'uA',    # alpha-L-iduronic acid

    'NGA':  'VB',    # N-acetyl-beta-D-galactosamine

    'NAG':  'YB',    # N-acetyl-beta-D-glucosamine

    'SLB':  'SB',    # 5-N-acetyl-beta-D-neuraminic acid

}

known_sugars = set(_residue_name_to_glycam_code.keys())

_glycam_prefix = {
    (0,):       '0',
    (2,):       '2',
    (3,):       '3',
    (4,):       '4',
    (6,):       '6',
    (2,3):      'Z',
    (2,4):      'Y',
    (2,6):      'X',
    (3,4):      'W',
    (3,6):      'V',
    (4,6):      'U',
    (2,3,4):    'T',
    (2,3,6):    'S',
    (2,4,6):    'R',
    (3,4,6):    'Q',
    (2,3,4,6):  'P',

}

missing_in_privateer= ['ALL', 'AFL', 'WOO', 'HSY', 'XYP', 'SHD', 'GIV', '4N2',
    'T6T', 'G6D', ]
# GLA appears twice in list



def find_glycan_template_name(residue):
    bonded_atoms = []
    neighbors = residue.neighbors
    if not len(neighbors):
        return "CCD_"+residue.name
    core_name = _residue_name_to_glycam_code[residue.name]
    for n in neighbors:
        bonds = residue.bonds_between(n)
        assert len(bonds) == 1
        for atom in bonds[0].atoms:
            if atom.residue == residue:
                bonded_atoms.append(atom)
    bonded_atom_numbers = tuple(sorted(
        [int(s) for atom in bonded_atoms for s in atom.name if s.isdigit() and s != '1']
    ))
    if not len(bonded_atom_numbers):
        bonded_atom_numbers = (0,)
    return 'GLYCAM_'+_glycam_prefix[bonded_atom_numbers]+core_name
