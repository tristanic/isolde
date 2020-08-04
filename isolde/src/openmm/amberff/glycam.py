# @Author: Tristan Croll <tic20>
# @Date:   15-Feb-2019
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 29-Jul-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright:2016-2019 Tristan Croll

glycam_suffix_to_ccd_name = {
    'TA':   None,   # alpha-D-talopyranose - not in CCD
    'TB':   None,   # beta-D-talopyranose - not in CCD
    'NA':   'AFD',  # alpha-D-allopyranose
    'NB':   'ALL',  # beta-D-allopyranose
    'EA':   'SHD',  # alpha-D-altropyranose
    'EB':   None,   # beta-D-altropyranse. Not currently in CCD - does appear as part of disaccharides 9QZ and M13
    'FA':   'FCA',  # alpha-D-fucopyranose
    'FB':   'FCB',  # beta-D-fucopyranose
    'LA':   'GLA',  # alpha-D-galactopyranose
    'LB':   'GAL',  # beta-D-galactopyranose
    'GA':   'GLC',  # alpha-D-glucopyranose
    'GB':   'BGC',  # beta-D-glucopyranose
    'KA':   None,   # alpha-D-gulopyranose. Doesn't appear in isolation in CCD
    'KB':   'GL0',  # beta-D-gulopyranose
    'MA':   'MAN',  # alpha-D-mannopyranose
    'MB':   'BMA',  # beta-D-mannopyranose
    'QA':   'G6D',  # alpha-D-quinovopyranose
    'QB':   None,   # beta-D-quinovopyranose. Not in CCD
    'HA':   'XXR',  # alpha-D-rhamnopyranose
    'HB':   None,   # beta-D-rhamnopyranose. Not in CCD
    'tA':   None,   # alpha-L-talopyranose. Not in CCD
    'tB':   None,   # beta-L-talopyranose. Not in CCD
    'nA':   None,   # alpha-L-allopyranose. Not in CCD
    'nB':   'WOO',  # beta-L-allopyranose
    'eA':   'Z6H',  # alpha-L-altropyranose
    'eB':   None,   # beta-L-altropyranose. Not in CCD
    'fA':   'FUC',  # alpha-L-fucopyranose
    'fB':   'FUL',  # beta-L-fucopyranose
    'lA':   'GXL',  # alpha-L-galactopyranose
    'lB':   'GIV',  # beta-L-galactopyranose
    'gA':   None,   # alpha-L-glucopyranose. Not in CCD
    'gB':   None,   # beta-L-glucopyranose. Not in CCD
    'kA':   'GUP',  # alpha-L-gulopyranose
    'kB':   None,   # beta-L-gulopyranose. Not in CCD
    'mA':   None,   # alpha-L-mannopyranose. Not in CCD
    'mB':   None,   # beta-L-mannopyranose. Not in CCD
    'qA':   None,   # alpha-L-quinovopyranose. Not in CCD
    'qB':   None,   # beta-L-quinovopyranose. Not in CCD
    'hA':   'RAM',  # alpha-L-rhamnopyranose.
    'hB':   'RM4',  # beta-L-rhamnopyranose
    'VA':   'A2G',  # N-acetyl-alpha-D-galactosamine
    'VB':   'NGA',  # N-acetyl-beta-D-galactosamine
    'YA':   'NDG',  # N-acetyl-alpha-D-glucosamine
    'YB':   'NAG',  # N-acetyl-beta-D-glucosamine
    'WA':   'BM3',  # N-acetyl-alpha-D-mannosamine
    'WB':   'BM7',  # N-acetyl-beta-D-mannosamine
    'vA':   None,   # N-acetyl-alpha-L-galactosamine. Not in CCD
    'vB':   None,   # N-acetyl-beta-L-galactosamine. Not in CCD
    'yA':   None,   # N-acetyl-alpha-L-glucosamine. Not in CCD
    'yB':   None,   # N-acetyl-beta-L-glucosamine. Not in CCD
    'wA':   None,   # N-acetyl-alpha-L-mannosamine. Not in CCD
    'wB':   None,   # N-acetyl-beta-L-mannosamine. Not in CCD
    'Bc':   None,   # Not even sure what this is. Not clearly described in GLYCAM code
    'bc':   None,
    'YN':   'PA1',  # alpha-D-glucosamine
    'Yn':   'GCS',  # beta-D-glucosamine
    'YNP':  'PA1',  # protonated alpha-D-glucosamine
    'YnP':  'GCS',  # protonated beta-D-glucosamine
    'OA':   'ADA',  # alpha-D-galacturonic acid
    'OB':   'GTR',  # beta-D-galacturonic acid
    'ZB':   'BDP',  # beta-D-glucuronic acid
    'ZA':   'GCU',  # alpha-D-glucuronic acid
    'UA':   None,   # alpha-D-iduronic acid
    'UB':   None,   # beta-D-iduronic acid
    'oA':   None,   # alpha-L-galacturonic acid
    'oB':   None,   # beta-L-galacturonic acid
    'zA':   None,   # alpha-L-glucuronic acid
    'zB':   None,   # beta-L-glucuronic acid
    'uA':   'IDR',  # alpha-L-iduronic acid
    'uB':   None,   # beta-L-iduronic acid
    'AA':  '64K',   # alpha-D-arabinopyranose
    'AB':   None,   # beta-D-arabinopyranose
    'DA':  'LDY',   # alpha-D-lyxopyranose
    'DB':   None,   # beta-D-lyxopyranose
    'RA':   None,   # alpha-D-ribopyranose
    'RB':  'RIP',   # beta-D-ribopyranose
    'XA':  'XYS',   # alpha-D-xylopyranose
    'XB':  'XYP',   # beta-D-xylopyranose
    'aA':  'ARA',   # alpha-L-arabinopyranose
    'aB':  'ARB',   # beta-L-arabinopyranose
    'dA':   None,   # alpha-L-lyxopyranose
    'dB':   None,   # beta-L-lyxopyranose
    'rA':   None,   # alpha-L-ribopyranose
    'rB':  '0MK',   # beta-L-ribopyranose
    'xA':  'HSY',   # alpha-L-xylopyranose
    'xB':  'LXC',   # beta-L-xylopyranose
    'AE':   'ABE',  # alpha-D-abequopyranose
    'Ae':   None,   # beta-D-abequopyranose
    'TV':   'TYV',  # alpha-D-tyvelopyranose
    'Tv':   None,   # beta-D-tyvelopyranose
    'tV':   None,   # alpha-L-tyvelopyranose
    'tv':   None,   # beta-L-tyvelopyranose
    'CA':   None,   # alpha-D-fructopyranose
    'CB':  'BDF',   # beta-D-fructopyranose
    'PA':   None,   # alpha-D-psicopyranose
    'PB':   None,   # beta-D-psicopyranose
    'BA':   None,   # alpha-D-sorbopyranose
    'BB':   None,   # beta-D-sorbopyranose
    'JA':  'T6T',   # alpha-D-tagatopyranose
    'JB':   None,   # beta-D-tagatopyranose
    'cA':   None,   # alpha-L-fructopyranose
    'cB':   None,   # beta-L-fructopyranose
    'pA':   None,   # alpha-L-psicopyranose
    'pB':   None,   # beta-L-psicopyranose
    'bA':   'SOE',  # alpha-L-sorbopyranose
    'bB':   None,   # beta-L-sorbopyranose
    'jA':   None,   # alpha-L-tagatopyranose
    'jB':   None,   # beta-L-tagatopyranose

    'SA':   'SIA',  # alpha-D-sialic acid
    'SB':   'SLB',  # beta-D-sialic acid
    'sA':   None,   # alpha-L-sialic acid
    'sB':   None,   # beta-L-sialic acid
    'GL':   'NGC',  # alpha-D-N-glycolylneuraminic acid
    'gL':   None,   # alpha-L-N-glycloylneuraminic acid

    'AD':   'BXY',  # alpha-D-arabinofuranose
    'AU':   'BXX',  # beta-D-arabinofuranose
    'DD':   None,   # alpha-D-lyxofuranose
    'DU':   None,   # beta-D-lyxofuranose
    'RD':   'RIB',  # alpha-D-ribofuranose
    'RU':   'BDR',  # beta-D-ribofuranose
    'XD':   None,   # alpha-D-xylofuranose
    'XU':   'XYZ',  # beta-D-xylofuranose
    'aD':   'AHR',  # alpha-L-arabinofuranose
    'aU':   'FUB',  # beta-L-arabinofuranose
    'dD':   None,   # alpha-L-lyxofuranose
    'dU':   None,   # beta-L-lyxofuranose
    'rD':   'Z6J',  # alpha-L-ribofuranose
    'rU':   '32O',  # beta-L-ribofuranose
    'xD':   None,   # alpha-L-xylofuranose
    'xU':   None,   # beta-L-xylofuranose
    'CD':   'Z9N',  # alpha-D-fructofuranose
    'CU':   'FRU',  # beta-D-fructofuranose
    'PD':   'PSV',  # alpha-D-psicofuranose
    'PU':   None,   # beta-D-psicofuranose
    'BD':   None,   # alpha-D-sorbofuranose
    'BU':   None,   # beta-D-sorbofuranose
    'JD':   None,   # alpha-D-tagatofuranose
    'JU':   None,   # beta-D-tagatofuranose
    'cD':   None,   # alpha-L-fructofuranose
    'cU':   'LFR',  # beta-L-fructofuranose
    'pD':   None,   # alpha-L-psicofuranose
    'pU':   None,   # beta-L-psicofuranose
    'bD':   None,   # alpha-L-sorbofuranose
    'bU':   None,   # beta-L-sorbofuranose
    'jD':   None,   # alpha-L-tagatofuranose
    'jU':   None,   # beta-L-tagatofuranose
}

residue_name_to_glycam_code = {val: key for key, val in glycam_suffix_to_ccd_name.items() if val is not None}

anchor_name_to_ccd_template = {
    "OLS": "SER_LL_DHG",
    "OLT": "THR_LL_DHG1",
    "NLN": "ASN_LL",
    "COLS": "SER_LEO2_DHG",
    "COLT": "THR_LEO2_DHG1",
    "CNLN": "ASN_LEO2",
    "NOLS": "SER_LSN3_DHG",
    "NOLT": "THR_LSN3_DHG1",
    "NNLN": "ASN_LSN3",
    # "MEX": "74C",
    # "ROH": "OH",
    # "SO3": "SO3",
    "ZOLS": "SER_LFZW_DHG",
    "ZOLT": "THR_LFZW_DHG1",

}


_special_reducing_termini = {
    'SLB':  '2',
}

_anchor_name_map = {
    'ASN'   : 'NLN',
    'SER'   : 'OLS',
    'THR'   : 'OLT',
    'HYP'   : 'OLP',
}

known_sugars = set(residue_name_to_glycam_code.keys())

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

glycam_prefix_to_description = {
    '0':    'non-reducing terminal',
    '2':    'O2-linked',
    '3':    'O3-linked',
    '4':    'O4-linked',
    '6':    'O6-linked',
    'Z':    '2,3-linked',
    'Y':    '2,4-linked',
    'X':    '2,6-linked',
    'W':    '3,4-linked',
    'V':    '3,6-linked',
    'U':    '4,6-linked',
    'T':    '2,3,4-linked',
    'S':    '2,3,6-linked',
    'R':    '2,4,6-linked',
    'Q':    '3,4,6-linked',
    'P':    '2,3,4,6-linked'
}

def find_glycan_template_name_and_link(residue):
    '''
    Work out the template name for a sugar residue, and return this with (if
    applicable) the amino acid residue it's bonded to.
    '''
    bonded_atoms = []
    link_res = None
    r = residue
    neighbors = r.neighbors
    from chimerax.atomic import Residue
    core_name = residue_name_to_glycam_code[r.name]
    reducing_terminal_number = _special_reducing_termini.get(r.name, '1')
    for n in neighbors:
        if n.polymer_type==Residue.PT_AMINO:
            link_res = n
        bonds = r.bonds_between(n)
        if not len(bonds) == 1:
            raise RuntimeError('Residue {}{} has multiple bonds to neighbor {}{}'.format(
                r.chain_id, r.number, n.chain_id, n.number
            ))
        for atom in bonds[0].atoms:
            if atom.residue == residue:
                bonded_atoms.append(atom)
    bonded_atom_numbers = tuple(sorted(
        [int(s) for atom in bonded_atoms for s in atom.name if s.isdigit() and s != reducing_terminal_number]
    ))
    if not len(bonded_atom_numbers):
        bonded_atom_numbers = (0,)
    try:
        prefix = _glycam_prefix[bonded_atom_numbers]
    except:
        raise RuntimeError('Residue {}{} has incorrect bonding: {}'.format(residue.chain_id, residue.number, bonded_atom_numbers))
    return ('GLYCAM_'+prefix+core_name, link_res)

def find_glycan_anchor_name(residue):
    base_name = _anchor_name_map.get(residue.name, None)
    if base_name is None:
        raise RuntimeError('Residue {} {}{} is linked to a glycan, but no '
            'parameters are available for this type of linkage.'.format(
            residue.name, residue.chain_id, residue.number
            ))
    n_term=False
    c_term=False
    N = residue.find_atom('N')
    h_count = 0
    for na in N.neighbors:
        if na.residue != residue:
            break
        if na.element.name=='H':
            h_count += 1
    if h_count == 3:
        n_term=True

    C = residue.find_atom('C')
    o_count = 0
    for ca in C.neighbors:
        if ca.residue != residue:
            break
        if ca.element.name == 'O':
            o_count += 1
    if o_count == 2:
        c_term=True
    if n_term:
        prefix='N'
    elif c_term:
        prefix='C'
    else:
        prefix=''
    return 'GLYCAM_'+prefix+base_name


def is_sugar_template(residue_template):
    '''
    Work through a GLYCAM template to decide if it's a sugar or not. For our
    purposes, it's enough to determine that it contains a 5- or 6-membered ring
    with 4 or 5 carbons and one oxygen.
    '''
    h_atoms = [atom for atom in residue_template.atoms if not atom.name.startswith('H')]
    kept_atoms = []
    current_atom = h_atoms[0]
    def find_ring(atom, traversed_atoms=[]):
        # if len(traversed_atoms):
            # print('Branching from: {}, traversed {}'.format(traversed_atoms[-1].name,
            #     [a.name for a in traversed_atoms]))
        while True:
            # print(atom.name)
            bonds = atom.bonds
            if len(traversed_atoms):
                last_atom = traversed_atoms[-1]
            else:
                last_atom = None
            traversed_atoms.append(atom)
            next_atoms = [a for b in bonds for a in (b.atom1, b.atom2) if not a.name.startswith('H') and a not in (atom, last_atom) ]
            # print("Next atoms: {}".format([a.name for a in next_atoms]))
            if len(set(next_atoms).intersection(set(traversed_atoms))):
                for i, t in enumerate(traversed_atoms):
                    if t in next_atoms:
                        break
                ring_atoms = traversed_atoms[i:]
                # print("Found ring: {}".format([a.name for a in ring_atoms]))
                if len(ring_atoms) in (5,6):
                    c_count = 0
                    o_count = 0
                    for a in ring_atoms:
                        if 'C' in a.name:
                            c_count += 1
                        elif 'O' in a.name:
                            o_count += 1
                    if c_count in (4, 5) and o_count == 1:
                        return True
                else:
                    # avoid getting stuck in an infinite loop
                    next_atoms = [a for a in next_atoms if a not in traversed_atoms]
            if len(next_atoms) == 0:
                return False
            while len(next_atoms) > 1:
                ta = traversed_atoms.copy()
                result = find_ring(next_atoms.pop(-1), ta)
                if result:
                    return result
            atom = next_atoms[0]
    return find_ring(current_atom)
