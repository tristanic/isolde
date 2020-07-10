# @Author: Tristan Croll <tic20>
# @Date:   15-Feb-2019
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 27-Apr-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright:2016-2019 Tristan Croll



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
    'A2G':  'VA',    # N-acetyl-alpha-D-galactosamine

    'NAG':  'YB',    # N-acetyl-beta-D-glucosamine

    'SLB':  'SB',    # 5-N-acetyl-beta-D-neuraminic acid
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
    core_name = _residue_name_to_glycam_code[r.name]
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
    n_term=True
    c_term=True
    N = residue.find_atom('N')
    for na in N.neighbors:
        if na.residue != residue:
            n_term=False
    C = residue.find_atom('C')
    for ca in C.neighbors:
        if ca.residue != residue:
            c_term = False
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
