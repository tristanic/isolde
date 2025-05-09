# @Author: Tristan Croll <tic20>
# @Date:   11-Jun-2019
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tic20
# @Last modified time: 04-Jan-2021
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll


_valid_chain_id_chars = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890'
def next_chain_id(model, prefix=''):
    current_chain_ids = model.residues.unique_chain_ids
    for id_char in _valid_chain_id_chars:
        if prefix+id_char not in current_chain_ids:
            return prefix+id_char
    if prefix=='':
        return next_chain_id(model, prefix='A')
    index = _valid_chain_id_chars.find(prefix)
    return next_chain_id(model, prefix= _valid_chain_id_chars[index+1])

def add_hydrogen_to_atom(atom, coord, name = None):
    '''
    Add a single hydrogen atom to the given heavy atom, at the given coordinate.
    Simple-minded tool, taking no notice of chemistry or geometry.
    '''
    from chimerax.atomic import Element
    r = atom.residue
    s = atom.structure
    if name is not None:
        if name in r.atoms.names:
            raise TypeError('This atom name is already taken!')
    else:
        existing_names = [a.name for a in atom.neighbors if a.element.name == 'H']
        if len(existing_names):
            last_digits = [int(n[-1]) for n in existing_names if n[-1].isdigit()]
            if len(last_digits):
                num = max(last_digits)+1
            else:
                num=1
        else:
            num=1

        name = 'H'+atom.name[1:]+str(num)

    from chimerax.atomic.struct_edit import add_atom
    na=add_atom(name, Element.get_element('H'),r,coord,bonded_to=atom)
    return na

def make_fmet_from_met(session, residue):
    err_string = "make_fmet_from_met() is only applicable to N-terminal methionine residues!"
    if residue.name != 'MET':
        raise TypeError(err_string)
    n = residue.find_atom('N')
    if n is None:
        raise TypeError('Missing N atom! Is your residue complete?')
    for nb in n.neighbors:
        if nb.residue != residue:
            raise TypeError(err_string)
        if nb.element.name == 'H':
            nb.delete()
    from chimerax.mmcif import find_template_residue
    import numpy
    fme = find_template_residue(session, 'FME')
    ref_atoms = ('C','CA','N')
    r_coords = numpy.array([residue.find_atom(a).coord for a in ref_atoms])
    t_coords = numpy.array([fme.find_atom(a).coord for a in ref_atoms])
    from chimerax.geometry import align_points
    tf, rms = align_points(t_coords, r_coords)
    from chimerax.atomic.struct_edit import add_atom
    atom_pairs = (('N','H'),('N','CN'),('CN','O1'),('CN','HCN'))
    for (bname, aname) in atom_pairs:
        ta = fme.find_atom(aname)
        add_atom(aname, ta.element, residue, tf*ta.coord, bonded_to=residue.find_atom(bname))
    residue.name='FME'    

def fix_amino_acid_protonation_state(residue):
    if residue.name not in ('GLU', 'ASP'):
        raise TypeError('This method is only applicable to GLU and ASP residues!')
    if residue.name == 'GLU':
        a = residue.find_atom('HE2')
        if a is not None:
            a.delete()
    elif residue.name == 'ASP':
        a = residue.find_atom('HD2')
        if a is not None:
            a.delete()


def set_his_protonation_state(residue, position='ND'):
    '''
    Add/remove hydrogens to a histidine residue to chane its protonation state.

    * Args:

        - residue: a :class:`Residue` object
        - position: one of 'ND', 'NE' or 'both'
    '''
    if residue.name != 'HIS':
        raise TypeError('This method is only applicable to histidine residues!')
    from chimerax.atomic.struct_edit import add_dihedral_atom
    he2 = residue.find_atom('HE2')
    hd1 = residue.find_atom('HD1')
    position = position.upper()
    if position == 'ND' and hd1 and not he2:
        return
    if position=='NE' and he2 and not hd1:
        return
    if position=='BOTH' and he2 and hd1:
        return
    if he2:
        he2.delete()
    if hd1:
        hd1.delete()
    bond_length = 1.01
    angle = 125.0
    dihedral = 180.0
    hd1_dihedral_atoms = [residue.find_atom(name) for name in ('ND1', 'CE1', 'NE2')]
    if not all(hd1_dihedral_atoms):
        raise TypeError('This residue is missing at least one necessary heavy atom!')
    if position in ('ND', 'BOTH'):
        atom = add_dihedral_atom('HD1', 'H', *hd1_dihedral_atoms, bond_length, angle, dihedral, bonded=True)
        nd1 = hd1_dihedral_atoms[0]
        atom.bfactor = nd1.bfactor
        atom.occupancy = nd1.occupancy
    if position in ('NE', 'BOTH'):
        atom = add_dihedral_atom('HE2', 'H', *hd1_dihedral_atoms[::-1], bond_length, angle, dihedral, bonded=True)
        ne2 = hd1_dihedral_atoms[-1]
        atom.bfactor = ne2.bfactor
        atom.occupancy = ne2.occupancy

def add_amino_acid_residue(model, resname, prev_res=None, next_res=None,
        chain_id=None, number=None, center=None, insertion_code=' ', add_b_factor=0,
        occupancy=1, phi=-135, psi=135):
    session = model.session
    if (not chain_id or not number or center is None) and (not prev_res and not next_res):
        raise TypeError('If no anchor residues are specified, chain ID, '
            'number and center must be provided!')
    if prev_res and next_res:
        raise TypeError('Cannot specify both previous and next residues!')
    other_atom = None
    insertion_point = None
    import numpy
    if prev_res:
        ref_res = prev_res
        b_factor = prev_res.atoms[numpy.in1d(prev_res.atoms.names, ('N','CA','C','O','CB'))].bfactors.mean()+add_b_factor
        pri = model.residues.index(prev_res)
        if pri > 0 and pri < len(model.residues)-1:
            insertion_point = model.residues[pri+1]
        catom = prev_res.find_atom('C')
        for n in catom.neighbors:
            if n.residue != prev_res:
                raise TypeError('This residue already has another bonded to its '
                    'C terminus!')
        if chain_id is not None:
            session.logger.warning('AddAA: chain_id argument is ignored when adding to an existing residue.')
        chain_id = prev_res.chain_id
        oxt = prev_res.find_atom('OXT')
        if oxt is not None:
            oxt.delete()
    elif next_res:
        ref_res = next_res
        b_factor = next_res.atoms[numpy.in1d(next_res.atoms.names, ('N','CA','C','O','CB'))].bfactors.mean()+add_b_factor
        insertion_point = next_res
        natom = next_res.find_atom('N')
        for n in natom.neighbors:
            if n.residue != next_res:
                raise TypeError('This residue already has another bonded to its '
                    'N terminus!')
        if chain_id is not None:
            session.logger.warning('AddAA: chain_id argument is ignored when adding to an existing residue.')
        chain_id = next_res.chain_id
        for hname in ('H2', 'H3'):
            h = next_res.find_atom(hname)
            if h is not None:
                h.delete()
            h = next_res.find_atom('H1')
            if h is not None:
                h.name='H'
            if next_res.name == 'PRO':
                h = next_res.find_atom('H')
                if h:
                    h.delete()
    else:
        b_factor = add_b_factor
        ref_res = None
    if number is None:
        if prev_res:
            number = prev_res.number + 1
        elif next_res:
            number = next_res.number - 1

    from chimerax import mmcif
    tmpl = mmcif.find_template_residue(session, resname)
    from .place_ligand import new_residue_from_template
    import numpy
    # delete extraneous atoms
    r = new_residue_from_template(model, tmpl, chain_id, [0,0,0], number,
            insert_code=insertion_code, b_factor=b_factor, precedes=insertion_point)
    r.atoms[numpy.in1d(r.atoms.names, ['OXT', 'HXT', 'H2', 'H1', 'HN1', 'HN2'])].delete()

    # Translate and rotate residue to (roughly) match the desired position
    if not next_res and not prev_res:
        r.atoms.coords += numpy.array(center) - r.atoms.coords.mean(axis=0)
        b_factor = max(add_b_factor, 5)
    else:
        from chimerax.geometry import align_points
        from chimerax.atomic.struct_edit import add_bond
        if prev_res:
            correct_o_position(prev_res, psi)
            add_bond(r.find_atom('N'), prev_res.find_atom('C'))
            n_pos = _find_next_N_position(prev_res)
            ca_pos = _find_next_CA_position(n_pos, prev_res)
            c_pos = _find_next_C_position(ca_pos, n_pos, prev_res, phi)
            target_coords = numpy.array([n_pos, ca_pos, c_pos])
            align_coords = numpy.array([r.find_atom(a).coord for a in ['N', 'CA', 'C']])
        elif next_res:
            add_bond(r.find_atom('C'), next_res.find_atom('N'))
            c_pos = _find_prev_C_position(next_res, psi)
            ca_pos = _find_prev_CA_position(c_pos, next_res)
            #o_pos = _find_prev_O_position(c_pos, next_res)
            n_pos = _find_prev_N_position(c_pos, ca_pos, next_res, psi)
            target_coords = numpy.array([c_pos, ca_pos, n_pos])
            align_coords = numpy.array([r.find_atom(a).coord for a in ['C', 'CA', 'N']])

        tf = align_points(align_coords, target_coords)[0]
        r.atoms.coords = tf*r.atoms.coords
        correct_o_position(r, psi)
    if r.name in ('GLU', 'ASP'):
        fix_amino_acid_protonation_state(r)
    if r.name == 'PRO':
        r.atoms[r.atoms.names=='H'].delete()

    r.atoms.bfactors = b_factor
    r.atoms.occupancies = occupancy

    if ref_res is not None:
        from . import copy_atom_style_from
        copy_atom_style_from(session, r.atoms, ref_res)
    model.atoms.selecteds=False
    r.atoms.selecteds=True
    return r

def correct_o_position(res, psi):
    n, ca, c, o = [res.find_atom(name) for name in ['N', 'CA', 'C', 'O']]
    from chimerax.geometry import rotation, dihedral
    d = dihedral(*[a.coord for a in (n, ca, c, o)])
    r = rotation(c.coord-ca.coord, psi+180-d, center=c.coord)
    from chimerax.atomic import Atoms
    Atoms([o]).transform(r)

def _find_next_N_position(prev_res):
    from chimerax.atomic.struct_edit import find_pt
    bond_length = 1.34
    angle = 120
    dihedral = 180
    a1 = prev_res.find_atom('C')
    a2 = prev_res.find_atom('CA')
    a3 = prev_res.find_atom('O')
    return find_pt(*[a.coord for a in [a1, a2, a3]], bond_length, angle, dihedral)

def _find_next_CA_position(n_pos, prev_res):
    from chimerax.atomic.struct_edit import find_pt
    bond_length = 1.48
    angle = 124
    omega = 180
    c = prev_res.find_atom('C')
    ca = prev_res.find_atom('CA')
    return find_pt(n_pos, *[a.coord for a in [c, ca]], bond_length, angle, omega)

def _find_next_C_position(ca_pos, n_pos, prev_res, phi):
    from chimerax.atomic.struct_edit import find_pt
    bond_length = 1.53
    angle = 120
    c = prev_res.find_atom('C')
    return find_pt(ca_pos, n_pos, c.coord, bond_length, angle, phi)


def _find_prev_C_position(next_res, psi):
    from chimerax.atomic.struct_edit import find_pt
    bond_length = 1.34
    angle = 120
    a1 = next_res.find_atom('N')
    a2 = next_res.find_atom('CA')
    a3 = next_res.find_atom('C')
    return find_pt(*[a.coord for a in [a1, a2, a3]], bond_length, angle, psi)

def _find_prev_CA_position(c_pos, next_res):
    from chimerax.atomic.struct_edit import find_pt
    bond_length = 1.53
    angle = 120
    omega = 180
    n = next_res.find_atom('N')
    ca = next_res.find_atom('CA')
    return find_pt(c_pos, *[a.coord for a in [n, ca]], bond_length, angle, omega)

def _find_prev_O_position(c_pos, next_res):
    from chimerax.atomic.struct_edit import find_pt
    bond_length = 1.22
    angle = 120
    dihedral = 0
    n = next_res.find_atom('N')
    ca = next_res.find_atom('CA')
    return find_pt(c_pos, *[a.coord for a in (n, ca)], bond_length, angle, dihedral)

def _find_prev_N_position(c_pos, ca_pos, next_res, psi):
    from chimerax.atomic.struct_edit import find_pt
    nn = next_res.find_atom('N')
    bond_length = 1.44
    angle = 120
    return find_pt(ca_pos, c_pos, nn.coord, bond_length, angle, psi)

def current_and_possible_disulfides(model, cutoff_distance=3.0):
    import numpy
    from chimerax.atomic import Atoms
    atoms = model.atoms
    cys_s = atoms[numpy.logical_and(atoms.residues.names=='CYS', atoms.names=='SG')]
    from chimerax.atom_search import AtomSearchTree
    tree = AtomSearchTree(cys_s)
    current_disulfides = set()
    possible_disulfides = set()
    ambiguous = set()
    for s in cys_s:
        nearby = tree.search(s.coord, cutoff_distance)
        # Search will always find the input atom, so want results with at least
        # (hopefully only!) two atoms.
        if len(nearby) > 1:
            if len(nearby) > 2:
                ambiguous.add(frozenset([a.residue for a in nearby]))
                continue
            others = [a for a in nearby if a != s]
            for a in others:
                if a in s.neighbors:
                    current_disulfides.add(frozenset((s.residue,a.residue)))
                else:
                    possible_disulfides.add(frozenset((s.residue,a.residue)))
    return current_disulfides, possible_disulfides, ambiguous

def create_disulfide(cys1, cys2):
    from chimerax.core.errors import UserError
    if cys1.structure != cys2.structure:
        raise UserError('Both cysteine residues must be in the same model!')
    structure = cys1.structure
    s1 = cys1.find_atom('SG')
    s2 = cys2.find_atom('SG')
    if s1 is None or s2 is None:
        raise UserError('Missing SG atom! Are both residues complete cysteines?')
    if s2 in s1.neighbors:
        raise UserError('These residues are already disulfide bonded!')
    for s in (s1, s2):
        for a in s.neighbors:
            if a.element.name == 'H':
                a.delete()
    from chimerax.atomic.struct_edit import add_bond
    add_bond(s1, s2)

_CYS_ALIGN_ATOMS=('CA', 'CB', 'SG')
def break_disulfide(cys1, cys2):
    from chimerax.core.errors import UserError
    from chimerax.atomic import Atoms
    s1 = cys1.find_atom('SG')
    s2 = cys2.find_atom('SG')
    if s1 is None or s2 is None:
        raise UserError('Missing SG atom! Are both residues complete cysteines?')
    if s2 not in s1.neighbors:
        raise UserError('These residues are not disulfide bonded!')
    has_hydrogens = ('H' in cys1.atoms.element_names)
    b = Atoms((s1, s2)).intra_bonds[0]
    b.delete()
    if has_hydrogens:
        from chimerax.build_structure import modify_atom
        modify_atom(s1, s1.element, 2, res_name = 'CYS', connect_back=False)
        modify_atom(s2, s2.element, 2, res_name = 'CYS', connect_back=False)
