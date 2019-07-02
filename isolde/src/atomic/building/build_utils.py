# @Author: Tristan Croll <tic20>
# @Date:   11-Jun-2019
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 11-Jun-2019
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll



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

    na = s.new_atom(name, Element.get_element('H'))
    na.coord = coord
    na.bfactor = atom.bfactor
    na.aniso_u6 = atom.aniso_u6
    na.draw_mode = atom.draw_mode
    r.add_atom(na)
    s.new_bond(atom, na)
    return na

def current_and_possible_disulfides(model, cutoff_distance=3.0):
    import numpy
    from chimerax.atomic import Atoms
    atoms = model.atoms
    cys_s = atoms[numpy.logical_and(atoms.residues.names=='CYS', atoms.names=='SG')]
    from chimerax.atomic.search import AtomSearchTree
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
    structure.new_bond(s1, s2)

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
        m = cys1.structure
        import numpy
        from chimerax.core.geometry import align_points, rotation, angle
        from chimerax.atomic import TmplResidue
        from math import copysign
        templ = TmplResidue.get_template('CYS')
        tatom = templ.find_atom('HG')
        align_coords = numpy.array([templ.find_atom(a).coord for a in _CYS_ALIGN_ATOMS])
        for (s, cys, s_prime) in ((s1, cys1, s2), (s2, cys2, s1)):
            target_coords = numpy.array([cys.find_atom(a).coord for a in _CYS_ALIGN_ATOMS])
            p, _ = align_points(align_coords, target_coords)
            # Adding the hydrogens straight from the template can leave them
            # essentially on top of each other. We need to check the H-SG-SG'
            # angle, and if it's too tight rotate the H around the CB-SG bond
            h = m.new_atom(tatom.name, tatom.element)
            h.coord = p*tatom.coord
            cys.add_atom(h)
            m.new_bond(s, h)
            ang = angle(h.coord, s.coord, s_prime.coord)
            if abs(ang) < 30:
                r_angle = copysign(90, ang)-ang
                r = rotation(target_coords[2]-target_coords[1], r_angle, center=target_coords[2])
                h.coord = r*h.coord
            h.draw_mode = h.STICK_STYLE
