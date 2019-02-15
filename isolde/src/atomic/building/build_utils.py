
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
            last_digits = [int(n[-1]) for n in names if n[-1].isdigit()]
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
