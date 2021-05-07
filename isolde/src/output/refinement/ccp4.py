def refmac_distance_restraint(atom1, atom2, distance=None, sigma=0.1):
    def atom_descriptor(atom):
        residue = atom.residue
        insertion_code = residue.insertion_code
        if insertion_code=='':
            insertion_code='.'
        return 'chain {} resi {} ins {} atom  {:<3}'.format(residue.chain_id, residue.number,
            insertion_code, atom.name)
    if distance is None:
        from chimerax.geometry import distance as gdist
        distance = gdist(atom1.coord, atom2.coord)
    restraint_string = ('exte dist first {} second {} value {:.5f} sigma {:.1f}').format(
        atom_descriptor(atom1), atom_descriptor(atom2), distance, sigma)
    return restraint_string

def refmac_distance_restraints(session, model, distance_cutoff=4.5, include_waters=False,
        file_name='RESTRAINTS.txt'):
    import numpy
    m = model
    from chimerax.atomic import AtomicStructures
    if isinstance(m, AtomicStructures):
        if len(m) != 1:
            from chimerax.core.errors import UserError
            raise UserError('Please specify a single atomic model!')
        m = m[0]
    residues = m.residues
    if not include_waters:
        residues = residues[residues.names !='HOH']
    atoms = residues.atoms[residues.atoms.element_names!='H']
    coords = atoms.coords
    from chimerax.geometry import find_close_points
    seen = set()
    with open(file_name, 'wt') as rfile:
        rfile.write(
            '# ISOLDE Restraints File\n'
            '# \n'
            '# Restraints to ISOLDE output geometry\n'
        )
        for i, atom in enumerate(atoms):
            query_coord = numpy.array([coords[i]])
            indices = find_close_points(query_coord, coords, distance_cutoff)[1]
            for ind in indices:
                atom2 = atoms[ind]
                if atom2 == atom:
                    continue
                # Do not include restraints for neighbors or 1-3 relationships
                if any([atom2==n or atom2 in n.neighbors for n in atom.neighbors]):
                    continue
                # Don't double-count
                pair = frozenset((atom, atom2))
                if pair in seen:
                    continue
                rfile.write(refmac_distance_restraint(atom, atom2)+'\n')
                seen.add(pair)
    import os
    session.logger.info(
        f'Top-out distance restraints file for REFMAC5 written to {file_name}. '
        f'This is essentially equivalent to a ProSMART restraints file, restraining '
        f'interatomic distances to their current values. Use it in the "External restraints" '
        f'section of a Refmac5 job in the CCP-EM GUI, or at the command line as: \n'
        f'refmac5 {{all other arguments}} < {file_name}'
    )


def register_ccp4_commands(logger):
    from chimerax.core.commands import (
        register, CmdDesc,
        FileNameArg, BoolArg, FloatArg
    )
    from chimerax.atomic import AtomicStructuresArg             

    cmd_desc = CmdDesc(
        synopsis = "Write REFMAC input to restrain a model to its current geometry",
        required = [
            ('model', AtomicStructuresArg),
        ],
        keyword = [
            ('distance_cutoff', FloatArg),
            ('include_waters', BoolArg),
            ('file_name', FileNameArg)
        ]
    )
    register('isolde write refmacRestraints', cmd_desc, refmac_distance_restraints, logger=logger)

