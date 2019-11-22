



def remove_nucleic_terminal_h(model):
    from chimerax.atomic import Residue
    nucleic = model.residues[model.residues.polymer_types==Residue.PT_NUCLEIC]
    if not len(nucleic):
        return

    p = nucleic.atoms[nucleic.atoms.names=='P']
    for phos in p:
        for n in phos.neighbors:
            if n.element.name=='H':
                n.delete()
