import numpy
from chimerax.core.atomic import Residue
def is_continuous_protein_chain(sel):
    '''
    Checks if the residues in a selection are all protein, and form a 
    continuous chain with no breaks. 
    NOTE: only one atom from each residue need be selected.
    '''
    try:
        p = find_polymer(sel)
    except:
        return False
    res = sel.unique_residues
    if not numpy.all(res.polymer_types == Residue.PT_AMINO):
        return False
        
    indices = sorted(p.indices(res))
    first = indices[0]
    # Check if the indices are continuous
    return all(a == b for a, b in enumerate(indices, indices[0]))    

def find_polymer(sel):
    '''
    Find the polymer containing the first residue in a selection.
    '''
    us = sel.unique_structures
    if len(us) != 1:
        raise TypeError('All atoms must be from the same structure!')
    m = us[0]
    res = sel.unique_residues
    polymers = m.polymers(
        missing_structure_treatment = m.PMS_NEVER_CONNECTS)
    first_index = -1
    r0 = res[0]
    for p in polymers:
        first_index = p.index(r0)
        if first_index != -1:
            return p
    raise IndexError('Polymer not found!')
            
def atom_residues(atoms):

    rt = {}
    for a in atoms:
        rt[a.residue] = 1
    rlist = list(rt.keys())
    return rlist
    
    
def molecule_from_atoms(m, atoms, name = None, bonds = True, ribbons = True, pseudobonds = True):

    from chimerax.core.atomic import AtomicStructure, Structure
    structure_class = AtomicStructure if isinstance(m, AtomicStructure) else Structure
    cm = structure_class(m.session, name = (name or m.name), auto_style = False)
    cm.ss_assigned = True
#    cm.color = m.color
    cm.display = m.display
#    cm.lineWidth = m.lineWidth
#    cm.pointSize = m.pointSize
#    cm.ballScale = m.ballScale

#    cm.pdbVersion = m.pdbVersion
#    if hasattr(m, 'pdbHeaders'):
#        cm.setAllPDBHeaders(m.pdbHeaders)
#    if hasattr(m, 'mmCIFHeaders'):
#        cm.mmCIFHeaders = m.mmCIFHeaders

    rmap = {}
    rlist = atom_residues(atoms)
    rorder = dict((r,i) for i,r in enumerate(m.residues))
    rlist.sort(key = lambda r: rorder[r])
    for r in rlist:
        cr = cm.new_residue(r.name, r.chain_id, r.number)
#        cr.isHet = r.isHet
        cr.is_helix = r.is_helix
        cr.is_strand = r.is_strand
        cr.ribbon_color = r.ribbon_color
#        cr.ribbonStyle = r.ribbonStyle
#        cr.ribbonDrawMode = r.ribbonDrawMode
        if ribbons:
            cr.ribbon_display = r.ribbon_display
        else:
            cr.ribbon_display = False
        rmap[r] = cr

    amap = {}
    for a in atoms:
        ca = cm.new_atom(a.name, a.element_name)
        ca.coord = a.coord
#        ca.altLoc = a.altLoc
        ca.color = a.color
        ca.draw_mode = a.draw_mode
        ca.display = a.display
        ca.bfactor = a.bfactor
        amap[a] = ca
        cr = rmap[a.residue]
        cr.add_atom(ca)
    
    if bonds:
        for b in atom_bonds(atoms):
            a1, a2 = b.atoms
            cb = cm.new_bond(amap[a1], amap[a2])
            cb.color = b.color
            cb.radius = b.radius
    #        cb.drawMode = b.drawMode
            cb.display = b.display
            cb.halfbond = b.halfbond
    
    if pseudobonds:
        for name, pbg in m.pbg_map.items():
            cpbgs = {}
            for pb in pbg.pseudobonds:
                a1, a2 = pb.atoms
                if a1 not in amap or a2 not in amap:
                    continue
                cpbg = cpbgs.get(name)
                if cpbg is None:
                    cpbgs[name] = cpbg = cm.pseudobond_group(name)
                cpb = cpbg.new_pseudobond(amap[a1],amap[a2])
                cpb.display = pb.display
                cpb.color = pb.color
                cpb.radius = pb.radius
                cpb.halfbond = pb.halfbond

    cm.new_atoms()

    return cm
