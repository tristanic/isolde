# Copyright 2017 Tristan Croll
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


import numpy
from simtk.openmm.app import modeller
from simtk.openmm import unit

def add_hydrogens_to_model(model):
    top = sh.openmm_topology_from_model(model)
    a = model.atoms
    b = model.bonds
    
    builder = modeller.Modeller(top, a.coords/10)
    
    # Add missing hydrogens to the model, using default settings
    builder.addHydrogens()
    
    # Now builder.topology has been extended to include the hydrogens,
    # and builder.positions has coordinates.

    # Hydrogens are added interleaved between heavy atoms
    index_offset = 0
    
    newtop = builder.topology
    positions = builder.positions / unit.angstrom
    
    heavy_atom_map = {} # key: index in new topology. Value: index in old 
                        #       topology (and in model)
    top_hydrogens = {}
    model_hydrogens = {}
    
    
    for i, atom in enumerate(newtop.atoms()):
        if atom.element.symbol == 'H':
            index_offset += 1
            top_hydrogens[i] = atom
            na = model.new_atom(atom.name, 'H')
            na.coord = positions[i]
            model_hydrogens[i] = na
        else:
            oi = i - index_offset
            heavy_atom_map[i] = oi
            a[oi].coord = positions[i]
            
    
    
    for atom0, atom1 in newtop.bonds():
        if atom1.element.symbol == 'H':
            m_atom0 = model.atoms[heavy_atom_map[atom0.index]]
            m_atom1 = model_hydrogens[atom1.index]
            res = m_atom0.residue
            res.add_atom(m_atom1)
            model.new_bond(m_atom0, m_atom1)
    
            
def openmm_topology_from_model(model):
    '''
    Take an AtomicStructure model from ChimeraX and return an OpenMM
    topology (e.g. for use with the OpenMM Modeller class).
    '''
    a = model.atoms
    b = model.bonds
    n = len(a)
    r = a.residues
    aname = a.names
    ename = a.element_names
    rname = r.names
    rnum = r.numbers
    cids = r.chain_ids
    from simtk.openmm.app import Topology, Element
    from simtk import unit
    top = Topology()
    cmap = {}
    rmap = {}
    atoms = {}
    for i in range(n):
        cid = cids[i]
        if not cid in cmap:
            cmap[cid] = top.addChain()   # OpenMM chains have no name
        rid = (rname[i], rnum[i], cid)
        if not rid in rmap:
            rmap[rid] = top.addResidue(rname[i], cmap[cid])
        element = Element.getBySymbol(ename[i])
        atoms[i] = top.addAtom(aname[i], element,rmap[rid])

    a1, a2 = b.atoms
    for i1, i2 in zip(a.indices(a1), a.indices(a2)):
        if -1 not in [i1, i2]:
            top.addBond(atoms[i1],  atoms[i2])
            
    return top
        
    
        
    
    
    
    

    
    
