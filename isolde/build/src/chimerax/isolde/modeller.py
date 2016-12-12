import numpy
from simtk.openmm.app import modeller
from simtk.openmm import unit
from . import sim_interface as sh

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
    
            
        
    
        
    
    
    
    

    
    
