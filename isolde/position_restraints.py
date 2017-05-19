import numpy
from chimerax.core.atomic import Residue, Atoms, AtomicStructure
from .restraints_base import Position_Restraint, Position_Restraints   

MAX_RESTRAINT_FORCE = 100.0 # kJ/mol/A

class Atom_Position_Restraints(Position_Restraints):
    '''
    For restraining of atoms to specific points in space.
    '''
    def __init__(self, atoms_or_list, include_hydrogens = False, create_target = False):
        '''
        Initialises restraint objects for the atoms in the given selection.
        By default, hydrogens are ignored.
        '''
        if type(atoms_or_list) == Atoms:
            if not include_hydrogens:
                atoms = atoms_or_list
                atoms = atoms.filter(atoms.element_names != 'H')
            restraints = []
            for a in atoms:
                restraints.append(Position_Restraint(a))
        else:
            restraints = atoms_or_list
        super().__init__(restraints)
        
        
        # Create a minimal AtomicStructure containing only atom information
        # to act as anchor points for Pseudobonds representing the active
        # restraints.
        from chimerax.core.commands.split import molecule_from_atoms 
        master_model = self.atoms.unique_structures[0]     
        session = master_model.session  
        from .util import molecule_from_atoms


        # In order to allow pseudobonds between atoms in different 
        # AtomicStructure objects, the Pseudobondgroup must be in the 
        # global manager.
        pbg = self._pseudobondgroup = session.pb_manager.get_group(master_model.name + 'restraint pseudobonds')
        if pbg not in session.models.list():
            master_model.add([pbg])
        
        if create_target:
            t = self._target_model = molecule_from_atoms(
                master_model, self.atoms, name = 'xyz restraint targets', 
                ribbons = False, bonds = False, pseudobonds = False)
            t.atoms.displays = False
            t.atoms.draw_modes = t.atoms.BALL_STYLE
            t.atoms.radii = t.atoms.radii * 0.5
            master_model.add([t])

        
            # Create a Pseudobond for each restrainable atom, connecting it
            # to its target. This will be automatically displayed if the 
            # target atom is displayed. 
            for r, ta in zip(self, t.atoms):
                r._target_atom = ta
                pb = r._pseudobond = pbg.new_pseudobond(r.atom, ta)
        
        else:
            ta = self._restraints[0]._target_atom
            if ta is not None:
                t = self._target_model = ta.structure
            else:
                t = self._target_model = None
        
        
        
        # Pseudobond radii will scale according to the force they're applying,
        # with set minimum and maximum radii.
        self.min_pb_radius = 0.05
        self.max_pb_radius = 0.25
        
        self._min_force = 0
        self._max_force = MAX_RESTRAINT_FORCE # from kJ/mol/A to kJ/mol/nm

        # Define a color map to indicate how well the restraint is satisfied
        from . import color
        self._cmap = color.standard_three_color_scale('GWR', self._min_force, self._max_force)

        
        self._pb_radii = numpy.array([self.min_pb_radius]*len(self), numpy.float32)
        self._pb_colors = numpy.zeros([4, len(self)], numpy.uint8)
    
    def update_pseudobonds(self, forces, indices):
        '''
        Update the radii and the colors of the restraint pseudobonds
        according to the forces being applied to their target atoms.
        '''
        if self._target_model is not None:
            forces = forces.astype(numpy.float32)
            maxf = self._max_force
            minf = self._min_force
            maxr = self.max_pb_radius
            minr = self.min_pb_radius
            radii = (forces - minf) / (maxf-minf) * (maxr-minr) + minr
            radii[radii < minr] = minr
            radii[radii > maxr] = maxr
            
            colors = self._cmap.get_colors(forces)
            
            self._pseudobondgroup.pseudobonds.radii[indices] = radii
            self._pseudobondgroup.pseudobonds.colors[indices] = colors
    
    
    
            
        
