import numpy
from chimerax.core.atomic import Residue, Atoms, AtomicStructure, Structure
from .restraints_base import Position_Restraint, Position_Restraints   
from .restraints_base import HIDE_ISOLDE

MAX_RESTRAINT_FORCE = 2000.0 # kJ/mol/A


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
        master_model = self.master_model = self.atoms.unique_structures[0]     
        session = self.session = master_model.session  

        self._master_display_change_handler = None
        self._master_delete_handler = None
        
        # In order to allow pseudobonds between atoms in different 
        # AtomicStructure objects, the Pseudobondgroup must be in the 
        # global manager.
        pbg = self._pseudobondgroup = session.pb_manager.get_group(master_model.name + 'restraint pseudobonds')
        if pbg not in session.models.list():
            master_model.add([pbg])
        
        if create_target:
            t = self._target_model = restraint_indicator_from_atoms(
                master_model, self.atoms, name = 'xyz restraint targets')
            # Atoms with hide == True will override display == True
            t.atoms.hides |= HIDE_ISOLDE
            #t.atoms.draw_modes = t.atoms.BALL_STYLE
            #t.atoms.radii = t.atoms.radii * 0.5
            master_model.add([t])

        
            # Create a Pseudobond for each restrainable atom, connecting it
            # to its target. This will be automatically displayed if the 
            # target atom is displayed. 
            for r, ta in zip(self, t.atoms):
                r._target_atom = ta
                pb = r._pseudobond = pbg.new_pseudobond(r.atom, ta)
                pb.display = False
            
            from chimerax.core.atomic import get_triggers
            t = get_triggers(self.session)
            self._master_display_change_handler = t.add_handler(
                'changes', self._update_target_display)
            self._master_delete_handler = session.triggers.add_handler(
                'remove models', self._check_if_master_deleted)
     
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
    
    @property
    def hidden_indicators(self):
        return self._target_model.atoms[self.spring_constants == 0]
    
    def _update_target_display(self, trigger_name, changes):
        if 'display changed' in changes.atom_reasons():
            if self.master_model in changes.modified_atoms().unique_structures:
                self.target_indicators.displays = self.atoms.displays
        #~ if 'hide changed' in changes.atom_reasons():
            #~ if self._target_model in changes.modified_atoms().unique_structures:
                #~ self.hidden_indicators.hides = True
    
    def _check_if_master_deleted(self, trigger_name, models):
        if self.master_model in models:
            self._remove_handlers()
            self.master_model = None
    
    def _remove_handlers(self):
        if self._master_display_change_handler is not None:
            from chimerax.core.atomic import get_triggers
            t = get_triggers(self.session)
            t.remove_handler(self._master_display_change_handler)
            self._master_display_change_handler = None
        if self._master_delete_handler is not None:
            self.session.triggers.remove_handler(self._master_delete_handler)
            self._master_delete_handler = None
    
    
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
            pb = self._pseudobondgroup.pseudobonds
            
            all_radii = pb.radii
            all_radii[indices] = radii
            pb.radii = all_radii
            all_colors = pb.colors
            all_colors[indices] = colors
            pb.colors = all_colors
    
    def cleanup(self):
        self._remove_handlers()
    
    
class PositionRestraintIndicator(Structure):
    '''
    For each position restraint, we want a neat way to display its target
    position with a bond back to the atom itself. Redrawing needs to be
    fast, since at times there will be many restraints in play at once.
    The simplest way to accomplish this is to work with the pre-existing
    infrastructure for fast handling of atoms and bonds, and build the 
    targets as fake Atoms. This allows us to use the existing Pseudobonds
    functionality to draw the bonds between atoms and their targets, giving
    us a lot of functionality for free. However, it does mean we need to
    subclass AtomicStructure, since we need to take all control over 
    visualisation away from the main ChimeraX interface.
    '''
    def __init__(self, session, name = None, auto_style = False):
        super().__init__(session, name = name, auto_style = auto_style)
        self._geometry_set = False
    
    def atomspec_has_atoms(self):
        return False
        
    #def atomspec_atoms(self):
    #    return None
    
    @property
    def ribbon_display_count(self):
        return 0
    
    def _update_level_of_detail(self, total_atoms):
        atoms = self.atoms  # optimzation, avoid making new numpy array from C++
        avis = atoms.visibles
        p = self._atoms_drawing
        if p is None:
            if avis.sum() == 0:
                return
            self._atoms_drawing = p = self.new_drawing('atoms')
            self._atoms_drawing.custom_x3d = self._custom_atom_x3d
            # Update level of detail of spheres
        self.set_atom_pin_geometry(self._atoms_drawing, total_atoms)
        
    def set_atom_pin_geometry(self, drawing, natoms = None):
        if natoms == 0:
            return
        ta = drawing.triangles
        if ta is None or not self._geometry_set:
            from .geometry import pin_drawing
            va, na, ta = pin_drawing(1.0, 0.1, 3.0)
            drawing.vertices = va
            drawing.normals = na
            drawing.triangles = ta
            self._geometry_set = True

def atom_residues(atoms):

    rt = {}
    for a in atoms:
        rt[a.residue] = 1
    rlist = list(rt.keys())
    return rlist
     
def restraint_indicator_from_atoms(m, atoms, name = None):

    cm = PositionRestraintIndicator(m.session, name = (name or m.name), auto_style = False)
    cm.ss_assigned = True
    cm.display = m.display

    rmap = {}
    rlist = atom_residues(atoms)
    rorder = dict((r,i) for i,r in enumerate(m.residues))
    rlist.sort(key = lambda r: rorder[r])
    for r in rlist:
        cr = cm.new_residue(r.name, r.chain_id, r.number)
        cr.is_helix = r.is_helix
        cr.is_strand = r.is_strand
        cr.ribbon_display = False
        rmap[r] = cr

    amap = {}
    for a in atoms:
        ca = cm.new_atom(a.name, a.element_name)
        ca.coord = a.coord
        ca.color = a.color
        ca.draw_mode = a.draw_mode
        ca.display = a.display
        ca.bfactor = a.bfactor
        amap[a] = ca
        cr = rmap[a.residue]
        cr.add_atom(ca)

    cm.new_atoms()
    
    return cm


#~ class CustomPropertyManager:
    #~ '''
    #~ A handler to maintain custom atomic properties for a model.
    #~ '''
    #~ def __init__(self, model, level):
        #~ '''
        #~ Create a container for custom per-atom or per-residue properties.
        #~ Args:
            #~ model:
                #~ A Structure or AtomicStructure
            #~ level:
                #~ 'atom' or 'residue'
        #~ '''
        #~ # A dict with direct 1:1 mapping between reference atom/residue 
        #~ # and index in the property value arrays. Only used to update the 
        #~ # index lookups when the number of atoms and/or residues within the 
        #~ # model changes.
        #~ self._dict_map = {}
        #~ self.model = model
        #~ self._level = level
        #~ self._properties = {}
    
    #~ def add_new_custom_property(self, name, dtype, default_value):
        #~ import numpy

   #~ def _model_changed(self, new_reference_list):
        #~ import numpy
        #~ new_array = numpy.empty(len(new_reference_list), self.dtype)
        #~ new_dict_map = {}
        #~ for i, a in enumerate(new_reference_list):
            #~ try:
                #~ new_array[i] = self.data[self._dict_map[a]]
            #~ except KeyError:
                #~ new_array[i] = self.default_value
            #~ new_dict_map[a] = i
        #~ self.data = data


#~ class CustomProperty:
    #~ '''
    #~ Associates a set of atoms (or residues) with an array of custom values
    #~ or objects.
    #~ '''
    #~ def __init__(self, reference_objects, dtype, default_value):
        #~ ''' Prepare the custom property array with default values. '''
        #~ self._ref_objects = reference_objects
        #~ self._ref_type = type(reference_objects[0])
        #~ self._plural_ref_type = type(reference_objects)
        #~ self.dtype = dtype
        
        #~ import numpy
        #~ data = self.data = numpy.array([default_value]*len(reference_objects), dtype)
        #~ for i, r in enumerate(reference_objects):
            #~ self._dict_map[r] = i
    
    #~ def __getitem__(self, key):
        #~ import numpy
        #~ if not len(self):
            #~ return None
        #~ if isinstance(key,(int, numpy.integer, slice, numpy.ndarray)):
            #~ return self.data[key]
        #~ if isinstance(key, self._ref_type):
            #~ index = self._ref_objects.index(key)
            #~ return self.data[index]
        #~ if isinstance(key, self._plural_ref_type):
            #~ indices = self._ref_objects.indices(key)
            #~ return self.data[indices]
        
        #~ raise TypeError('Unrecognised key!')
    
         
    
    

        
        
        
    
#~ class AtomicStructureWithCustomProperties(AtomicStructure):
    #~ def __init__(self, *args, **kwargs)
        #~ self.custom_atomic_properties = CustomPropertyManager()
        #~ super().__init__(*args, **kwargs)
            
    #~ def add_custom_atomic_property(self, key, dtype, default_val):
        #~ '''
        #~ Add a custom property to all atoms in this structure.
        #~ Args:
            #~ key:
                #~ A string, number or hashable object to be used for recall
                #~ of this property
            #~ dtype:
                #~ A numpy dtype indicating the type of data stored
            #~ default_val:
                #~ The default value to be associated with new atoms. The 
                #~ initial array will be populated with this value.
        #~ '''
        
