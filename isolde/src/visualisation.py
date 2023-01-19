# @Author: Tristan Croll <tic20>
# @Date:   20-Apr-2018
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tic20
# @Last modified time: 23-May-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright:2016-2019 Tristan Croll



from enum import IntEnum
class map_styles(IntEnum):
    mesh_square             = 0
    mesh_triangle           = 1
    solid_t20               = 2
    solid_t40               = 3
    solid_t60               = 4
    solid_t80               = 5
    solid_opaque            = 6

human_readable_map_display_styles = {
    map_styles.mesh_square: "Mesh (squares)",
    map_styles.mesh_triangle: "Mesh (triangles)",
    map_styles.solid_t20: "Solid (20% opacity)",
    map_styles.solid_t40: "Solid (40% opacity)",
    map_styles.solid_t60: "Solid (60% opacity)",
    map_styles.solid_t80: "Solid (80% opacity)",
    map_styles.solid_opaque: "Solid (opaque)"
    }

# array of settings to apply to the map depending on the chosen
# representation.
map_style_settings = {
    map_styles.mesh_square: {'style': 'mesh', 'square_mesh': True, 'transparency':0, 'flip_normals': True},
    map_styles.mesh_triangle: {'style': 'mesh', 'square_mesh': False, 'transparency':0, 'flip_normals': True},
    map_styles.solid_t20: {'style': 'surface', 'transparency': 0.8},
    map_styles.solid_t40: {'style': 'surface', 'transparency': 0.6},
    map_styles.solid_t60: {'style': 'surface', 'transparency': 0.4},
    map_styles.solid_t80: {'style': 'surface', 'transparency': 0.2},
    map_styles.solid_opaque: {'style': 'surface', 'transparency': 0.0}
    }

def default_atom_visualisation(model):
    session = model.session
    from chimerax.core.objects import Objects
    from chimerax.std_commands import color
    from chimerax.atomic import Atom
    atoms = model.atoms
    ao = Objects(atoms=atoms)
    color.color(session, ao, color='bychain')
    color.color(session, ao, color='byhetero')
    model.bonds.radii = 0.2
    model.bonds.displays=True
    atoms.draw_modes = Atom.STICK_STYLE
    atoms.displays = True
    from chimerax.clipper.util import nonpolar_hydrogens
    atoms[nonpolar_hydrogens(atoms)].displays=False
    model.residues.ribbon_displays = True

class VisualisationStateMgr:
    '''
    Use to modify the visualisation of a structure, storing the current visualisation state 
    for later reversion.
    '''
    IGNORED_BASE_COLOR = (0.75,0.75,0.75)
    def __init__(self, session, isolde):
        self.session = session
        self.isolde = isolde
        isolde.triggers.add_handler(isolde.SIMULATION_STARTED, self._sim_start_cb)
        isolde.triggers.add_handler(isolde.SIMULATION_TERMINATED, self._sim_end_cb)
                
    def _sim_start_cb(self, *_):
        m = self._current_model = self.isolde.selected_model
        self.store_original_visualisation()
        self.prepare_sim_visualisation()
    
    def _sim_end_cb(self, trigger_name, reason):
        m = self.isolde.selected_model
        from chimerax.isolde.openmm.openmm_interface import SimHandler
        if reason == SimHandler.REASON_MODEL_DELETED or m != self._current_model:
            return
        self.revert_visualisation()


    def store_original_visualisation(self):
        '''
        Store the current visualisation state of the model, so it can be
        reverted once the simulation is done.
        '''
        m = self._current_model
        atoms = m.atoms
        bonds = m.bonds
        residues = m.residues
        self.spotlight_mode = None
        from chimerax.clipper.symmetry import get_symmetry_handler
        sym = self.symmetry_handler = get_symmetry_handler(m)
        if sym:
            self.spotlight_mode = sym.spotlight_mode
        self._original_atom_states = (
            atoms.colors,
            atoms.draw_modes,
            atoms.displays,
            atoms.radii
            )
        self._original_bond_states = (
            bonds.displays,
            bonds.radii,
        )
        self._original_residue_states = (
            residues.ribbon_displays,
            residues.ribbon_colors
        )

    def revert_visualisation(self):
        '''
        Return the model visualisation to the way it was before the simulation
        began.
        '''
        if not hasattr(self, '_original_atom_states'):
            from chimerax.core.errors import UserError
            raise UserError('No stored state to revert!')
        m = self._current_model
        atoms = m.atoms
        bonds = m.bonds
        residues = m.residues
        # If atoms are added or deleted while a simulation is running, the
        # stored settings will become invalid. In that case, just revert to a
        # default visualisation.
        try:
            atoms.colors, atoms.draw_modes, atoms.displays, atoms.radii = \
                self._original_atom_states
            bonds.displays, bonds.radii = \
                self._original_bond_states
            residues.ribbon_displays, residues.ribbon_colors = \
                self._original_residue_states
        except:
            default_atom_visualisation(self.model)

        sym = self.symmetry_handler
        if sym:
            sym.spotlight_mode = self.spotlight_mode

    def prepare_sim_visualisation(self):
        m = self._current_model
        sc = self.isolde.sim_manager.sim_construct
        m.residues.ribbon_displays = False
        fixed_bonds = sc.fixed_atoms.intra_bonds
        fixed_bonds.radii *= self.isolde.params.fixed_bond_radius_ratio
        from .constants import control
        sc.all_atoms.hides &= ~control.HIDE_ISOLDE
        ea = sc.excluded_atoms
        if ea is not None:
            import numpy
            colors = ea.colors.astype(float)/255
            colors[:,:3] = (colors[:,:3]*0.3+numpy.array(self.IGNORED_BASE_COLOR)*0.7)
            colors = (colors*255).astype(numpy.uint8)
            ea.colors = colors
