from chimerax.core.settings import Settings
from .constants import defaults
from .validation.constants import rotarama_defaults
from .restraints.constants import restraint_color_defaults
class _IsoldeBasicSettings(Settings):
    AUTO_SAVE = {
        'num_selection_padding_residues':   defaults.SELECTION_SEQUENCE_PADDING,
        'soft_shell_cutoff_distance':       defaults.SOFT_SHELL_CUTOFF,
        'hard_shell_cutoff_distance':       defaults.HARD_SHELL_CUTOFF,
        'fixed_bond_radius_ratio':          defaults.FIXED_BOND_RADIUS_RATIO,
        'hide_surroundings_during_sim':     defaults.HIDE_SURROUNDINGS_DURING_SIM,
        #'remask_maps_during_sim':           defaults.REMASK_MAPS_DURING_SIM,
        'map_mask_radius':                  defaults.STANDARD_MAP_MASK_RADIUS,
        'map_oversampling_rate':            defaults.MAP_SHANNON_RATE,
        'selection_outline_width':          defaults.SELECTION_OUTLINE_WIDTH,

        'sim_fidelity_mode':                defaults.SIM_FIDELITY_MODE,
        'trajectory_smoothing':             defaults.TRAJECTORY_SMOOTHING,
        'smoothing_alpha':                  defaults.SMOOTHING_ALPHA,

        'phenix_base_path':                 None,
    }

class _IsoldeColorSettings(Settings):
    AUTO_SAVE = {
        'rama_color_scale': [
            rotarama_defaults.MAX_FAVORED_COLOR,
            rotarama_defaults.ALLOWED_COLOR,
            rotarama_defaults.OUTLIER_COLOR,
            rotarama_defaults.NA_COLOR],
        'rota_color_scale': [
            rotarama_defaults.MAX_FAVORED_COLOR,
            rotarama_defaults.ALLOWED_COLOR,
            rotarama_defaults.OUTLIER_COLOR],
        'proper_dihedral_restraint_color_scale': [
            restraint_color_defaults.PROPER_DIHEDRAL_RESTRAINT_SATISFIED_COLOR,
            restraint_color_defaults.PROPER_DIHEDRAL_RESTRAINT_STRAINED_COLOR,
            restraint_color_defaults.PROPER_DIHEDRAL_RESTRAINT_SEVERE_COLOR],
        'distance_restraint_colors': [
            restraint_color_defaults.DISTANCE_RESTRAINT_BOND_COLOR,
            restraint_color_defaults.DISTANCE_RESTRAINT_TARGET_COLOR],
        'adaptive_dihedral_restraints_color_scale': [
            restraint_color_defaults.ADAPTIVE_DIHEDRAL_RESTRAINT_SATISFIED_COLOR,
            restraint_color_defaults.ADAPTIVE_DIHEDRAL_RESTRAINT_STRAINED_COLOR,
            restraint_color_defaults.ADAPTIVE_DIHEDRAL_RESTRAINT_SEVERE_COLOR],
        'adaptive_distance_restraints_color_scale': [
            restraint_color_defaults.ADAPTIVE_DISTANCE_RESTRAINT_SATISFIED_COLOR,
            restraint_color_defaults.ADAPTIVE_DISTANCE_RESTRAINT_TOO_CLOSE_COLOR,
            restraint_color_defaults.ADAPTIVE_DISTANCE_RESTRAINT_TOO_FAR_COLOR],
        'position_restraint_colors': [
            restraint_color_defaults.POSITION_RESTRAINT_PIN_COLOR,
            restraint_color_defaults.POSITION_RESTRAINT_BOND_COLOR],
        'tugging_arrow_color':  restraint_color_defaults.TUGGING_ARROW_COLOR,    

    }

class _IsoldeAdvancedSettings(Settings):
    AUTO_SAVE = {
        'openmm_default_platform':          defaults.OPENMM_DEFAULT_PLATFORM,
        'openmm_gpu_device_index':          defaults.DEVICE_INDEX,

    }

def register_settings_options(session):
    from chimerax.ui.options import (
        ColorOption, BooleanOption, IntOption, FloatOption,
    )

basic_settings = None # set during bundle initialisation
color_settings = None # set during bundle initialisation
advanced_settings = None # set during bundle initialisation